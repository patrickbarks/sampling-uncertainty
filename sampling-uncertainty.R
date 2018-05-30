
###### 1. Preliminaries

## load libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(grid)
library(gridExtra)
library(popbio)
library(popdemo)
library(rstan)

## set options for rstan library, and compile relevent models
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## source required functions
source('sampling-uncertainty-functions.R')

## load COMPADRE db
load('dat/COMPADRE_v.X.X.X.RData')

## convert COMPADRE db to tidy tibble
compadre_tidy <- compadre$metadata %>% 
  as_tibble() %>% 
  mutate(matA = lapply(compadre$mat, function(x) x$matA),
         matU = lapply(compadre$mat, function(x) x$matU),
         matF = lapply(compadre$mat, function(x) x$matF),
         matC = lapply(compadre$mat, function(x) x$matC),
         matrixClass = compadre$matrixClass) %>% 
  mutate(SpeciesPop = as.integer(as.factor(paste(SpeciesAuthor, MatrixPopulation))))

## fix typo in A matrix for Eriogonum longifolium (3.420 should be 0.342)
satterthwaite_fix <- which(
  compadre_tidy$SpeciesAuthor == 'Eriogonum_longifolium_var._gnaphalifolium_2' &
    compadre_tidy$MatrixPopulation == 'Unburned' &
    compadre_tidy$MatrixStartYear == 1991
)

compadre_tidy$matA[[satterthwaite_fix]][5,5] <- 0.342





###### 2. Subset COMPADRE to datasets of interest, and calculate life expectancy
# (l0) and variance(log lambda)

## initial subset: herbaceous perennial, unmanipulated, wild, dimension > 2,
# repro/survival processes divided, periodicity == 1, no NA
compadre_tidy_sub <- compadre_tidy %>% 
  filter(OrganismType == 'Herbaceous perennial',
         MatrixSplit == 'Divided',
         MatrixFec == 'Yes',
         MatrixTreatment == 'Unmanipulated',
         MatrixCaptivity == 'W',
         MatrixDimension > 2,
         AnnualPeriodicity == 1,
         SurvivalIssue <= 1.01) %>% 
  mutate(NA_matU = map_lgl(matU, ~ any(is.na(.x)))) %>% 
  mutate(NA_matA = map_lgl(matA, ~ any(is.na(.x)))) %>% 
  filter(NA_matU == FALSE, NA_matU == FALSE)

## subset to populations with at least 3 annual matrices
# matA must be ergodic, matU/matF must not be constant over all years, lambda != 0
dat_replicated <- compadre_tidy_sub %>% 
  filter(MatrixComposite == 'Individual',
         MatrixEndYear - MatrixStartYear == 1) %>%
  mutate(MatrixYears = paste(MatrixStartYear, MatrixEndYear)) %>%
  mutate(lambda = map_dbl(matA, GetLambda)) %>% 
  mutate(ergodic = map_lgl(matA, popdemo::is.matrix_ergodic)) %>% 
  filter(lambda != 0) %>% 
  group_by(SpeciesPop) %>%
  mutate(n_year = length(unique(MatrixYears[ergodic == TRUE])),
         matU_equal = CheckMatsEqual(matU),
         matF_equal = CheckMatsEqual(matF)) %>% 
  ungroup() %>%
  filter(matU_equal == FALSE, matF_equal == FALSE) %>% 
  filter(n_year >= 3)

## calculate mean matU for each population and estimate l0, and ensure perennial
# based on mean matU (survivorship to age 3 > 0)
dat_l0 <- dat_replicated %>% 
  mutate(startLife = map_int(matrixClass,
                             ~ min(which(.x$MatrixClassOrganized == 'active')))) %>% 
  group_by(SpeciesPop, SpeciesAuthor, MatrixPopulation) %>%
  summarize(matU_mean = MeanMat(matU),
            startLife = unique(startLife)) %>% 
  ungroup() %>% 
  mutate(l0 = map2_dbl(matU_mean, startLife, GetLifeExpect)) %>% 
  mutate(perennial = map2_lgl(matU_mean, startLife, CheckPerennial)) %>% 
  filter(perennial == TRUE)

## calculate variance(log lambda) and join with life expectancy, and add col for
# matrix dimension
dat_full <- dat_replicated %>% 
  group_by(SpeciesAccepted, SpeciesAuthor, MatrixPopulation) %>% 
  summarize(var_log_lambda = var(log(lambda[ergodic == TRUE]), na.rm = TRUE),
            n_year = unique(n_year)) %>% 
  ungroup() %>% 
  left_join(dat_l0, by = c('SpeciesAuthor', 'MatrixPopulation')) %>%
  filter(!is.na(l0)) %>% 
  mutate(mat_dim = map_int(matU_mean, ~ nrow(.x))) %>% 
  mutate(mat_dim = ifelse(mat_dim > 7, '8+', mat_dim)) %>% 
  mutate(mat_dim = paste('matrix_dim =', mat_dim)) %>% 
  mutate(spp_int = as.integer(as.factor(SpeciesAccepted)))

# ensure relationship between l0 and var(log lambda) not just function of matrix dim
p0 <- ggplot(dat_full, aes(l0, var_log_lambda)) +
  geom_point(shape = 1, size = 3, stroke = 0.7, alpha = 0.9) +
  scale_x_log10() +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ mat_dim, nrow = 3) +
  geom_smooth(method = 'lm') +
  xlab('Life expectancy (years)') +
  ylab(expression(paste('Variance(log ', lambda, ')'))) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.x = element_text(margin = margin(.4, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .4, 0, 0, unit = 'cm')))

# save to file
ggsave('img/fig_0.png', p0, height = 8, width = 10, units = 'in', dpi = 300)





###### 3. Model relationship between l0 and var(log lambda), assuming zero
# measurement error

# compile stan model
stan_regress_hier <- stan_model('stan/regress_hier.stan')

# arrange data for stan
dat_stan <- list(N = nrow(dat_full),
                 N_spp = length(unique(dat_full$spp_int)),
                 spp = dat_full$spp_int,
                 x = log10(dat_full$l0),
                 y = log10(dat_full$var_log_lambda))

# fit stan model
stan_fit <- sampling(
  stan_regress_hier,
  data = dat_stan,
  warmup = 2000,
  iter = 4000,
  thin = 2,
  chains = 2
)

# model diagnostics
# library(shinystan)
# launch_shinystan(stan_fit)

# extract posterior samples for intercept and slope
mu_alpha <- rstan::extract(stan_fit, 'mu_alpha')$mu_alpha
mu_beta <- rstan::extract(stan_fit, 'mu_beta')$mu_beta

# extract posterior samples for best fit line and 95% credible interval
pred_x <- seq(min(dat_stan$x), max(dat_stan$x), length.out = 50)

df_pred <- tibble(mu_alpha, mu_beta, pred_x = list(pred_x)) %>% 
  mutate(pred = pmap(list(mu_alpha, mu_beta, pred_x), ~ ..1 + ..2 * ..3)) %>% 
  dplyr::select(pred_x, pred) %>% 
  unnest() %>% 
  mutate(pred_x = 10^pred_x, pred = 10^pred) %>% 
  group_by(pred_x) %>% 
  summarize(pred_med = quantile(pred, 0.500),
            pred_low = quantile(pred, 0.025),
            pred_upp = quantile(pred, 0.975))

# plot
tt1 <- theme(panel.grid = element_blank(),
            text = element_text(size = 18),
            axis.title = element_text(size = 20),
            panel.background = element_rect(fill = 'grey92'),
            axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
            axis.title.y = element_text(margin = margin(0, .2, 0, 0, unit = 'cm')))

p1 <- ggplot(dat_full, aes(x = pred_x)) +
  geom_line(data = df_pred, aes(y = pred_med), lwd = 1.6, col = 'darkblue') +
  geom_ribbon(data = df_pred, aes(ymin = pred_low, ymax = pred_upp), alpha = 0.2) +
  geom_point(aes(l0, var_log_lambda), shape = 1, size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10(breaks = 10^(-5:2), labels = LabelFn) +
  coord_cartesian(xlim = c(1.05, 450), ylim = c(0.000008, 3.5)) +
  xlab('Life expectancy (years)') +
  ylab(expression(paste('Variance(log ', lambda, ')'))) +
  tt1

# save to file
ggsave('img/fig_1.png', p1, height = 8, width = 10, units = 'in', dpi = 300)





##### 4. Model sampling uncertainty in transition rates and derived parameters
# for population 'B' from Kiviniemi (Plant Ecology, 2002)

# read stage sample size data
kiviniemi_sample_sizes <- read_csv('dat/kiviniemi_sample_sizes.csv') %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup()

# focal species and population
focal_spp <- 'Agrimonia_eupatoria'
focal_pop <- 'B'

# set group factor levels and labels
stage_names <- c('Seedling', 'Juvenile', 'Vegetative', 'Reproductive')
gr_mean <- "1993-98 (Mean)"
gr_pooled <- "1993-98 (Pooled)"
gr_levels <- c("1993", "1994", "1995", "1996", "1997", gr_mean, gr_pooled)
gr_labels <- c("1993-94", "1994-95", "1995-96", "1996-97", "1997-98", gr_mean, gr_pooled)

# subset compadre to focal species and population
kiviniemi_full <- compadre_tidy %>% 
  filter(SpeciesAuthor == focal_spp) %>% 
  filter(MatrixPopulation == focal_pop) %>% 
  mutate(Group = ifelse(MatrixComposite == 'Mean', gr_mean, MatrixStartYear))

# subset to Individual mats and join sample sizes
kiviniemi_individ <- kiviniemi_full %>% 
  filter(MatrixComposite == 'Individual') %>% 
  mutate(ergodic = map_lgl(matA, popdemo::is.matrix_ergodic)) %>% 
  left_join(kiviniemi_sample_sizes, by = c('MatrixPopulation', 'MatrixStartYear'))

# which Individual mats ergodic?
gr_ind_ergodic <- kiviniemi_individ$Group[kiviniemi_individ$ergodic == TRUE]

## matrix properties
# possible transitions in given life cycle
posU <- mean(kiviniemi_individ$matU) > 0
posF <- mean(kiviniemi_individ$matF) > 0
posA <- mean(kiviniemi_individ$matA) > 0

posU_ind <- which(posU, arr.ind = T) %>% as_tibble() %>% mutate(posU = TRUE)
posF_ind <- which(posF, arr.ind = T) %>% as_tibble() %>% mutate(posF = TRUE)
posA_ind <- which(posA, arr.ind = T) %>% as_tibble() %>% mutate(posA = TRUE)

# matrix dimensions, and row/col indices
mat_dim <- nrow(posU)
mat_ind <- which(!is.na(posU), arr.ind = T) %>% as_tibble()

# generate pooled matrices
pooled <- kiviniemi_individ %>% 
  mutate(countU = map2(N, matU, Tr2Count)) %>% 
  mutate(countF = map2(N, matF, Tr2Count)) %>% 
  group_by(MatrixPopulation) %>% 
  summarize(countU = list(Reduce('+', countU)),
            countF = list(Reduce('+', countF)),
            N = list(Reduce('+', N))) %>% 
  mutate(matU = map2(N, countU, Count2Tr)) %>% 
  mutate(matF = map2(N, countF, Count2Tr)) %>% 
  mutate(matA = map2(matU, matF, ~ .x + .y)) %>% 
  mutate(Group = gr_pooled)

# sampling distribution for individual matrices
individ_sim <- kiviniemi_individ %>% 
  mutate(simU = map2(matU, N, ~ SimMatUWrapper(.x, posU, .y, 1000))) %>% 
  mutate(simF = map2(matF, N, ~ SimMatFWrapper(.x, posF, .y, 1000))) %>% 
  mutate(simA = map2(simF, simU, ~ map2(.x, .y, `+`))) %>% 
  dplyr::select(Group, simU, simF, simA) %>% 
  unnest() %>% 
  group_by(Group) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  dplyr::select(Group, rep, simU, simF, simA)

individ_tr <- individ_sim %>% 
  group_by(Group, rep) %>% 
  do(Mat2Tr(mat_ind, .$simA)) %>% 
  ungroup() %>% 
  left_join(posA_ind, by = c('row', 'col')) %>% 
  mutate(val = ifelse(is.na(posA), NA, val))

# sampling distribution for mean matrices
mean_tr <- individ_tr %>% 
  group_by(rep, row, col) %>% 
  summarize(val = mean(val), posA = unique(posA)) %>% 
  ungroup() %>% 
  mutate(Group = gr_mean) %>% 
  dplyr::select(Group, rep, row, col, val, posA)

mean_sim <- mean_tr %>% 
  left_join(posU_ind, by = c('col', 'row')) %>% 
  left_join(posF_ind, by = c('col', 'row')) %>% 
  mutate(tr_U = ifelse(is.na(posU), 0, val)) %>% 
  mutate(tr_F = ifelse(is.na(posF), 0, val)) %>% 
  arrange(Group, rep, col, row) %>% 
  group_by(Group, rep) %>% 
  do(Tr2Mat(.$tr_U, .$tr_F, mat_dim)) %>% 
  ungroup()

# sampling distribution for pooled matrices
pooled_sim <- pooled %>% 
  mutate(simU = map2(matU, N, ~ SimMatUWrapper(.x, posU, .y, 1000))) %>% 
  mutate(simF = map2(matF, N, ~ SimMatFWrapper(.x, posF, .y, 1000))) %>% 
  mutate(simA = map2(simF, simU, ~ map2(.x, .y, `+`))) %>% 
  dplyr::select(Group, simU, simF, simA) %>% 
  unnest() %>% 
  group_by(Group) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  dplyr::select(Group, rep, simU, simF, simA)

pooled_tr <- pooled_sim %>% 
  group_by(Group, rep) %>% 
  do(Mat2Tr(mat_ind, .$simA)) %>% 
  ungroup() %>% 
  left_join(posA_ind, by = c('row', 'col')) %>% 
  mutate(val = ifelse(is.na(posA), NA, val))

# point estimates for all mats
point_tr <- kiviniemi_full %>% 
  dplyr::select(Group, matA) %>% 
  rbind(tibble(Group = gr_pooled, matA = pooled$matA)) %>% 
  group_by(Group) %>% 
  do(Mat2Tr(mat_ind, .$matA)) %>% 
  ungroup() %>% 
  left_join(posA_ind, by = c('row', 'col')) %>% 
  mutate(val = ifelse(is.na(posA), NA, val)) %>% 
  mutate(Group = factor(Group, levels = gr_levels, labels = gr_labels)) %>%
  mutate(group_int = as.numeric(Group)) %>% 
  mutate(col = factor(col, labels = stage_names)) %>% 
  mutate(row = factor(row, labels = stage_names))

point_derived <- kiviniemi_full %>% 
  dplyr::select(Group, matU, matA) %>% 
  rbind(tibble(Group = gr_pooled, matU = pooled$matU, matA = pooled$matA)) %>%
  mutate(lambda = map_dbl(matA, GetLambda)) %>% 
  mutate(l0 = map_dbl(matU, GetLifeExpect)) %>% 
  mutate(Group = factor(Group, levels = gr_levels, labels = gr_labels))

point_var_lambda <- kiviniemi_individ %>%
  filter(ergodic == TRUE) %>% 
  mutate(lambda = map_dbl(matA, lambda)) %>% 
  summarize(var_lambda = var(log(lambda)), N = n())

# join individual, mean, and pooled data
tr_full <- rbind.data.frame(individ_tr, mean_tr, pooled_tr)

mat_derived_full <- rbind.data.frame(individ_sim, mean_sim, pooled_sim) %>% 
  mutate(lambda = map_dbl(simA, popbio::lambda)) %>% 
  mutate(l0 = map_dbl(simU, GetLifeExpect))

## calculate confidence intervals and organize data for plotting
tr_plot <- tr_full %>% 
  group_by(Group, row, col) %>% 
  summarize(med = quantile(val, 0.500, na.rm = T),
            low90 = quantile(val, 0.050, na.rm = T),
            upp90 = quantile(val, 0.950, na.rm = T),
            low99 = quantile(val, 0.005, na.rm = T),
            upp99 = quantile(val, 0.995, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Group = factor(Group, levels = gr_levels, labels = gr_labels)) %>%
  mutate(group_int = as.numeric(Group)) %>% 
  mutate(col = factor(col, labels = stage_names)) %>% 
  mutate(row = factor(row, labels = stage_names))

l0_plot <- mat_derived_full %>% 
  group_by(Group) %>% 
  summarize(l0_med = quantile(l0, 0.500),
            l0_low90 = quantile(l0, 0.050),
            l0_upp90 = quantile(l0, 0.950),
            l0_low99 = quantile(l0, 0.005),
            l0_upp99 = quantile(l0, 0.995)) %>% 
  mutate(Group = factor(Group, levels = gr_levels, labels = gr_labels))

lambda_plot <- mat_derived_full %>% 
  group_by(Group) %>% 
  summarize(med = quantile(log(lambda), 0.500),
            low90 = quantile(log(lambda), 0.050),
            upp90 = quantile(log(lambda), 0.950),
            low99 = quantile(log(lambda), 0.005),
            upp99 = quantile(log(lambda), 0.995)) %>% 
  mutate(Group = factor(Group, levels = gr_levels, labels = gr_labels))

var_lambda <- mat_derived_full %>% 
  filter(Group %in% gr_ind_ergodic) %>% 
  group_by(Group) %>% 
  mutate(rep = sample(rep)) %>% 
  ungroup() %>% 
  group_by(rep) %>% 
  summarize(var_lambda = var(log(lambda)), N = n()) %>% 
  mutate(var_lambda_chisq = map2_dbl(N, var_lambda, SampleVarSim)) %>% 
  summarize(var_l_med = quantile(var_lambda, 0.500),
            var_l_low90 = quantile(var_lambda, 0.050),
            var_l_upp90 = quantile(var_lambda, 0.950),
            var_l_low99 = quantile(var_lambda, 0.005),
            var_l_upp99 = quantile(var_lambda, 0.995),
            var_l_med_chisq = quantile(var_lambda_chisq, 0.500),
            var_l_low90_chisq = quantile(var_lambda_chisq, 0.050),
            var_l_upp90_chisq = quantile(var_lambda_chisq, 0.950),
            var_l_low99_chisq = quantile(var_lambda_chisq, 0.005),
            var_l_upp99_chisq = quantile(var_lambda_chisq, 0.995)) %>% 
  mutate(Group = '1993-1997') %>% 
  mutate(group_int = 1)

## additional data frames to help with plotting
# year labels
x_labs <- tr_plot %>% 
  dplyr::select(group_int, Group) %>% 
  unique() %>% 
  arrange(group_int)

# df to 'x' out stage-transitions that are not possible
stage_display <- posA_ind %>% 
  mutate(display = '') %>% 
  right_join(expand.grid(col = 1:mat_dim, row = 1:mat_dim), by = c('col', 'row')) %>% 
  mutate(display = ifelse(is.na(display), 'x', display)) %>% 
  mutate(col = factor(col, labels = stage_names)) %>% 
  mutate(row = factor(row, labels = stage_names)) %>% 
  mutate(group_int = mean(x_labs$group_int))

## plot stage-transition matrix
tt2 <- theme(panel.background = element_rect(fill = 'grey93'),
              panel.grid = element_blank(),
              strip.text = element_text(size = 16),
              plot.title = element_text(size = 17, hjust = 0, vjust = 0),
              axis.title = element_text(size = 19),
              axis.text.x = element_text(size = 14.5, angle = 60, hjust = 1),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(margin = margin(0, .4, 0, 0, unit = 'cm')),
              axis.title.y.right = element_text(hjust = 0, vjust = -0.9),
              plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), 'pt'))

p2 <- ggplot(tr_plot, aes(x = group_int)) +
  geom_text(data = stage_display, aes(y = 0.5, label = display),
            size = 4.5, hjust = 0.5, vjust = 0.5) +
  geom_linerange(aes(ymin = low90, ymax = upp90), size = 1.4, col = 'grey10') +
  geom_linerange(aes(ymin = low99, ymax = upp99), size = 0.6, col = 'grey10') +
  geom_point(data = point_tr, aes(y = val), shape = 1, size = 3.5, stroke = 0.7) +
  scale_x_continuous(breaks = x_labs$group_int, labels = x_labs$Group, expand = c(0.08, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.1, 0),
                     sec.axis = dup_axis(breaks = NULL, labels = NULL, name = 'Stage, time t+1')) +
  facet_grid(row ~ col) +
  xlab('Year(s)') + ylab('Transition rate') + ggtitle('Stage, time t') +
  tt2

## plot derived parameters
tt3 <- theme(panel.background = element_rect(fill = 'grey93'),
              panel.grid = element_blank(),
              strip.text = element_text(size = 17, hjust = 0),
              plot.title = element_text(size = 18, vjust = 0),
              axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(margin = margin(0, .5, 0, 0, unit = 'cm')))

p3a <- ggplot(l0_plot, aes(x = Group)) +
  geom_linerange(aes(ymin = l0_low90, ymax = l0_upp90), size = 1.4, col = 'grey10') +
  geom_linerange(aes(ymin = l0_low99, ymax = l0_upp99), size = 0.6, col = 'grey10') +
  geom_point(data = point_derived, aes(y = l0), shape = 1, size = 3.5, stroke = 0.7) +
  scale_y_log10() +
  xlab(' ') + ylab(NULL) + ggtitle('Life expectancy (years)') +
  tt3

p3b <- ggplot(lambda_plot, aes(x = Group)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.3) +
  geom_linerange(aes(ymin = low90, ymax = upp90), size = 1.4, col = 'grey10') +
  geom_linerange(aes(ymin = low99, ymax = upp99), size = 0.6, col = 'grey10') +
  geom_point(data = point_derived, aes(y = log(lambda)), shape = 1, size = 3.5, stroke = 0.7) +
  xlab('Year(s)') + ylab(NULL) + ggtitle(expression(paste('log ', lambda))) +
  tt3

p3c <- ggplot(var_lambda, aes(x = group_int)) +
  geom_linerange(aes(ymin = var_l_low90, ymax = var_l_upp90), size = 1.4) +
  geom_linerange(aes(ymin = var_l_low99, ymax = var_l_upp99), size = 0.6) +
  geom_point(data = point_var_lambda, aes(x = 1, y = var_lambda), shape = 1, size = 3.5, stroke = 0.7) +
  coord_cartesian(ylim = c(c(0.0002, 0.03))) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = 1, labels = var_lambda$Group) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1), labels = LabelFn) +
  xlab(NULL) + ylab(NULL) + ggtitle(expression(paste('Variance(log ', lambda, ')'))) +
  tt3

# arrange all panels
g3a <- ggplotGrob(p3a)
g3b <- ggplotGrob(p3b)
g3c <- ggplotGrob(p3c)

g3 <- arrangeGrob(cbind(g3a, g3b, g3c, size = 'first'))

## save plots to file
ggsave('img/fig_2.png', p2, height = 10, width = 10, units = 'in', dpi = 300)
ggsave('img/fig_3.png', g3, height = 4.8, width = 10, units = 'in', dpi = 300)





##### 5. Add error bars from the Kiviniemi population to initial plot

# arrange data
l0_plot_pooled <- l0_plot %>% 
  filter(Group == '1993-98 (Pooled)')

kiviniemi_plot <- var_lambda %>%
  dplyr::select(Group, group_int, var_l_low90, var_l_upp90, var_l_low99, var_l_upp99) %>% 
  mutate(SpeciesAuthor = focal_spp) %>% 
  mutate(MatrixPopulation = focal_pop) %>% 
  left_join(dat_full, by = c('SpeciesAuthor', 'MatrixPopulation')) %>%
  mutate(l0_low90 = l0_plot_pooled$l0_low90) %>% 
  mutate(l0_upp90 = l0_plot_pooled$l0_upp90) %>% 
  mutate(l0_low99 = l0_plot_pooled$l0_low99) %>% 
  mutate(l0_upp99 = l0_plot_pooled$l0_upp99)

# plot
p4 <- ggplot(kiviniemi_plot, aes(l0, var_log_lambda)) +
  geom_point(data = dat_full, shape = 1, size = 3, stroke = 0.7, alpha = 0.3) +
  geom_point(shape = 1, size = 3) +
  geom_linerange(aes(ymin = var_l_low90, ymax = var_l_upp90), size = 1.3) +
  geom_linerange(aes(ymin = var_l_low99, ymax = var_l_upp99), size = 0.5) +
  geom_errorbarh(aes(xmin = l0_low90, xmax = l0_upp90), size = 1.3, height = 0) +
  geom_errorbarh(aes(xmin = l0_low99, xmax = l0_upp99), size = 0.5, height = 0) +
  scale_x_log10() +
  scale_y_log10(breaks = 10^(-5:2), labels = LabelFn) +
  coord_cartesian(xlim = c(1.05, 450), ylim = c(0.000008, 3.5)) +
  xlab('Life expectancy (years)') +
  ylab(expression(paste('Variance(log ', lambda, ')'))) +
  tt1

# save to file
ggsave('img/fig_4.png', p4, height = 8, width = 10, units = 'in', dpi = 300)





##### 6. Redo initial analysis with simulated measurement error in l0 and
# var(log lambda)

# simulate error in l0 and var(log lambda)
dat_full_sim <- dat_full %>% 
  mutate(spp_int = as.integer(as.factor(SpeciesAuthor))) %>%
  mutate(dim_int = as.integer(as.factor(mat_dim))) %>%
  mutate(log_x = log10(l0),
         log_y = log10(var_log_lambda)) %>%
  mutate(log_x_se = runif(n(), 0.02, 0.2)) %>%
  mutate(log_y_se = runif(n(), 0.08, 1)) %>%
  mutate(log_x_low = log_x - 1.96 * log_x_se, log_x_upp = log_x + 1.96 * log_x_se) %>% 
  mutate(log_y_low = log_y - 1.96 * log_y_se, log_y_upp = log_y + 1.96 * log_y_se) %>% 
  mutate(x_low = 10^log_x_low, x_upp = 10^log_x_upp) %>% 
  mutate(y_low = 10^log_y_low, y_upp = 10^log_y_upp)

# compile stan model
stan_regress_hier_error <- stan_model('stan/regress_hier_error.stan')

# arrange data for stan
dat_stan <- list(N = nrow(dat_full_sim),
                 N_spp = length(unique(dat_full_sim$spp_int)),
                 spp = dat_full_sim$spp_int,
                 x_mean = dat_full_sim$log_x,
                 x_se = dat_full_sim$log_x_se,
                 y_mean = dat_full_sim$log_y,
                 y_se = dat_full_sim$log_y_se)

# fit stan model
stan_fit_error <- sampling(
  stan_regress_hier_error,
  data = dat_stan,
  warmup = 2000,
  iter = 4000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 12)
)

# model diagnostics
# library(shinystan)
# launch_shinystan(stan_fit)

# posterior samples for intercept and slope
mu_alpha_error <- rstan::extract(stan_fit_error, 'mu_alpha')$mu_alpha
mu_beta_error <- rstan::extract(stan_fit_error, 'mu_beta')$mu_beta

# posterior samples for best fit line and 95% credible interval
pred_error <- tibble(mu_alpha_error, mu_beta_error, pred_x = list(pred_x)) %>% 
  mutate(pred = pmap(list(mu_alpha_error, mu_beta_error, pred_x), ~ ..1 + ..2 * ..3)) %>% 
  dplyr::select(pred_x, pred) %>% 
  unnest() %>% 
  mutate(pred_x = 10^pred_x, pred = 10^pred) %>% 
  group_by(pred_x) %>% 
  summarize(pred_med = quantile(pred, 0.500),
            pred_low = quantile(pred, 0.025),
            pred_upp = quantile(pred, 0.975))

# plot
p5 <- ggplot(dat_full_sim, aes(x = l0, y = var_log_lambda)) +
  geom_point(shape = 1, size = 3, stroke = 0.7, alpha = 0.25) +
  geom_linerange(aes(ymin = y_low, ymax = y_upp), alpha = 0.25) +
  geom_errorbarh(aes(xmin = x_low, xmax = x_upp), height = 0, alpha = 0.25) +
  geom_line(aes(x = pred_x, y = pred_med),
            data = pred_error, inherit.aes = F, lwd = 1.6, col = 'darkblue') +
  geom_ribbon(aes(x = pred_x, ymin = pred_low, ymax = pred_upp),
              data = pred_error, inherit.aes = F, alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10(breaks = 10^(-5:2), labels = LabelFn) +
  coord_cartesian(xlim = c(1.05, 450), ylim = c(0.000008, 3.5)) +
  xlab('Life expectancy (years)') +
  ylab(expression(paste('Variance(log ', lambda, ')'))) +
  tt1

# save to file
ggsave('img/fig_5.png', p5, height = 8, width = 10, units = 'in', dpi = 300)


