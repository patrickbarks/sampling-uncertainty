

# simulate from sampling-distribution of stage-specific transition rates
# using Dirichlet distribution, conjugate to multinomial
SimStageMatU <- function(x, vital_ind, n) {
  colsum <- sum(x)
  if(length(vital_ind) == 0 | is.na(n)) {
    out <- x
  } else {
    if(colsum > 1) x <- x / colsum
    mortality <- 1 - sum(x)
    rates <- c(x[vital_ind], mortality)
    frequencies <- round(rates * n)
    rates_sim <- MCMCpack::rdirichlet(1, frequencies + 1)
    out <- vector(mode = 'numeric', length = length(x))
    out[vital_ind] <- rates_sim[1:length(vital_ind)]
  }
  return(out)
}

# simulate from sampling-distribution of stage-specific fecundity rates
# using Gamma distribution, conjugate to Poisson
SimStageMatF <- function(x, vital_ind, n) {
  if(length(vital_ind) == 0 | is.na(n)) {
    out <- x
  } else {
    rates <- x[vital_ind]
    y <- round(rates * n)
    rates_sim <- sapply(y, function(y) rgamma(1, shape = 1 + y, rate = 0 + n))
    out <- vector(mode = 'numeric', length = length(x))
    out[vital_ind] <- rates_sim
  }
  return(out)
}

# generate simulated matU by applying SimStageMatU to each column of matU
SimMatU <- function(matU, posU, N) {
  if(class(matU) == 'list') matU <- matU[[1]]
  if(class(N) == 'list') N <- unlist(N)
  
  simU <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
  
  for(i in 1:ncol(matU)) {
    vital_ind <- which(posU[,i] > 0)
    simU[,i] <- SimStageMatU(matU[,i], vital_ind, N[i])
  }
  return(simU)
}

# generate simulated matF by applying SimStageMatF to each column of matF
SimMatF <- function(matF, posF, N) {
  if(class(matF) == 'list') matF <- matF[[1]]
  if(class(N) == 'list') N <- unlist(N)
  
  simF <- matrix(0, nrow = nrow(matF), ncol = ncol(matF))
  
  for(i in 1:ncol(matF)) {
    vital_ind <- which(posF[,i] > 0)
    simF[,i] <- SimStageMatF(matF[,i], vital_ind, N[i])
  }
  return(simF)
}

# wrapper function to simulate nsim matU (i.e. apply SimMatU fn nsim times)
SimMatUWrapper <- function(matU, posU, N, nsim) {
  return(replicate(nsim, SimMatU(matU, posU, N), simplify = FALSE))
}

# wrapper function to simulate nsim matF (i.e. apply SimMatF fn nsim times)
SimMatFWrapper <- function(matF, posF, N, nsim) {
  return(replicate(nsim, SimMatF(matF, posF, N), simplify = FALSE))
}

# calculate life expectancy from matU (return NA if error or l0 < 0)
GetLifeExpect <- function(matU, start = 1) {
  dim <- nrow(matU)
  N <- try(solve(diag(dim) - matU), silent = TRUE)
  if(class(N) == 'try-error') {
    l0 <- NA_real_
  } else {
    l0 <- colSums(N)[start]
  }
  l0 <- ifelse(l0 < 0, NA_real_, l0)
  return(l0)
}

# calculate lambda from matA (return NA if error)
GetLambda <- function(matA) {
  lambda <- try(popbio::lambda(matA), silent = TRUE)
  if(class(lambda) == 'try-error') {
    lambda <- NA_real_
  }
  return(lambda)
}

# extract element-wise transition rates from matrix
Mat2Tr <- function(mat_indices, mat) {
  mat_indices$val <- unlist(c(mat))
  return(mat_indices)
}

# construct matU, matF, and matA, given transition rates and matrix dimensions
Tr2Mat <- function(tr_U, tr_F, mat_dim) {
  simU <- matrix(tr_U, nrow = mat_dim, ncol = mat_dim)
  simF <- matrix(tr_F, nrow = mat_dim, ncol = mat_dim)
  simA <- simU + simF
  return(tibble(simU = list(simU), simF = list(simF), simA = list(simA)))
}

# transform transition rates to integer counts, given stage sample-sizes N
Tr2Count <- function(N, matU) {
  return(sapply(seq_along(N), function(i) round(N[i] * matU[,i])))
}

# transform transition counts to transition rates, given stage sample-sizes N
Count2Tr <- function(N, countU) {
  return(sapply(seq_along(N), function(i) countU[,i] / N[i]))
}

# calculate element-wise mean matrix given list of matrices
MeanMat <- function(mat_list) {
  len <- length(mat_list)
  if (len > 1) {
    matMean <- Reduce('+', mat_list) / len
  } else {
    matMean <- mat_list[[1]]
  }
  return(list(matMean))
}

# check whether all matrices in list are equal
CheckMatsEqual <- function(mat_list) {
  len <- length(mat_list)
  if (len > 1) {
    all_equal <- all(sapply(mat_list, function(x) all(x == mat_list[[1]])))
  } else {
    all_equal <- TRUE
  }
  return(all_equal)
}

# check whether life cycle perennial (lx[3] > 0)
CheckPerennial <- function(matU, startLife, N = 3) {
  matUtemp <- matU
  lx <- vector(mode = 'numeric', length = N)
  
  for (i in 1:N) {
    lx[i] <- sum(matUtemp[,startLife])
    matUtemp <- matUtemp %*% matU
  }
  
  return(lx[N] > 0)
}

# simulate from sampling distribution of sample mean
SampleMeanSim <- function(y) {
  if(all(is.na(y))) {
    return(NA_real_)
  } else {
    y_mean <- mean(y)
    y_se <- sd(y) / sqrt(length(y))
    return(rnorm(1, y_mean, y_se))
  }
}

# simulate from sampling distribution of sample variance
SampleVarSim <- function(N, y) {
  return((N - 1) * y / rchisq(1, N - 1))
}

# logit transform
logit <- function(p) {
  log(p / (1 - p))
}

# inverse logit transform
logit_inv <- function(a) {
  exp(a) / (exp(a) + 1)
}

# format plot labels using decimal notation rather than scientific notation
LabelFn <- function(x) {
  return(formatC(x, format = 'fg'))
}








