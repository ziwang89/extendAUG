#' @title Main augmented estimate function for extendAUG
#'
#' @param m size of validation set (i.e., 100, 200, 400, & 800)
#' @param x n*i matrix (predictors with dim n: size of full data times i: num of predictors)
#' @param s surrogate phenotypes
#' @param y true phenotypes (with NA)
#' @param val_x value of x
#' @param val_s value of s
#' @param val_y value of y
#'
#' @return estimates, variance, and ci of three models (M1 valid, M2 native & M4 augment)
#' @import matlib MASS parallel
#' @importFrom stats glm qnorm rbinom
#' @importFrom utils stack
#' @export


augmented_est <- function(m, x, s, y, val_x, val_s, val_y){

  n = dim(x)[1]
  nr = dim(x)[2]

  # the ratio of the validation set size to the full data set size
  ## the validation ratio
  rho <- m/n
  val_x <- as.matrix(val_x)
  x <- as.matrix(x)


  ## Model(1): validation data only
  # glm for validation set only with true disease
  glm1 <- glm(val_y~ val_x[,-1], family="binomial")
  ## Point estimate for Model(1)
  beta_hat <- glm1$coefficients
  var_beta_hat <- summary(glm1)$coefficients[,2]**2

  ## bridge model
  # glm for validation set only with surrogate
  glm2 <- glm(val_s~val_x[,-1],family="binomial")
  gamma_hat <- glm2$coefficients
  var_gamma_hat <- summary(glm2)$coefficients[,2]**2

  ## Model(2): naive approach
  # glm for all surrogate
  glm3 <- glm(s ~ x[,-1],family="binomial")
  # Point estimate for Model(2)
  gamma_bar <- glm3$coefficients
  var_gamma_bar <- summary(glm3)$coefficients[,2]**2

  id <- sample(1:n, m, replace = FALSE)

  ## Preparation for Model(4) the augmented estimator
  ####################################################
  # to get the augmented estimator
  Z1 <- beta_hat %*% t(x)
  Z2 <- gamma_hat %*% t(x)
  Z3 <- gamma_bar %*% t(x)
  val_Z1 <- Z1[id]
  val_Z2 <- Z2[id]
  val_Z3 <- Z3[id]

  info_matrix1 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix1[i,j] <- sum( val_x[,i]^2 * exp(val_Z1) / (1+exp(val_Z1))^2 ) /m
      }else
        info_matrix1[i,j] <- sum( val_x[,i] * val_x[,j] * exp(val_Z1) / (1+exp(val_Z1))^2) /m
      info_matrix1[j,i] <- info_matrix1[i,j]
    }
  }

  info_matrix2 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix2[i,j] <- sum(x[,i]^2 * exp(Z3) / (1+exp(Z3))^2 ) /n
      }else
        info_matrix2[i,j] <- sum(x[,i] * x[,j] * exp(Z3) / (1+exp(Z3))^2) /n
      info_matrix2[j,i] <- info_matrix2[i,j]
    }
  }

  p1 <- matrix(0, nrow = nr, ncol =nr)
  p2 <- matrix(0, nrow = nr, ncol =nr)
  p3 <- matrix(0, nrow = nr, ncol =nr)

  for (i in 1:m){

    phi1 <- (val_y[i] - (exp(val_Z1[i])/(1+exp(val_Z1[i])))) %*% as.numeric(val_x[i,])
    phi1 <- inv(info_matrix1) %*% as.vector(phi1)
    p1 <- p1 + phi1 %*% t(phi1)

    phi2 <- (val_s[i] - (exp(val_Z3[i])/(1+exp(val_Z3[i])))) %*% as.numeric(val_x[i,])
    phi2 <- inv(info_matrix2) %*% as.vector(phi2)
    p2 <- p2 + phi2 %*% t(phi2)

    p3 <- p3 + phi1 %*% t(phi2)

  }

  sigma <- p1 / m
  sigma_star <- (1-rho) * p2 / m
  omega <- (1-rho) * p3 / m

  ## Point estimate for Model(4)
  beta_aug = t(beta_hat - omega %*% inv(sigma_star) %*% (gamma_hat - gamma_bar))
  var_beta_aug = (sigma - omega %*% inv(sigma_star) %*% t(omega)) /m
  var_beta_aug = diag(var_beta_aug)
  # #########################

  #estimators
  est <- rbind(beta_hat, gamma_bar, beta_aug)  # combine estimators: M1: GLM1 --- M2: GLM3 ---- model4

  #confidence interval
  # betagold.ci <- c(t(beta_gold + qnorm(0.975)*sqrt(var_beta_gold)%*%t(c(-1,1))))
  betahat.ci <- c(t(beta_hat + qnorm(0.975)*sqrt(var_beta_hat)%*%t(c(-1,1))))
  gammabar.ci <- c(t(gamma_bar + qnorm(0.975)*sqrt(var_gamma_bar)%*%t(c(-1,1))))
  betaaug.ci <- c(t(c(beta_aug) + qnorm(0.975)*sqrt(var_beta_aug)%*%t(c(-1,1))))

  ci <- rbind(betahat.ci, gammabar.ci, betaaug.ci) # combine CIs: M1: GLM1 --- M2: GLM3 ---- model4

  # variance
  var <- rbind(var_beta_hat, var_gamma_bar, var_beta_aug) # combine Vars: M1: GLM1 --- M2: GLM3 ---- model4

  # diff <- var_beta_hat - var_beta_aug
  return(list(beta_aug = est, ci_aug = ci, var_aug =var))
}
