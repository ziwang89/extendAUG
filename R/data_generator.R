#' @title Data generator function for extendAUG
#'
#' @param n size of full data: integer, greater than 1000 (i.e., 2,000, 5,000, 10,000)
#' @param m size of validation set (i.e., 100, 200, 400, & 800)
#'
#' @return validation matrix of x, y, s, and their values (for augmented_est function)
#' @import matlib MASS mice robustHD parallel
#' @importFrom mice complete
#' @importFrom wakefield r_sample_binary
#' @importFrom robustHD standardize
#' @importFrom stats qnorm rbinom
#' @importFrom utils stack
#' @export


data_generator <- function(n = 2000, m = 100){
  result_list <- list()

  se1 = 0.90 # sensitivity
  sp1 = 0.95 # specificity

  b0 = -1.0 ## lower the value is, rare the disease is   exp(-1) = 0.3679 --> 36.79% paitents would have the disease
  b1 = 0.5
  b2 = 1

  beta <- c(b0,b1,b2)
  # beta <- c(b0,b1)
  l = length(beta)

  print(paste("data size is ", n,"; betas are ", b0,", ", b1, ", ", b2,"; m is ", m, sep = ""))

  # to generate exposures
  ## binary exposure
  x1 = r_sample_binary(n, x = 1:2, prob = NULL, name = "Binary")      # --- Family history of breast or ovarian cancer Y/N
  ## conti exposure
  x2 = sample(size = n, 12:80, replace = TRUE)     # ---- age

  # to standardize conti exposure
  x2_temp <- standardize(x2)
  # to combine exposures
  x_temp = cbind(x1,x2_temp)

  z = beta %*% t(cbind(1,x_temp))
  pr = 1/(1+exp(-z))     # pass through an inv-logit function


  # to generate  outcome variable
  y = rbinom(n,1,pr)     # bernoulli response variable

  alpha1 = se1        # sensitivity pr(s=1|y=1) for non-exposure group
  alpha2 = sp1

  pr_s = alpha1*(y==1) + (1-alpha2)*(y==0)
  # to generate  naive outcome variable
  s = rbinom(n,1,pr_s)

  ## step I: missing
  # missing 10% data in age > 50
  ind_x2_g50 <- which(x2 >50)
  lenNA_g50 <- length(ind_x2_g50)*0.1
  ind_foo_g50 <- sample(ind_x2_g50, size = lenNA_g50)
  x2[ind_foo_g50] <- NA

  # missing 5% data age <=  50
  ind_x2_l50 <- which(x2 <= 50)
  lenNA_l50 <- length(ind_x2_l50)*0.05
  ind_foo_l50 <- sample(ind_x2_l50, size = lenNA_l50)
  x2[ind_foo_l50] <- NA
  # mimic done
  # ---------------------------

  ## step II: Imputing Missing Data
  x_temp = cbind(x1,x2)
  #summary(x_temp)
  # Imputing Missing Data
  tempDate <- mice(x_temp, printFlag = FALSE)
  #summary(tempDate)
  completedData <- complete(tempDate,1)
  #summary(completedData)
  x_temp <- completedData

  id <- sample(1:n, m, replace = FALSE)
  x <- cbind(1, x_temp)

  #nr <- dim(x)[2]
  x[,3] <- standardize(x[,3])

  val_x = x[id,]          # validation set of x
  val_y = y[id]           # validation set of y
  val_s = s[id]           # validation set of s

  result_list <- list(x = x,s = s,y = y, val_x = val_x, val_s = val_s, val_y = val_y)
  return(result_list)
}
