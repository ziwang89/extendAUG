#' @title Run one-time simulation for extendAUG
#'
#' @param Nsim number of iterations: integer,1
#' @param n size of full data: integer (i.e., 2,000,5,000,10,000)
#' @param x_beta coefficients of the covariates ("1" or "2":  1 -- coefficient of x1; 2 -- coefficient of x2)
#'
#' @return Augmented estimator of three models for four size of validation set (m = 100,200, 400, &800)
#' @import matlib MASS wakefield mice robustHD mice
#' @importFrom mice complete
#' @importFrom wakefield r_sample_binary
#' @importFrom robustHD standardize
#' @export

one_time_simulation <- function(Nsim = 1, n = 2000 , x_beta = 2){

  niters = Nsim

  se1 = 0.90 # sensitivity
  sp1 = 0.95 # specificity

  bias = array(dim = c(24, 4))
  se = array(dim = c(24, 4))
  index = 0

  ind_m = 1
  est_list = list()


  for (ind_m in c(1,2,3,4)){ # for loop for different m value: size of validation set
    m_list = c(100,200,400,800)
    m = m_list[ind_m]

    b0 = -1.0 ## lower the value is, rare the disease is   exp(-1) = 0.3679 --> 36.79% paitents would have the disease
    b1 = 0.5
    b2 = 1

    beta <- c(b0,b1,b2)
    l = length(beta)

    index = index + 1

    print(paste("data size is ", n,"; betas are ", b0,", ", b1, ", ", b2,"; m is ", m, sep = ""))

    final_result_EST = array(dim=c(niters*3,l))
    final_result_VAR = array(dim=c(niters*3,l))
    final_result_CI = array(dim=c(niters*3,l*2))

    i = niters

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
    tempDate <- mice(x_temp, print = FALSE)
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

    # step II: apply Augmentation method

    # assign the results to large matrix
    augmented_result <- augmented_est(m, x, s, y, val_x, val_s, val_y)


    # final result: EST
    final_result_EST[(i*3 -2):(i*3-1),] = augmented_result$beta_aug[c(1,2),]
    final_result_EST[(i*3 -0),] = augmented_result$beta_aug[3,]

    # plot the results
    #j = 2 # coeeficent of x1
    j = x_beta + 1 # coeeficent of x2

    model1 <- NA;
    model2 <- NA;
    model4 <- NA;

    model1[i] <- final_result_EST[i*3-2,j]
    model2[i] <- final_result_EST[i*3-1,j]
    model4[i] <- final_result_EST[i*3-0,j]


    est <- data.frame(model1, model2, model4)
    names(est) = c("Model (1)", "Model (2)", "Model (4)")

    est_list[ind_m] = list(est)
  }

  names(est_list) <- c("m_100", "m_200", "m_400", "m_800")
  return(est_list)
}


