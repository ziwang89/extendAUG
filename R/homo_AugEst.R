#' @title Run extendAUG simulation (Replication and Extension of the Augmented Estimation)
#'
#' @param Nsim number of iterations: integer
#' @param n size of full data: integer (i.e., 2,000,5,000,10,000)
#' @param x_beta coefficients of the covariates ("1" or "2":  1 -- coefficient of x1; 2 -- coefficient of x2)
#'
#' @return Augmented estimator of three models for four size of validation set (m = 100,200, 400, &800)
#' @import matlib MASS wakefield mice robustHD mice
#' @importFrom mice mice complete
#' @importFrom wakefield r_sample_binary
#' @importFrom robustHD standardize
#' @importFrom stats glm
#' @export
#'
homo_AugEst <- function(n = 2000, Nsim = 10, x_beta = 2, parallel_ifrun = FALSE){


  ### -----------

  if (is.infinite(n) | n%%1 !=0){
    stop("Error: n (Size of Full Data) should be a finite integer")
  }else if (is.infinite(Nsim) | Nsim%%1!=0){
    stop("Error: Nsim (Number of Simulation) should be a finite integer")
  }else if (!(x_beta %in% c(1,2))){
    stop("Error: Not Valid x_beta Input (x_beta is only 1 or 2)")
  }
  else
  {
    if (parallel_ifrun == FALSE){
      run_one <- function(Nsim = 1 , n, x_beta){
        tryCatch(
          {
            # run one_time_simulation
            out = one_time_simulation(Nsim , n, x_beta)

            return(out)

          },error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })
      }

      rep_result = replicate(n = Nsim , run_one(Nsim = 1, n, x_beta))
      sim_m100 = list()
      sim_m200 = list()
      sim_m400 = list()
      sim_m800 = list()

      for (i in 1:Nsim) {
        sim_m100[[i]] = rep_result[,i]$m_100
        sim_m200[[i]] = rep_result[,i]$m_200
        sim_m400[[i]] = rep_result[,i]$m_400
        sim_m800[[i]] = rep_result[,i]$m_800
      }
      homo_result = list(sim_100 =  stack(unlist(sim_m100)),
                         sim_200 =  stack(unlist(sim_m200)),
                         sim_400 =  stack(unlist(sim_m400)),
                         sim_800 =  stack(unlist(sim_m800)))

      return(homo_result)

    }else{

      run_para <- function(Nsim , n, x_beta){ ### parallelization
        tryCatch(
          {
            # run one_time_simulation
            rep_AugEst <- function(n = n, Nsim = Nsim, x_beta = x_beta){

              sim_m100 = list()
              sim_m200 = list()
              sim_m400 = list()
              sim_m800 = list()

              RNGkind("L'Ecuyer-CMRG")
              set.seed (670)


              if (Sys.info()[1] == "Windows"){
                cl = makeCluster(detectCores()-2)
                clusterSetRNGStream(cl, iseed=2021)
                sim_results = parLapply(cl, 1:Nsim, one_time_simulation, n=n , x_beta = x_beta)
                stopCluster(cl)
              } else {
                sim_results = mclapply (1:Nsim,  mc.cores=16, one_time_simulation, n = n, x_beta = x_beta)
              }

              for (i in 1:Nsim) {
                sim_m100[[i]] = sim_results[[i]]$m_100
                sim_m200[[i]] = sim_results[[i]]$m_200
                sim_m400[[i]] = sim_results[[i]]$m_400
                sim_m800[[i]] = sim_results[[i]]$m_800
              }
              retern_results = list(sim_100 =  stack(unlist(sim_m100)),
                                    sim_200 =  stack(unlist(sim_m200)),
                                    sim_400 =  stack(unlist(sim_m400)),
                                    sim_800 =  stack(unlist(sim_m800)))

              return(retern_results)
            }
            out =  rep_AugEst(n = n, Nsim = Nsim, x_beta = x_beta)
            return(out)
          },error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })


      }

      homo_result = run_para(Nsim , n, x_beta)
      return(homo_result)

    }
  }


}

