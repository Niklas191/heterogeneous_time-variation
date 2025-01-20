# ==============================================================================
# Capturing heterogeneous time-variation in covariate effects in non-proportional hazard regression models
# Niklas Hagemann, Thomas Kneib and Kathrin Moellenhoff
#
# Code for the simulation study
# ==============================================================================

# Packages and functions  ======================================================

library(pammtools) # package for PEM data transformation
library(mgcv) # package for fitting GAMs
library(foreach) # package for parallel computation
library(doParallel) # package for parallel computation
library(dplyr) # package for data manipulation
library(pec) # package for IBS
library(survival) # package for handling survival datasets
library(CoxFlexBoost) # package for data generation

insertSource("pec_bug_fix.R", package = "pammtools") # This manually solves a bug which can occur when using pec::pec

hetero <- function(t, g, x2){
  if(g == 1){
    effect <- 0.5
  } else if(g == 2){
    effect <- 1
  } else if(g == 3){
    effect <- 1.5
  } else if(g == 4){
    effect <- 2
  }
  return(effect*x2)
}

timeeff <- function(t, x2){
  COEF <- COEFs[1,]
  effect <- COEF[1] + COEF[2]*t + COEF[3]*(t^2) + COEF[4]*(t^3) + COEF[5]*(t^4)
  return((0.5+effect)*x2)
}

FRE <- function(t, g, x2){
  COEF <- COEFs[g,]
  effect <- COEF[1] + COEF[2]*t + COEF[3]*(t^2) + COEF[4]*(t^3) + COEF[5]*(t^4)
  return((0.5+effect)*x2)
}

lambda_fct1 <- function(t, x){
  t <- 10*t
  lambda_value <- exp(0.3*t - 0.3*x[1] + FRE(t, g = x[3], x2 = x[2]))
  return(lambda_value)
}

lambda_fct2 <- function(t, x){
  t <- 10*t
  lambda_value <- exp(0.3*t - 0.3*x[1] + 0.5 * hetero(t, g = x[3], x2 = x[2]) + 0.5 * timeeff(t, x2 = x[2]))
  return(lambda_value)
}

lambda_fct3 <- function(t, x){
  t <- 10*t
  lambda_value <- exp(0.3*t - 0.3*x[1] + hetero(t, g = x[3], x2 = x[2]))
  return(lambda_value)
}

lambda_fct4 <- function(t, x){
  t <- 10*t
  lambda_value <- exp(0.3*t - 0.3*x[1] + timeeff(t, x2 = x[2]))
  return(lambda_value)
}

cens_fct <- function(time, rate = 0.4){
  censor_time <- rexp(n = length(time), rate = rate)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  return(cbind(t_obs, event))
}


# Setup  =======================================================================
# HERE THE USER SPECIFIES WHICH SIMULATION SCENARIO TO USE AND OTHER DETAILS!

n_seed <- 1000 # Number of simulation repetitions
n_total <- 400 # Total number of observations
lambda_fct <- lambda_fct1 # Which scenario (I)-(IV) should be simulated?
n_cores <- 16 # How many cores to use for the parallel computation
global_seed <- 1234567

# Seeds  =======================================================================

set.seed(global_seed) 
seeds <- sample(1:1000000, n_seed)
print(seeds)

# Parallel computation =========================================================

my.cluster = parallel::makeCluster(n_cores, type="FORK", outfile="") # This is for Linux

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Simulation ===================================================================

restult_list <- foreach(i = 1:n_seed, .packages = c("pammtools", "mgcv", "dplyr")) %dopar% {

  simu_seed <- seeds[i]
  set.seed(simu_seed)
  
  COEFs <- matrix(c(0, 1.8, -0.48, 0.032, 0,
                    2, -3.0667, 1.3733, -0.2133, 0.0107,
                    0, -0.87, 0.5447, -0.0816, 0.0037,
                    1, 0.6267, -0.3873, 0.0581, -0.0027), 
                  nrow = 4, 
                  ncol = 5,
                  byrow = TRUE)

  x1 <- runif(n_total, 0, 1)
  x2 <- runif(n_total, 0, 1)
  g <- rep(1:4, each = (n_total/4))
  
  x_mat <- as.matrix(cbind(x1, x2, g))
  x_df <- as.data.frame(x_mat)
  
  while(!data_done) {
    skip_to_next <- FALSE
    # Sometimes the numerical inversion fails. Just trying it again solves this. 
    tryCatch({
      sim_df <- rSurvTime(lambda = lambda_fct4, x = x_mat, cens_fct = cens_fct, upper = 1)
    }, 
    error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next){
      next
    } else{
      data_done <- TRUE
    }
  }

  sim_df$time <- round(sim_df$time, 6) # Needed for numerical stability.
  sim_df$time[sim_df$time <= 1e-06] <- 1e-06 # Needed for PAMMs to be well defined. Also, it is more realistic that the patient does not die during the examination but at least a few minutes after. 
  
  names(sim_df) <- c("time", "event", "x1", "x2", "g")
  sim_df$g <- as.factor(sim_df$g) 

  prop_cen <- 1 - mean(sim_df$event)
  
  pamm_sim_df <- sim_df %>% as.data.frame() %>% as_ped(Surv(time = time, event = event)~., cut = unique(sim_df$time), id = "id")

  pamm_sim_df$tstart <- round(pamm_sim_df$tstart, 6) # Needed for numerical stability.
  pamm_sim_df$tend <- round(pamm_sim_df$tend, 6) # Needed for numerical stability.

  model_fs <- pamm(formula = ped_status ~ 
                     s(tend, bs = "ps", k = 13, m = c(3, 1)) +  
                     s(tend, g, bs = "fs", xt  = list(bs = "ps"), k = 13, m = c(3, 1), by = x2) +
                     x1, 
                   data = pamm_sim_df,
                   family = poisson(), 
                   offset = offset, 
                   method = "REML")

  model_add <- pamm(formula = ped_status ~ 
                      s(tend, bs = "ps", k = 13, m = c(3, 1)) +  
                      s(tend, bs = "ps", k = 13, m = c(3, 1), by = x2) +
                      s(tend, g, bs = "re", by = x2) +
                      x1, 
                    data = pamm_sim_df,
                    family = poisson(), 
                    offset = offset, 
                    method = "REML")
  
  model_re <- pamm(formula = ped_status ~ 
                     s(tend, bs = "ps", k = 13, m = c(3, 1)) +  
                     s(tend, g, bs = "re", by = x2) +
                     x1, 
                   data = pamm_sim_df,
                   family = poisson(), 
                   offset = offset, 
                   method = "REML")
  
  model_ps <- pamm(formula = ped_status ~ 
                     s(tend, bs = "ps", k = 13, m = c(3, 1)) +  
                     s(tend, bs = "ps", k = 13, m = c(3, 1), by = x2) +
                     x1, 
                   data = pamm_sim_df,
                   family = poisson(), 
                   offset = offset, 
                   method = "REML")
  
  aic <- cbind(c("FS", "ADD", "PS", "RE"), 
               c(AIC(model_fs), AIC(model_add), AIC(model_ps), AIC(model_re)))
  bic <- cbind(c("FS", "ADD", "PS", "RE"), 
               c(BIC(model_fs), BIC(model_add), BIC(model_ps), BIC(model_re)))
  logL <- cbind(c("FS", "ADD", "PS", "RE"), 
                   c(logLik(model_fs), logLik(model_add), logLik(model_ps), logLik(model_re)))

  pec_fs = pec(model_fs, 
               Surv(time, event) ~ 1, 
               data = sim_df, 
               times = seq(0.001, 1, by = 0.001),
               start = 0.001, 
               exact = FALSE,
               verbose = FALSE)
  
  pec_add = pec(model_add, 
                Surv(time, event) ~ 1, 
                data = sim_df, 
                times = seq(0.001, 1, by = 0.001),
                start = 0.001, 
                exact = FALSE,
                verbose = FALSE)
  
  pec_ps = pec(model_ps, 
               Surv(time, event) ~ 1, 
               data = sim_df, 
               times = seq(0.001, 1, by = 0.001),
               start = 0.001, 
               exact = FALSE,
               verbose = FALSE)
  
  pec_re = pec(model_re, 
               Surv(time, event) ~ 1, 
               data = sim_df, 
               times = seq(0.001, 1, by = 0.001),
               start = 0.001, 
               exact = FALSE,
               verbose = FALSE)
  

  CRPS <- cbind(c("FS", "ADD", "PS", "RE"), c(
               pec::crps(pec_fs)[2,1], 
               pec::crps(pec_add)[2,1], 
               pec::crps(pec_ps)[2,1], 
               pec::crps(pec_re)[2,1]))
  
  simulation_date <- Sys.Date()
  
  save(model_fs, 
       model_add,
       model_ps,
       model_re,
       model_linear,
       model_smooth,
       prop_cen, 
       aic, 
       bic,
       logL,
       PEC, 
       CRPS,
       simu_seed,
       simulation_date, file = paste0(save_path, "_results_", i, ".Rdata"))
}

