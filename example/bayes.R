# Vancomycin Dose-Adjustment Regimen - Bayesian Optimised Dosing Regimen
# -----------------------------------------------------------------------------
# Takes induction dataset and changes the dosage after the time when the first
#   blood sample would be taken. Uses empirical bayesian estimation to 
#   determine patient pharmacokinetic characteristics. This is then used to
#   determine the best dose for the patient.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
  setwd(paste0(getwd(), "/example"))

# Load package libraries
  library(dplyr)  # Split and rearrange data - required for mrgsolve
  library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # Graphical package

# Source regimen
  source("regimen.R")
  
# Read .rds files for population, input and induction datasets
  # pop_tb <- readr::read_rds("output/vanc_population.rds")
  # input_tb <- readr::read_rds("output/vanc_induction.rds")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# As there are variable sample times, the latest blood sample is marked using
#   the dose_num object, which specifies which dose is currently the focus.
# This is no more reliable than using sample_time
  bayes_fn <- function(induction_tb) {
  # Change ID2 column name to ID and determine final time
    induction_tb <- dplyr::rename(induction_tb, ID = ID2)
    final_time <- tail(induction_tb, 1) %>% dplyr::pull(time)
  # Determine when to take sample according to exerpt from guidelines:
  # "For patients with normal renal function, the first TDM (trough) should 
  #  occur just prior (within one hour) to the fourth dose or on day 3, 
  #  whichever occurs earlier. In patients with GFR between 20-39 mL/min check 
  #  trough level before (within one hour) of the third dose."
    subj_crcl <- unique(induction_tb$CRCL)
    if (subj_crcl <= 39) {
      if (subj_crcl < 20) {
        stop("Function not designed for patients with CrCL < 20 mL/min.")
      }
      dose_num <- 3
    } else {
      dose_num <- 4
    }
    bayesin_tb <- induction_tb
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Begin loop for successive dose adjustment
    sample_all <- c(NULL)
  # Repeat the loop until dose has been optimised for the entirety of patients
  #   treatment over the course of 168 hours
    repeat {
    # Determine sample time and sample concentration based on the designated
    #   sample number assigned initially, or after dose adjustment.
      sample_time <- dplyr::filter(bayesin_tb, evid != 0) %>%
        dplyr::slice(dose_num) %>%
        dplyr::pull(time)
    # Collect all sample times in one variable
      sample_all <- c(sample_all, sample_time)
    # Determine previous dose and frequency
      last_amt <- dplyr::filter(bayesin_tb, time == sample_time) %>%
        dplyr::pull(amt)
      last_frq <- dplyr::filter(bayesin_tb, evid != 0 & time >= sample_time) %>%
        dplyr::pull(time) %>% diff() %>% unique()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Empirical Bayesian Estimation
    # Loop until successful minimisation
    # Define run number (no runs completed)
      run_num <- 0
    # Create sample dataset, which provides no info on DV beyond the current
    #   sample
      sample_tb <- bayesin_tb %>%
        dplyr::filter(evid == 1) %>%
        dplyr::mutate(DV = dplyr::if_else(time %in% sample_all, DV, NA_real_))
    # Repeat the following until empirical Bayesian estimation is successful
      repeat {
      # Extract omega distribution data from model
        ETABSV <- mrgsolve::omat(mod, make = TRUE) %>% diag()
        n_eta <- length(ETABSV)
        mod_etas <- 1:n_eta
      # Initial estimates for Bayes parameters
        init_par <- exp(double(n_eta))
        if (run_num > 0) {
          init_par <- init_par*exp(runif(n_eta, min = -0.01, max = 0.01))
        }
      # Previous dependent variable values
        prev_DV <- dplyr::filter(sample_tb, time %in% sample_all) %>% 
          dplyr::pull(DV)
      # Define bayesian estimation function
        bayes_estimate_fn <- function(par) {
        # Describe parameters to be optimised within mrgsolve data set
          ETA <- log(par)
          ETA_df <- as.data.frame(matrix(
            c(unique(sample_tb$ID), ETA), 
            nrow = 1
          ))
          names(ETA_df) <- c("ID", paste0("ETA", 1:length(ETA)))
        # Define mrgsolve dataset
          estim_tb <- sample_tb %>%
            dplyr::select(ID, time, evid, amt, cmt, rate, WT, DIAL, CRCL) %>%
        # Define ETAs as estimated
            dplyr::inner_join(ETA_df, by = "ID") %>%  # merge in ETA data
        # Run data through mrgsolve, with idata using initial tumour size
            mrgsolve::data_set(x = mod, data = .) %>%
            mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
            mrgsolve::mrgsim() %>%
            tibble::as_tibble() %>% 
        # Ensure IPRED has finite values
            dplyr::mutate(IPRE = dplyr::if_else(
              !is.finite(IPRE) | IPRE < .Machine$double.eps,
              .Machine$double.eps,
              IPRE
            ))
        # Define yhat
          yhat <- dplyr::filter(estim_tb, time %in% sample_all) %>%
            dplyr::pull(IPRE)
        # Posterior log-likelihood
        # Error model: IPRE*(1 + ERR_PRO) + ERR_ADD
        # Can be simplified to: IPRE + W*ERR
        # Where W = sqrt(pow(IPRE*ERR_PRO, 2) + pow(ERR_ADD, 2))
          loglikpost_sd <- sqrt((yhat*mod$ERR_PRO)^2 + mod$ERR_ADD^2)
          loglikpost <- dnorm(prev_DV, mean = yhat, sd = loglikpost_sd, log = T)
        # Prior log-likelihood
          loglikprior <- dnorm(ETA, mean = 0, sd = sqrt(ETABSV), log = T)
        # Return objective function value to be minimised
          return(-1*sum(loglikpost, loglikprior))
        }  # end bayes_estimate_fn
      # Run bayes_estimate_fn using optim()
        bayes_estimate <- try(optim(init_par, bayes_estimate_fn, method = "L-BFGS-B",
          lower = rep(0.001, times = n_eta), upper = rep(1000, times = n_eta),
          control = list(
            parscale = init_par, fnscale = bayes_estimate_fn(init_par),
            factr = 1e12, pgtol = 1e-8
          )
        ))
      # Stop function if there is an error with Bayesian estimation
        if (class(bayes_estimate) == "try-error") browser()  # error catch
        minimised <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
        if (bayes_estimate$message == minimised) break
        run_num <- run_num + 1
      }  # end loop as successful minimisation has occurred
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate concentrations according to new Bayes estimates
      estim_ETA_df <- as.data.frame(matrix(
        c(unique(sample_tb$ID), bayes_estimate$par), 
        ncol = length(bayes_estimate$par) + 1
      ))
      names(estim_ETA_df) <- c("ID", paste0("ETA", 1:length(bayes_estimate$par)))
    # Create input dataset for simulation
      input_sim_tb <- sample_tb %>%
        dplyr::select(ID, time, evid, amt, cmt, rate, WT, DIAL, CRCL) %>%
    # Add in values for ETA values that correspond with EBEs
        dplyr::inner_join(estim_ETA_df, by = "ID")  # merge in ETA data
    # Simulate new concentrations
      output_sim_tb <- input_sim_tb %>%
        dplyr::filter(time <= sample_time) %>%
        mrgsolve::data_set(x = mod, data = .) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optimise dose for the individual using individual Bayes predicted 
    #   concentrations (and compartment amounts) at time of last sampling.
    #   Only optimise if the AUC0-24 for the last day is outside of the target
    #   range.
    # Determine AUC over the last 24 hours then define when next sample occurs
      last_bayes_auc <- output_sim_tb %>%
        dplyr::filter(time %in% c(sample_time - 24, sample_time) & evid == 1) %>%
        dplyr::pull(AUC) %>% diff()
      sample_next <- dplyr::if_else(last_bayes_auc >= 800, 24, 72)
      next_time <- sample_time + sample_next
      if (next_time > final_time) next_time <- final_time
    # Define trough target and upper bound, begin if statement
      auc_lower <- 400
      auc_target <- 500
      auc_upper <- 550
      if (last_bayes_auc < auc_lower | last_bayes_auc >= auc_upper | last_amt == 0) {
      # Modify model code ready for simulation
        optim_mod <- mrgsolve::init(mod, list(
          CENT = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(CENT),
          PERI = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(PERI),
          AUC = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(AUC)
        ))
      # Set initial dose and error estimates
        init_par <- c(1000, 0.01)
      # Subset input dataset so only future concentrations are predicted
        test_times <- seq(sample_time, final_time, by = last_frq)
        input_optim_tb <- input_sim_tb %>% 
          dplyr::filter(time >= sample_time & time <= next_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% test_times, 1, 0),
            rate = dplyr::if_else(evid != 0, -2, 0)
          )
      # Find the doses that maximise the likelihood of trough concentrations
      #   being the target concentration
        optimise_dose_fn <- function(par) {
        # Add fitted parameters to the input data set, then simulate 
        #   concentration-time profiles with fitted doses
          output_optim_tb <- input_optim_tb %>%
            dplyr::mutate(amt = dplyr::if_else(
              evid == 1, par[1], amt
            )) %>%
            mrgsolve::data_set(x = optim_mod, data = .) %>%
            mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
            mrgsolve::mrgsim(start = sample_time, end = next_time) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(IPRE = dplyr::if_else(
              !is.finite(IPRE) | IPRE < .Machine$double.eps,
              .Machine$double.eps,
              IPRE
            ))
        # Define yhat and the residual
          yhat <- output_optim_tb %>%
            dplyr::filter(time %in% c(next_time - 24, next_time) & evid == 1) %>%
            dplyr::pull(AUC) %>% diff()
          res <- dnorm(auc_target, yhat, yhat*par[2], log = T)
        # Objective function value to be minimised
          return(-1*sum(res))
        }
        optimise_dose <- try(optim(init_par, optimise_dose_fn, method = "L-BFGS-B",
          lower = c(0.0001, 0.0001),
          upper = c(5000, Inf),
          control = list(parscale = init_par, factr = 1e7)
        ))
        if (class(optimise_dose) == "try-error") browser()  # error catch
      # Administer the individual the optimised dose
        exact_amt <- optimise_dose$par[1]
        if (exact_amt < 1) exact_amt <- 0
        dose_amt <- ceiling(exact_amt/250)*250
      } else {
      # Give previous dose
        dose_amt <- last_amt
      }  # end if at target auc
    # Check to see if dose change is required
      if (dose_amt == last_amt) {
        input_final_tb <- bayesin_tb %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      } else {
    # If dose dose need changing, adjust dataset to dosing frequency
    # This is mainly here in case variable dosing frequency is tested, as it
    # stands this works well, but may not work without it.
      # Define new dose times
        dose_times <- seq(sample_time, final_time, by = last_frq)
      # Create dataset that is adjusted for these new dose times
        adjusted_tb <- dplyr::filter(bayesin_tb, time >= sample_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% dose_times, 1, 0),
            amt = dplyr::if_else(evid == 1, dose_amt, 0),
            rate = dplyr::if_else(amt != 0, -2, 0)
          )
      # Create input dataset using these new dose times
        input_final_tb <- dplyr::filter(bayesin_tb, time < sample_time) %>%
          dplyr::bind_rows(adjusted_tb) %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      }
    # Simulate to represent time passing since last sample
      bayesout_tb <- mod %>%
        mrgsolve::data_set(data = input_final_tb) %>%
        mrgsolve::idata_set(data = pop_tb) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>%
        mrgsolve::mrgsim() %>%
        tibble::as_tibble()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Determine when the next blood sample should be taken. If dose was held
    #   next blood sample should be in 24 hours. Otherwise, next sample should
    #   be in 72 hours.
    # This is done above as it impacts dosing optimisation
      # sample_next <- dplyr::if_else(auc24_val >= 800, 24, 72)
      dose_num <- dose_num + sample_next/last_frq
      if (dose_num >= dim(dplyr::filter(bayesout_tb, evid != 0))[1]) break
      bayesin_tb <- bayesout_tb
    }  # end repeat
    dplyr::select(bayesout_tb, -ID)
  }  # brackets closing "bayes_fn"
  
  tictoc::tic()
  final_tb <- output_tb %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, bayes_fn))  # create new list column using bayes_fn
  tictoc::toc() 
    
  readr::write_rds(output_tb, path = "output/vanc_regimen_bayes_fix.rds")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot patient data
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Set palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  
# Plot median patient vancomycin concentrations with 90% confidence intervals
  plot_tb <- dplyr::select(output_tb, ID, bayes) %>%
    tidyr::unnest()
  p <- NULL
  p <- ggplot(data = plot_tb)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = "red", size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = "red", size = 1, alpha = 0.25)
  p <- p + labs(x = "Time (hours)", y = "Vancomycin Concentration (mg/L)")
  p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p <- p + facet_wrap(~CRCL)
  p
  
# Plot % patients with AUC > 400 mg.h/L
  auc_tb <- dplyr::filter(plot_tb, time %in% (0:7*24)) %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      dplyr::mutate(data, dAUC = c(0, diff(AUC)))
    })) %>% tidyr::unnest() %>%
    dplyr::filter(time != 0) %>%
    dplyr::mutate(time = time/24)
  
  auc_target_bycrcl <- auc_tb %>%
    dplyr::group_by(time, CRCL) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      tibble::tibble(
        gt400 = sum(data$dAUC > 400)/dim(data)[1],
        gt700 = sum(data$dAUC > 700)/dim(data)[1],
        gt400lt700 = sum(data$dAUC >= 400 & data$dAUC <= 700)/dim(data)[1]
      )
    })) %>% tidyr::unnest()
  
  p <- NULL
  p <- ggplot(data = auc_target_bycrcl)
  p <- p + ggtitle("Probability of Target Attainment - Model-Based TDM")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p <- p + facet_wrap(~CRCL)
  p
  