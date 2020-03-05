# server.R script for Bayesian Forecasting Web-Application
# Reactive objects (i.e., those dependent on widget input) are written here
# ------------------------------------------------------------------------------
## Server: all reactive values, output and ui elements
  server <- function(input, output, session) {
    
# Prepare reactive environment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set intial conditions for data input UI
    
  # Define initial values for reactiveValues object
  r <- reactiveValues(
  # These two values define the initial state of the data input UI
  # n - number of rows
  # df - values within those rows
  # Dose Data Input UI
    dose_n = 2,
    dose_df = data.frame(
      "Date" = c(Sys.Date() - 1, Sys.Date() - 1),
      "Time" = c("09:00", "21:00"),
      "Dose" = c(2000, 1500),
      "Freq" = c(2, 2),
      "Repeat" = c(0, 1)
    ),
  # Plasma Sample Input UI
    conc_n = 1,
    conc_df = data.frame(
      "Date" = c(Sys.Date()),
      "Time" = c("20:00"),
      "Conc" = c(23)
    )
  )  # r.reactiveValues
    
# Define data input UI
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data input UI is "reactive" therefore it is coded in the server
# Function used to define the UI is located in global.R
# R.doseInput and R.concInput call this function, providing the reactive
# components from the reactiveValues object `r`.
  
  # Define dose input
    R.doseInput <- reactive({
      dose.input.ui(r$dose_n, r$dose_df, "dose")
    })  # R.doseInput
    
  # Render the UI for the dose input reactive function
    output$doseUI <- renderUI(R.doseInput())
    
  # Define concentration input
    R.concInput <- reactive({
      conc.input.ui(r$conc_n, r$conc_df)
    })  # R.concInput
    
  # Render the UI for the concentration input reactive function
    output$concUI <- renderUI(R.concInput())
    
# Observe logic for data input UI
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Handles the add, remove and save buttons for the data input UI
# Allows for rows to be added or removed depending on the amount of doses and
# plasma concentration samples the user would like to use.
    
  # Use observeEvent to determine if user wants to add or remove a row of boxes
  # Also use to save current state of input.boxes for use
  # It then saves those values to our reactiveValue object (r)
  # Function dose.input.update and conc.input.update are defined in global.R
  # These functions save users current input, must be adjusted according to
  #   any adjustments made to dose.input.ui functions
    
  # Dose data add, remove and save rows
    observeEvent(input$doseAdd, {
      r$dose_n <- r$dose_n + 1
      r$dose_df <- dose.input.update(input, r$dose_n, "dose")
    # If using dateInput, must provide a value
      if (is.na(r$dose_df$Date[r$dose_n])) {
        r$dose_df <- rbind(r$dose_df, r$dose_df[r$dose_n - 1, ])
      }
    })  # observeEvent.doseAdd

    observeEvent(input$doseRem, {
      if (r$dose_n > 1) {
        r$dose_n <- r$dose_n - 1
      }  # if
      r$dose_df <- dose.input.update(input, r$dose_n, "dose")
    })  # observeEvent.doseRem

    observeEvent(input$doseSave, {
      r$dose_df <- dose.input.update(input, r$dose_n, "dose")
      print(r$dose_df)
    })  # observeEvent.doseSave
    
  # Concentration data add, remove and save rows
    observeEvent(input$concAdd, {
      r$conc_n <- r$conc_n + 1
      r$conc_df <- conc.input.update(input, r$conc_n)
    # If using dateInput, must provide a value
      if (is.na(r$conc_df$Date[r$conc_n])) {
        r$conc_df <- rbind(r$conc_df, r$conc_df[r$conc_n - 1, ])
      }
    })  # observeEvent.concAdd

    observeEvent(input$concRem, {
      if (r$conc_n > 1) {
        r$conc_n <- r$conc_n - 1
      }  # if
      r$conc_df <- conc.input.update(input, r$conc_n)
    })  # observeEvent.concRem

    observeEvent(input$concSave, {
      r$conc_df <- conc.input.update(input, r$conc_n)
      print(r$conc_df)
    })  # observeEvent.concSave
    
# Transform data input into long-format
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Combines and prepares data provided from the user to create a dataset that
# can be interpreted by `mrgsolve`, with the resulting output being used for
# empirical Bayesian estimation.
    
  # Define initial dose time (finds earliest dose time)
    R.input.start <- reactive({
      r$dose_df %>%
        mutate(DateTime = paste0(Date, " ", Time, ":00")) %>%  # dplyr::mutate
        pull(DateTime) %>%  # dplyr::pull
        as.POSIXlt() %>% 
        min()
    })  # R.input.start
    
  # Create long-format dose dataset
    R.input.dose <- reactive({
      r$dose_df %>%
        mutate(DateTime = paste0(Date, " ", Time, ":00")) %>%  # dplyr::mutate
        mutate(time = as.numeric(as.POSIXlt(DateTime) - R.input.start(), "hours")) %>%
        mutate(evid = 1, rate = -2, cmt = 1, ii = 24/Freq, addl = Repeat, DV = 0) %>%
        select(time, amt = Dose, evid, rate, cmt, ii, addl, DV)  # dplyr::select
    })  # R.input.dose
    
  # Create long-format concentration dataset
    R.input.conc <- reactive({
      r$conc_df %>%
        mutate(DateTime = paste0(Date, " ", Time, ":00")) %>%  # dplyr::mutate
        mutate(time = as.numeric(as.POSIXlt(DateTime) - R.input.start(), "hours")) %>%
        mutate(amt = 0, evid = 0, rate = 0, cmt = 1, ii = 0, addl = 0) %>%
        select(time, amt, evid, rate, cmt, ii, addl, DV = Conc)  # dplyr::select
    })  # R.input.conc
    
# Run empirical Bayesian estimation step
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Takes the user input and determines the EBEs for the patient data submitted
    
  # Determine empirical Bayesian estimates
    R.output.ebe <- reactive({
    # Extract PPV from OMEGA block
      PPV_vec <- diag(omat(mod, make = TRUE))
    # Define initial parameters
      init_par <- exp(double(length(PPV_vec)))
    # Bind and order dose and concentration data
      input_data <- bind_rows(R.input.dose(), R.input.conc()) %>%
        arrange(time, desc(amt)) %>%
        mutate(WT = input$covWt, CRCL = input$covCrcl, DIAL = as.numeric(input$covDial))
    # Optimise objective function value of bayes.estimate function in global.R
      bayes_estimate <- try(optim(init_par, bayes.estimate, method = "L-BFGS-B",
        input_df = input_data, PPV = PPV_vec,
        lower = rep(0.001, times = length(PPV_vec)), 
        upper = rep(1000, times = length(PPV_vec)),
        control = list(
          parscale = init_par, 
          fnscale = abs(bayes.estimate(init_par, input_data, PPV_vec)),
          factr = 1e12, pgtol = 1e-8
        )
      ))  # try.optim.bayes.estimate
      if (bayes_estimate$message == optim_success) {
        return(bayes_estimate)
      } else {
        return(list(
          par = double(length(PPV_vec))
        ))
      }
    })  # R.output.ebe
    
# Define future dosing input UI
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Using the same function as was used for dose data input, create single line
# of dosing UI that is populated with the last set of doses input by the user.
# Date and time values should be equal to the next dose that would be given 
# based on the users last set of doses.
    
  # Define future dose input
    R.futInput <- reactive({
      latest_dose <- r$dose_df[r$dose_n, ] %>%
        mutate(DateTime = as.POSIXct(paste0(Date, " ", Time, ":00"))) %>%
  # Determine time until next dose in seconds
        mutate(NextDoseSecs = 60*60*(Repeat + 1)*24/Freq) %>%
        mutate(NewDateTime = DateTime + NextDoseSecs) %>%
        mutate(NewDate = as.Date(NewDateTime, "%Y-%m-%d", tz = Sys.timezone())) %>%
        mutate(NewTime = format(NewDateTime, "%H:%M")) %>%
        select(Date = NewDate, Time = NewTime, Dose, Freq, Repeat)
      dose.input.ui(1, latest_dose, "fut")
    })
    
  # Render the UI for the dose input reactive function
    output$futUI <- renderUI(R.futInput())
    
# Create reactive plot for Bayesian forecasting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Using the EBEs display a plot of the individual predicted concentrations prior 
# to the sample, and the forecasted concentrations dependent on future dosing.
# Future concentrations are simulated to take up one thirds of the plot. This
# can be altered by changing how `sim_prop` and `sim_freq` are determined. Care
# should be taken with `sim_freq` in longer time courses, may result in large
# simulations, and therefore slower run-time of the application.
    
  # Create long-format future dose dataset
    R.input.future <- reactive({
    # Read in values from future dose input boxes
      fut_df <- dose.input.update(input, 1, "fut")
    # If UI has not been rendered, set values to initial values
      if (dim(fut_df)[1] == 0) {
        fut_df <- r$dose_df[r$dose_n, ]
      }
    # Reuse dose.input.update function to extract input data
      fut_df %>%
    # Then use same process as done for historic dosing and concentration data
        mutate(DateTime = paste0(Date, " ", Time, ":00")) %>%  # dplyr::mutate
        mutate(time = as.numeric(as.POSIXlt(DateTime) - R.input.start(), "hours")) %>%
        mutate(evid = 1, rate = -2, cmt = 1, ii = 24/Freq, addl = Repeat, DV = 0) %>%
        select(time, amt = Dose, evid, rate, cmt, ii, addl, DV)  # dplyr::select
    })
    
  # Create long-format simulation dataset
    R.input.sim <- reactive({
    # Extract optimised ETA values
      ETA <- log(R.output.ebe()$par)
      ETA_df <- as.data.frame(matrix(c(2, ETA), ncol = length(ETA) + 1))
      names(ETA_df) <- c("ID", paste0("ETA", 1:length(ETA)))
    # Extract future dose time and frequency
      dose_time <- unique(R.input.future()$time)
      dose_freq <- unique(R.input.future()$ii)
    # Define proportion and frequency of simulated concentrations for forecasting
      sim_prop <- 3  # original values take up 1/sim_prop of the plot
      sim_freq <- dose_freq/6  # sim_freq times between each dose
    # Create simulation times
      conc_df <- data.frame(
        time = seq(0, sim_prop*dose_time, by = sim_freq),
        amt = 0, evid = 0, rate = 0, cmt = 1, ii = 0, addl = 0, DV = 0
      )
    # Bind and order input data
      ipre_df <- bind_rows(R.input.dose(), R.input.conc(), R.input.future(), conc_df) %>%
        arrange(time, desc(amt)) %>%
    # Add covariate data, ID = 2 is for individual predicted data
        mutate(ID = 2, WT = input$covWt, CRCL = input$covCrcl) %>%
        mutate(DIAL = as.numeric(input$covDial)) %>%
    # Define ETAs as estimated
        inner_join(ETA_df, by = "ID") %>%  # merge in ETA data
    # Add flag to differentiate forecasted data from past data
        mutate(FLAG = if_else(time >= dose_time & evid == 0, 1, 0)) %>%
    # Assign observed data to a new column so it can be plotted
        mutate(OBS = if_else(DV != 0, DV, NA_real_))  # numeric NA
    # New ID to identify individual predicted concentrations and ETA values
      pred_df <- ipre_df %>%  
        mutate(ID = 1) %>%
        mutate_at(paste0("ETA", 1:length(ETA)), function(x) 0) %>%
    # Filter out forecasted values
        filter(time <= dose_time)
    # Bind output for simultation using mrgsolve
      bind_rows(pred_df, ipre_df)  
    })
    
  # Simulate concentration based on past and future concentrations
    R.plot.data <- reactive({
    # Simulate both predicted and individual predicted concentrations
      R.input.sim() %>%
        mrgsolve::data_set(x = mod, data = .) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt, addl, ii, FLAG, OBS) %>% 
        mrgsolve::mrgsim() %>%
        as.data.frame() %>%
    # Change class of plotting identifiers to factor
        mutate(ID = factor(ID), FLAG = factor(FLAG))
    })
    
  # Render plot to UI    
  # Lines show predictions
  # - colour denotes prediction type
  # - linetype denotes historic vs. future predictions
  # Dots show plasma concentration samples
  # Dashed green horizontal lines denote target concentrations (SA Health)
    output$plot <- renderPlot({
    # Define plotdata
      plot_df <- R.plot.data()
      levels(plot_df$ID) <- c("Population", "Individual")
    # Define palette
      palette1 <- with(cbPalette, c("grey", "red"))
    # Create plot
      p <- NULL
      p <- ggplot(data = plot_df)
      p <- p + geom_line(aes(x = time, y = IPRE, linetype = FLAG,
        colour = ID), size = 0.8)
      p <- p + geom_point(aes(x = time, y = OBS), shape = 1, size = 2, 
        colour = cbPalette$blue)
      p <- p + geom_hline(yintercept = c(15, 20), linetype = "dashed", 
        colour = cbPalette$green)
      p <- p + scale_x_continuous("Time (hours)", breaks = x.breaks.fn)
      p <- p + scale_y_continuous("Concentration (mg/L)")
      p <- p + scale_colour_manual("Prediction", values = palette1)
      p <- p + scale_linetype_manual(NULL, values = c("solid", "dashed"))
      p <- p + guides(linetype = FALSE)
      p
    })
    
# Calculate AUC for UI output
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Using the simulated data that is plotted, the AUC for certain time-frames of
# the data can be provided to the user. The AUC is calculated for the 24 hours
# prior to the forecasted dose, and for the 24 hours at the end of the 
# forecasted doses.
    
  # Calculate AUC for historic and future doses
    R.output.auc <- reactive({
    # Define the final historic and final future dose
      sim_df <- filter(R.plot.data(), ID == 2 & evid == 0)
      hst_time <- pull(R.input.future(), time)
      fut_time <- with(R.input.future(), time + ii*(1+addl))
    # Calculate AUC
      hst_auc <- filter(sim_df, time %in% c(hst_time - 24, hst_time)) %>%
        pull(AUC) %>%
        diff()
      fut_auc <- filter(sim_df, time %in% c(fut_time - 24, fut_time)) %>%
        pull(AUC) %>%
        diff()
    # Return output
      list(hst = hst_auc, fut = fut_auc)
    })
    
  # Render text for valueBox AUC values in the UI
    output$hauc <- renderText({
      hauc <- signif(R.output.auc()$hst, 3)
      return(paste(hauc, "mg.h/L"))
    })
    
    output$fauc <- renderText({
      fauc <- signif(R.output.auc()$fut, 3)
      return(paste(fauc, "mg.h/L"))
    })
    
# Define session behaviour and error messages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reactive objects handling errors and debug
    
  # Close the R session when browser closes
    session$onSessionEnded(function() {
      stopApp()
    })

  # Open console for R session
    observe(label = "console", {
      if(input$console != 0) {
        options(browserNLdisabled = TRUE)
        isolate(browser())
      }
    })
  
  }  #server