# global.R script for Bayesian Forecasting Web-Application
# Objects that are not reactive are written here
# -----------------------------------------------------------------------------
# Application built using:
# R version 3.4.4 (2018-03-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# Attached packages:
# - shiny 1.3.2
# - shinydashboard 0.7.1
# - plyr 1.8.4
# - dplyr 0.8.1
# - mrgsolve 0.9.1
# - ggplot2 3.1.1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load package libraries
  library(shinydashboard)  # shiny UI package
  library(plyr)  # iterative functions, ldply, llply, ddply
  library(dplyr)  # data manipulation, load after plyr
  library(mrgsolve)  # popPK model simulation, requires dplyr
  library(ggplot2)  # graphical visualisation of data

# Define options
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Set colourblind palette
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
  
# Define successful convergence of optimisation message
  optim_success <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

# Define model objects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Source models
  source("model.R")
  
# Define functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Data Input UI
# This is the function used to create the rendered data input UI
# This could be directly contained with a renderUI function, but is separated
# for clarity.
# The function takes the values for n and df and creates UI accordingly.
# llply is used to create a list() of UI elements so that you add additional
# elements under each previous set of UI elements
# inputId = paste("text", i) to ensure each input is unique
# label = ifelse(...) is used to give a title to the first set of input.boxes
# value = df[i, #] is used to index the saved value in values$df
# Dose data input UI
  dose.input.ui <- function(n, df, name) {
    if(n > 0) {
      llply(seq_len(n), function(i) {
        fluidRow(
          div(class = "MyClass",  #1
            dateInput(paste0(name, "Date", i),
              ifelse(i == 1, "Date", NA),
              df[i, 1]
            )  # input$doseDate-i
          ),  # div
          div(class = "MyClass",  #2
            textInput(paste0(name, "Time", i),
              ifelse(i == 1, "Time", NA),
              df[i, 2]
            )  # input$doseTime-i
          ),  # div
          div(class = "MyClass",  #3
            numericInput(paste0(name, "Dose", i),
              ifelse(i == 1, "Dose", NA),
              df[i, 3]
            )  # input$doseDose-i
          ),  # div
          div(class = "MyClass",  #4
            numericInput(paste0(name, "Freq", i),
              ifelse(i == 1, "Frequency", NA),
              df[i, 4]
            )  # input$doseFreq-i
          ),  # div
          div(class = "MyClass",  #5
            numericInput(paste0(name, "Rep", i),
              ifelse(i == 1, "Repeated Dose", NA),
              df[i, 5]
            )  # input$doseRep-i
          ),  # div
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"  # make inline
          )),  # tags$head
          tags$head(tags$style(type = "text/css", 
            ".MyClass {max-width: 110px}"  # set width
          ))  # tags$head
        )  # fluidRow
      })  # llply
    }  # if
  }  # dose.input.ui.function

# Concentration data input function
  conc.input.ui <- function(n, df) {
    if(n > 0) {
      llply(seq_len(n), function(i) {
        fluidRow(
          div(class = "MyClass",  #1
            dateInput(paste0("concDate", i),
              ifelse(i == 1, "Date", NA),
              df[i, 1]
            )  # input$concDate-i
          ),  # div
          div(class = "MyClass",  #2
            textInput(paste0("concTime", i),
              ifelse(i == 1, "Time", NA),
              df[i, 2]
            )  # input$concTime-i
          ),  # div
          div(class = "MyClass",  #3
            numericInput(paste0("concConc", i),
              ifelse(i == 1, "Concentration", NA),
              df[i, 3]
            )  # input$concConc-i
          ),  # div
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"  # make inline
          )),  # tags$head
          tags$head(tags$style(type = "text/css", 
            ".MyClass {max-width: 110px}"  # set width
          ))  # tags$head
        )  # fluidRow
      })  # llply
    }  # if
  }  # conc.input.ui.function
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Data Input Update
# Used with  observeEvent to determine if user wants to add or remove a row of 
# boxes. Also used to save current state of input.boxes for use.
# Each set either increases or decreases the number of boxes to render (n)
# ldply to create a data.frame with a number of rows (designated by n) with
# the values of that data.frame being the input of the user
# It then saves those values to our reactiveValue object (r)
# To determine what the input of the user is you need to refer to the input
# object like you normally would for accessing values from the UI
# However our input names are input$#i (e.g. input$doseDate1, input$doseDate2 etc.)
# To access these variable input names we use get() which allows us to use
# i (a sequential number determined by ldply) to access each input
# Dose data input update
  dose.input.update <- function(input, n, name) {
    ldply(seq_len(n), function(i) {
      data.frame(
        "Date" = get("input")[[paste0(name, "Date", i)]],
        "Time" = get("input")[[paste0(name, "Time", i)]],
        "Dose" = get("input")[[paste0(name, "Dose", i)]],
        "Freq" = get("input")[[paste0(name, "Freq", i)]],
        "Repeat" = get("input")[[paste0(name, "Rep", i)]])
    })  # ldply
  }  # dose.input.update.function
  
# Concentration data input update
  conc.input.update <- function(input, n) {
    ldply(seq_len(n), function(i) {
      data.frame(
        "Date" = get("input")[[paste0("concDate", i)]],
        "Time" = get("input")[[paste0("concTime", i)]],
        "Conc" = get("input")[[paste0("concConc", i)]])
    })  # ldply
  }  # conc.input.update.function
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Empirical Bayesian Estimation
# Designed to take estimates for ETAs and return objective function value.
# Also takes an input `mrgsolve` data.frame and PPV values from the OMEGA block
  bayes.estimate <- function(par, input_df, PPV) {
  # Describe parameters to be optimised within mrgsolve data set
    ETA <- log(par)
    ETA_df <- as.data.frame(matrix(c(1, ETA), ncol = length(ETA) + 1))
    names(ETA_df) <- c("ID", paste0("ETA", 1:length(ETA)))
  # Create ETA columns for each ETA
    output_df <- input_df %>%
      mutate(ID = 1) %>%
      inner_join(ETA_df, by = "ID") %>%  # merge in ETA data
  # Run data through mrgsolve
      mrgsolve::data_set(x = mod, data = .) %>%
      mrgsolve::carry_out(amt, evid, rate, cmt, addl, ii) %>% 
      mrgsolve::mrgsim() %>%
      as.data.frame() %>%
  # Ensure IPRED has finite values
      mutate(IPRE = dplyr::if_else(
        !is.finite(IPRE) | IPRE < .Machine$double.eps,
        .Machine$double.eps,
        IPRE
      ))
  # Define DV and DVhat
    DV <- filter(input_df, evid == 0) %>% pull(DV)
    DVhat <- filter(output_df, evid == 0) %>% pull(IPRE)
  # Posterior log-likelihood
  # Error model: IPRE*(1 + ERR_PRO) + ERR_ADD
  # Can be simplified to: IPRE + W*ERR
  # Where W = sqrt(pow(IPRE*ERR_PRO, 2) + pow(ERR_ADD, 2))
    loglikpost_sd <- sqrt((DVhat*mod$ERR_PRO)^2 + mod$ERR_ADD^2)
    loglikpost <- dnorm(DV, mean = DVhat, sd = loglikpost_sd, log = T)
  # Prior log-likelihood
    loglikprior <- dnorm(ETA, mean = 0, sd = sqrt(PPV), log = T)
  # Return objective function value to be minimised
    return(-1*sum(loglikpost, loglikprior))
  }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Plot Auto-Axis Breaks Functions
# Function to automatically set the breaks for the dose forecast plot, so that
# they are always divisible by 24.
  x.breaks.fn <- function(x) {
  # Extract upper and lower limit
    upper <- x[[2]]
    lower <- x[[1]]
  # Define how breaks are determined
  # Breaks from zero to the upper value divided by 24 rounded up
  # Then multiple these numbers by 24, to get 0, 24, 48, ..., ceiling(upper/24)
    24*(0:(ceiling(upper/24)))
  }
  