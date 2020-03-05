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
    dose_n = 3,
    dose_df = data.frame(
      "Date" = c(Sys.Date(), Sys.Date(), Sys.Date()),
      "Time" = c("09:00", "09:00", "09:00"),
      "Dose" = c(2000, 2000, 2000),
      "Freq" = c(2, 2, 2),
      "Repeat" = c(0, 0, 0)
    ),
  # Plasma Sample Input UI
    conc_n = 1,
    conc_df = data.frame(
      "Date" = c(Sys.Date()),
      "Time" = c("09:00"),
      "Conc" = c(0)
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
      dose.input.ui(r$dose_n, r$dose_df)
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
      r$dose_df <- dose.input.update(input, r$dose_n)
    # If using dateInput, must provide a value
      if (is.na(r$dose_df$Date[r$dose_n])) {
        r$dose_df <- rbind(r$dose_df, r$dose_df[r$dose_n - 1, ])
      }
    })  # observeEvent.doseAdd

    observeEvent(input$doseRem, {
      if (r$dose_n > 1) {
        r$dose_n <- r$dose_n - 1
      }  # if
      r$dose_df <- dose.input.update(input, r$dose_n)
    })  # observeEvent.doseRem

    observeEvent(input$doseSave, {
      r$dose_df <- dose.input.update(input, r$dose_n)
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