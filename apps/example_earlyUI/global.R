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
  dose.input.ui <- function(n, df) {
    if(n > 0) {
      llply(seq_len(n), function(i) {
        fluidRow(
          div(class = "MyClass",  #1
            dateInput(paste0("doseDate", i),
              ifelse(i == 1, "Date", NA),
              df[i, 1]
            )  # input$doseDate-i
          ),  # div
          div(class = "MyClass",  #2
            textInput(paste0("doseTime", i),
              ifelse(i == 1, "Time", NA),
              df[i, 2]
            )  # input$doseTime-i
          ),  # div
          div(class = "MyClass",  #3
            numericInput(paste0("doseDose", i),
              ifelse(i == 1, "Dose", NA),
              df[i, 3]
            )  # input$doseDose-i
          ),  # div
          div(class = "MyClass",  #4
            numericInput(paste0("doseFreq", i),
              ifelse(i == 1, "Frequency", NA),
              df[i, 4]
            )  # input$doseFreq-i
          ),  # div
          div(class = "MyClass",  #5
            numericInput(paste0("doseRep", i),
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
  dose.input.update <- function(input, n) {
    ldply(seq_len(n), function(i) {
      data.frame(
        "Date" = get("input")[[paste0("doseDate", i)]],
        "Time" = get("input")[[paste0("doseTime", i)]],
        "Dose" = get("input")[[paste0("doseDose", i)]],
        "Freq" = get("input")[[paste0("doseFreq", i)]],
        "Repeat" = get("input")[[paste0("doseRep", i)]])
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