# ui.R script for Bayesian Forecasting Web-Application
# The user-interface and widget input for the Shiny application is defined here
# Sends user-defined input to server.R, calls created output from server.R
# ------------------------------------------------------------------------------

# Header: More can be done here, but just simple for now
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  header <- dashboardHeader(
  	titleWidth = 200,
    title = "Vancomycin Dosing"
  )  #dashboardHeader
  
# Sidebar: Can contain tabs, can also contain inputs if desired
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determines how many tabs exist in the body of the ui
  sidebar <- dashboardSidebar(
  	width = 200, #Width of sidebar the same as width of header
    
  # Sidebar options
    sidebarMenu(
  		menuItem("Patient Information", tabName = "patient-tab", 
  		  icon = icon("child")
  		),  # menuItem.patient-tab
  		menuItem("Bayesian Forecasting", tabName = "bayes-tab", 
  		  icon = icon("prescription-bottle")
  		),  # menuItem.patient-tab
  		br(),
  		div(style = "padding-left: 30px", actionButton("console", "Debug Console"))
    )  # sidebarMenu
    
  )  # dashboardSidebar
  
# Body: Main content of each tab, as determined by sidebar
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  body <- dashboardBody(
    tabItems(
      
  ### Patient Information Tab
  		tabItem(tabName = "patient-tab",
  		    
      # Dosing Input
		    box(title = "Dose Data Input", width = 12, status = "primary",
		      uiOutput("doseUI"),
          fluidRow(
            actionButton("doseAdd", "Add"),
            actionButton("doseRem", "Remove"),
            actionButton("doseSave", "Save")
          ),  # fluidRow
		      align = "center"
		    ),  # box.dose-input
  		  
  		  fluidRow(
  	  # Covariate Input
    			box(title = "Patient Data Input", width = 6, status = "primary",
    			  numericInput("covWt", "Weight (kg)",
    			    value = 75, min = 7.5, max = 750
    			  ),  # numericInput.patWt
    			  numericInput("covCrcl", "Creatinine Clearance (mL/min)",
    			    value = 100, min = 0, max = 200
    			  ),  # numericInput.patCrcl
    			  radioButtons("covDial", "Haemodialysis Status", inline = T,
    			    choices = list(
    			      "Yes" = 1,
    			      "No" = 0
    			    ),  # choices.list
    			    selected = 0
    			  )  # radioButtons.patDial
    		  ),  # box.patient-input
  		    
  	  # Plasma Sample Input
  		    box(title = "Sample Data Input", width = 6, status = "primary",
  		      uiOutput("concUI"),
            fluidRow(
              actionButton("concAdd", "Add"),
              actionButton("concRem", "Remove"),
              actionButton("concSave", "Save")
            ),  # fluidRow
		        align = "center"
  		    )  # box.concentrations-input
  		  )  # fluidRow
  		),  # patient-tab
      
  ### Bayesian Forecasting Tab
  		tabItem(tabName = "bayes-tab",
  		  
  		  box(title = "Predicted Concentration-Time Profile", width = 12,
  		    plotOutput("plot")
  		  ),  # box.plot-output
  		  
		    box(title = "Dose Planner", width = 12,
		      uiOutput("futUI"),
		      align = "center"
		    ),  # box.future-input
  		  
  		  fluidRow(
    		  valueBox(subtitle = "AUC (Last 24 Hours)", width = 4, 
            color = "light-blue", value = uiOutput("hauc"), 
    		    icon = icon("backward")
    		  ),  # valueBox.hauc-output
    		  valueBox(subtitle = "AUC (Final Forecasted Dose)", width = 4, 
            color = "purple", value = uiOutput("fauc"),
    		    icon = icon("forward")
    		  ),  # valueBox.fauc-output
    		  valueBox(subtitle = "Target AUC", width = 4, 
    		    value = "400 - 700", color = "red", 
    		    icon = icon("bullseye")
    		  )  # valueBox.targetAUC
  		  )  # fluidRow
  		  
  		)  # bayes-tab
      
    )  # tabItem
  )  # dashboardBody

# UI end
  dashboardPage(header, sidebar, body)
  