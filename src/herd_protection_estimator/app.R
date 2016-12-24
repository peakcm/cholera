#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
require(deSolve)
require(ggplot2)
require(shiny)

source_https <- function(url, ...) {
  # load package
  require(RCurl)
  
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}


#### Load libraries and functions ####
source_https("https://raw.githubusercontent.com/peakcm/cholera/master/src/calculate_VE.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/calculate_Re.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/Seasonality.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/prob_outbreak_fcn.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/SIRV_model.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/Run_SIRV_model.R",
             "https://raw.githubusercontent.com/peakcm/cholera/master/src/revaccination.R")

# setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
# source("src/calculate_Re.R")
# source("src/calculate_VE.R")
# source("src/Seasonality.R")
# source("src/prob_outbreak_fcn.R")
# source("src/SIRV_model.R")
# source("src/Run_SIRV_model.R")
# source("src/revaccination.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(tags$style(HTML("
        .selectize-input, .selectize-dropdown {
                            font-size: 50%;
                            }
                            "))),
   # Application title
   titlePanel("Herd Immunity Estimator"),
  
  plotOutput("R_t"),
  
  textOutput("lose_herd"),
  textOutput("DHI"),
  textOutput("vax_consumed"),
  
  submitButton("Submit"),
  
  
  hr(),
  
 # Sidebar with a slider input for number of bins 
 fluidRow(
    column(3,
       sliderInput("R_0",
                   "Basic Reproductive Number:",
                   min = 0.5,
                   max = 3,
                   value = 1.25,
                   step = 0.05),
       radioButtons("VE_shape", 
                    "Vaccine:",
                    choices = c("Shanchol", "Dukoral", "Perfect"),
                    selected = "Shanchol"),
      radioButtons("outcome_of_interest",
                   "Outcome of Interest:",
                   choices=c("Effective Reproductive Number" = "Re", "Probability of Outbreak" = "prob"),
                   selected = "Re"),
      sliderInput("outbreak_size",
                  "Minimum size of 'Outbreak':",
                  min = 1, 
                  max = 100,
                  value = 10)),
    column(3, offset = 1,
       sliderInput("mig",
                   "Average Residence Time (years):",
                   min = 1,
                   max = 40,
                   value = 20,
                   step = 1),
       sliderInput("birth",
                   "Life Expectancy (years):",
                   min = 40,
                   max = 100,
                   value = 70,
                   step = 10),
       sliderInput("beta_amp",
                   "Amplitude of Seasonality:",
                   min = 0,
                   max = 1,
                   value = 0.2,
                   step = 0.1)),
    column(3, offset = 1,
           sliderInput("VaxCov",
                       "Initial Vaccine Coverage:",
                       min = 0,
                       max = 1,
                       value = 1,
                       step = 0.01),
         sliderInput("vac_routine_count",
                     "Number of Routine Vaccine Courses per Day:",
                     min = 0,
                     max = 200,
                     value = 0,
                     step = 10),
         sliderInput("vac_mass_freq",
                 "Frequency of Mass Vaccination (years):",
                 min = 0,
                 max = 10,
                 value = 0,
                 step = 1),
        sliderInput("vac_mass_frac",
                   "Fraction Vaccinated during Mass Campaigns:",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.01),
       # sliderInput("vac_birth_frac",
       #             "Fraction Vaccinated (Birth):",
       #             min = 0,
       #             max = 1,
       #             value = 1,
       #             step = 0.01),
       # sliderInput("vac_mig_frac",
       #             "Fraction Vaccinated (Immigrants):",
       #             min = 0,
       #             max = 1,
       #             value = 1,
       #             step = 0.01),
       sliderInput("vac_max",
                   "Maximum Vaccine Courses Allocated:",
                   min = 0,
                   max = 5e5,
                   value = 5e5,
                   step = 1e4)
    )
 )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
    times <- seq(0,356*5)
    # Calculate elements of VE
    max_V_months = 48
    V_comps_per_month = 0.5 # choose from 0.5, 1, 1.5, etc.
    n.comps.V = max_V_months*V_comps_per_month
    N = 1e5
    
    model <- reactive({
      
      VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = input$VE_shape, bound = TRUE, max_V_months = max_V_months)

      params <- list(beta=input$R_0/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                     beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                     beta_amp = input$beta_amp,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                     beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                     gamma=1/2,                     # Duration of disease
                     sigma=1/1.4,                   # Incubation period
                     birth_death_rate=1/(input$birth*365), # Average birth and death rate
                     nat_wane=1/(365*10),         # Rate of natural immunity waning
                     mig_rates_constant = TRUE,      # TRUE if migration rates are constant
                     mig_in= 1/(input$mig*365),             # Rate of immigration
                     mig_out=1/(input$mig*365),             # Rate of emigration
                     foreign_infection=0.00,        # Proportion of immigrants who are infected
                     n.comps.V=n.comps.V,           # Number of V compartments
                     VE=VE,                         # Vaccine efficacy over time
                     V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                     vac_routine_count = input$vac_routine_count,         # Number of courses given each day
                     
                     vac_mass_freq = round(input$vac_mass_freq*365),                  # Days between re-vaccination campaigns
                     vac_mass_frac = input$vac_mass_frac,
                     vac_birth_frac = 0,            # Fraction of babies vaccinated
                     vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
                     vac_max = input$vac_max,                 # Maximum number of vaccines to be given
                     vac_recip = c("routine_S", "mass_all"),  # Recipients of vaccination ("all", "S", "migrant", "birth")
                     vac_stopper = 1e10          # Don't vaccinate after this day
      )
      inits = rep(0, 7+params$n.comps.V)
      inits[1] = N*(1-input$VaxCov) # initially susceptible
      inits[2] = N*(input$VaxCov) # initially vaccinated
      inits[params$n.comps.V+3] = 0 # initially infected
      inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
      
      #### Test function ####
      test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
    })
    
    model_no_vax <- reactive({
      
      VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = input$VE_shape, bound = TRUE, max_V_months = max_V_months)
      
      params <- list(beta=input$R_0/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                     beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                     beta_amp = input$beta_amp,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                     beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                     gamma=1/2,                     # Duration of disease
                     sigma=1/1.4,                   # Incubation period
                     birth_death_rate=1/(input$birth*365), # Average birth and death rate
                     nat_wane=1/(365*10),         # Rate of natural immunity waning
                     mig_rates_constant = TRUE,      # TRUE if migration rates are constant
                     mig_in= 1/(input$mig*365),             # Rate of immigration
                     mig_out=1/(input$mig*365),             # Rate of emigration
                     foreign_infection=0.00,        # Proportion of immigrants who are infected
                     n.comps.V=n.comps.V,           # Number of V compartments
                     VE=VE,                         # Vaccine efficacy over time
                     V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                     
                     vac_routine_count = 0,          # Fraction of the current S pop that is vaccinated on each day
                     
                     vac_mass_freq = 0,                  # Days between re-vaccination campaigns
                     vac_mass_frac = 0,
                     vac_birth_frac = 0,            # Fraction of babies vaccinated
                     vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
                     vac_max = 0,                 # Maximum number of vaccines to be given
                     vac_recip = c("routine_S", "mass_S"),  # Recipients of vaccination ("all", "S", "migrant", "birth")
                     vac_stopper = 1e10          # Don't vaccinate after this day
      )
      inits = rep(0, 7+params$n.comps.V)
      inits[1] = N # initially susceptible
      inits[2] = 0 # initially vaccinated
      inits[params$n.comps.V+3] = 0 # initially infected
      inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
      
      #### Test function ####
      test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
    })
  
  
    output$R_t <- renderPlot({
      if (input$outcome_of_interest %in% "Re"){
        ggplot(model(), aes(x = time/365, y = Re)) + geom_line()  + 
          theme_bw() + xlab("Years since Initial Vaccination") + ylab("Effective Reproductive Number") + scale_x_continuous(breaks = seq(0, 5, 1)) +
          ylim(0,max(3, input$R_0)) + geom_hline(yintercept=1, col="red") +
          geom_line(data = model_no_vax(), aes(x = time/365, y = Re), col = "grey", lty = "dashed")
      } else if (input$outcome_of_interest %in% "prob"){
        ggplot(model(), aes(x = time/365, y = prob_outbreak_fcn(Re, input$outbreak_size))) + geom_line()  + 
          theme_bw() + xlab("Years Since Initial Vaccination") + ylab(paste("Probability of >", input$outbreak_size, "cases")) + scale_x_continuous(breaks = seq(0, 5, 1)) + ylim(0,1) +
          geom_line(data = model_no_vax(), aes(x = time/365, y = prob_outbreak_fcn(Re, input$outbreak_size)), col = "grey", lty = "dashed")
      }
     })
   
   output$lose_herd <- renderText({
     paste("Herd Immunity is first lost after: ", round(which(model()$Re > 1)[1]/365,2), " years.")
   })  
   
   output$DHI <- renderText({
     paste("Herd Immunity gained by Vaccination: ", round(sum(round(model()$Re,3) < 1)/365, 2) - round(sum(round(model_no_vax()$Re,3) < 1)/365, 2), "years.")
   })
   
   output$vax_consumed <- renderText({
     paste("The number of vaccine courses consumed is: ", round(max(model()$Vax)), "(", round(100*max(model()$Vax)/input$vac_max), "% of the", input$vac_max, "allotted )")
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

