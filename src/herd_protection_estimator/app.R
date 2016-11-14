#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
source("src/revaccination.R")
require(ggplot2)

require(shiny)

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
  
  hr(),
  
 # Sidebar with a slider input for number of bins 
 fluidRow(
    column(3,
      # submitButton("Submit"),
       sliderInput("R_0",
                   "Basic Reproductive Number:",
                   min = 0.5,
                   max = 3,
                   value = 1.25,
                   step = 0.05),
       sliderInput("VaxCov",
                   "Initial Vaccine Coverage:",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.01),
       radioButtons("VE_shape", 
                    "Vaccine:",
                    choices = c("Dukoral", "Shanchol"),
                    selected = "Dukoral"),
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
         checkboxGroupInput("vac_recip",
                            "Ongoing Vaccine Recipients:", 
                            choices = c("Births" = "birth", "Migrants" = "migrant", "Routine for Susceptibles" = "S", "Routine for All" = "all"),
                            selected = c("all")),
         sliderInput("vac_routine_freq",
                 "Frequency of Routine Revaccination (years):",
                 min = 0,
                 max = 10,
                 value = 0,
                 step = 1),
       sliderInput("vac_routine_frac",
                   "Fraction Vaccinated (Routinely):",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.01),
       sliderInput("vac_birth_frac",
                   "Fraction Vaccinated (Birth):",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.01),
       sliderInput("vac_mig_frac",
                   "Fraction Vaccinated (Immigrants):",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.01),
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
                     mig_in= 1/(input$mig*365),             # Rate of immigration
                     mig_out=1/(input$mig*365),             # Rate of emigration
                     foreign_infection=0.00,        # Proportion of immigrants who are infected
                     n.comps.V=n.comps.V,           # Number of V compartments
                     VE=VE,                         # Vaccine efficacy over time
                     V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                     
                     vac_routine_freq = round(input$vac_routine_freq*365),                  # Days between re-vaccination campaigns
                     vac_routine_frac = input$vac_routine_frac,
                     vac_birth_frac = input$vac_birth_frac,            # Fraction of babies vaccinated
                     vac_mig_frac = input$vac_mig_frac,            # Fraction of immigrants vaccinated upon arrival
                     vac_max = input$vac_max,                 # Maximum number of vaccines to be given
                     vac_recip = input$vac_recip  # Recipients of vaccination ("all", "S", "migrant", "birth")
      )
      inits = rep(0, 7+params$n.comps.V)
      inits[1] = 100000*(1-input$VaxCov) # initially susceptible
      inits[2] = 100000*(input$VaxCov) # initially vaccinated
      inits[params$n.comps.V+3] = 0 # initially infected
      inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
      
      #### Test function ####
      test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
    })
  
  
    output$R_t <- renderPlot({
      if (input$outcome_of_interest %in% "Re"){
        ggplot(model(), aes(x = time/365, y = Re)) + geom_line()  + 
          theme_bw() + xlab("Years since Vaccination") + ylab("R Effective") + scale_x_continuous(breaks = seq(0, 5, 1)) +
          ylim(0,max(3, input$R_0)) + geom_hline(yintercept=1, col="red")
      } else if (input$outcome_of_interest %in% "prob"){
        ggplot(model(), aes(x = time/365, y = prob_outbreak_fcn(Re, input$outbreak_size))) + geom_line()  + 
          theme_bw() + xlab("Years since Vaccination") + ylab(paste("Probability of >", input$outbreak_size, "cases")) + scale_x_continuous(breaks = seq(0, 5, 1)) + ylim(0,1)
      }
     })
   
   output$lose_herd <- renderText({
     paste("Herd Immunity is first lost after: ", round(which(model()$Re > 1)[1]/365,2), " years.")
   })  
   
   output$DHI <- renderText({
     paste("The Cumulative Duration of Herd Immunity is: ", round(sum(round(model()$Re,3) < 1)/365, 2), " years.")
   })
   
   output$vax_consumed <- renderText({
     paste("The number of vaccine courses consumed is: ", round(max(model()$Vax)), "(", round(100*max(model()$Vax)/input$vac_max), "% of the", input$vac_max, "allocated )")
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

