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
require(ggplot2)

require(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Herd Immunity Estimator"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("R_0",
                     "Basic Reproductive Number:",
                     min = 0.5,
                     max = 3,
                     value = 1.25),
         sliderInput("VaxCov",
                     "Vaccine Coverage:",
                     min = 0,
                     max = 1,
                     value = 1),
         radioButtons("VE_shape", 
                      "Vaccine:",
                      choices = c("Shanchol", "Dukoral"),
                      selected = "Shanchol"),
         sliderInput("mig",
                     "Migration Rate (1/years):",
                     min = 0,
                     max = 1/2,
                     value = 0),
         sliderInput("birth",
                     "Birth/Death Rate (1/years):",
                     min = 0,
                     max = 1/40,
                     value = 0),
         sliderInput("beta_amp",
                     "Amplitude of Seasonality:",
                     min = 0,
                     max = 1,
                     value = 0.2)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("R_t")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$R_t <- renderPlot({
     
     times <- seq(0,356*5)
     
     # Calculate elements of VE
     max_V_months = 48
     V_comps_per_month = .5
     
     n.comps.V = max_V_months*V_comps_per_month
     # VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE)
     VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = input$VE_shape, bound = TRUE, max_V_months = max_V_months)
     # VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect", bound = TRUE, max_V_months = max_V_months)
     
     params <- list(beta=input$R_0/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                    beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                    beta_amp = input$beta_amp,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                    beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                    gamma=1/2,                     # Duration of disease
                    sigma=1/1.4,                   # Incubation period
                    birth_death_rate=input$birth/365, # Average birth and death rate
                    nat_wane=1/(365*10),         # Rate of natural immunity waning
                    mig_in= input$mig/365,             # Rate of immigration
                    mig_out=input$mig/365,             # Rate of emigration
                    foreign_infection=0.00,        # Proportion of immigrants who are infected
                    n.comps.V=n.comps.V,           # Number of V compartments
                    VE=VE,                         # Vaccine efficacy over time
                    V_step=V_comps_per_month/30.5  # Average time in each vaccine compartment is one month
     )
     inits = rep(0, 6+params$n.comps.V)
     inits[1] = 100000*(1-input$VaxCov) # initially susceptible
     inits[2] = 100000*(input$VaxCov) # initially vaccinated
     inits[params$n.comps.V+3] = 0 # initially infected
     
     #### Test function ####
     test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
     
     ggplot(test.run, aes(x = time/365, y = Re)) + geom_line() + geom_hline(yintercept=1, col="red") + 
       theme_bw() + xlab("Years") + ylab("R Effective") + scale_x_continuous(breaks = seq(0, 5, 1)) +
       ylim(0,max(3, input$R_0))
     })
}

# Run the application 
shinyApp(ui = ui, server = server)

