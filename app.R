library(shiny)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(shinydashboard)
library(deSolve)

# Define UI for application that draws a histogram
ui <- navbarPage(theme = bslib::bs_theme(bootswatch = "flatly"),
                 title = "Misclassification Model",
                 #To use shinydashboard functions with a navbarpage layout
                 header = tagList(
                   useShinydashboard(),
                   withMathJax()
                 ),
                 tabPanel("Model Assumptions",
                          h5(strong("Transmission Model Assumptions")),
                          h6("1) number of primary cases occurring in each infection cluster is modelled as a geometric distribution with mean \\(R_p\\)"),
                          h6("2) each primary case generates secondary cases according to a negative binomial distribution with mean \\(R_s\\)), and dispersion parameter \\(K_s\\)"),
                          h6("3) Each secondary cases, continues to generate later-generation secondary cases with the same \\(R_s\\) and \\(K_s\\) parameters."),
                          br(),
                          h5(strong("Surveillance Model Assumptions")),
                          h6("1) number of primary cases occurring in each infection cluster is modelled as a geometric distribution with mean \\(R_p\\)"),
                          h6("2) each primary case generates secondary cases according to a negative binomial distribution with mean \\(R_s\\)), and dispersion parameter \\(K_s\\)"),
                          h6("3) Each secondary cases, continues to generate later-generation secondary cases with the same \\(R_s\\) and \\(K_s\\) parameters."),
                          
                          br(),
                 ),
                 tabPanel("Primary Cases",
                          h2("Primary Case Classification"),
                          br(),
                          
                          fluidRow(
                            column(2, offset=0.75, 
                                   sliderInput("NumPrimaryCasesCluster",
                                               label=div(style="text-align:center",
                                                         "Number of Primary Cases in a Cluster (\\(R_p\\))"),
                                               min = 1,
                                               max = 2,
                                               value = 1.5,
                                               step = 0.2)),
                            column(2, offset=0.75,
                                   sliderInput("PObs",
                                               label=div(style="text-align:center",
                                                         "Independent Probability of Observing Each Infection(\\(P_{obs}\\))"),
                                               min=0,
                                               max=1,
                                               value=0.5,
                                               step = 0.1)),
                          ),
                          br(),
                          p(strong("Model Focus:"), "Our focus is on ‘weak’ or ‘subcritical’ transmission, which we define as occurring when the effective reproduction number, 
                          (\\(R_s\\)), is <1, meaning sustained transmission is not possible. stochastic and often sporadic nature of subcritical transmission, high-quality
                          surveillance data is often challenging to obtain. Case misclassification may result from unobserved cases, or from transmission chain entanglement", style = "font-size:12px"),
                          
                          mainPanel((fluidRow(                              
                            column(3, valueBoxOutput("Cp_p", width=NULL)),
                            column(3, valueBoxOutput("Cp_s", width=NULL)),
                            column(3, valueBoxOutput("Cp_o", width=NULL)))
                          ),
                          )
                 ),
                 
                 
                 
                 tabPanel("Secondary Cases",
                          h2("Secondary Case Classification"),
                          h5("Introduction of disease may not lead to an outbreak depending on the stochasticity of transmission.
                             We refer to the probability that an outbreak does not occur as the probability of self-limited spread, (\\(P_{e}\\))."),
                          h5(strong("Select values on the sliders below to see the impact of these variables on the probability of secondary case classification.")),
                          uiOutput("warning_text2"),
                          br(),
                          fluidRow(
                            column(2, offset=0.75, 
                                   sliderInput("R",
                                               label=div(style="text-align:center",
                                                         "Reproductive number (R)"),
                                               min = 0.1,
                                               max = 10,
                                               value = 2.5,
                                               step = 0.5)),
                            column(2, offset=0.75,
                                   shinyWidgets::sliderTextInput('k', 
                                                                 label=div(style="text-align:center",
                                                                           'Dispersion parameter (k)'),
                                                                 choices=c(0.1, 0.2, 0.5, 5),
                                                                 selected=0.2, grid=T)),
                            column(2, offset=0.75,
                                   sliderInput("Cth",
                                               label=div(style="text-align:center",
                                                         "Outbreak size threshold (\\(C_{th}\\))*"),
                                               min = 1,
                                               max = 20,
                                               value = 2.5)),
                            column(2, offset=0.75,
                                   sliderInput("Z",
                                               label=div(style="text-align:center",
                                                         "Number of disease introductions (\\(\\zeta\\))**"),
                                               min=1,
                                               max=5,
                                               value=1))),
                          
                          br(),
                          fluidRow(
                            column(5, uiOutput('pe_box', width=3)),
                            column(5, uiOutput('pe_box_outbreak',width=3))),
                 ),
)



# Define server logic required to draw a histogram
server <- function(input, output)  {
  output$Cp_p <- renderValueBox({
    valueBox(
      value = paste0(round(input$numsusceptibles*input$avgcontacts, digits=2)),
      "Classifier Probability: True primary classified as true primary",
      color = "green")
  })
  output$Cp_s <- renderValueBox({
    valueBox(
      value = paste0(round(input$disprev*input$probinfection, digits=2)),
      "Classifier Probability: True primary classified as secondary",
      color = "purple")
  })
  output$Cp_o <- renderValueBox({
    valueBox(
      value = paste0(round(input$numsusceptibles*input$avgcontacts*input$disprev*input$probinfection, digits=2)),
      "Classifier Probability: True primary is unobserved",
      color="light-blue")
  })
  
  clust_size <- reactive({
    input$Z:input$Cth
  })
  
  log_s_ij_arr <- reactive({
    exp(lgamma((clust_size()-input$Z)+input$k*clust_size()) - lgamma((clust_size()-input$Z)+1) - lgamma(input$k*clust_size()) +
          input$k * clust_size() * log(input$k/(input$R+input$k)) +
          (clust_size()-input$Z) * log(input$R/(input$R+input$k)))
  })
  
  s_ij_arr <- reactive({
    log_s_ij_arr()/clust_size()*input$Z
  })
  pe <- reactive({
    sum(s_ij_arr())*100
  })
  
  pe_outbreak <- reactive({
    100-pe()  
  })
  
  output$pe_box <- renderValueBox({
    valueBox(value=paste0(signif(pe(),digits=4),"%"),
             "Probability of self-limited spread",
             icon=icon("times"),
             color="green",
             width=3)
  })
  
  output$pe_box_outbreak <- renderValueBox({
    valueBox(value=paste0(signif(pe_outbreak(),digits=4),"%"),
             "Probability of outbreak occurring",
             icon=icon("check"),
             color="red",
             width=3)
  })
  
  SEIR_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -(beta/gamma) * (S/1000) * I
      dE <- (beta/gamma) * (S/1000) * I - (1/gamma) * E
      dI <-  (1/gamma) * E - (1/delta)*I
      dR <-  (1/delta) * I
      return(list(c(dS, dE, dI, dR)))
    })
  }
  
  SEIR_values_1 <- reactive({
    req(input$time_values, input$beta, input$gamma, input$delta)
    ode(y = c(S = 1000, E=1, I = 0, R = 0),
        times = seq(input$time_values[1], input$time_values[2]),
        func = SEIR_equations,
        parms = c(beta = input$beta, gamma = input$gamma, delta=input$delta))
  })
  
  output$SEIR_plot <- renderPlot({
    val <- as.data.frame(SEIR_values_1())
    with(val, {
      plot(time, S, type = "l", col = "blue",
           xlab = "Time period (days)", ylab = "Number of individuals",
           ylim = c(0,1000), lwd=4)
      lines(time, E, col="orange", lwd=4)
      lines(time, I, col = "red", lwd=4)
      lines(time, R, col = "green", lwd=4)
    })
    legend("topright", legend=c("Susceptible", "Exposed", "Infectious", "Recovered"),
           col = c("blue", "orange", "red", "green"), lwd =4, bty = "n")
  })
  
  
  output$warning_text1 <- renderUI({
    HTML(paste0("<f><font color = red><i>", "These models are intended for exploratory analysis only. Please refer to the", 
                "</f></font></i>", "<a href=https://doi.org/10.1101/2021.07.05.21260043>"," corresponding publication ", "</a>", 
                "<f><font color = red><i>", "for description of modeling limitations. 
        Please involve an expert in computational modeling for any policy-decision making that is influenced by this web interface.", "</f></font></i>"))
  })
  
  output$warning_text2 <- renderUI({
    HTML(paste0("<f><font color = red><i>", "These models are intended for exploratory analysis only. Please refer to the", 
                "</f></font></i>", "<a href=https://doi.org/10.1101/2021.07.05.21260043>"," corresponding publication ", "</a>", 
                "<f><font color = red><i>", "for description of modeling limitations. 
        Please involve an expert in computational modeling for any policy-decision making that is influenced by this web interface.", "</f></font></i>"))
  })
  
  output$warning_text3 <- renderUI({
    HTML(paste0("<f><font color = red><i>", "These models are intended for exploratory analysis only. Please refer to the", 
                "</f></font></i>", "<a href=https://doi.org/10.1101/2021.07.05.21260043>"," corresponding publication ", "</a>", 
                "<f><font color = red><i>", "for description of modeling limitations. 
        Please involve an expert in computational modeling for any policy-decision making that is influenced by this web interface.", "</f></font></i>"))
  })
}


# Run the application
shinyApp(ui = ui, server = server)
