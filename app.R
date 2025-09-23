library(shiny)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(shinydashboard)
library(deSolve)

# Define UI for application that draws a histogram
ui <- navbarPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  title = "Misclassification Model",
  #To use shinydashboard functions with a navbarpage layout
  header = tagList(
    useShinydashboard(),
    withMathJax()
  ),
  tabPanel(
    "Model Assumptions",
    h5(strong("Transmission Model Assumptions")),
    h6("1) Each primary case generates secondary cases according to a negative binomial distribution with mean \\(R_s\\), and dispersion parameter \\(K_s\\)."),
    h6("2) Each secondary case continues to generate later-generation secondary cases with the same \\(R_s\\) and \\(K_s\\) parameters."),
    h6("3) A 'true' chain consists of all cases that can be directly linked to a true primary case."),
    h6("4) If one or more of these chains occur in overlapping time and space, they become entangled and form a 'true cluster'."),
    h6("5) The number of true primary cases in each cluster is modelled as a geometric distribution with mean \\(R_p\\)."),
    br(),
    h5(strong("Surveillance Model Assumptions")),
    h6("1) Each infection is observed with independent probability \\(p_{obs}\\)"),
    h6("2) The surveillance system observes a cluster (a 'reported cluster') and assigns one case as primary and all others as secondary."),
    h6("3) True chains do not span multiple reported clusters."),
    br(),
  ),
  tabPanel(
    "Primary Cases",
    h2("Primary Case Classification"),
    br(),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "NumPrimaryCasesCluster",
          label=div(style="text-align:center",
                    "Number of primary cases per cluster (\\(R_p\\))"),
          min = 1,
          max = 2,
          value = 1.5,
          step = 0.2
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "PObs",
          label=div(style="text-align:center",
                    "Probability of observing each case (\\(P_{obs}\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
    ),
    br(),
    p(strong("Model Focus:"), "Our focus is on ‘weak’ or ‘subcritical’ transmission, which we define as occurring when the effective reproduction number, 
    (\\(R_s\\)), is <1, meaning sustained transmission is not possible. Due to the stochastic and sporadic nature of subcritical transmission, high-quality
    surveillance data is often challenging to obtain. Case misclassification may result from unobserved cases, or from transmission chain entanglement", style = "font-size:12px"),
    h3("Classifier probabilities:"),
    mainPanel(
      (
        fluidRow(
          column(3, uiOutput("Cp_p", width = 3)),
          column(3, uiOutput("Cp_s", width = 3)),
          column(3, uiOutput("Cp_o", width = 3))
        )
      ),
    )
  ),
  tabPanel(
    "Secondary Cases",
    h2("Secondary Case Classification"),
    h5(strong("Select values on the sliders below to see the impact of these variables on the probability of secondary case classification.")),
    br(),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "Rp",
          label=div(style="text-align:center",
                    "Number of primary cases per cluster (\\(R_p\\))"),
          min = 1,
          max = 2,
          value = 1.5,
          step = 0.2
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Rs",
          label=div(style="text-align:center",
                    "Reproductive number (\\(R_s\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      )
      ),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "ks",
          label=div(style="text-align:center",
                    "Dispersion (\\(k_s\\))"),
          min = 0,
          max = 2,
          value = 1.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "PObs",
          label=div(style="text-align:center",
                    "Probability of observing each case (\\(P_{obs}\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      
      br(),
      h3("Classifier probabilities:"),
      fluidRow(
        column(3, uiOutput("Cs_p", width = 3)),
        column(3, uiOutput("Cs_s", width = 3)),
        column(3, uiOutput("Cs_o", width = 3))
      )
    ),
  ),
  tabPanel(
    "Case classification",
    h2("Case classifier performance"),
    h5(strong("Select values on the sliders below to see the impact of these variables on the performance of case classification.")),
    # uiOutput("warning_text2"),
    br(),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "Rp",
          label=div(style="text-align:center",
                    "Number of primary cases per cluster (\\(R_p\\))"),
          min = 1,
          max = 2,
          value = 1.5,
          step = 0.2
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Rs",
          label=div(style="text-align:center",
                    "Reproductive number (\\(R_s\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      )
    ),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "ks",
          label=div(style="text-align:center",
                    "Dispersion (\\(k_s\\))"),
          min = 0,
          max = 2,
          value = 1.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "PObs",
          label=div(style="text-align:center",
                    "Probability of observing each case (\\(P_{obs}\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      
      br(),
      h3("Classification probabilities:"),
      fluidRow(
        column(3, uiOutput("Pp_p", width = 3)),
        column(3, uiOutput("Ps_s", width = 3))
        ),
      fluidRow(
        column(6, uiOutput("theta", width = 6))
      )
    ),
  ),
  tabPanel(
    "Rs bias",
    h2("Bias in Rs inference"),
    h5(strong("Select values on the sliders below to see the impact of these variables on the performance of R estimation.")),
    # uiOutput("warning_text2"),
    br(),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "Rp",
          label=div(style="text-align:center",
                    "Number of primary cases per cluster (\\(R_p\\))"),
          min = 1,
          max = 2,
          value = 1.5,
          step = 0.2
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Rs",
          label=div(style="text-align:center",
                    "Reproductive number (\\(R_s\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      )
    ),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "ks",
          label=div(style="text-align:center",
                    "Dispersion (\\(k_s\\))"),
          min = 0,
          max = 2,
          value = 1.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "PObs",
          label=div(style="text-align:center",
                    "Probability of observing each case (\\(P_{obs}\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      
      br(),
      fluidRow(
        column(4, uiOutput("Robs", width = 4)),
        column(4, uiOutput("delta", width = 4))
      )
    ),
  ),
  tabPanel(
    "Odds ratio",
    h2("Bias in observed odds ratio"),
    h5(strong("Select values on the sliders below to see the impact of these variables on the performance of trait classification.")),
    # uiOutput("warning_text2"),
    br(),
    fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "Rp",
          label=div(style="text-align:center",
                    "Number of primary cases per cluster (\\(R_p\\))"),
          min = 1,
          max = 2,
          value = 1.5,
          step = 0.2
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Rs",
          label=div(style="text-align:center",
                    "Reproductive number (\\(R_s\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "ks",
          label=div(style="text-align:center",
                    "Dispersion (\\(k_s\\))"),
          min = 0,
          max = 2,
          value = 1.5,
          step = 0.1
        )
      )
    ),
      fluidRow(
      column(
        3, offset = 0.75,
        sliderInput(
          "PObs",
          label=div(style="text-align:center",
                    "Probability of observing each case (\\(P_{obs}\\))"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Ip",
          label = div(style = "text-align:center",
                      "Fraction of primary cases that are trait positive"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      column(
        3, offset = 0.75,
        sliderInput(
          "Is",
          label = div(style = "text-align:center",
                      "Fraction of secondary cases that are trait positive"),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        )
      ),
      
      br(),
      fluidRow(
        column(3, uiOutput("Qprime_pp", width = 3)),
        column(3, uiOutput("Qprime_sp", width = 3)),
        column(3, uiOutput("OR", width = 3))
      )
    ),
  )
)


# Define server logic required to draw a histogram
server <- function(input, output)  {
  Cp_p_eval <- reactive({
    input$PObs / (1 - input$PObs + input$PObs * input$NumPrimaryCasesCluster)
  })
  Cp_s_eval <- reactive({
    input$PObs - Cp_p_eval()
  })
  output$Cp_p <- renderValueBox({
    valueBox(
      value = paste0(round(Cp_p_eval(), digits=2)),
      "True primary classified as true primary",
      color = "green")
  })
  output$Cp_s <- renderValueBox({
    valueBox(
      value = paste0(round(Cp_s_eval(), digits=2)),
      "True primary classified as secondary",
      color = "purple")
  })
  output$Cp_o <- renderValueBox({
    valueBox(
      value = paste0(round(1 - input$PObs, digits=2)),
      "True primary is unobserved",
      color="light-blue")
  })
  
  lx_y <- function(x, y, ks, Rs) {
    exp(
      lgamma(y + ks*x) - 
        lgamma(y + 1) - 
        lgamma(ks*x) + 
        ks * x * log(ks/(Rs + ks)) + 
        y * log(Rs/(Rs + ks))
    )
  }
  l_i_n_c <- function(i, n, ks, Rs) {
    (i / n) * lx_y(n, n-i, ks, Rs)
  }
  p_i <- function(i, R_p) {
    ((R_p - 1)/R_p) ^ (i - 1) / R_p
  }
  
  Cs_p_fn <- function(Rp, Rs, pobs, ks, jmax = 200, nmax = 200) {
    n_vals <- 2:max(nmax, jmax)
    i_vals <- 1:(jmax - 1)
    
    l_mat <- outer(i_vals, n_vals, Vectorize(function(i, n) l_i_n_c(i, n, ks, Rs)))
    
    p_vec <- p_i(i_vals, Rp)
    mat <- sweep(l_mat, 1, p_vec, `*`)
    
    cs_n <- apply(mat, 1, function(row) rev(cumsum(rev(row))))
    cs_n <- t(cs_n)
    
    rj_vals <- sapply(2:jmax, function(j) {
      i_range <- 1:(j - 1)
      n_index <- j - min(n_vals) + 1
      sum(cs_n[i_range, n_index]) * (1 - Rs) / (Rp * Rs)
    })
    
    j_vals <- 2:jmax
    sum(rj_vals * pobs * (1 - pobs) ^ (j_vals - 1))
  }
  Cs_p_eval <- reactive({
    Cs_p_fn(input$Rp, input$Rs, input$PObs, input$ks, 200)
  })
  Cs_s_eval <- reactive({
    input$PObs - Cs_p_eval()
  })
  output$Cs_p <- renderValueBox({
    valueBox(
      value = paste0(round(Cs_p_eval(), digits=2)),
      "True secondary classified as true primary",
      color = "green")
  })
  output$Cs_s <- renderValueBox({
    valueBox(
      value = paste0(round(Cs_s_eval(), digits=2)),
      "True secondary classified as secondary",
      color = "purple")
  })
  output$Cs_o <- renderValueBox({
    valueBox(
      value = paste0(round(1 - input$PObs, digits=2)),
      "True secondary is unobserved",
      color="light-blue")
  })
  Pp_p_eval <- reactive({
    Cp_p_eval() * (1 - input$Rs)/(Cp_p_eval() * (1 - input$Rs) + Cs_p_eval() * input$Rs)
  })
  Ps_s_eval <- reactive({
    Cs_s_eval() * input$Rs / (Cs_s_eval() * input$Rs + Cp_s_eval() * (1 - input$Rs))
  })
  theta_eval <- reactive({
    (Cs_s_eval() * input$Rs + Cp_p_eval() * (1 - input$Rs)) / input$PObs
  })
  output$Pp_p <- renderValueBox({
    valueBox(
      value = paste0(round(Pp_p_eval(), digits = 2)),
      "Probability that a primary classification is a true primary case",
      color = "green"
    )
  })
  output$Ps_s <- renderValueBox({
    valueBox(
      value = paste0(round(Ps_s_eval(), digits = 2)),
      "Probability that a secondary classification is a true secondary case",
      color = "purple"
    )
  })
  output$theta <- renderValueBox({
    valueBox(
      value = paste0(round(theta_eval(), digits = 2)),
      "Accuracy that an observed case has the correct primary vs secondary assignment",
      color = "light-blue"
    )
  })
  Robs_eval <- reactive({
    (Cp_s_eval() * (1 - input$Rs) + input$Rs * Cs_s_eval()) / input$PObs
  })
  delta_eval <- reactive({
    Robs_eval() - input$Rs
  })
  output$Robs <- renderValueBox({
    valueBox(
      value = paste0(round(Robs_eval(), digits = 2)),
      "Observed reproduction number",
      color = "green"
    )
  })
  output$delta <- renderValueBox({
    valueBox(
      value = paste0(round(delta_eval(), digits = 2)),
      "Bias of R estimate",
      color = "purple"
    )
  })
  Qpm_eval <- reactive({
    (1 - input$Rs) * (1 - input$Ip)
  })
  Qpp_eval <- reactive({
    (1 - input$Rs) * input$Ip
  })
  Qsm_eval <- reactive({
    input$Rs * (1 - input$Is)
  })
  Qsp_eval <- reactive({
    input$Rs * input$Is
  })
  Qprime_pm_eval <- reactive({
    (Qpm_eval() * Cp_p_eval() + Qsm_eval() * Cs_p_eval()) / input$PObs
  })
  Qprime_pp_eval <- reactive({
    (Qpp_eval() * Cp_p_eval() + Qsp_eval() * Cs_p_eval()) / input$PObs
  })
  Qprime_sm_eval <- reactive({
    (Qpm_eval() * Cp_s_eval() + Qsm_eval() * Cs_s_eval()) / input$PObs
  })
  Qprime_sp_eval <- reactive({
    (Qpp_eval() * Cp_s_eval() + Qsp_eval() * Cs_s_eval()) / input$PObs
  })
  OR_eval <- reactive({
    Qprime_pp_eval() * Qprime_sm_eval() / (Qprime_pm_eval() * Qprime_sp_eval())
  })
  output$Qprime_pp <- renderValueBox({
    valueBox(
      value = paste0(round(Qprime_pp_eval(), digits = 2)),
      "Probability that a randomly observed case is primary and trait positive",
      color = "green"
    )
  })
  output$Qprime_sp <- renderValueBox({
    valueBox(
      value = paste0(round(Qprime_sp_eval(), digits = 2)),
      "Probability that a randomly observed case is secondary and trait positive",
      color = "purple"
    )
  })
  output$OR <- renderValueBox({
    valueBox(
      value = paste0(round(OR_eval(), digits = 2)),
      "Observed odds ratio that a trait positive case is primary",
      color = "light-blue"
    )
  })
  
  ##
  
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
