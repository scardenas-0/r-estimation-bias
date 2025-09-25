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
      )
    ),
    
    br(),
    p(strong("Model Focus:"), "Our focus is on ‘weak’ or ‘subcritical’ transmission, which we define as occurring when the effective reproduction number, 
    (\\(R_s\\)), is <1, meaning sustained transmission is not possible. Due to the stochastic and sporadic nature of subcritical transmission, high-quality
    surveillance data is often challenging to obtain. Case misclassification may result from unobserved cases, or from transmission chain entanglement", style = "font-size:12px"),
    h3("Classifier probabilities:"),
    fluidRow(
      column(9, actionButton("go", "Run calculations"))
    ),
    br(),
    fluidRow(
      column(3, uiOutput("Cp_p", width = 3)),
      column(3, uiOutput("Cp_s", width = 3)),
      column(3, uiOutput("Cp_o", width = 3))
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
        column(9, actionButton("go", "Run calculations"))
      ),
      br(),
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
    uiOutput("warning_text1"),
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
        column(9, actionButton("go", "Run calculations"))
      ),
      br(),
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
        column(9, actionButton("go", "Run calculations"))
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
        column(9, actionButton("go", "Run calculations"))
      ),
      br(),
      fluidRow(
        column(3, uiOutput("Qprime_pp", width = 3)),
        column(3, uiOutput("Qprime_sp", width = 3)),
        column(3, uiOutput("OR", width = 3))
      )
    ),
  ),
  tabPanel(
    "Odds ratio data",
    h2("Bias in observed odds ratio from data"),
    h5(strong("Select values on the sliders below to see the impact of these variables on the performance of trait classification.")),
    uiOutput("warning_text2"),
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
          "ks",
          label=div(style="text-align:center",
                    "Dispersion (\\(k_s\\))"),
          min = 0,
          max = 2,
          value = 0.5,
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
    ),
    fluidRow(
      column(
        3, offset = 0.75,
        numericInput(
          "TotalCases",
          label = div(style = "text-align:center",
                      "Total number of cases"),
          value = 829,
          step = 1
        )
      ),
      column(
        3, offset = 0.75,
        numericInput(
          "SecCases",
          label = div(style = "text-align:center",
                      "Number of secondary cases"),
          value = 533,
          step = 1
        )
      ),
      column(
        3, offset = 0.75,
        numericInput(
          "ppCases",
          label = div(style = "text-align:center",
                      "Number of trait-positive primary cases"),
          value = 161,
          step = 1
        )
      ),
      column(
        3, offset = 0.75,
        numericInput(
          "spCases",
          label = div(style = "text-align:center",
                      "Number of trait-positive secondary cases"),
          value = 238,
          step = 1
        )
      )
    ),
    br(),
    fluidRow(
      column(9, actionButton("go", "Run calculations"))
    ),
    br(),
    fluidRow(
      column(3, uiOutput("Prop_pp", width = 3)),
      column(3, uiOutput("Prop_sp", width = 3)),
      column(3, uiOutput("OR_inf", width = 3))
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output)  {
  Cp_p_eval <- eventReactive(input$go, {
    input$PObs / (1 - input$PObs + input$PObs * input$NumPrimaryCasesCluster)
  })
  Cp_s_eval <- eventReactive(input$go, {
    input$PObs - Cp_p_eval()
  })
  Cx_o_eval <- eventReactive(input$go, {
    1 - input$PObs
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
      value = paste0(round(Cx_o_eval(), digits=2)),
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
  Cs_p_eval <- eventReactive(input$go, {
    Cs_p_fn(input$Rp, input$Rs, input$PObs, input$ks, 200)
  })
  Cs_s_eval <- eventReactive(input$go, {
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
      value = paste0(round(Cx_o_eval(), digits=2)),
      "True secondary is unobserved",
      color="light-blue")
  })
  Pp_p_eval <- eventReactive(input$go, {
    Cp_p_eval() * (1 - input$Rs)/(Cp_p_eval() * (1 - input$Rs) + Cs_p_eval() * input$Rs)
  })
  Ps_s_eval <- eventReactive(input$go, {
    Cs_s_eval() * input$Rs / (Cs_s_eval() * input$Rs + Cp_s_eval() * (1 - input$Rs))
  })
  theta_eval <- eventReactive(input$go, {
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
  Robs_eval <- eventReactive(input$go, {
    (Cp_s_eval() * (1 - input$Rs) + input$Rs * Cs_s_eval()) / input$PObs
  })
  delta_eval <- eventReactive(input$go, {
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
  Qprime_pp_eval <- eventReactive(input$go, {
    (Qpp_eval() * Cp_p_eval() + Qsp_eval() * Cs_p_eval()) / input$PObs
  })
  Qprime_sm_eval <- reactive({
    (Qpm_eval() * Cp_s_eval() + Qsm_eval() * Cs_s_eval()) / input$PObs
  })
  Qprime_sp_eval <- eventReactive(input$go, {
    (Qpp_eval() * Cp_s_eval() + Qsp_eval() * Cs_s_eval()) / input$PObs
  })
  OR_eval <- eventReactive(input$go, {
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
  snCases <- reactive({
    (input$SecCases - input$spCases) / input$TotalCases
  })
  pnCases <- reactive({
    (input$TotalCases - input$SecCases - input$ppCases) / input$TotalCases
  })
  spCases <- reactive({
    input$spCases / input$TotalCases
  })
  ppCases <- reactive({
    input$ppCases / input$TotalCases
  })
  Cs_p_eval2 <- reactive({
    Cs_p_fn(input$Rp, input$SecCases / input$TotalCases, input$PObs, input$ks, 200)
  })
  Cs_s_eval2 <- reactive({
    input$PObs - Cs_p_eval2()
  })
  data_OR_denom <- reactive({
    Cp_s_eval() * Cs_p_eval2() - Cp_p_eval() * Cs_s_eval2()
  })
  Q_pn_data_eval <- reactive({
    input$PObs * (snCases() * Cs_p_eval2() - pnCases() * Cs_s_eval2()) / data_OR_denom()
  })
  Q_pp_data_eval <- reactive({
    input$PObs * (spCases() * Cs_p_eval2() - ppCases() * Cs_s_eval2()) / data_OR_denom()
  })
  Q_sn_data_eval <- reactive({
    -input$PObs * (snCases() * Cp_p_eval() - pnCases() * Cp_s_eval()) / data_OR_denom()
  })
  Q_sp_data_eval <- reactive({
    -input$PObs * (spCases() * Cp_p_eval() - ppCases() * Cp_s_eval()) / data_OR_denom()
  })
  Prop_pp_eval <- eventReactive(input$go, {
    Q_pp_data_eval() / (Q_pn_data_eval() + Q_pp_data_eval())
  })
  Prop_sp_eval <- eventReactive(input$go, {
    Q_sp_data_eval() / (Q_sn_data_eval() + Q_sp_data_eval())
  })
  OR_inf_eval <- eventReactive(input$go, {
    (Q_pp_data_eval() * Q_sn_data_eval()) / (Q_pn_data_eval() * Q_sp_data_eval())
  })
  output$Prop_pp <- renderValueBox({
    valueBox(
      value = paste0(round(
        Prop_pp_eval(), digits = 2
      )),
      "Inferred proportion of primary cases that are positive",
      color = "green"
    )
  })
  output$Prop_sp <- renderValueBox({
    valueBox(
      value = paste0(round(
        Prop_sp_eval(),
        digits = 2
      )),
      "Inferred proportion of secondary cases that are positive",
      color = "purple"
    )
  })
  output$OR_inf <- renderValueBox({
    valueBox(
      value = paste0(round(
        OR_inf_eval(), digits = 2
      )),
      "Inferred odds ratio for positive cases being primary",
      color = "light-blue"
    )
  })
  
  
  output$warning_text1 <- renderUI({
    HTML(paste0("<f><font color = red><i>", "These models are intended for exploratory analysis only. Please refer to the", 
                # "</f></font></i>", "<a href=https://doi.org/10.1101/2021.07.05.21260043>",
                " corresponding publication ", "</a>", 
                "<f><font color = red><i>", "for description of modeling limitations. 
        Please involve an expert in computational modeling for any policy-decision making that is influenced by this web interface.", "</f></font></i>"))
  })
  
  output$warning_text2 <- renderUI({
    HTML(paste0("<f><font color = red><i>", "Ensure that R_p is not too large compared to your data...", "</f></font></i>"))
  })
  
}


# Run the application
shinyApp(ui = ui, server = server)
