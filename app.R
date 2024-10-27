library(shiny)
library(ggplot2)
library(dplyr)
library(RANN)

# UI definition
ui <- fluidPage(
  titlePanel("Population Stability Index (PSI) Demo"),
  
  sidebarLayout(
    sidebarPanel(
      # Separate sample sizes for each distribution
      h4("Sample Sizes"),
      numericInput("sample_size_p", 
                   "Reference Sample Size:", 
                   value = 10000,
                   min = 1000,
                   max = 100000),
      numericInput("sample_size_q", 
                   "Comparison Sample Size:", 
                   value = 10000,
                   min = 1000,
                   max = 100000),
      
      # Reference distribution parameters
      h4("Reference Distribution (P)"),
      numericInput("alpha_p", 
                   withMathJax("\\(\\alpha_P\\):"), 
                   value = 2,
                   min = 0.1,
                   max = 10,
                   step = 0.1),
      numericInput("beta_p", 
                   withMathJax("\\(\\beta_P\\):"), 
                   value = 100,
                   min = 0.1,
                   max = 100,
                   step = 0.1),
      
      # Comparison distribution parameters
      h4("Comparison Distribution (Q)"),
      numericInput("alpha_q", 
                   withMathJax("\\(\\alpha_Q\\):"), 
                   value = 2,
                   min = 2,
                   max = 10,
                   step = 0.1),
      numericInput("beta_q", 
                   withMathJax("\\(\\beta_Q\\):"), 
                   value = 100,
                   min = 0.1,
                   max = 100,
                   step = 0.1),
      
      # Binning controls
      h4("Binning Options"),
      selectInput("binning_method",
                  "Binning Method:",
                  choices = c("Equal Width" = "equal_width",
                              "Quantile Based" = "quantile"),
                  selected = "equal_width"),
      sliderInput("n_bins",
                  "Number of Bins:",
                  min = 5,
                  max = 1000,
                  value = 100),
      
      # Add checkbox for true density curves
      checkboxInput("show_true_density", 
                    "Show True Density Curves", 
                    value = TRUE),
      
      # Resample button
      actionButton("resample", "Resample Data", 
                   class = "btn-primary",
                   style = "margin-top: 20px; width: 100%")
    ),
    
    mainPanel(
      plotOutput("density_plot"),
      h4("PSI Values"),
      verbatimTextOutput("psi_results")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  # Create reactive value for random seed
  rv <- reactiveValues(seed = 123)
  
  # Update seed when resample button is clicked
  observeEvent(input$resample, {
    rv$seed <- sample.int(1e5, 1)
  })
  
  # Function to get breaks based on method
  get_breaks <- function(data_p, data_q, n_bins, method = "equal_width") {
    combined_data <- c(data_p, data_q)
    
    if (method == "equal_width") {
      # Equal width bins
      breaks <- seq(min(combined_data), max(combined_data), length.out = n_bins + 1)
    } else {
      # Quantile-based bins
      probs <- seq(0, 1, length.out = n_bins + 1)
      breaks <- quantile(combined_data, probs = probs, type = 1)
    }
    
    return(unique(breaks))
  }
  
  # Function to calculate histogram data with custom breaks
  calculate_histogram <- function(data, breaks) {
    hist_data <- hist(data, breaks = breaks, plot = FALSE)
    # Create data frame with break points and densities
    n_bins <- length(hist_data$density)
    
    # For steps, we need n_bins + 1 x-coordinates (all break points)
    # and n_bins densities repeated
    data.frame(
      x = breaks,
      density = c(hist_data$density, hist_data$density[n_bins])
    )
  }
  
  # Function to calculate PSI using histograms
  calculate_PSI <- function(data_p, data_q, bins, method = "equal_width") {
    breaks <- get_breaks(data_p, data_q, bins, method)
    
    p_counts <- hist(data_p, breaks = breaks, plot = FALSE)$counts
    q_counts <- hist(data_q, breaks = breaks, plot = FALSE)$counts
    
    p_prop <- p_counts / sum(p_counts)
    q_prop <- q_counts / sum(q_counts)
    
    epsilon <- 1e-10
    p_prop <- pmax(p_prop, epsilon)
    q_prop <- pmax(q_prop, epsilon)
    
    PSI <- sum((p_prop - q_prop) * log(p_prop / q_prop))
    return(list(PSI = PSI, breaks = breaks))
  }
  
  # KNN-based KL divergence estimation
  kl_divergence_knn <- function(x, y, k = 5) {
    n <- length(x)
    m <- length(y)
    
    nn_x <- nn2(matrix(x, ncol = 1), matrix(x, ncol = 1), k = k + 1)
    dist_x <- nn_x$nn.dists[, k + 1]
    
    nn_y <- nn2(matrix(y, ncol = 1), matrix(x, ncol = 1), k = k)
    dist_y <- nn_y$nn.dists[, k]
    
    kl_estimate <- mean(log(dist_y / dist_x)) + log(m / (n - 1))
    return(kl_estimate)
  }
  
  # Function to calculate true PSI using numerical integration
  calculate_true_PSI <- function(alpha_p, beta_p, alpha_q, beta_q) {
    integrand <- function(x) {
      p <- dbeta(x, alpha_p, beta_p)
      q <- dbeta(x, alpha_q, beta_q)
      
      I <- (p - q) * log(p / q)
      I <- ifelse(is.nan(I), 0, I)
      I
    }
    
    integrate(integrand, lower = 0, upper = 1)$value
  }
  
  # Reactive values for the simulated data
  simulated_data <- reactive({
    set.seed(rv$seed)
    data_p <- rbeta(input$sample_size_p, input$alpha_p, input$beta_p)
    data_q <- rbeta(input$sample_size_q, input$alpha_q, input$beta_q)
    list(p = data_p, q = data_q)
  }) %>% 
    bindEvent(rv$seed, input$sample_size_p, input$sample_size_q, 
              input$alpha_p, input$beta_p, input$alpha_q, input$beta_q)
  
  # Reactive for plotting data
  plot_data <- reactive({
    data <- simulated_data()
    breaks <- get_breaks(data$p, data$q, input$n_bins, input$binning_method)
    
    # Calculate histogram data for both distributions
    hist_p <- calculate_histogram(data$p, breaks)
    hist_q <- calculate_histogram(data$q, breaks)
    
    # Combine histogram data
    hist_data <- rbind(
      data.frame(hist_p, group = "Reference"),
      data.frame(hist_q, group = "Comparison")
    )
    
    # Generate x values for true density curves
    x <- seq(0, 1, length.out = 2001)
    density_p <- dbeta(x, input$alpha_p, input$beta_p)
    density_q <- dbeta(x, input$alpha_q, input$beta_q)
    
    true_densities <- data.frame(
      x = rep(x, 2),
      density = c(density_p, density_q),
      group = rep(c("Reference", "Comparison"), each = length(x))
    )
    
    list(hist_data = hist_data, 
         densities = true_densities, 
         breaks = breaks)
  })
  
  # Density plot with histograms
  output$density_plot <- renderPlot({
    plot_data <- plot_data()
    
    p <- ggplot() +
      # Add step plots for histograms
      geom_step(data = plot_data$hist_data,
                aes(x = x, y = density, color = group))
    
    # Add true density curves only if checkbox is checked
    if (input$show_true_density) {
      p <- p + geom_line(data = plot_data$densities,
                         aes(x = x, y = density, color = group),
                         linewidth = 1)
    }
    
    p <- p + labs(title = "Density Comparison",
                  subtitle = sprintf("Reference n=%d, Comparison n=%d\nBinning: %s", 
                                     input$sample_size_p, input$sample_size_q,
                                     ifelse(input$binning_method == "equal_width", 
                                            "Equal Width", "Quantile Based")),
                  x = "Value",
                  y = "Density") +
      theme_minimal() +
      scale_color_manual(values = c("Reference" = "blue", "Comparison" = "red"))
    
    # Add rug plot for quantile-based binning
    if (input$binning_method == "quantile") {
      p <- p + geom_rug(data = data.frame(x = plot_data$breaks),
                        aes(x = x),
                        color = "black",
                        alpha = 0.5)
    }
    
    p
  })
  
  # Calculate and display PSI values
  output$psi_results <- renderText({
    data <- simulated_data()
    
    # Calculate different PSI estimates
    hist_psi_result <- calculate_PSI(data$p, data$q, input$n_bins, input$binning_method)
    hist_psi <- hist_psi_result$PSI
    
    # Calculate KNN-based PSI
    D_KL_pq <- kl_divergence_knn(data$p, data$q)
    D_KL_qp <- kl_divergence_knn(data$q, data$p)
    knn_psi <- D_KL_pq + D_KL_qp
    
    # Calculate true PSI
    true_psi <- calculate_true_PSI(input$alpha_p, input$beta_p, 
                                   input$alpha_q, input$beta_q)
    
    paste0("Histogram-based PSI (", input$n_bins, " bins, ", 
           ifelse(input$binning_method == "equal_width", "equal width", "quantile-based"),
           "): ", round(hist_psi, 4),
           "\nk-NN based PSI: ", round(knn_psi, 4),
           "\nTrue PSI: ", round(true_psi, 4))
  })
}

# Run the app
shinyApp(ui = ui, server = server)