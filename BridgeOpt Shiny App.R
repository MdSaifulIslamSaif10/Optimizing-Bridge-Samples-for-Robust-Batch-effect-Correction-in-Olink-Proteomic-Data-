library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(MASS)
library(tidyr)
library(OlinkAnalyze)
library(stringr)
library(dplyr)
library(broom)
library(Matrix)
# Define simulation functions
create_cov_matrix_from_distribution <- function(num_assays,
                                                var_range = c(0.1, 4.0),
                                                rho_range = c(0.1, 0.8)) {
  variances <- runif(num_assays, min = var_range[1], max = var_range[2])
  sigma_vec <- sqrt(variances)
  D <- diag(sigma_vec)
  
  rho_mat <- matrix(0, num_assays, num_assays)
  for (i in 1:num_assays) {
    for (j in i:num_assays) {
      if (i == j) {
        rho_mat[i, j] <- 1
      } else {
        rho_val <- runif(1, min = rho_range[1], max = rho_range[2])
        rho_mat[i, j] <- rho_val
        rho_mat[j, i] <- rho_val
      }
    }
  }
  
  cov_matrix <- D %*% rho_mat %*% D
  
  if (!all(eigen(cov_matrix, only.values = TRUE)$values > 0)) {
    rho_mat <- rho_mat * 0.9
    diag(rho_mat) <- 1
    cov_matrix <- D %*% rho_mat %*% D
    if (!all(eigen(cov_matrix, only.values = TRUE)$values > 0)) {
      cov_matrix <- as.matrix(nearPD(cov_matrix)$mat)
    }
  }
  
  return(cov_matrix)
}
generate_batch_data <- function(num_samples, num_assays, effect_size, similar_assays,
                                case_proportion = 0.5, var_range = c(0.1, 4.0),
                                rho_range = c(0.1, 0.6), mean_range = c(2, 20), batch_id = NULL, batch_type = NULL, study_design) {
  if (study_design == "Study Design II: cases and controls collected in separate batches") {
    prefix <- ifelse(batch_type == "cases", "A", "C")
    sample_ids <- paste0(prefix, batch_id, "_", 1:num_samples)
    num_cases <- if (batch_type == "cases") num_samples else 0
    num_controls <- num_samples - num_cases
  } else {
    num_cases <- round(num_samples * case_proportion)
    num_controls <- num_samples - num_cases
    if (num_cases < 1 || num_controls < 1) stop("Invalid case/control counts.")
    case_ids <- paste0("A", batch_id, "_C", 1:num_cases)
    control_ids <- paste0("A", batch_id, "_N", 1:num_controls)
    sample_ids <- c(case_ids, control_ids)
  }
  
  batch_data <- data.frame(
    SampleID = rep(sample_ids, each = num_assays),
    OlinkID = rep(paste0("OID", sprintf("%05d", 1:num_assays)), num_samples),
    UniProt = rep(paste0("P", sprintf("%04d", 1:num_assays)), num_samples),
    Assay = rep(paste0("Assay_", 1:num_assays), num_samples),
    MissingFreq = runif(num_samples * num_assays, 0.001, 0.002),
    QC_Warning = sample(c("Pass", "Warning"), num_samples * num_assays, replace = TRUE, prob = c(1, 0)),
    LOD = 1
  )
  
  different_assays <- num_assays - similar_assays
  npx_disease_means_similar <- seq(mean_range[1], mean_range[2], length.out = similar_assays)
  npx_disease_means_different <- seq(mean_range[1], mean_range[2], length.out = different_assays)
  npx_disease_means <- c(npx_disease_means_similar, npx_disease_means_different)
  npx_control_means <- c(npx_disease_means_similar, npx_disease_means_different + effect_size)
  
  cov_matrix <- create_cov_matrix_from_distribution(num_assays = num_assays,
                                                    var_range = var_range,
                                                    rho_range = rho_range)
  
  npx_data <- rbind(
    if (num_cases > 0) mvrnorm(n = num_cases, mu = npx_disease_means, Sigma = cov_matrix),
    if (num_controls > 0) mvrnorm(n = num_controls, mu = npx_control_means, Sigma = cov_matrix)
  )
  
  npx_long <- as.data.frame(npx_data) %>%
    setNames(paste0("Assay_", 1:num_assays)) %>%
    mutate(SampleID = sample_ids) %>%
    tidyr::pivot_longer(cols = starts_with("Assay_"), names_to = "Assay", values_to = "NPX")
  
  batch_data <- batch_data %>%
    dplyr::left_join(npx_long, by = c("SampleID", "Assay")) %>%
    mutate(
      Panel_Version = "v.1201",
      PlateID = paste0("Example_Data_", batch_id, "_CAM.csv"),
      Subject = rep(paste0("ID", sample(1:num_samples, num_samples, replace = TRUE)), each = num_assays),
      Treatment = rep(c(rep("Treated", num_cases), rep("Untreated", num_controls)), each = num_assays),
      Site = rep(sample(c("Site_A", "Site_B", "Site_C", "Site_D"), num_samples, replace = TRUE), each = num_assays),
      Time = rep(sample(c("Baseline", "Week.6", "Week.12"), num_samples, replace = TRUE), each = num_assays),
      Project = paste0("data", batch_id),
      Panel = "Olink Cardiometabolic"
    )
  
  return(batch_data)
}
select_bridge_samples <- function(batch_data, bridge_size, sample_missing_freq = 0.1) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  bridge_samples <- batch_data %>%
    dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
    olink_bridgeselector(sampleMissingFreq = sample_missing_freq, n = bridge_size) %>%
    dplyr::select(SampleID)
  
  return(bridge_samples)
}
add_batch_effects <- function(batch_data, num_assays, bridge_samples = NULL,
                              mean_additive = 1.5, sd_additive = 0.75,
                              mean_multiplicative = 1.2, sd_multiplicative = 0.2,
                              noise_mean = 0, noise_sd = 0.5) {
  original_samples <- if (!is.null(bridge_samples)) {
    batch_data %>% filter(!(SampleID %in% bridge_samples$SampleID))
  } else {
    batch_data
  }
  
  ref_stats <- original_samples %>%
    group_by(Assay) %>%
    dplyr::summarise(
      estimated_mean = mean(NPX, na.rm = TRUE),
      .groups = "drop"
    )
  
  batch_effects <- ref_stats %>%
    mutate(
      additive_effects = rnorm(n(), mean = mean_additive, sd = sd_additive),
      multiplicative_effects = rnorm(n(), mean = mean_multiplicative, sd = sd_multiplicative)
    )
  
  batch_data_updated <- batch_data %>%
    dplyr::left_join(batch_effects, by = "Assay") %>%
    mutate(
      NPX = (NPX - estimated_mean) * multiplicative_effects + additive_effects + estimated_mean,
      Project = "data2"
    ) %>%
    dplyr::select(-additive_effects, -multiplicative_effects, -estimated_mean)
  
  batch_data_updated <- batch_data_updated %>%
    group_by(Assay) %>%
    mutate(
      NPX = NPX + rnorm(n(), mean = noise_mean, sd = noise_sd)
    ) %>%
    ungroup()
  
  return(batch_data_updated)
}
combine_batches <- function(batch1_data, batch2_data, remove_bridge = FALSE, overlap_samples = NULL) {
  combined_data <- dplyr::bind_rows(batch1_data, batch2_data)
  
  if (remove_bridge && !is.null(overlap_samples)) {
    combined_data <- combined_data %>%
      dplyr::filter(!(SampleID %in% overlap_samples & Project == "data2"))
  }
  
  return(combined_data)
}
perform_statistical_analysis <- function(data, similar_assays, alpha = 0.05) {
  t_test_results <- data %>%
    dplyr::group_by(Assay) %>%
    do(tidy(t.test(NPX ~ Treatment, data = .))) %>%
    dplyr::ungroup() %>%
    mutate(
      AssayNumber = as.numeric(str_extract(Assay, "\\d+")),
      truth = ifelse(AssayNumber <= similar_assays, "null", "alternative"),
      significant = p.value < alpha
    )
  
  power <- t_test_results %>%
    dplyr::filter(truth == "alternative") %>%
    summarise(power = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(power)
  
  typeI_error <- t_test_results %>%
    dplyr::filter(truth == "null") %>%
    summarise(rate = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(rate)
  
  fdr <- t_test_results %>%
    dplyr::filter(significant) %>%
    summarise(fdr = ifelse(n() == 0, 0, mean(truth == "null", na.rm = TRUE))) %>%
    dplyr::pull(fdr)
  
  return(list(power = power, typeI_error = typeI_error, fdr = fdr))
}
normalize_data <- function(batch1_data, batch2_data, overlap_samples, project_names = c("data1", "data2")) {
  overlap_samples_list <- list("DF1" = overlap_samples$SampleID, "DF2" = overlap_samples$SampleID)
  
  npx_br_data <- olink_normalization_bridge(
    project_1_df = batch1_data %>% as_tibble(),
    project_2_df = batch2_data %>% as_tibble(),
    bridge_samples = overlap_samples_list,
    project_1_name = project_names[1],
    project_2_name = project_names[2],
    project_ref_name = project_names[1]
  )
  
  npx_br_data_no_overlap <- npx_br_data %>%
    dplyr::filter(!(SampleID %in% overlap_samples$SampleID & Project == project_names[2]))
  
  return(npx_br_data_no_overlap)
}
# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "Batch Effect Simulation App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Simulation", tabName = "simulation", icon = icon("calculator")),
      menuItem("User Guide", tabName = "user_guide", icon = icon("info-circle"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "simulation",
              fluidRow(
                box(
                  title = "Simulation Parameters",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 4,
                  selectInput("study_design", "Study Design",
                              choices = c("Study Design I: cases and controls mixed within batches",
                                          "Study Design II: cases and controls collected in separate batches"),
                              selected = "Study Design I: cases and controls mixed within batches"),
                  conditionalPanel(
                    condition = "input.study_design == 'Study Design II: cases and controls collected in separate batches'",
                    numericInput("num_samples_batch1", "Number of Samples (Batch 1 - Cases)", value = 100, min = 10, max = 200),
                    numericInput("num_samples_batch2", "Number of Samples (Batch 2 - Controls)", value = 100, min = 10, max = 200)
                  ),
                  conditionalPanel(
                    condition = "input.study_design == 'Study Design I: cases and controls mixed within batches'",
                    numericInput("num_samples_batch1", "Number of Samples (Batch 1)", value = 100, min = 10, max = 200),
                    numericInput("num_samples_batch2", "Number of Samples (Batch 2)", value = 100, min = 10, max = 200),
                    numericInput("case_proportion_batch1", "Case Proportion (Batch 1)", value = 0.7, min = 0, max = 1, step = 0.1),
                    numericInput("case_proportion_batch2", "Case Proportion (Batch 2)", value = 0.3, min = 0, max = 1, step = 0.1)
                  ),
                  numericInput("num_assays", "Number of Assays", value = 100, min = 50, max = 200),
                  numericInput("similar_assays", "Number of Similar Assays", value = 50, min = 10, max = 100),
                  numericInput("n_iterations", "Number of Iterations", value = 40, min = 1, max = 500),
                  textInput("effect_sizes", "Effect Sizes (comma-separated)", value = "0.1,0.3,0.5,0.7,0.9"),
                  textInput("bridge_sizes", "Bridge Sizes (comma-separated)", value = "2,4,6,8,10,12,14,16,18,20,22,24,26,28,30"),
                  numericInput("mean_minimum", "NPX Mean Range Minimum", value = 2, min = 1, max = 10),
                  numericInput("mean_maximum", "NPX Mean Range Maximum", value = 20, min = 10, max = 50),
                  numericInput("var_minimum", "Data Generation Parameter: Variance Minimum", value = 0.1, min = 0.01, max = 1.0),
                  numericInput("var_maximum", "Data Generation Parameter: Variance Maximum", value = 4.0, min = 1.0, max = 10.0),
                  numericInput("rho_minimum", "Data Generation Parameter: Correlation Minimum", value = 0.1, min = 0.0, max = 0.5),
                  numericInput("rho_maximum", "Data Generation Parameter: Correlation Maximum", value = 0.6, min = 0.5, max = 1.0),
                  numericInput("alpha", "Significance Level (alpha)", value = 0.05, min = 0.01, max = 0.1),
                  selectInput("batch_strength", "Batch Effect Strength", choices = c("Weak", "Standard", "Strong"), selected = "Standard"),
                  numericInput("mean_additive", "Batch Effect Creation Parameter: Additive Effect Mean", value = 0.5, min = -5.0, max = 5.0),
                  numericInput("sd_additive", "Batch Effect Creation Parameter: Additive Effect SD", value = 1.2, min = 0.1, max = 2.0),
                  numericInput("mean_multiplicative", "Batch Effect Creation Parameter: Multiplicative Effect Mean", value = 1.2, min = 0.5, max = 2.0),
                  numericInput("sd_multiplicative", "Batch Effect Creation Parameter: Multiplicative Effect SD", value = 0.2, min = 0.05, max = 0.5),
                  numericInput("noise_mean", "Batch Effect Creation Parameter: Noise Mean", value = 0, min = -1.0, max = 1.0),
                  numericInput("noise_sd", "Batch Effect Creation Parameter: Noise SD", value = 0.5, min = 0.1, max = 1.0),
                  actionButton("run_simulation", "Run Simulation", class = "btn-primary")
                ),
                box(
                  title = "Results",
                  status = "success",
                  solidHeader = TRUE,
                  width = 8,
                  tabsetPanel(
                    tabPanel("Detailed Results",
                             DTOutput("results_table")
                    ),
                    tabPanel("Optimal Bridge Sizes",
                             DTOutput("optimal_table")
                    ),
                    tabPanel("Power Plot",
                             plotOutput("power_plot")
                    ),
                    tabPanel("Type I Error Plot",
                             plotOutput("typeI_error_plot")
                    ),
                    tabPanel("FDR Plot",
                             plotOutput("fdr_plot")
                    )
                  )
                )
              )
      ),
      tabItem(tabName = "user_guide",
              fluidRow(
                box(
                  title = "User Guide",
                  status = "info",
                  solidHeader = TRUE,
                  width = 12,
                  p("Welcome to the BridgeOpt Shiny application! This easy-to-use tool is designed to help you optimize the selection of bridge samples for correcting batch effects in Olink proteomic studies. BridgeOpt provides a hands-on way to simulate different scenarios and find the optimal number of bridge samples in Olink proteomics data. Let’s walk through how to use it step by step!"),
                  h3("Overview"),
                  p("The BridgeOpt Shiny app offers three main features accessible from the left sidebar menu. The Simulation Panel lets users configure experiment parameters, run simulations, and view results. The user guide panel provides background information and usage tips, while the Video Demonstration offers a visual walkthrough of the app."),
                  h3("Simulation Panel"),
                  p("The Simulation Panel is the main workspace in this Shinny application. It is where users set up, run, and visualize the results of their simulations. It’s broken down into a few key areas to help users customize and run their simulations, which are described below."),
                  h4("1. Study Design Selection"),
                  tags$ul(
                    tags$li(strong("Design I: Mixed Batches"), " In this setup, both cases (e.g., individuals with a disease or condition) and controls (e.g., healthy individuals) are included within each batch. The user can choose the proportion of cases and controls in each batch, for example, 70% cases in one batch and 30% in another."),
                    tags$li(strong("Design II: Separate Batches:"), " In this design, one batch contains only cases, and the other only the control group. The user can pick the design that matches their study setup, and the app will adjust the options accordingly.")
                  ),
                  h4("2. Parameter Specification"),
                  p("This section lets users fine-tune their simulation to reflect their real-world experiment. They can adjust several settings to see how they impact optimal bridge sample selection. Here’s what each option does:"),
                  tags$ul(
                    tags$li("Number of Samples per Batch: Users can specify how many samples they want to simulate in each batch. For example, setting both to 100 means it will simulate 100 samples in Batch 1 and 100 in Batch 2."),
                    tags$li("Total Number of Assays: This is the number of proteins or assays users want to simulate. For instance, 92 assays mimic a typical Olink Target 96 panel."),
                    tags$li("Number of Non-differential Assays: This parameter specifies how many assays are assumed to show no true difference between cases and controls. For example, setting this value to 50 in a total of 100 assays means that half of the assays are truly differential and half are truly non-differential. These assays serve as a reference for evaluating the simulation performances using a two-sample t-test."),
                    tags$li("Effect Sizes: Users can specify expected biological differences between cases and controls as a list of values (e.g., 0.1, 0.3, 0.5). Each value represents the magnitude of a true effect. Testing several effect sizes allows users to observe how performance changes with varying signal strength."),
                    tags$li("Bridge Sizes: Users can input different bridge sizes (e.g., 2, 4, 6, 8, 10) to find the optimal bridge sample size. The app tests each value to identify the optimal number of bridge samples. It is important to note that the app will not test any value beyond the user input bridge sizes. We recommend testing at least up to 30 bridge sizes."),
                    tags$li("Batch Effect Strength: This option controls the strength of batch effects between batches. Users can choose among “Weak,” “Standard,” or “Strong” settings, which represent small, moderate, and large batch-related shifts, respectively. Each level corresponds to predefined parameter values for additive effects, multiplicative effects, and random noise. These parameters can also be manually adjusted for greater flexibility. The “Standard” setting is the default and recommended option, as it is designed to be more conservative than the batch effects typically encountered in real-world datasets."),
                    tags$li("Significance Level (α): Users define the statistical significance threshold (typically 0.05). This determines how strictly the simulation evaluates whether observed differences are statistically meaningful.")
                  ),
                  h4("3. Additional Optional Parameters"),
                  tags$ul(
                    tags$li("NPX Mean Range: Users can define the average protein expression levels used during data generation in the simulation. By default, this range is set between 2 and 20 NPX, which reflects the mean values typically observed in real Olink datasets."),
                    tags$li("Variance Range: User can input the range of variability in the simulated data (default: 0.1 to 4.0), allowing the model to represent different degrees of dispersion across proteins or assays."),
                    tags$li("Correlation Range: User can determine the minimum and maximum correlation values between assays (default: 0.1 to 0.6), capturing how closely related the protein measurements are. These variance and correlation settings are used to construct the covariance matrix for generating data from a multivariate normal distribution.")
                  ),
                  h3("Running the Simulation"),
                  p("After defining all parameters, users can initiate the simulation by clicking the “Run Simulation” button. The progress bar provides real-time feedback, showing the current iteration and which combination of effect size and bridge size is being processed. The app runs multiple iterations (based on the number of iterations) for each effect size and bridge size combination. This repetition is required for stable averages of power, type I error, and FDR. We recommend specifying at least 30 iterations to run simulations."),
                  h3("Results Display"),
                  p("Once the simulation is complete, the application automatically presents the outcomes across several result tabs:"),
                  tags$ol(
                    tags$li("Detailed Results: A comprehensive summary table showing performance metrics, including statistical power, Type I error, and false discovery rate (FDR) for every tested combination of effect size and bridge size, both before and after normalization in each iteration."),
                    tags$li("Optimal Bridge Sizes: This table provides the optimal bridge size recommendations for each effect size based on simulated power and Type-I-error."),
                    tags$li("Visualization of Simulation Results: The simulation results are presented under three scenarios: (1) gold standard, (2) uncorrected data, and (3) normalized data after applying Olink normalization correction. These scenarios are compared across three plots displayed in the app: the power plot, Type I Error plot, and FDR plot. The Power plot shows the ability to identify real effects across varying effect sizes and bridge sample numbers, where higher power values indicate stronger detection of truly differential proteins. The Type I Error plot reflects the rate of false positives, which should remain close to the predefined significance level (e.g., α = 0.05). The FDR Plot illustrates the proportion of false discoveries among the significant results, where lower values indicate higher confidence in the detected effects. In all plots, blue lines represent different bridge sizes, with larger bridge sizes shown in brighter blue. The red dashed line indicates the gold-standard scenario, while the orange dashed line represents the “no correction” scenario.")
                  )
                )
              )
      )
    )
  )
)
# Server Definition
server <- function(input, output, session) {
  simulation_results <- reactiveVal(NULL)
  
  observeEvent(input$batch_strength, {
    if (input$batch_strength == "Weak") {
      updateNumericInput(session, "mean_additive", value = 0.5)
      updateNumericInput(session, "sd_additive", value = 0.5)
      updateNumericInput(session, "mean_multiplicative", value = 1.1)
      updateNumericInput(session, "sd_multiplicative", value = 0.1)
      updateNumericInput(session, "noise_mean", value = 0)
      updateNumericInput(session, "noise_sd", value = 0.3)
    } else if (input$batch_strength == "Standard") {
      updateNumericInput(session, "mean_additive", value = 1.0)
      updateNumericInput(session, "sd_additive", value = 1.2)
      updateNumericInput(session, "mean_multiplicative", value = 1.2)
      updateNumericInput(session, "sd_multiplicative", value = 0.2)
      updateNumericInput(session, "noise_mean", value = 0)
      updateNumericInput(session, "noise_sd", value = 0.5)
    } else if (input$batch_strength == "Strong") {
      updateNumericInput(session, "mean_additive", value = 1.5)
      updateNumericInput(session, "sd_additive", value = 1.75)
      updateNumericInput(session, "mean_multiplicative", value = 1.3)
      updateNumericInput(session, "sd_multiplicative", value = 0.3)
      updateNumericInput(session, "noise_mean", value = 0)
      updateNumericInput(session, "noise_sd", value = 0.7)
    }
  })
  
  observeEvent(input$run_simulation, {
    withProgress(message = "Running simulation...", value = 0, {
      params <- list(
        study_design = input$study_design,
        effect_sizes = as.numeric(unlist(strsplit(input$effect_sizes, ","))),
        bridge_sizes = as.numeric(unlist(strsplit(input$bridge_sizes, ","))),
        n_iterations = input$n_iterations,
        num_samples_batch1 = input$num_samples_batch1,
        num_samples_batch2 = input$num_samples_batch2,
        case_proportion_batch1 = if (input$study_design == "Study Design I: cases and controls mixed within batches") input$case_proportion_batch1 else 1.0,
        case_proportion_batch2 = if (input$study_design == "Study Design I: cases and controls mixed within batches") input$case_proportion_batch2 else 0.0,
        num_assays = input$num_assays,
        similar_assays = input$similar_assays,
        mean_range = c(input$mean_minimum, input$mean_maximum),
        var_range = c(input$var_minimum, input$var_maximum),
        rho_range = c(input$rho_minimum, input$rho_maximum),
        alpha = input$alpha,
        mean_additive = input$mean_additive,
        sd_additive = input$sd_additive,
        mean_multiplicative = input$mean_multiplicative,
        sd_multiplicative = input$sd_multiplicative,
        noise_mean = input$noise_mean,
        noise_sd = input$noise_sd,
        sample_missing_freq = 0.1
      )
      
      if (any(is.na(params$effect_sizes)) || any(is.na(params$bridge_sizes))) {
        showNotification("Please enter valid comma-separated numbers for effect sizes and bridge sizes.", type = "error")
        return()
      }
      if (params$similar_assays >= params$num_assays) {
        showNotification("Number of similar assays must be less than total assays.", type = "error")
        return()
      }
      if (max(params$bridge_sizes) > params$num_samples_batch1) {
        showNotification("Bridge sizes cannot exceed number of samples in batch 1.", type = "error")
        return()
      }
      if (params$study_design == "Study Design I: cases and controls mixed within batches") {
        if (params$case_proportion_batch1 < 0 || params$case_proportion_batch1 > 1 ||
            params$case_proportion_batch2 < 0 || params$case_proportion_batch2 > 1) {
          showNotification("Case proportions must be between 0 and 1.", type = "error")
          return()
        }
      }
      
      total_steps <- length(params$effect_sizes) * length(params$bridge_sizes) * params$n_iterations
      current_step <- 0
      
      final_results_df <- data.frame()
      optimal_bridge_df <- data.frame()
      
      epsilon_power <- 0.01
      epsilon_typeI <- 0.01
      
      for (eff in params$effect_sizes) {
        cat("Effect Size:", eff, "\n")
        
        power_results_list <- list()
        
        for (bridge_size in params$bridge_sizes) {
          cat("Bridge Sample Size:", bridge_size, "\n")
          
          results_list <- list(
            power_Gold_Standard = numeric(),
            power_with_batch_effect = numeric(),
            power_after_batch_effect_removal = numeric(),
            typeI_error_Gold_Standard = numeric(),
            typeI_error_with_batch_effect = numeric(),
            typeI_error_after_batch_effect_removal = numeric(),
            FDR_Gold_Standard = numeric(),
            FDR_with_batch_effect = numeric(),
            FDR_after_batch_effect_removal = numeric()
          )
          
          for (iteration in 1:params$n_iterations) {
            current_step <- current_step + 1
            setProgress(value = current_step / total_steps,
                        detail = paste("Iteration:", iteration, "Effect Size:", eff, "Bridge Size:", bridge_size))
            
            set.seed(123 + iteration)
            
            batch1_data <- generate_batch_data(
              num_samples = params$num_samples_batch1,
              num_assays = params$num_assays,
              effect_size = eff,
              similar_assays = params$similar_assays,
              case_proportion = params$case_proportion_batch1,
              var_range = params$var_range,
              rho_range = params$rho_range,
              mean_range = params$mean_range,
              batch_id = "1",
              batch_type = if (params$study_design == "Study Design II: cases and controls collected in separate batches") "cases" else NULL,
              study_design = params$study_design
            )
            batch2_data <- generate_batch_data(
              num_samples = params$num_samples_batch2,
              num_assays = params$num_assays,
              effect_size = eff,
              similar_assays = params$similar_assays,
              case_proportion = params$case_proportion_batch2,
              var_range = params$var_range,
              rho_range = params$rho_range,
              mean_range = params$mean_range,
              batch_id = "2",
              batch_type = if (params$study_design == "Study Design II: cases and controls collected in separate batches") "controls" else NULL,
              study_design = params$study_design
            )
            
            combined_data_gold <- combine_batches(batch1_data, batch2_data)
            stats_gold <- perform_statistical_analysis(combined_data_gold, params$similar_assays, params$alpha)
            
            bridge_samples_all <- batch1_data %>%
              dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
              mutate(Project = "data2")
            
            batch2_data_with_bridge <- dplyr::bind_rows(batch2_data, bridge_samples_all)
            batch2_data_updated <- add_batch_effects(
              batch_data = batch2_data_with_bridge,
              num_assays = params$num_assays,
              bridge_samples = bridge_samples_all,
              mean_additive = params$mean_additive,
              sd_additive = params$sd_additive,
              mean_multiplicative = params$mean_multiplicative,
              sd_multiplicative = params$sd_multiplicative,
              noise_mean = params$noise_mean,
              noise_sd = params$noise_sd
            )
            
            bridge_samples_subset <- select_bridge_samples(batch1_data, bridge_size, params$sample_missing_freq)
            overlap_samples <- bridge_samples_subset
            
            batch2_data_for_norm <- batch2_data_updated %>%
              dplyr::filter(!(SampleID %in% bridge_samples_all$SampleID) | SampleID %in% overlap_samples$SampleID)
            
            combined_data_batch <- combine_batches(
              batch1_data, batch2_data_updated, remove_bridge = TRUE, overlap_samples = bridge_samples_all$SampleID
            )
            stats_batch <- perform_statistical_analysis(combined_data_batch, params$similar_assays, params$alpha)
            
            normalized_data <- normalize_data(
              batch1_data = batch1_data,
              batch2_data = batch2_data_for_norm,
              overlap_samples = overlap_samples
            )
            stats_normalized <- perform_statistical_analysis(normalized_data, params$similar_assays, params$alpha)
            
            results_list$power_Gold_Standard <- c(results_list$power_Gold_Standard, stats_gold$power)
            results_list$power_with_batch_effect <- c(results_list$power_with_batch_effect, stats_batch$power)
            results_list$power_after_batch_effect_removal <- c(results_list$power_after_batch_effect_removal, stats_normalized$power)
            results_list$typeI_error_Gold_Standard <- c(results_list$typeI_error_Gold_Standard, stats_gold$typeI_error)
            results_list$typeI_error_with_batch_effect <- c(results_list$typeI_error_with_batch_effect, stats_batch$typeI_error)
            results_list$typeI_error_after_batch_effect_removal <- c(results_list$typeI_error_after_batch_effect_removal, stats_normalized$typeI_error)
            results_list$FDR_Gold_Standard <- c(results_list$FDR_Gold_Standard, stats_gold$fdr)
            results_list$FDR_with_batch_effect <- c(results_list$FDR_with_batch_effect, stats_batch$fdr)
            results_list$FDR_after_batch_effect_removal <- c(results_list$FDR_after_batch_effect_removal, stats_normalized$fdr)
          }
          
          power_results_df <- data.frame(
            Bridge_Size = bridge_size,
            Effect_Size = eff,
            power_Gold_Standard = mean(results_list$power_Gold_Standard, na.rm = TRUE),
            power_with_batch_effect = mean(results_list$power_with_batch_effect, na.rm = TRUE),
            power_after_batch_effect_removal = mean(results_list$power_after_batch_effect_removal, na.rm = TRUE),
            typeI_error_Gold_Standard = mean(results_list$typeI_error_Gold_Standard, na.rm = TRUE),
            typeI_error_with_batch_effect = mean(results_list$typeI_error_with_batch_effect, na.rm = TRUE),
            typeI_error_after_batch_effect_removal = mean(results_list$typeI_error_after_batch_effect_removal, na.rm = TRUE),
            FDR_Gold_Standard = mean(results_list$FDR_Gold_Standard, na.rm = TRUE),
            FDR_with_batch_effect = mean(results_list$FDR_with_batch_effect, na.rm = TRUE),
            FDR_after_batch_effect_removal = mean(results_list$FDR_after_batch_effect_removal, na.rm = TRUE)
          )
          
          power_results_list[[as.character(bridge_size)]] <- power_results_df
        }
        
        power_results_df <- do.call(rbind, power_results_list)
        
        power_results_df_tmp <- power_results_df %>%
          mutate(
            power_diff = abs(power_after_batch_effect_removal - power_Gold_Standard),
            typeI_error_diff = abs(typeI_error_after_batch_effect_removal - typeI_error_Gold_Standard)
          )
        
        min_power_diff <- min(power_results_df_tmp$power_diff, na.rm = TRUE)
        min_typeI_error_diff <- min(power_results_df_tmp$typeI_error_diff, na.rm = TRUE)
        
        optimal_row_power <- power_results_df_tmp %>%
          dplyr::filter(power_diff <= min_power_diff + epsilon_power) %>%
          dplyr::slice(which.min(Bridge_Size))
        
        optimal_row_typeI <- power_results_df_tmp %>%
          dplyr::filter(typeI_error_diff <= min_typeI_error_diff + epsilon_typeI) %>%
          dplyr::slice(which.min(Bridge_Size))
        
        optimal_bridge_df <- rbind(
          optimal_bridge_df,
          data.frame(
            Effect_Size = eff,
            Optimal_Bridge_Size_power = optimal_row_power$Bridge_Size,
            Power_After = optimal_row_power$power_after_batch_effect_removal,
            Power_Gold = optimal_row_power$power_Gold_Standard,
            Power_Diff = optimal_row_power$power_diff,
            Min_Power_Diff = min_power_diff,
            Power_Diff_Difference = optimal_row_power$power_diff - min_power_diff,
            Optimal_Bridge_Size_typeI = optimal_row_typeI$Bridge_Size,
            TypeI_Error_After = optimal_row_typeI$typeI_error_after_batch_effect_removal,
            TypeI_Error_Gold = optimal_row_typeI$typeI_error_Gold_Standard,
            TypeI_Error_Diff = optimal_row_typeI$typeI_error_diff,
            Min_TypeI_Error_Diff = min_typeI_error_diff,
            TypeI_Error_Diff_Difference = optimal_row_typeI$typeI_error_diff - min_typeI_error_diff
          )
        )
        
        final_results_df <- rbind(
          final_results_df,
          power_results_df_tmp
        )
      }
      
      results <- list(final_results = final_results_df, optimal_bridge = optimal_bridge_df)
      if (is.null(results)) return()
      
      setProgress(1.0, detail = "Simulation complete!")
      simulation_results(results)
    })
  })
  
  output$results_table <- renderDT({
    req(simulation_results())
    data <- simulation_results()$final_results
    
    numeric_cols <- sapply(data, is.numeric)
    numeric_col_indices <- which(numeric_cols) - 1
    
    datatable(
      data,
      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
      rownames = FALSE
    ) %>% formatRound(columns = numeric_col_indices, digits = 4)
  })
  
  output$optimal_table <- renderDT({
    req(simulation_results())
    data <- simulation_results()$optimal_bridge %>%
      dplyr::select(Effect_Size, Optimal_Bridge_Size_power, Optimal_Bridge_Size_typeI)
    
    numeric_cols <- sapply(data, is.numeric)
    numeric_col_indices <- which(numeric_cols) - 1
    
    datatable(
      data,
      options = list(pageLength = 10, autoWidth = TRUE),
      rownames = FALSE
    ) %>% formatRound(columns = numeric_col_indices, digits = 4)
  })
  
  output$power_plot <- renderPlot({
    req(simulation_results())
    results <- simulation_results()$final_results
    
    gold_standard <- results %>%
      group_by(Effect_Size) %>%
      summarize(power_Gold_Standard = mean(power_Gold_Standard, na.rm = TRUE),
                .groups = "drop")
    
    batch_effect <- results %>%
      group_by(Effect_Size) %>%
      summarize(power_with_batch_effect = mean(power_with_batch_effect, na.rm = TRUE),
                .groups = "drop")
    
    ggplot(results, aes(x = Effect_Size, y = power_after_batch_effect_removal,
                        color = as.numeric(Bridge_Size), group = Bridge_Size)) +
      geom_line(size = 1.5) +
      geom_point(size = 3) +
      scale_color_gradientn(colors = c("#A6CEE3", "#1F78B4", "#08306B"), name = "Bridge Size") +
      geom_line(data = gold_standard, aes(x = Effect_Size, y = power_Gold_Standard,
                                          linetype = "Gold Standard"),
                color = "red", size = 2, inherit.aes = FALSE) +
      geom_point(data = gold_standard, aes(x = Effect_Size, y = power_Gold_Standard),
                 color = "red", size = 4, inherit.aes = FALSE) +
      geom_line(data = batch_effect, aes(x = Effect_Size, y = power_with_batch_effect,
                                         linetype = "No Batch Correction"),
                color = "orange", size = 2, inherit.aes = FALSE) +
      geom_point(data = batch_effect, aes(x = Effect_Size, y = power_with_batch_effect),
                 color = "orange", size = 4, inherit.aes = FALSE) +
      labs(title = "Average Power vs Effect Size",
           x = "Effect Size",
           y = "Power",
           color = "Bridge Size",
           linetype = NULL) +
      guides(color = guide_colorbar(order = 1),
             linetype = guide_legend(title = "Reference Lines",
                                     override.aes = list(color = c("red", "orange")),
                                     order = 2)) +
      theme_classic(base_size = 20) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 35),
        axis.title.x = element_text(face = "bold", size = 35),
        axis.title.y = element_text(face = "bold", size = 35),
        axis.text.x = element_text(face = "bold", size = 35),
        axis.text.y = element_text(face = "bold", size = 35),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(face = "bold", size = 25)
      ) +
      scale_y_continuous(limits = c(0, 1))
  })
  
  output$typeI_error_plot <- renderPlot({
    req(simulation_results())
    results <- simulation_results()$final_results
    
    gold_standard <- results %>%
      group_by(Effect_Size) %>%
      summarize(typeI_error_Gold_Standard = mean(typeI_error_Gold_Standard, na.rm = TRUE),
                .groups = "drop")
    
    batch_effect <- results %>%
      group_by(Effect_Size) %>%
      summarize(typeI_error_with_batch_effect = mean(typeI_error_with_batch_effect, na.rm = TRUE),
                .groups = "drop")
    
    ggplot(results, aes(x = Effect_Size, y = typeI_error_after_batch_effect_removal,
                        color = as.numeric(Bridge_Size), group = Bridge_Size)) +
      geom_line(size = 1.5) +
      geom_point(size = 3) +
      scale_color_gradientn(colors = c("#A6CEE3", "#1F78B4", "#08306B"), name = "Bridge Size") +
      geom_line(data = gold_standard, aes(x = Effect_Size, y = typeI_error_Gold_Standard,
                                          linetype = "Gold Standard"),
                color = "red", size = 2, inherit.aes = FALSE) +
      geom_point(data = gold_standard, aes(x = Effect_Size, y = typeI_error_Gold_Standard),
                 color = "red", size = 4, inherit.aes = FALSE) +
      geom_line(data = batch_effect, aes(x = Effect_Size, y = typeI_error_with_batch_effect,
                                         linetype = "No Batch Correction"),
                color = "orange", size = 2, inherit.aes = FALSE) +
      geom_point(data = batch_effect, aes(x = Effect_Size, y = typeI_error_with_batch_effect),
                 color = "orange", size = 4, inherit.aes = FALSE) +
      labs(title = "Average Type I Error vs Effect Size",
           x = "Effect Size",
           y = "Type I Error",
           color = "Bridge Size",
           linetype = NULL) +
      guides(color = guide_colorbar(order = 1),
             linetype = guide_legend(title = "Reference Lines",
                                     override.aes = list(color = c("red", "orange")),
                                     order = 2)) +
      theme_classic(base_size = 20) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 35),
        axis.title.x = element_text(face = "bold", size = 35),
        axis.title.y = element_text(face = "bold", size = 35),
        axis.text.x = element_text(face = "bold", size = 35),
        axis.text.y = element_text(face = "bold", size = 35),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(face = "bold", size = 25)
      ) +
      scale_y_continuous(limits = c(0, 1))
  })
  
  output$fdr_plot <- renderPlot({
    req(simulation_results())
    results <- simulation_results()$final_results
    
    gold_standard <- results %>%
      group_by(Effect_Size) %>%
      summarize(FDR_Gold_Standard = mean(FDR_Gold_Standard, na.rm = TRUE),
                .groups = "drop")
    
    batch_effect <- results %>%
      group_by(Effect_Size) %>%
      summarize(FDR_with_batch_effect = mean(FDR_with_batch_effect, na.rm = TRUE),
                .groups = "drop")
    
    ggplot(results, aes(x = Effect_Size, y = FDR_after_batch_effect_removal,
                        color = as.numeric(Bridge_Size), group = Bridge_Size)) +
      geom_line(size = 1.5) +
      geom_point(size = 3) +
      scale_color_gradientn(colors = c("#A6CEE3", "#1F78B4", "#08306B"), name = "Bridge Size") +
      geom_line(data = gold_standard, aes(x = Effect_Size, y = FDR_Gold_Standard,
                                          linetype = "Gold Standard"),
                color = "red", size = 2, inherit.aes = FALSE) +
      geom_point(data = gold_standard, aes(x = Effect_Size, y = FDR_Gold_Standard),
                 color = "red", size = 4, inherit.aes = FALSE) +
      geom_line(data = batch_effect, aes(x = Effect_Size, y = FDR_with_batch_effect,
                                         linetype = "No Batch Correction"),
                color = "orange", size = 2, inherit.aes = FALSE) +
      geom_point(data = batch_effect, aes(x = Effect_Size, y = FDR_with_batch_effect),
                 color = "orange", size = 4, inherit.aes = FALSE) +
      labs(title = "Average FDR vs Effect Size",
           x = "Effect Size",
           y = "FDR",
           color = "Bridge Size",
           linetype = NULL) +
      guides(color = guide_colorbar(order = 1),
             linetype = guide_legend(title = "Reference Lines",
                                     override.aes = list(color = c("red", "orange")),
                                     order = 2)) +
      theme_classic(base_size = 20) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 35),
        axis.title.x = element_text(face = "bold", size = 35),
        axis.title.y = element_text(face = "bold", size = 35),
        axis.text.x = element_text(face = "bold", size = 35),
        axis.text.y = element_text(face = "bold", size = 35),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(face = "bold", size = 25)
      ) +
      scale_y_continuous(limits = c(0, 0.6))
  })
}
# Run the Shiny App
shinyApp(ui = ui, server = server)