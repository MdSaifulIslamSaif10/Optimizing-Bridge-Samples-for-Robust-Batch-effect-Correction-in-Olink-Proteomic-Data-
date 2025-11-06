# Load required libraries
library(MASS)
library(tidyr)
library(OlinkAnalyze)
library(stringr)
library(dplyr)
library(broom)
library(truncnorm)
library(Matrix)
library(ggplot2)

# Function to create covariance matrix from distribution
create_cov_matrix_from_distribution <- function(num_assays,
                                                var_range = c(0.1, 4.0),
                                                rho_range = c(0.1, 0.7)) {
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

# Function 1: Generate Batch Data
generate_batch_data <- function(num_samples, num_assays, effect_size, similar_assays,
                                batch_type = "cases", batch_id = "1", noise_sd = 0.5) {
  # Generate sample IDs
  prefix <- ifelse(batch_type == "cases", "A", "C")
  sample_ids <- paste0(prefix, batch_id, "_", 1:num_samples)
  
  # Create metadata
  batch_data <- data.frame(
    SampleID = rep(sample_ids, each = num_assays),
    OlinkID = rep(paste0("OID", sprintf("%05d", 1:num_assays)), num_samples),
    UniProt = rep(paste0("P", sprintf("%04d", 1:num_assays)), num_samples),
    Assay = rep(paste0("Assay_", 1:num_assays), num_samples),
    MissingFreq = runif(num_samples * num_assays, 0.001, 0.002),
    QC_Warning = sample(c("Pass", "Warning"), num_samples * num_assays, replace = TRUE, prob = c(1, 0)),
    LOD = 1
  )
  
  # Define means for cases and controls (for NPX)
  different_assays <- num_assays - similar_assays
  npx_disease_means_similar <- seq(2, 20, length.out = similar_assays)
  npx_disease_means_different <- seq(2, 20, length.out = different_assays)
  npx_disease_means <- c(npx_disease_means_similar, npx_disease_means_different)
  npx_means <- if (batch_type == "controls") {
    c(npx_disease_means_similar, npx_disease_means_different + effect_size)
  } else {
    npx_disease_means
  }
  
  # Covariance matrix
  cov_matrix <- create_cov_matrix_from_distribution(num_assays = num_assays,
                                                    var_range = c(0.1, 4.0),
                                                    rho_range = c(0.1, 0.7))
  
  # Generate NPX data
  npx_data <- mvrnorm(n = num_samples, mu = npx_means, Sigma = cov_matrix)
  
  # Convert to long format
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
      Treatment = ifelse(batch_type == "cases", "Treated", "Untreated"),
      Site = rep(sample(c("Site_A", "Site_B", "Site_C", "Site_D"), num_samples, replace = TRUE), each = num_assays),
      Time = rep(sample(c("Baseline", "Week.6", "Week.12"), num_samples, replace = TRUE), each = num_assays),
      Project = paste0("data", batch_id),
      Panel = "Olink Cardiometabolic"
    )
  
  return(batch_data)
}

# Function 2: Select Bridge Samples
select_bridge_samples <- function(batch_data, bridge_size, sample_missing_freq = 0.1) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  bridge_samples <- batch_data %>%
    dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
    olink_bridgeselector(sampleMissingFreq = sample_missing_freq, n = bridge_size) %>%
    dplyr::select(SampleID)
  
  return(bridge_samples)
}

# Function 3: Add Batch Effects
add_batch_effects <- function(batch_data, num_assays, bridge_samples = NULL,
                              mean_additive = 0.5, sd_additive = 0.4,
                              mean_multiplicative = 1.1, sd_multiplicative = 0.1,
                              noise_mean = 0, noise_sd = 0.4) {
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
      multiplicative_effects = rnorm(
        n(), mean = mean_multiplicative, sd = sd_multiplicative
      )
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

# Function 4: Combine Batches
combine_batches <- function(batch1_data, batch2_data, remove_bridge = FALSE, overlap_samples = NULL) {
  combined_data <- dplyr::bind_rows(batch1_data, batch2_data)
  
  if (remove_bridge && !is.null(overlap_samples)) {
    combined_data <- combined_data %>%
      dplyr::filter(!(SampleID %in% overlap_samples & Project == "data2"))
  }
  
  return(combined_data)
}

# Function 5: Perform Statistical Analysis
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
  
  power_binary <- t_test_results %>%
    dplyr::filter(truth == "alternative") %>%
    summarise(power = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(power)
  
  typeI_error_binary <- t_test_results %>%
    dplyr::filter(truth == "null") %>%
    summarise(rate = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(rate)
  
  fdr_binary <- t_test_results %>%
    dplyr::filter(significant) %>%
    summarise(fdr = ifelse(n() == 0, 0, mean(truth == "null", na.rm = TRUE))) %>%
    dplyr::pull(fdr)
  
  return(list(
    power_binary = power_binary,
    typeI_error_binary = typeI_error_binary,
    fdr_binary = fdr_binary
  ))
}

# Function 6: Compute DSC
compute_dsc <- function(data) {
  data_wide <- data %>%
    dplyr::select(SampleID, Assay, NPX, Project) %>%
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  
  data_matrix <- as.matrix(data_wide %>% dplyr::select_if(is.numeric))
  rownames(data_matrix) <- data_wide$SampleID
  data_matrix_t <- t(data_matrix)
  batch_labels <- data_wide$Project
  names(batch_labels) <- data_wide$SampleID
  
  source("https://raw.githubusercontent.com/vittoriofortino84/COPS/master/R/evaluation_metrics.R")
  dsc <- DSC(data_matrix_t, batch_labels)
  
  return(dsc)
}

# Function 7: Normalize Data Using Bridge Samples
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

# Main Function for the Simulation with Optimization
run_simulation <- function(params) {
  final_results_df <- data.frame()
  optimal_bridge_df <- data.frame()
  
  epsilon_power <- 0.01
  epsilon_typeI <- 0.01
  
  for (eff in params$effect_sizes) {
    cat("Effect Size:", eff, "\n")
    
    for (iteration in 1:params$n_iterations) {
      cat("Iteration:", iteration, "\n")
      
      # Set seed for reproducibility within iteration
      set.seed(123 + iteration)
      
      # Generate datasets: cases in batch1, controls in batch2
      batch1_data <- generate_batch_data(
        num_samples = params$num_samples_batch1,
        num_assays = params$num_assays,
        effect_size = eff,
        similar_assays = params$similar_assays,
        batch_type = "cases",
        batch_id = "1",
        noise_sd = params$noise_sd
      )
      batch2_data <- generate_batch_data(
        num_samples = params$num_samples_batch2,
        num_assays = params$num_assays,
        effect_size = eff,
        similar_assays = params$similar_assays,
        batch_type = "controls",
        batch_id = "2",
        noise_sd = params$noise_sd
      )
      
      # CHANGE 1: Use all non-control batch1_data samples as bridge samples for batch effect creation
      bridge_samples_all <- batch1_data %>%
        dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
        mutate(Project = "data2")
      
      # CHANGE 2: Combine batch2_data with all batch1_data samples and apply batch effect
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
      
      # Compute gold standard once per dataset
      combined_data_gold <- combine_batches(batch1_data, batch2_data)
      stats_gold <- perform_statistical_analysis(combined_data_gold, params$similar_assays)
      dsc_gold <- compute_dsc(combined_data_gold)
      
      power_results_list <- list()
      
      for (bridge_size in params$bridge_sizes) {
        cat("Bridge Sample Size:", bridge_size, "\n")
        
        results_list <- list(
          power_Gold_Standard_binary = numeric(),
          power_with_batch_effect_binary = numeric(),
          power_after_batch_effect_removal_binary = numeric(),
          typeI_error_Gold_Standard_binary = numeric(),
          typeI_error_with_batch_effect_binary = numeric(),
          typeI_error_after_batch_effect_removal_binary = numeric(),
          FDR_Gold_Standard_binary = numeric(),
          FDR_with_batch_effect_binary = numeric(),
          FDR_after_batch_effect_removal_binary = numeric(),
          DSC_Gold_Standard = numeric(),
          DSC_with_batch_effect = numeric(),
          DSC_after_batch_effect_removal = numeric()
        )
        
        # CHANGE 3: Select bridge samples from original batch1_data
        bridge_samples_subset <- select_bridge_samples(batch1_data, bridge_size)
        overlap_samples <- bridge_samples_subset
        
        # CHANGE 4: Create batch2 data for normalization (excluding batch1 samples except selected bridge samples)
        batch2_data_for_norm <- batch2_data_updated %>%
          dplyr::filter(!(SampleID %in% bridge_samples_all$SampleID) | SampleID %in% overlap_samples$SampleID)
        
        # CHANGE 5: Combine batches, removing all batch1 samples for analysis
        combined_data_batch <- combine_batches(
          batch1_data, batch2_data_updated, remove_bridge = TRUE, overlap_samples = bridge_samples_all$SampleID
        )
        stats_batch <- perform_statistical_analysis(combined_data_batch, params$similar_assays)
        dsc_batch <- compute_dsc(combined_data_batch)
        
        # CHANGE 6: Normalize data using selected bridge samples
        normalized_data <- normalize_data(
          batch1_data = batch1_data,
          batch2_data = batch2_data_for_norm,
          overlap_samples = overlap_samples
        )
        stats_normalized <- perform_statistical_analysis(normalized_data, params$similar_assays)
        dsc_normalized <- compute_dsc(normalized_data)
        
        # Store results
        results_list$power_Gold_Standard_binary <- stats_gold$power_binary
        results_list$power_with_batch_effect_binary <- stats_batch$power_binary
        results_list$power_after_batch_effect_removal_binary <- stats_normalized$power_binary
        results_list$typeI_error_Gold_Standard_binary <- stats_gold$typeI_error_binary
        results_list$typeI_error_with_batch_effect_binary <- stats_batch$typeI_error_binary
        results_list$typeI_error_after_batch_effect_removal_binary <- stats_normalized$typeI_error_binary
        results_list$FDR_Gold_Standard_binary <- stats_gold$fdr_binary
        results_list$FDR_with_batch_effect_binary <- stats_batch$fdr_binary
        results_list$FDR_after_batch_effect_removal_binary <- stats_normalized$fdr_binary
        results_list$DSC_Gold_Standard <- dsc_gold
        results_list$DSC_with_batch_effect <- dsc_batch
        results_list$DSC_after_batch_effect_removal <- dsc_normalized
        
        # Store results for this bridge size
        power_results_df <- data.frame(
          Bridge_Size = bridge_size,
          Effect_Size = eff,
          Iteration = iteration,
          power_Gold_Standard_binary = mean(results_list$power_Gold_Standard_binary, na.rm = TRUE),
          power_with_batch_effect_binary = mean(results_list$power_with_batch_effect_binary, na.rm = TRUE),
          power_after_batch_effect_removal_binary = mean(results_list$power_after_batch_effect_removal_binary, na.rm = TRUE),
          typeI_error_Gold_Standard_binary = mean(results_list$typeI_error_Gold_Standard_binary, na.rm = TRUE),
          typeI_error_with_batch_effect_binary = mean(results_list$typeI_error_with_batch_effect_binary, na.rm = TRUE),
          typeI_error_after_batch_effect_removal_binary = mean(results_list$typeI_error_after_batch_effect_removal_binary, na.rm = TRUE),
          FDR_Gold_Standard_binary = mean(results_list$FDR_Gold_Standard_binary, na.rm = TRUE),
          FDR_with_batch_effect_binary = mean(results_list$FDR_with_batch_effect_binary, na.rm = TRUE),
          FDR_after_batch_effect_removal_binary = mean(results_list$FDR_after_batch_effect_removal_binary, na.rm = TRUE),
          DSC_Gold_Standard = mean(results_list$DSC_Gold_Standard, na.rm = TRUE),
          DSC_with_batch_effect = mean(results_list$DSC_with_batch_effect, na.rm = TRUE),
          DSC_after_batch_effect_removal = mean(results_list$DSC_after_batch_effect_removal, na.rm = TRUE)
        )
        
        power_results_list[[as.character(bridge_size)]] <- power_results_df
      }
      
      power_results_df <- do.call(rbind, power_results_list)
      
      power_results_df_tmp <- power_results_df %>%
        mutate(
          power_diff = abs(power_after_batch_effect_removal_binary - power_Gold_Standard_binary),
          typeI_error_diff = abs(typeI_error_after_batch_effect_removal_binary - typeI_error_Gold_Standard_binary)
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
          Iteration = iteration,
          Optimal_Bridge_Size_power = optimal_row_power$Bridge_Size,
          Power_After_binary = optimal_row_power$power_after_batch_effect_removal_binary,
          Power_Gold_binary = optimal_row_power$power_Gold_Standard_binary,
          Power_Diff = optimal_row_power$power_diff,
          Min_Power_Diff = min_power_diff,
          Power_Diff_Difference = optimal_row_power$power_diff - min_power_diff,
          Optimal_Bridge_Size_typeI = optimal_row_typeI$Bridge_Size,
          TypeI_Error_After_binary = optimal_row_typeI$typeI_error_after_batch_effect_removal_binary,
          TypeI_Error_Gold_binary = optimal_row_typeI$typeI_error_Gold_Standard_binary,
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
  }
  
  return(list(
    final_results = final_results_df,
    optimal_bridge = optimal_bridge_df
  ))
}

# Example Run with Parameters
params <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.5,
  sd_additive = 0.75,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results <- run_simulation(params)
print("Results for all bridge sizes:")
print(results$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results$optimal_bridge)




# Example Run with Parameters
params_50_samples <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 50,
  num_samples_batch2 = 50,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_50_samples <- run_simulation(params_50_samples)
print("Results for all bridge sizes:")
print(results_50_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_50_samples$optimal_bridge)







# Example Run with Parameters
params_100_samples <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_100_samples <- run_simulation(params_100_samples)
print("Results for all bridge sizes:")
print(results_100_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_100_samples$optimal_bridge)




# Example Run with Parameters
params_200_samples <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 200,
  num_samples_batch2 = 200,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_200_samples <- run_simulation(params_200_samples)
print("Results for all bridge sizes:")
print(results_200_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_200_samples$optimal_bridge)











# Example Run with Parameters
params_500_samples <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 500,
  num_samples_batch2 = 500,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_500_samples <- run_simulation(params_500_samples)
print("Results for all bridge sizes:")
print(results_500_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_500_samples$optimal_bridge)











# Example Run with Parameters
params_1000_samples <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 1000,
  num_samples_batch2 = 1000,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_1000_samples <- run_simulation(params_1000_samples)
print("Results for all bridge sizes:")
print(results_1000_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_1000_samples$optimal_bridge)























# Example Run with Parameters
params_weak<- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.25,
  sd_additive = 0.5,
  mean_multiplicative = 1.1,
  sd_multiplicative = 0.1,
  noise_mean = 0,
  noise_sd = 0.3
)

# Run the simulation
results_weak <- run_simulation(params_weak)
print("Results for all bridge sizes:")
print(results_weak$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_weak$optimal_bridge)











# Example Run with Parameters
params_standard <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_standard <- run_simulation(params_standard)
print("Results for all bridge sizes:")
print(results_standard$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_standard$optimal_bridge)















# Example Run with Parameters
params_strong <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 2.0,
  mean_multiplicative = 1.3,
  sd_multiplicative = 0.3,
  noise_mean = 0,
  noise_sd = 0.7
)

# Run the simulation
results_strong<- run_simulation(params_strong)
print("Results for all bridge sizes:")
print(results_strong$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_strong$optimal_bridge)



























# Example Run with Parameters
params_300_proteins <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 300,
  similar_assays = 150,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_300_proteins <- run_simulation(params_300_proteins)
print("Results for all bridge sizes:")
print(results_300_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_300_proteins$optimal_bridge)








# Example Run with Parameters
params_500_proteins <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 500,
  similar_assays = 250,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_500_proteins <- run_simulation(params_300_proteins)
print("Results for all bridge sizes:")
print(results_500_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_500_proteins$optimal_bridge)






# Example Run with Parameters
params_1000_proteins <- list(
  effect_sizes = c(0.1, 0.3, 0.5, 0.7, 0.9),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  num_assays = 1000,
  similar_assays = 500,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_1000_proteins <- run_simulation(params_1000_proteins)
print("Results for all bridge sizes:")
print(results_1000_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_1000_proteins$optimal_bridge)