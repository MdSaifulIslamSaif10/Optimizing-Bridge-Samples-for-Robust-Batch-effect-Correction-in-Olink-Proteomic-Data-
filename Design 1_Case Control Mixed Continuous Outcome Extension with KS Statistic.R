library(MASS)
library(tidyr)
library(OlinkAnalyze)
library(stringr)
library(dplyr)
library(broom)
library(sva)
library(truncnorm)
library(Matrix)
library(stats)
library(ggplot2)



# Function to create covariance matrix from distribution
create_cov_matrix_from_distribution <- function(num_assays, var_range = c(0.1, 4.0), rho_range = c(0.1, 0.6)) {
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

# Function 1: Generate Batch Data with Mandatory KS Statistic Control
generate_batch_data <- function(num_samples, num_assays, effect_size, similar_assays,
                                case_proportion = 0.5, batch_id = "1", true_betas = NULL,
                                noise_sd = 0.1, ks_target = 0.9, batch1_baseline_y = NULL) {
  num_cases <- round(num_samples * case_proportion)
  num_controls <- num_samples - num_cases
  if (num_cases < 1 || num_controls < 1) stop("Invalid case/control counts.")
  
  case_ids <- paste0("A", batch_id, "_C", 1:num_cases)
  control_ids <- paste0("A", batch_id, "_N", 1:num_controls)
  sample_ids <- c(case_ids, control_ids)
  
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
  npx_disease_means_similar <- seq(2, 20, length.out = similar_assays)
  npx_disease_means_different <- seq(2, 20, length.out = different_assays)
  npx_disease_means <- c(npx_disease_means_similar, npx_disease_means_different)
  npx_control_means <- c(npx_disease_means_similar, npx_disease_means_different + effect_size)
  
  cov_matrix <- create_cov_matrix_from_distribution(num_assays = num_assays, var_range = c(0.1, 4.0), rho_range = c(0.1, 0.6))
  npx_data <- rbind(
    if (num_cases > 0) mvrnorm(n = num_cases, mu = npx_disease_means, Sigma = cov_matrix),
    if (num_controls > 0) mvrnorm(n = num_controls, mu = npx_control_means, Sigma = cov_matrix)
  )
  
  npx_long <- as.data.frame(npx_data) %>%
    setNames(paste0("Assay_", 1:num_assays)) %>%
    mutate(SampleID = sample_ids) %>%
    tidyr::pivot_longer(cols = starts_with("Assay_"), names_to = "Assay", values_to = "NPX")
  
  # Baseline generation with KS control
  if (batch_id == "1") {
    baseline_y <- rnorm(num_samples, mean = 1, sd = 1)
  } else if (batch_id == "2") {
    if (is.null(batch1_baseline_y)) stop("batch1_baseline_y must be provided for Batch 2 KS control.")
    optimize_ks_diff <- function(params, y1, target_ks) {
      mean_diff <- params[1]
      sd_diff <- params[2]
      if (sd_diff <= 0) return(1e6)
      y2 <- rnorm(length(y1), mean = mean_diff, sd = sd_diff)
      ks_stat <- suppressWarnings(ks.test(y1, y2)$statistic)
      return(abs(ks_stat - target_ks))
    }
    set.seed(123 + as.numeric(batch_id))
    init_params <- c(mean_diff = 3, sd_diff = 1.5)
    opt_result <- optim(init_params, optimize_ks_diff, y1 = batch1_baseline_y, target_ks = ks_target,
                        method = "L-BFGS-B", lower = c(0.5, 0.5), upper = c(5, 2))
    baseline_y <- rnorm(num_samples, mean = opt_result$par[1], sd = opt_result$par[2])
  }
  
  if (is.null(true_betas)) {
    true_betas <- numeric(num_assays)
    true_betas[1:similar_assays] <- 0
    true_betas[(similar_assays + 1):num_assays] <- effect_size
  }
  
  batch_data <- batch_data %>%
    dplyr::left_join(npx_long, by = c("SampleID", "Assay")) %>%
    group_by(SampleID) %>%
    mutate(y_continuous = rep(baseline_y[match(unique(SampleID), sample_ids)], each = num_assays) +
             true_betas[as.numeric(str_extract(Assay, "\\d+"))] * NPX + rnorm(n(), mean = 0, sd = 0.05)) %>%
    ungroup() %>%
    mutate(Panel_Version = "v.1201",
           PlateID = paste0("Example_Data_", batch_id, "_CAM.csv"),
           Subject = rep(paste0("ID", sample(1:num_samples, num_samples, replace = TRUE)), each = num_assays),
           Treatment = rep(c(rep("Treated", num_cases), rep("Untreated", num_controls)), each = num_assays),
           Site = rep(sample(c("Site_A", "Site_B", "Site_C", "Site_D"), num_samples, replace = TRUE), each = num_assays),
           Time = rep(sample(c("Baseline", "Week.6", "Week.12"), num_samples, replace = TRUE), each = num_assays),
           Project = paste0("data", batch_id),
           Panel = "Olink Cardiometabolic")
  
  attr(batch_data, "true_betas") <- true_betas
  return(list(data = batch_data, baseline_y = baseline_y))
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
    dplyr::summarise(estimated_mean = mean(NPX, na.rm = TRUE), .groups = "drop")
  
  batch_effects <- ref_stats %>%
    mutate(additive_effects = rnorm(n(), mean = mean_additive, sd = sd_additive),
           multiplicative_effects = rnorm(n(), mean = mean_multiplicative, sd = sd_multiplicative))
  
  batch_data_updated <- batch_data %>%
    dplyr::left_join(batch_effects, by = "Assay") %>%
    mutate(NPX = (NPX - estimated_mean) * multiplicative_effects + additive_effects + estimated_mean,
           Project = "data2") %>%
    dplyr::select(-additive_effects, -multiplicative_effects, -estimated_mean)
  
  batch_data_updated <- batch_data_updated %>%
    group_by(Assay) %>%
    mutate(NPX = NPX + rnorm(n(), mean = noise_mean, sd = noise_sd)) %>%
    ungroup()
  
  return(batch_data_updated)
}

# Function 4: Combine Batches
combine_batches <- function(batch1_data, batch2_data, remove_bridge = FALSE, overlap_samples = NULL) {
  combined_data <- dplyr::bind_rows(batch1_data, batch2_data)
  if (remove_bridge && !is.null(overlap_samples)) {
    combined_data <- combined_data %>% dplyr::filter(!(SampleID %in% overlap_samples & Project == "data2"))
  }
  return(combined_data)
}

# Function 5: Perform Statistical Analysis
perform_statistical_analysis <- function(data, similar_assays, true_betas, alpha = 0.05) {
  lm_results <- data %>%
    dplyr::group_by(Assay) %>%
    do(tidy(lm(y_continuous ~ NPX, data = .))) %>%
    dplyr::filter(term == "NPX") %>%
    dplyr::ungroup() %>%
    mutate(AssayNumber = as.numeric(str_extract(Assay, "\\d+")),
           true_beta = true_betas[AssayNumber],
           significant = p.value < alpha,
           detected_beta = ifelse(significant, estimate, 0))
  
  power_continuous <- lm_results %>% filter(true_beta != 0) %>% summarise(power = mean(significant, na.rm = TRUE)) %>% dplyr::pull(power)
  typeI_error_continuous <- lm_results %>% filter(true_beta == 0) %>% summarise(rate = mean(significant, na.rm = TRUE)) %>% dplyr::pull(rate)
  fdr_continuous <- lm_results %>% dplyr::filter(significant) %>% summarise(fdr = ifelse(n() == 0, 0, mean(true_beta == 0, na.rm = TRUE))) %>% dplyr::pull(fdr)
  
  return(list(power_continuous = power_continuous, typeI_error_continuous = typeI_error_continuous, fdr_continuous = fdr_continuous))
}

# Function 6: Compute DSC
compute_dsc <- function(data) {
  data_wide <- data %>% dplyr::select(SampleID, Assay, NPX, Project) %>% tidyr::pivot_wider(names_from = Assay, values_from = NPX)
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
  overlap_samples_list <- list("DF1" = overlap_samples, "DF2" = overlap_samples)
  npx_br_data <- olink_normalization_bridge(project_1_df = batch1_data %>% as_tibble(),
                                            project_2_df = batch2_data %>% as_tibble(),
                                            bridge_samples = overlap_samples_list,
                                            project_1_name = project_names[1],
                                            project_2_name = project_names[2],
                                            project_ref_name = project_names[1])
  npx_br_data_no_overlap <- npx_br_data %>% dplyr::filter(!(SampleID %in% overlap_samples & Project == project_names[2]))
  return(npx_br_data_no_overlap)
}

# Function 8: Apply ComBat Correction
combat_correction <- function(data, batch_col = "Project", ref_batch = "data1") {
  wide_data <- data %>% dplyr::select(SampleID, Project, Treatment, Assay, NPX) %>% tidyr::pivot_wider(id_cols = c(SampleID, Project, Treatment),
                                                                                                       names_from = Assay, values_from = NPX)
  expr_mat <- wide_data %>% dplyr::select(-SampleID, -Project, -Treatment) %>% as.matrix()
  expr_mat <- t(expr_mat)
  colnames(expr_mat) <- wide_data$SampleID
  pheno <- wide_data %>% dplyr::select(SampleID, Project)
  pheno <- pheno[match(colnames(expr_mat), pheno$SampleID), ]
  batch_labels <- factor(pheno[[batch_col]])
  combat_expr <- ComBat(dat = expr_mat, batch = batch_labels, mod = model.matrix(~1, data = pheno), ref.batch = ref_batch)
  combat_wide <- as.data.frame(t(combat_expr))
  combat_wide$SampleID <- rownames(combat_wide)
  combat_long <- combat_wide %>% tidyr::pivot_longer(cols = -SampleID, names_to = "Assay", values_to = "NPX")
  final_corrected <- data %>% dplyr::select(-NPX) %>% dplyr::left_join(combat_long, by = c("SampleID", "Assay"))
  return(final_corrected)
}

# Main Function to Orchestrate the Simulation
run_simulation <- function(params) {
  final_results_df <- data.frame()
  optimal_bridge_df <- data.frame()
  
  epsilon_power <- 0.01
  epsilon_typeI <- 0.01
  
  for (eff in params$effect_sizes) {
    cat("Effect Size:", eff, "\n")
    
    for (iteration in 1:params$n_iterations) {
      cat("Iteration:", iteration, "\n")
      
      set.seed(123 + iteration)
      
      batch1_result <- generate_batch_data(
        num_samples = params$num_samples_batch1,
        num_assays = params$num_assays,
        effect_size = eff,
        similar_assays = params$similar_assays,
        case_proportion = params$case_proportion_batch1,
        batch_id = "1",
        true_betas = NULL,
        noise_sd = params$noise_sd,
        ks_target = params$ks_statistic
      )
      batch1_data <- batch1_result$data
      
      batch2_result <- generate_batch_data(
        num_samples = params$num_samples_batch2,
        num_assays = params$num_assays,
        effect_size = eff,
        similar_assays = params$similar_assays,
        case_proportion = params$case_proportion_batch2,
        batch_id = "2",
        true_betas = attr(batch1_data, "true_betas"),
        noise_sd = params$noise_sd,
        ks_target = params$ks_statistic,
        batch1_baseline_y = batch1_result$baseline_y
      )
      batch2_data <- batch2_result$data
      
      # Use all non-control batch1_data samples as bridge samples for batch effect creation
      bridge_samples_all <- batch1_data %>%
        dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
        mutate(Project = "data2")
      
      # Combine batch2_data with all batch1_data samples and apply batch effect once
      batch2_data_with_bridge <- dplyr::bind_rows(batch2_data, bridge_samples_all)
      batch2_data_updated <- add_batch_effects(
        batch_data = batch2_data_with_bridge,
        num_assays = params$num_assays,
        mean_additive = params$mean_additive,
        sd_additive = params$sd_additive,
        mean_multiplicative = params$mean_multiplicative,
        sd_multiplicative = params$sd_multiplicative,
        noise_mean = params$noise_mean,
        noise_sd = params$noise_sd
      )
      
      # Compute gold standard once per dataset
      combined_data_gold <- combine_batches(batch1_data, batch2_data)
      stats_gold <- perform_statistical_analysis(combined_data_gold, params$similar_assays, attr(batch1_data, "true_betas"))
      dsc_gold <- compute_dsc(combined_data_gold)
      
      power_results_list <- list()
      
      for (bridge_size in params$bridge_sizes) {
        cat("Bridge Sample Size:", bridge_size, "\n")
        
        results_list <- list(
          power_Gold_Standard_continuous = numeric(),
          power_with_batch_effect_continuous = numeric(),
          power_after_batch_effect_removal_continuous = numeric(),
          power_COMBAT_after_batch_effect_removal_continuous = numeric(),
          typeI_error_Gold_Standard_continuous = numeric(),
          typeI_error_with_batch_effect_continuous = numeric(),
          typeI_error_after_batch_effect_removal_continuous = numeric(),
          typeI_error_COMBAT_after_batch_effect_removal_continuous = numeric(),
          FDR_Gold_Standard_continuous = numeric(),
          FDR_with_batch_effect_continuous = numeric(),
          FDR_after_batch_effect_removal_continuous = numeric(),
          FDR_COMBAT_after_batch_effect_removal_continuous = numeric(),
          DSC_Gold_Standard = numeric(),
          DSC_with_batch_effect = numeric(),
          DSC_after_batch_effect_removal = numeric(),
          DSC_COMBAT_after_batch_effect_removal = numeric()
        )
        
        # Select bridge samples from original batch1_data
        bridge_samples_subset <- select_bridge_samples(batch1_data, bridge_size)
        overlap_samples <- bridge_samples_subset$SampleID
        
        # Create batch2 data for normalization (excluding batch1 samples except selected bridge samples)
        batch2_data_for_norm <- batch2_data_updated %>%
          dplyr::filter(!(SampleID %in% bridge_samples_all$SampleID) | SampleID %in% overlap_samples)
        
        # Combine batches, removing all batch1 samples for analysis
        combined_data_batch <- combine_batches(
          batch1_data, batch2_data_updated, remove_bridge = TRUE, overlap_samples = bridge_samples_all$SampleID
        )
        stats_batch <- perform_statistical_analysis(combined_data_batch, params$similar_assays, attr(batch1_data, "true_betas"))
        dsc_batch <- compute_dsc(combined_data_batch)
        
        # Normalize data using selected bridge samples
        normalized_data <- normalize_data(
          batch1_data = batch1_data,
          batch2_data = batch2_data_for_norm,
          overlap_samples = overlap_samples
        )
        stats_normalized <- perform_statistical_analysis(normalized_data, params$similar_assays, attr(batch1_data, "true_betas"))
        dsc_normalized <- compute_dsc(normalized_data)
        
        # Apply ComBat correction
        combat_data <- combat_correction(combined_data_batch)
        stats_combat <- perform_statistical_analysis(combat_data, params$similar_assays, attr(batch1_data, "true_betas"))
        dsc_combat <- compute_dsc(combat_data)
        
        # Store results
        results_list$power_Gold_Standard_continuous <- stats_gold$power_continuous
        results_list$power_with_batch_effect_continuous <- stats_batch$power_continuous
        results_list$power_after_batch_effect_removal_continuous <- stats_normalized$power_continuous
        results_list$power_COMBAT_after_batch_effect_removal_continuous <- stats_combat$power_continuous
        results_list$typeI_error_Gold_Standard_continuous <- stats_gold$typeI_error_continuous
        results_list$typeI_error_with_batch_effect_continuous <- stats_batch$typeI_error_continuous
        results_list$typeI_error_after_batch_effect_removal_continuous <- stats_normalized$typeI_error_continuous
        results_list$typeI_error_COMBAT_after_batch_effect_removal_continuous <- stats_combat$typeI_error_continuous
        results_list$FDR_Gold_Standard_continuous <- stats_gold$fdr_continuous
        results_list$FDR_with_batch_effect_continuous <- stats_batch$fdr_continuous
        results_list$FDR_after_batch_effect_removal_continuous <- stats_normalized$fdr_continuous
        results_list$FDR_COMBAT_after_batch_effect_removal_continuous <- stats_combat$fdr_continuous
        results_list$DSC_Gold_Standard <- dsc_gold
        results_list$DSC_with_batch_effect <- dsc_batch
        results_list$DSC_after_batch_effect_removal <- dsc_normalized
        results_list$DSC_COMBAT_after_batch_effect_removal <- dsc_combat
        
        # Store results for this bridge size
        power_results_df <- data.frame(
          Bridge_Size = bridge_size,
          KS_Statistic = params$ks_statistic,
          Effect_Size = eff,
          Iteration = iteration,
          power_Gold_Standard_continuous = mean(results_list$power_Gold_Standard_continuous, na.rm = TRUE),
          power_with_batch_effect_continuous = mean(results_list$power_with_batch_effect_continuous, na.rm = TRUE),
          power_after_batch_effect_removal_continuous = mean(results_list$power_after_batch_effect_removal_continuous, na.rm = TRUE),
          power_COMBAT_after_batch_effect_removal_continuous = mean(results_list$power_COMBAT_after_batch_effect_removal_continuous, na.rm = TRUE),
          typeI_error_Gold_Standard_continuous = mean(results_list$typeI_error_Gold_Standard_continuous, na.rm = TRUE),
          typeI_error_with_batch_effect_continuous = mean(results_list$typeI_error_with_batch_effect_continuous, na.rm = TRUE),
          typeI_error_after_batch_effect_removal_continuous = mean(results_list$typeI_error_after_batch_effect_removal_continuous, na.rm = TRUE),
          typeI_error_COMBAT_after_batch_effect_removal_continuous = mean(results_list$typeI_error_COMBAT_after_batch_effect_removal_continuous, na.rm = TRUE),
          FDR_Gold_Standard_continuous = mean(results_list$FDR_Gold_Standard_continuous, na.rm = TRUE),
          FDR_with_batch_effect_continuous = mean(results_list$FDR_with_batch_effect_continuous, na.rm = TRUE),
          FDR_after_batch_effect_removal_continuous = mean(results_list$FDR_after_batch_effect_removal_continuous, na.rm = TRUE),
          FDR_COMBAT_after_batch_effect_removal_continuous = mean(results_list$FDR_COMBAT_after_batch_effect_removal_continuous, na.rm = TRUE),
          DSC_Gold_Standard = mean(results_list$DSC_Gold_Standard, na.rm = TRUE),
          DSC_with_batch_effect = mean(results_list$DSC_with_batch_effect, na.rm = TRUE),
          DSC_after_batch_effect_removal = mean(results_list$DSC_after_batch_effect_removal, na.rm = TRUE),
          DSC_COMBAT_after_batch_effect_removal = mean(results_list$DSC_COMBAT_after_batch_effect_removal, na.rm = TRUE)
        )
        
        power_results_list[[as.character(bridge_size)]] <- power_results_df
      }
      
      power_results_df <- do.call(rbind, power_results_list)
      
      power_results_df_tmp <- power_results_df %>%
        mutate(
          power_diff = abs(power_after_batch_effect_removal_continuous - power_Gold_Standard_continuous),
          typeI_error_diff = abs(typeI_error_after_batch_effect_removal_continuous - typeI_error_Gold_Standard_continuous)
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
          KS_Statistic = params$ks_statistic,
          Effect_Size = eff,
          Iteration = iteration,
          Optimal_Bridge_Size_power = optimal_row_power$Bridge_Size,
          Power_After_continuous = optimal_row_power$power_after_batch_effect_removal_continuous,
          Power_Gold_continuous = optimal_row_power$power_Gold_Standard_continuous,
          Power_Diff = optimal_row_power$power_diff,
          Min_Power_Diff = min_power_diff,
          Power_Diff_Difference = optimal_row_power$power_diff - min_power_diff,
          Optimal_Bridge_Size_typeI = optimal_row_typeI$Bridge_Size,
          TypeI_Error_After_continuous = optimal_row_typeI$typeI_error_after_batch_effect_removal_continuous,
          TypeI_Error_Gold_continuous = optimal_row_typeI$typeI_error_Gold_Standard_continuous,
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
  effect_sizes = c(0.05, 0.1, 0.4),
  bridge_sizes = c(2, 4, 6, 8, 10, 12, 14, 16, 18),
  n_iterations = 3,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.5,
  case_proportion_batch2 = 0.5,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1,
  sd_additive = 0.75,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5,
  noise_sd_y = 0.05,
  ks_statistic = 0.6
)

# Run the simulation
results <- run_simulation(params)















# Example Run with Parameters
params_50_50 <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.5,
  case_proportion_batch2 = 0.5,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_50_50<- run_simulation(params_50_50)


# Print results
print("Results for all bridge sizes:")
print(results_50_50$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_50_50$optimal_bridge)












# Example Run with Parameters
params_60_40 <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.6,
  case_proportion_batch2 = 0.4,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_60_40<- run_simulation(params_60_40)


# Print results
print("Results for all bridge sizes:")
print(results_60_40$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_60_40$optimal_bridge)
















# Example Run with Parameters
params_70_30 <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_70_30<- run_simulation(params_70_30)


# Print results
print("Results for all bridge sizes:")
print(results_70_30$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_70_30$optimal_bridge)


















# Example Run with Parameters
params_80_20 <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.8,
  case_proportion_batch2 = 0.2,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_80_20<- run_simulation(params_80_20)


# Print results
print("Results for all bridge sizes:")
print(results_80_20$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_80_20$optimal_bridge)

















# Example Run with Parameters
params_90_10 <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.9,
  case_proportion_batch2 = 0.1,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_90_10<- run_simulation(params_90_10)


# Print results
print("Results for all bridge sizes:")
print(results_90_10$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_90_10$optimal_bridge)





























# Example Run with Parameters
params_50_samples <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 50,
  num_samples_batch2 = 50,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_50_samples<- run_simulation(params_50_samples)


# Print results
print("Results for all bridge sizes:")
print(results_50_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_50_samples$optimal_bridge)










# Example Run with Parameters
params_100_samples <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_100_samples<- run_simulation(params_100_samples)


# Print results
print("Results for all bridge sizes:")
print(results_100_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_100_samples$optimal_bridge)




















# Example Run with Parameters
params_200_samples <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 200,
  num_samples_batch2 = 200,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_200_samples<- run_simulation(params_200_samples)


# Print results
print("Results for all bridge sizes:")
print(results_200_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_200_samples$optimal_bridge)









# Example Run with Parameters
params_500_samples <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 500,
  num_samples_batch2 = 500,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_500_samples<- run_simulation(params_500_samples)


# Print results
print("Results for all bridge sizes:")
print(results_500_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_500_samples$optimal_bridge)








# Example Run with Parameters
params_1000_samples <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 1000,
  num_samples_batch2 = 1000,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_1000_samples<- run_simulation(params_1000_samples)


# Print results
print("Results for all bridge sizes:")
print(results_1000_samples$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_1000_samples$optimal_bridge)
















# Example Run with Parameters
params_weak <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 0.5,
  mean_multiplicative = 1.1,
  sd_multiplicative = 0.1,
  noise_mean = 0,
  noise_sd = 0.3
)

# Run the simulation
results_weak<- run_simulation(params_weak)


# Print results
print("Results for all bridge sizes:")
print(results_weak$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_weak$optimal_bridge)














# Example Run with Parameters
params_standard <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_standard<- run_simulation(params_standard)


# Print results
print("Results for all bridge sizes:")
print(results_standard$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_standard$optimal_bridge)







# Example Run with Parameters
params_strong <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.5,
  sd_additive = 1.75,
  mean_multiplicative = 1.3,
  sd_multiplicative = 0.3,
  noise_mean = 0,
  noise_sd = 0.7
)

# Run the simulation
results_strong<- run_simulation(params_strong)


# Print results
print("Results for all bridge sizes:")
print(results_strong$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_strong$optimal_bridge)





















# Example Run with Parameters
params_300_proteins <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 300,
  similar_assays = 150,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation
results_300_proteins= run_simulation(params_300_proteins)


# Print results
print("Results for all bridge sizes:")
print(results_300_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_300_proteins$optimal_bridge)











# Example Run with Parameters
params_500_proteins <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 500,
  similar_assays = 250,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)
# Run the simulation
results_500_proteins= run_simulation(params_500_proteins)


# Print results
print("Results for all bridge sizes:")
print(results_500_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_500_proteins$optimal_bridge)












# Example Run with Parameters
params_1000_proteins <- list(
  effect_sizes = c(0.1,0.3,0.5,0.7,0.9),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 40,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 1000,
  similar_assays = 500,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)
# Run the simulation
results_1000_proteins= run_simulation(params_1000_proteins)


# Print results
print("Results for all bridge sizes:")
print(results_1000_proteins$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_1000_proteins$optimal_bridge)







# Example Run with Parameters
params_ks_0<- list(
  effect_sizes = c(0.05,0.1,0.2,0.3,0.4),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5,
  noise_sd_y = 0.05,
  ks_statistic = 0
)

# Run the simulation
results_ks_0= run_simulation(params_ks_0)


# Print results
print("Results for all bridge sizes:")
print(params_ks_0$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(params_ks_0$optimal_bridge)







# Example Run with Parameters
params_ks_0.3<- list(
  effect_sizes = c(0.05,0.1,0.2,0.3,0.4),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5,
  noise_sd_y = 0.05,
  ks_statistic = 0.3
)

# Run the simulation
results_ks_0.3= run_simulation(params_ks_0.3)


# Print results
print("Results for all bridge sizes:")
print(results_ks_0.3$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_ks_0.3$optimal_bridge)








# Example Run with Parameters
params_ks_0.6<- list(
  effect_sizes = c(0.05,0.1,0.2,0.3,0.4),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5,
  noise_sd_y = 0.05,
  ks_statistic = 0.6
)

# Run the simulation
results_ks_0.6= run_simulation(params_ks_0.6)


# Print results
print("Results for all bridge sizes:")
print(results_ks_0.6$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_ks_0.6$optimal_bridge)













# Example Run with Parameters
params_ks_0.9<- list(
  effect_sizes = c(0.05,0.1,0.2,0.3,0.4),
  bridge_sizes = c(2,4,6, 8, 10, 12,14, 16, 18, 20, 22,24, 26,28,30),
  n_iterations = 100,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 50,
  alpha = 0.05,
  mean_additive = 1.0,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5,
  noise_sd_y = 0.05,
  ks_statistic = 0.9
)

# Run the simulation
results_ks_0.9= run_simulation(params_ks_0.9)


# Print results
print("Results for all bridge sizes:")
print(results_ks_0.9$final_results)
print("Optimal Bridge Sizes for power and Type I error:")
print(results_ks_0.9$optimal_bridge)

