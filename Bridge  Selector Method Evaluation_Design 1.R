library(MASS)
library(tidyr)
library(OlinkAnalyze)
library(stringr)
library(dplyr)
library(broom)
library(truncnorm)
library(Matrix)
library(stats)
library(ggplot2)
library(tibble)

# Function to create covariance matrix from distribution
create_cov_matrix_from_distribution <- function(num_assays,
                                                var_range = c(0.1, 4.0),
                                                rho_range = c(0.1, 0.6)) {
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
                                case_proportion = 0.5, batch_id = "1") {
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
  
  cov_matrix <- create_cov_matrix_from_distribution(num_assays = num_assays,
                                                    var_range = c(0.1, 4.0),
                                                    rho_range = c(0.1, 0.6))
  
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

# Function 2: Select Bridge Samples (Olink-based)
select_bridge_samples_olink <- function(batch_data, bridge_size, sample_missing_freq = 0.1) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  unique_samples <- batch_data %>% dplyr::distinct(SampleID) %>% pull(SampleID)
  if (bridge_size > length(unique_samples)) stop("Bridge size exceeds available samples.")
  
  bridge_samples <- batch_data %>%
    dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
    olink_bridgeselector(sampleMissingFreq = sample_missing_freq, n = bridge_size) %>%
    dplyr::select(SampleID)
  
  return(bridge_samples)
}

# Function 2a: Select Bridge Samples (Random)
select_bridge_samples_random <- function(batch_data, bridge_size) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  unique_samples <- batch_data %>% dplyr::distinct(SampleID)
  if (bridge_size > nrow(unique_samples)) stop("Bridge size exceeds available samples.")
  
  bridge_samples <- unique_samples %>%
    dplyr::sample_n(bridge_size, replace = FALSE) %>%
    dplyr::select(SampleID)
  
  return(bridge_samples)
}

# Function 2c: Select Bridge Samples (Quadrant-only)
select_bridge_samples_quadrant <- function(batch_data, bridge_size) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  unique_samples <- batch_data %>% dplyr::distinct(SampleID)
  if (bridge_size > nrow(unique_samples)) stop("Bridge size exceeds available samples.")
  
  npx_wide <- batch_data %>%
    dplyr::select(SampleID, Assay, NPX) %>%
    tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
    column_to_rownames("SampleID")
  
  pca_result <- prcomp(npx_wide, scale. = TRUE)
  pca_scores <- as.data.frame(pca_result$x[, 1:2])
  pca_scores$SampleID <- rownames(pca_scores)
  colnames(pca_scores)[1:2] <- c("PC1", "PC2")
  
  quadrant_samples <- pca_scores %>%
    dplyr::filter(PC1 < 0, PC2 > 0) %>%
    dplyr::slice_sample(n = bridge_size, replace = FALSE)
  
  if (nrow(quadrant_samples) < bridge_size) {
    if (nrow(quadrant_samples) == 0) {
      message("Quadrant method: No samples found in quadrant PC1 < 0, PC2 > 0. Falling back to random sampling.")
    } else {
      message(sprintf("Quadrant method: Only %d samples found in quadrant PC1 < 0, PC2 > 0. Supplementing with %d random samples.",
                      nrow(quadrant_samples), bridge_size - nrow(quadrant_samples)))
    }
    remaining <- bridge_size - nrow(quadrant_samples)
    available_samples <- pca_scores %>%
      dplyr::filter(!SampleID %in% quadrant_samples$SampleID)
    additional_samples <- available_samples %>%
      dplyr::slice_sample(n = min(remaining, nrow(available_samples)), replace = FALSE) %>%
      pull(SampleID)
    bridge_samples <- unique(c(quadrant_samples$SampleID, additional_samples))[1:bridge_size]
  } else {
    bridge_samples <- quadrant_samples$SampleID
  }
  
  return(data.frame(SampleID = bridge_samples))
}

# Function 2d: Select Bridge Samples (Coverage-Dispersion with Outlier Removal)
select_bridge_samples_coverage <- function(batch_data, bridge_size) {
  if (bridge_size < 1) stop("Bridge size must be at least 1.")
  
  unique_samples <- batch_data %>% dplyr::distinct(SampleID)
  if (bridge_size > nrow(unique_samples)) stop("Bridge size exceeds available samples.")
  
  # Step 1: Outlier detection based on Olink's qc_outliers logic
  df <- batch_data %>%
    dplyr::filter(!str_detect(SampleID, "CONTROL")) %>%
    dplyr::filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"))
  
  qc_outliers <- df %>%
    dplyr::group_by(Panel, SampleID) %>%
    dplyr::mutate(IQR = IQR(NPX, na.rm = TRUE), sample_median = median(NPX, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(SampleID, Panel, IQR, sample_median) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Panel) %>%
    dplyr::mutate(
      median_low = mean(sample_median, na.rm = TRUE) - 3 * sd(sample_median, na.rm = TRUE),
      median_high = mean(sample_median, na.rm = TRUE) + 3 * sd(sample_median, na.rm = TRUE),
      iqr_low = mean(IQR, na.rm = TRUE) - 3 * sd(IQR, na.rm = TRUE),
      iqr_high = mean(IQR, na.rm = TRUE) + 3 * sd(IQR, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Outlier = dplyr::if_else(sample_median < median_high & sample_median > median_low & IQR > iqr_low & IQR < iqr_high, 0, 1)) %>%
    dplyr::select(SampleID, Panel, Outlier)
  
  # Step 2: Remove outlier samples
  df_filtered <- df %>%
    dplyr::left_join(qc_outliers, by = c("SampleID", "Panel")) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Outliers = sum(Outlier)) %>%
    dplyr::filter(Outliers == 0) %>%
    dplyr::ungroup()
  
  # Step 3: Proceed with Coverage-Dispersion on remaining samples
  npx_wide <- df_filtered %>%
    dplyr::select(SampleID, Assay, NPX) %>%
    tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
    column_to_rownames("SampleID")
  
  pca_result <- prcomp(npx_wide, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$SampleID <- rownames(pca_df)
  colnames(pca_df)[1:2] <- c("PC1", "PC2")
  
  greedy_res <- greedy_CoverageDispersion(pca_df, k = bridge_size)
  bridge_samples <- pca_df$SampleID[greedy_res$indices]
  
  return(data.frame(SampleID = bridge_samples))
}

# Function 3: Add Batch Effects
add_batch_effects <- function(batch_data, num_assays, bridge_samples = NULL,
                              mean_additive = 0.5, sd_additive = 0.4,
                              mean_multiplicative = 1.1, sd_multiplicative = 0.1,
                              noise_mean = 0, noise_sd = 0.4) {
  original_samples <- if (!is.null(bridge_samples)) {
    batch_data %>% dplyr::filter(!(SampleID %in% bridge_samples$SampleID))
  } else {
    batch_data
  }
  
  ref_stats <- original_samples %>%
    dplyr::group_by(Assay) %>%
    dplyr::summarise(
      estimated_mean = mean(NPX, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(Assay, paste0("Assay_", 1:num_assays)))
  
  additive_effects_df <- data.frame(
    additive_effects_raw = rnorm(n = nrow(ref_stats), mean = mean_additive, sd = sd_additive)
  ) %>%
    mutate(
      abs_additive = abs(additive_effects_raw),
      order = row_number()
    ) %>%
    arrange(abs_additive) %>%
    mutate(
      assay_index = 1:n()
    ) %>%
    arrange(order) %>%
    dplyr::select(additive_effects_raw, assay_index)
  
  batch_effects <- ref_stats %>%
    mutate(
      assay_index = 1:n(),
      multiplicative_effects = rtruncnorm(
        n(), a = 0, b = Inf, mean = mean_multiplicative, sd = sd_multiplicative
      )
    ) %>%
    left_join(additive_effects_df, by = "assay_index") %>%
    mutate(
      additive_effects = additive_effects_raw
    ) %>%
    dplyr::select(Assay, estimated_mean, additive_effects, multiplicative_effects)
  
  batch_data_updated <- batch_data %>%
    dplyr::left_join(batch_effects, by = "Assay") %>%
    dplyr::mutate(
      NPX = (NPX - estimated_mean) * multiplicative_effects + additive_effects + estimated_mean,
      Project = "data2"
    ) %>%
    dplyr::select(-additive_effects, -multiplicative_effects, -estimated_mean)
  
  batch_data_updated <- batch_data_updated %>%
    dplyr::group_by(Assay) %>%
    dplyr::mutate(
      NPX = NPX + rnorm(n(), mean = noise_mean, sd = noise_sd)
    ) %>%
    dplyr::ungroup()
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
    dplyr::do(tidy(t.test(NPX ~ Treatment, data = .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      AssayNumber = as.numeric(stringr::str_extract(Assay, "\\d+")),
      truth = ifelse(AssayNumber <= similar_assays, "null", "alternative"),
      significant = p.value < alpha
    )
  
  power <- t_test_results %>%
    dplyr::filter(truth == "alternative") %>%
    dplyr::summarise(power = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(power)
  
  typeI_error <- t_test_results %>%
    dplyr::filter(truth == "null") %>%
    dplyr::summarise(rate = mean(significant, na.rm = TRUE)) %>%
    dplyr::pull(rate)
  
  fdr <- t_test_results %>%
    dplyr::filter(significant) %>%
    dplyr::summarise(fdr = ifelse(n() == 0, 0, mean(truth == "null", na.rm = TRUE))) %>%
    dplyr::pull(fdr)
  
  return(list(power = power, typeI_error = typeI_error, fdr = fdr))
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
  overlap_samples_list <- list("DF1" = overlap_samples, "DF2" = overlap_samples)
  
  npx_br_data <- olink_normalization_bridge(
    project_1_df = batch1_data %>% as_tibble(),
    project_2_df = batch2_data %>% as_tibble(),
    bridge_samples = overlap_samples_list,
    project_1_name = project_names[1],
    project_2_name = project_names[2],
    project_ref_name = project_names[1]
  )
  
  npx_br_data_no_overlap <- npx_br_data %>%
    dplyr::filter(!(SampleID %in% overlap_samples & Project == project_names[2]))
  
  return(npx_br_data_no_overlap)
}

# Function 9: Compute Coverage-Dispersion
compute_CoverageDispersion <- function(pca_df, bridge_ids) {
  bridge <- pca_df %>% dplyr::filter(SampleID %in% bridge_ids)
  non <- pca_df %>% dplyr::filter(!SampleID %in% bridge_ids)
  
  if (nrow(bridge) == 0 | nrow(non) == 0) return(NA_real_)
  
  D <- as.matrix(dist(rbind(non[, c("PC1", "PC2")], bridge[, c("PC1", "PC2")])))
  nb <- nrow(non)
  k <- nrow(bridge)
  nearest <- apply(D[1:nb, (nb + 1):(nb + k)], 1, min)
  
  SSCE <- sum(nearest^2)
  centre <- colMeans(pca_df[, c("PC1", "PC2")])
  SSTC <- sum(rowSums((non[, c("PC1", "PC2")] - centre)^2))
  
  1 - SSCE / SSTC
}

# Function 10: Greedy Coverage-Dispersion Selection
greedy_CoverageDispersion <- function(pca_df, k) {
  n <- nrow(pca_df)
  coords <- as.matrix(pca_df[, c("PC1", "PC2")])
  
  D <- as.matrix(dist(coords))
  
  x_bar <- colMeans(coords)
  SSTC <- sum(rowSums((coords - matrix(x_bar, nrow = n, ncol = 2, byrow = TRUE))^2))
  
  B <- integer(0)
  d_prev <- rep(NA_real_, n)
  CD_path <- numeric(k)
  
  for (t in seq_len(k)) {
    best_idx <- NA_integer_
    best_SSCE <- Inf
    
    for (cand in setdiff(seq_len(n), B)) {
      d_new <- ifelse(is.na(d_prev), D[, cand], pmin(d_prev, D[, cand]))
      SSCE <- sum(d_new^2)
      if (SSCE < best_SSCE) {
        best_SSCE <- SSCE
        best_idx <- cand
      }
    }
    
    B <- c(B, best_idx)
    d_prev <- ifelse(is.na(d_prev), D[, best_idx], pmin(d_prev, D[, best_idx]))
    CD_path[t] <- 1 - best_SSCE / SSTC
  }
  
  list(indices = B, CD_final = CD_path[k], CD_curve = CD_path)
}

# Modified Main Function to Orchestrate the Simulation
run_simulation <- function(params, effect_sizes = seq(0.1, 0.5, by = 0.1)) {
  all_results <- list()
  
  for (effect_size in effect_sizes) {
    cat("Running simulation for effect size:", effect_size, "\n")
    params$effect_size <- effect_size
    final_results_df <- data.frame()
    
    w1 <- 0.5
    w2 <- 0.5
    selection_methods <- c("Olink", "Random", "Quadrant", "Coverage_R2")
    bridge_size <- params$bridge_size
    
    power_diff_per_iteration <- list()
    typeI_error_diff_per_iteration <- list()
    CD_value_per_iteration <- list()
    loss_score_per_iteration <- list()
    
    for (iteration in 1:params$n_iterations) {
      cat("Iteration:", iteration, "\n")
      
      batch1_data <- generate_batch_data(
        num_samples = params$num_samples_batch1,
        num_assays = params$num_assays,
        effect_size = effect_size,
        similar_assays = params$similar_assays,
        case_proportion = params$case_proportion_batch1,
        batch_id = "1"
      )
      batch2_data <- generate_batch_data(
        num_samples = params$num_samples_batch2,
        num_assays = params$num_assays,
        effect_size = effect_size,
        similar_assays = params$similar_assays,
        case_proportion = params$case_proportion_batch2,
        batch_id = "2"
      )
      
      npx_wide <- batch1_data %>%
        dplyr::select(SampleID, Assay, NPX) %>%
        tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
        column_to_rownames("SampleID")
      pca_result <- prcomp(npx_wide, scale. = TRUE)
      pca_df <- as.data.frame(pca_result$x[, 1:2])
      pca_df$SampleID <- rownames(pca_df)
      colnames(pca_df)[1:2] <- c("PC1", "PC2")
      
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
      
      combined_data_gold <- combine_batches(batch1_data, batch2_data)
      stats_gold <- perform_statistical_analysis(combined_data_gold, params$similar_assays, params$alpha)
      dsc_gold <- compute_dsc(combined_data_gold)
      
      for (method in selection_methods) {
        cat("Selection Method:", method, "\n")
        
        if (method == "Olink") {
          bridge_samples <- select_bridge_samples_olink(batch1_data, bridge_size)
        } else if (method == "Random") {
          bridge_samples <- select_bridge_samples_random(batch1_data, bridge_size)
        } else if (method == "Quadrant") {
          bridge_samples <- select_bridge_samples_quadrant(batch1_data, bridge_size)
        } else if (method == "Coverage_R2") {
          bridge_samples <- select_bridge_samples_coverage(batch1_data, bridge_size)
        }
        bridge_ids <- unique(bridge_samples$SampleID)
        cd_value <- compute_CoverageDispersion(pca_df, bridge_ids)
        
        batch2_data_for_norm <- batch2_data_updated %>%
          dplyr::filter(!(SampleID %in% bridge_samples_all$SampleID) | SampleID %in% bridge_ids)
        
        combined_data_batch <- combine_batches(
          batch1_data, batch2_data_updated, remove_bridge = TRUE, overlap_samples = bridge_samples_all$SampleID
        )
        stats_batch <- perform_statistical_analysis(combined_data_batch, params$similar_assays, params$alpha)
        dsc_batch <- compute_dsc(combined_data_batch)
        
        normalized_data <- normalize_data(
          batch1_data = batch1_data,
          batch2_data = batch2_data_for_norm,
          overlap_samples = bridge_ids
        )
        stats_normalized <- perform_statistical_analysis(normalized_data, params$similar_assays, params$alpha)
        dsc_normalized <- compute_dsc(normalized_data)
        
        loss_score <- w1 * abs(stats_normalized$power - stats_gold$power) +
          w2 * abs(stats_normalized$typeI_error - stats_gold$typeI_error)
        
        power_diff <- stats_gold$power - stats_normalized$power
        typeI_error_diff <- stats_gold$typeI_error - stats_normalized$typeI_error
        
        if (!method %in% names(power_diff_per_iteration)) {
          power_diff_per_iteration[[method]] <- numeric()
          typeI_error_diff_per_iteration[[method]] <- numeric()
          CD_value_per_iteration[[method]] <- numeric()
          loss_score_per_iteration[[method]] <- numeric()
        }
        power_diff_per_iteration[[method]] <- c(power_diff_per_iteration[[method]], power_diff)
        typeI_error_diff_per_iteration[[method]] <- c(typeI_error_diff_per_iteration[[method]], typeI_error_diff)
        CD_value_per_iteration[[method]] <- c(CD_value_per_iteration[[method]], cd_value)
        loss_score_per_iteration[[method]] <- c(loss_score_per_iteration[[method]], loss_score)
        
        results_df <- data.frame(
          Effect_Size = effect_size,
          Bridge_Size = bridge_size,
          Selection_Method = method,
          Iteration = iteration,
          power_Gold_Standard = stats_gold$power,
          power_with_batch_effect = stats_batch$power,
          power_after_batch_effect_removal = stats_normalized$power,
          typeI_error_Gold_Standard = stats_gold$typeI_error,
          typeI_error_with_batch_effect = stats_batch$typeI_error,
          typeI_error_after_batch_effect_removal = stats_normalized$typeI_error,
          FDR_Gold_Standard = stats_gold$fdr,
          FDR_with_batch_effect = stats_batch$fdr,
          FDR_after_batch_effect_removal = stats_normalized$fdr,
          DSC_Gold_Standard = dsc_gold,
          DSC_with_batch_effect = dsc_batch,
          DSC_after_batch_effect_removal = dsc_normalized,
          CD_value = cd_value,
          power_diff = power_diff,
          typeI_error_diff = typeI_error_diff,
          loss_score = loss_score
        )
        
        final_results_df <- rbind(final_results_df, results_df)
      }
    }
    
    final_results_df <- final_results_df %>%
      group_by(Selection_Method, Effect_Size) %>%
      mutate(
        mean_power_diff = mean(power_diff, na.rm = TRUE),
        mean_typeI_error_diff = mean(typeI_error_diff, na.rm = TRUE),
        mean_CD_value = mean(CD_value, na.rm = TRUE),
        mean_loss_score = mean(loss_score, na.rm = TRUE)
      ) %>%
      ungroup()
    
    attr(final_results_df, "power_diff_per_iteration") <- power_diff_per_iteration
    attr(final_results_df, "typeI_error_diff_per_iteration") <- typeI_error_diff_per_iteration
    attr(final_results_df, "CD_value_per_iteration") <- CD_value_per_iteration
    attr(final_results_df, "loss_score_per_iteration") <- loss_score_per_iteration
    
    all_results[[as.character(effect_size)]] <- final_results_df
  }
  
  combined_results <- do.call(rbind, all_results)
  return(combined_results)
}

# Function to Plot Boxplots
plot_boxplot <- function(data, metric, title) {
  metric_sym <- rlang::sym(metric)
  
  desired_order <- c("Coverage_R2", "Olink", "Random", "Quadrant")
  
  method_labels <- c("Coverage_R2" = "Coverage",
                     "Olink" = "Olink",
                     "Random" = "Random",
                     "Quadrant" = "Quadrant")
  
  data <- data %>%
    mutate(Selection_Method = factor(Selection_Method, levels = desired_order))
  
  stats_tbl <- data %>%
    group_by(Selection_Method, Effect_Size) %>%
    summarise(
      mean_val = mean(!!metric_sym, na.rm = TRUE),
      med_val = median(!!metric_sym, na.rm = TRUE),
      sd_val = sd(!!metric_sym, na.rm = TRUE),
      y_max = max(!!metric_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_lab = y_max + 0.05 * diff(range(data[[metric]], na.rm = TRUE)),
      label = sprintf("μ = %.3f\nσ = %.3f", mean_val, sd_val)
    )
  
  y_axis_label <- switch(metric,
                         "power_diff" = "Power Difference",
                         "typeI_error_diff" = "Type I Error Difference",
                         "CD_value" = "Dispersion-Coverage Score",
                         "loss_score" = "Loss Score",
                         metric)
  
  ggplot(data, aes(x = Selection_Method, y = !!metric_sym)) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue", width = .48) +
    geom_jitter(width = .12, size = 2, alpha = .6) +
    geom_point(data = stats_tbl, aes(y = med_val),
               shape = 23, size = 4, fill = "white", colour = "black") +
    geom_point(data = stats_tbl, aes(y = mean_val),
               shape = 21, size = 4, fill = "red", colour = "black") +
    geom_label(data = stats_tbl, aes(y = y_lab, label = label),
               size = 6.5, fontface = "bold",
               label.padding = unit(0.15, "lines"),
               label.size = 0, fill = "white", alpha = .8) +
    scale_x_discrete(drop = FALSE,
                     limits = desired_order,
                     labels = method_labels) +
    facet_wrap(~ Effect_Size, labeller = label_both) +
    labs(title = title,
         x = "Feature Selection Method",
         y = y_axis_label) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 25),
      axis.title = element_text(face = "bold", size = 25),
      axis.text = element_text(angle = 45, hjust = 1,
                               face = "bold", size = 25),
      plot.caption = element_text(size = 25, hjust = 1),
      strip.text = element_text(size = 20, face = "bold")
    )
}

# Example Run with Parameters
params <- list(
  bridge_size = 8,
  n_iterations = 50,
  num_samples_batch1 = 100,
  num_samples_batch2 = 100,
  case_proportion_batch1 = 0.7,
  case_proportion_batch2 = 0.3,
  num_assays = 100,
  similar_assays = 40,
  alpha = 0.05,
  mean_additive = 0.5,
  sd_additive = 1.2,
  mean_multiplicative = 1.2,
  sd_multiplicative = 0.2,
  noise_mean = 0,
  noise_sd = 0.5
)

# Run the simulation for multiple effect sizes
effect_sizes <- seq(0.1, 0.9, by = 0.2)
results <- run_simulation(params, effect_sizes)
print("Simulation Results:")
print(results)

# New Function to Plot Power vs Effect Size
plot_power_vs_effect_size <- function(data) {
  avg_data <- data %>%
    group_by(Selection_Method, Effect_Size) %>%
    summarise(
      Gold_Power = mean(power_Gold_Standard, na.rm = TRUE),
      Method_Power = mean(power_after_batch_effect_removal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(Gold_Power, Method_Power),
      names_to = "Type",
      values_to = "Power"
    ) %>%
    mutate(
      Method = if_else(Type == "Gold_Power", "Gold Standard", Selection_Method),
      Method = factor(Method, levels = c("Gold Standard", "Coverage_R2", "Olink", "Random", "Quadrant"))
    )
  
  ggplot(avg_data, aes(x = Effect_Size, y = Power, color = Method, group = Method)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Gold Standard" = "red", "Coverage_R2" = "blue",
                                  "Olink" = "green", "Random" = "purple", "Quadrant" = "orange")) +
    labs(title = "Power vs Effect Size", x = "Effect Size", y = "Power") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 10),
      legend.position = "top"
    )
}

# New Function to Plot Type I Error vs Effect Size
plot_typeI_error_vs_effect_size <- function(data) {
  avg_data <- data %>%
    group_by(Selection_Method, Effect_Size) %>%
    summarise(
      Gold_TypeI = mean(typeI_error_Gold_Standard, na.rm = TRUE),
      Method_TypeI = mean(typeI_error_after_batch_effect_removal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(Gold_TypeI, Method_TypeI),
      names_to = "Type",
      values_to = "TypeI_Error"
    ) %>%
    mutate(
      Method = if_else(Type == "Gold_TypeI", "Gold Standard", Selection_Method),
      Method = factor(Method, levels = c("Gold Standard", "Coverage_R2", "Olink", "Random", "Quadrant"))
    )
  
  ggplot(avg_data, aes(x = Effect_Size, y = TypeI_Error, color = Method, group = Method)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Gold Standard" = "red", "Coverage_R2" = "blue",
                                  "Olink" = "green", "Random" = "purple", "Quadrant" = "orange")) +
    labs(title = "Type I Error vs Effect Size", x = "Effect Size", y = "Type I Error") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 10),
      legend.position = "top"
    )
}

# Generate new line plots
plot_power_vs_effect_size(results)
plot_typeI_error_vs_effect_size(results)