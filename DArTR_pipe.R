#!/usr/bin/env Rscript

################################################################################
#### FUll DArT & SNMF Pipeline -- created 21/5/25 by JB #####
#### FINAL VERSION with Diagnostics -- 25/8/25 #####
################################################################################

################################################################################
# --- ANALYSIS CONTROL PANEL ---
################################################################################
DATASET_FILTER <- "all_UK_samples|combined_samples"  # combined_samples or all_UK_samples 
TRUNCATION_FILTER <- "t100"
PARAMS_FILTER <- "global_strict|strict2|strict1|medium|lenient|very_lenient" # global_strict|strict2|strict1|medium|lenient|very_lenient


## Either use all UK pops for corehunter, or use a specified subset (set to false and change list)
USE_ALL_POPS_FOR_COREHUNTER <- TRUE 
COREHUNTER_UK_POPS <- c("Whitewell_(WW)", "Balnaguard_Glen_(BG)", "Clashindarroch_(CD)", "")


# --- Paths and Settings ---
BASE_DIR <- "/ssd0/jbaker/stacks/final_param_results"
LOG_FILE <- file.path(BASE_DIR, "dart_pipeline_progress.log")
POP_COORDS_PATH <- "/ssd0/jbaker/stacks/Master_popmap.csv" 

POP_MAP_UK_PATH <- "/ssd0/jbaker/stacks/all_UK_popmap"
COMBINED_SAMPLES_POP_MAP_PATH <- "/ssd0/jbaker/stacks/all_sample_pop_map"
MASTER_MAP_PATH <- "/ssd0/jbaker/stacks/combined_samples/master_genotype_to_targetid_mapping.csv"

# ============================================================================
# === LIBRARY AND FUNCTION SETUP
# ============================================================================

cat("--- Verifying package installations ---\n")
required_packages <- c(
    # Data Manipulation & Core Libraries
    "dplyr", "readr", "tidyr", 
    
    # Genetics & VCF Handling
    "vcfR", "dartRverse", "poppr", "pegas", "hierfstat",
    "SNPRelate", "LEA", "plink", "HardyWeinberg",
    
    # Plotting & Visualization
    "ggplot2", "patchwork", "ragg", "scales", "RColorBrewer",
    "ggraph", "mapmixture","networkD3",
    
    # Spatial & Mapping
    "rnaturalearth", "rnaturalearthdata", "sf",
    
    # Utilities
    "devtools"
)
missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE, quietly = TRUE)]
if (length(missing_packages) > 0) stop(paste("ERROR: Missing packages:", paste(missing_packages, collapse = ", ")))

# <<< DIAGNOSTIC FEATURE >>> - Centralized logging to file and console
log_file_path <- LOG_FILE
cat("=== DArT Pipeline Log -- Run Started:", as.character(Sys.time()), "===\n", file = log_file_path, append = FALSE)
log_message <- function(..., level = "INFO") {
    msg <- sprintf(...)
    formatted_msg <- paste0("[", level, " ", Sys.time(), "] ", msg, "\n")
    cat(formatted_msg)
    cat(formatted_msg, file = log_file_path, append = TRUE)
}

create_entropy_plot <- function(ce_summary_df, best_K, run_name) {
    # Create the ggplot
    entropy_plot <- ggplot(ce_summary_df, aes(x = K, y = mean_ce)) +
        geom_line(color = "gray50", linewidth = 0.8) +
        geom_point(size = 3, color = "steelblue") +
        # Add a separate, larger red point for the best K
        geom_point(data = . %>% filter(K == best_K), color = "red", size = 5, shape = 18) +
        # <<< NEW: Add text labels for diagnostics >>>
        geom_text(aes(label = round(mean_ce, 2)), vjust = -1.5, size = 2.5, color = "black") +
        labs(
            title = paste("Mean Cross-entropy for", run_name),
            subtitle = paste("Best K (based on lowest mean score) =", best_K),
            x = "Number of Ancestral Populations (K)",
            y = "Mean Cross-entropy"
        ) +
        theme_minimal(base_size = 10) +
        # Add a bit of vertical space for the text labels
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

    return(entropy_plot)
}

# <<< FINAL PLOTTING FUNCTION  >>>
generate_snmf_map <- function(admixture_data_for_plot, pop_coords_df, map_boundary, pie_size_val, plot_title_text) {
    
    k_val <- ncol(admixture_data_for_plot) - 2 # Infer K

    # --- Data Alignment ---
    sites_in_admix <- sort(unique(admixture_data_for_plot$pop))
    coords_for_plot <- pop_coords_df %>%
        filter(Site %in% sites_in_admix) %>%
        distinct(Site, .keep_all = TRUE) %>%
        .[match(sites_in_admix, .$Site), ]

    admixture_renamed <- admixture_data_for_plot %>%
        select(pop, id, everything()) %>%
        rename(Site = pop, ID = id)

    # --- Plotting ---
    map_plot <- tryCatch({
        mapmixture(
            admixture_df = admixture_renamed,
            coords_df = coords_for_plot,
            pie_size = pie_size_val,
            boundary = map_boundary,
            cluster_cols = get_contrasting_colors(k_val),
            plot_title = plot_title_text,
            plot_title_size = 12
        )
    }, error = function(e) {
        log_message("   - ERROR: mapmixture() failed for K=%d. Error: %s", k_val, as.character(e), level = "ERROR")
        return(NULL)
    })
    
    if(!is.null(map_plot)) {
        map_plot <- map_plot + theme(
            # --- AESTHETIC TWEAKS ---
            plot.margin = unit(c(0, 0, 0, 0), "pt"), # Reduced whitespace
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.3, "lines"), # Made legend colours smaller
            legend.spacing.y = unit(0.05, "lines"),
            legend.box.margin = margin(0, 0, 0, 0, "pt")
        )
    }
    
    return(map_plot)
}


# --- HELPER FUNCTION FOR COLORS ---
get_contrasting_colors <- function(K) {
    # (This function is unchanged)
    if (is.na(K) || K < 1) return(c("grey"))
    if (K <= 8) {
        return(RColorBrewer::brewer.pal(max(3, K), "Dark2")[1:K])
    } else {
        return(scales::hue_pal()(K))
    }
}

# # --- FINAL ROBUST FUNCTION TO CALCULATE SHARED PRIVATE ALLELES (V2) ---
# calculate_shared_private_alleles <- function(genlight_obj, min_n = 5) {
    
#     log_message("         - Starting shared private allele calculation (min n = %d)...", min_n)
    
#     # 1. Get counts for all populations in the original object
#     pop_counts <- table(pop(genlight_obj))
    
#     # 2. Identify the names of populations to KEEP (those with >= min_n individuals)
#     pops_to_keep <- names(pop_counts[pop_counts >= min_n])
    
#     # 3. Identify the names of populations to DROP (those with < min_n individuals)
#     pops_to_drop <- names(pop_counts[pop_counts < min_n])
    
#     if (length(pops_to_keep) < 5) {
#         log_message("           - WARNING: Fewer than two populations with n >= %d. Skipping calculation.", min_n, level="WARN")
#         return(data.frame()) # Return an empty data frame
#     }
    
#     # 4. Use gl.drop.pop with the list of populations to drop
#     log_message("           - Filtering to %d populations with n >= %d.", length(pops_to_keep), min_n)
#     gl_filt <- gl.drop.pop(genlight_obj, pop.list = pops_to_drop, verbose = 0)
#     # --- WORKAROUND END ---
    
#     # --- Step 2: Manually create a presence/absence matrix (version-agnostic method) ---
#     geno_matrix <- as.matrix(gl_filt) # Convert genlight to a 0/1/2 matrix
    
#     pa_list <- lapply(pops_to_keep, function(pop_name) {
#         inds_in_pop <- indNames(gl_filt)[pop(gl_filt) == pop_name]
#         pop_geno_matrix <- geno_matrix[inds_in_pop, , drop = FALSE]
#         allele_present <- colSums(pop_geno_matrix, na.rm = TRUE) > 0
#         return(as.integer(allele_present))
#     })
    
#     pa_matrix <- do.call(cbind, pa_list)
#     colnames(pa_matrix) <- pops_to_keep
    
#     # --- Step 3 & 4: Loop through pairs and count alleles private to that pair ---
#     pop_pairs <- combn(pops_to_keep, 2, simplify = FALSE)
    
#     results <- lapply(pop_pairs, function(pair) {
#         pop1 <- pair[1]
#         pop2 <- pair[2]
#         other_pops <- setdiff(pops_to_keep, pair)
        
#         # An allele is "shared private" if it is present in pop1 AND pop2...
#         alleles_in_pair <- pa_matrix[, pop1] == 1 & pa_matrix[, pop2] == 1
        
#         # ...AND it is absent from all other populations being considered.
#         presence_in_others <- rowSums(pa_matrix[, other_pops, drop=FALSE])
#         alleles_absent_elsewhere <- presence_in_others == 0
        
#         shared_private_count <- sum(alleles_in_pair & alleles_absent_elsewhere)
#         data.frame(pop1 = pop1, pop2 = pop2, shared_pa_count = shared_private_count)
#     })
    
#     # Step 5: Combine results
#     final_df <- dplyr::bind_rows(results) %>%
#         dplyr::filter(shared_pa_count > 0)
        
#     log_message("         - Found %d population pairs with at least one shared private allele.", nrow(final_df))
#     return(final_df)
# }

# ==============================================================================
# --- MAIN PROCESSING LOOP --
# ==============================================================================

#    --- Phase 1: VCF Loading and Preparation ---

vcf_files_all <- list.files(path = BASE_DIR, pattern = "populations\\.snps\\.vcf$", recursive = TRUE, full.names = TRUE)
vcf_files <- vcf_files_all
if (DATASET_FILTER != "") vcf_files <- vcf_files[grepl(DATASET_FILTER, vcf_files)]
if (TRUNCATION_FILTER != "") vcf_files <- vcf_files[grepl(TRUNCATION_FILTER, vcf_files)]
if (PARAMS_FILTER != "") vcf_files <- vcf_files[grepl(PARAMS_FILTER, vcf_files)]
vcf_files <- vcf_files[!grepl("/F1/", vcf_files)]
log_message("Processing %d VCF files after filtering.", length(vcf_files))

for (i in seq_along(vcf_files)) {
  vcf_path <- vcf_files[i]
  run_name <- basename(dirname(vcf_path))
  log_message("================== STARTING RUN %d/%d: '%s' ==================", i, length(vcf_files), run_name)
   # Define is_uk_run here so it's available for all phases
  is_uk_run <- grepl("all_UK_samples", run_name)

  tryCatch({
    output_dir <- dirname(vcf_path)
    
    # --- Phase 1: VCF to Genlight ---
    log_message("--- Phase 1: Loading and Preparing Genlight Object ---")
    vcf <- read.vcfR(vcf_path)
    genlight_obj_raw <- vcfR2genlight(vcf)
    log_message("  - VCF loaded. Initial object has %d individuals and %d loci.", nInd(genlight_obj_raw), nLoc(genlight_obj_raw))
    
    initial_locus_count <- nLoc(genlight_obj_raw)
    ploidy_levels <- ploidy(genlight_obj_raw)
    if (!all(ploidy_levels == 2)) {
      n_removed <- sum(ploidy_levels != 2)
      log_message("  - Found and removed %d non-diploid individuals.", n_removed)
      genlight_obj <- new("genlight", gen = genlight_obj_raw@gen[ploidy_levels == 2], ind.names = indNames(genlight_obj_raw)[ploidy_levels == 2], loc.names = locNames(genlight_obj_raw), ploidy = as.integer(ploidy_levels[ploidy_levels == 2]))
    } else {
      log_message("  - All individuals are diploid.")
      genlight_obj <- genlight_obj_raw
    }
    genlight_obj <- dartR.base::gl.compliance.check(genlight_obj, verbose = 0)

#    --- Phase 1.5: Read Depth and Call Rate Analysis ---

log_message("--- Phase 1.5: Analyzing Read Depth ---")

# Extract the Depth ("DP") field from the VCF
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

if (!is.null(dp)) {
    # Calculate depth per individual
    ind_depth <- colMeans(dp, na.rm = TRUE)
    # Calculate depth per locus (SNP)
    loc_depth <- rowMeans(dp, na.rm = TRUE)

    # --- Create Visual Output ---
    # Plot 1: Mean read depth per individual (existing)
    p_ind_depth <- ggplot(data.frame(depth = ind_depth), aes(x = depth)) +
        geom_density(fill = "steelblue", alpha = 0.7) +
        labs(title = "Read Depth per Individual", x = "Mean Depth", y = "Density") +
        theme_minimal()

    # Plot 2: Mean read depth per locus (existing)
    p_loc_depth <- ggplot(data.frame(depth = loc_depth), aes(x = depth)) +
        geom_density(fill = "darkorange", alpha = 0.7) +
        labs(title = "Read Depth per Locus (SNP)", x = "Mean Depth", y = "Density") +
        theme_minimal()

    # Plot 3: Individual Call Rate (Coverage) ---
    log_message("     - Calculating individual call rates (coverage)...")
    gt <- extract.gt(vcf, element = "GT")
    ind_call_rate <- colMeans(!is.na(gt))
    
    p_ind_callrate <- ggplot(data.frame(call_rate = ind_call_rate), aes(x = call_rate)) +
        geom_histogram(fill = "darkgreen", alpha = 0.7, bins=30) +
        scale_x_continuous(labels = scales::percent) +
        labs(title = "Call Rate per Individual (Coverage)", x = "Proportion of Genotyped Loci", y = "Count of Individuals") +
        theme_minimal()

    # Combine all three plots and save
    depth_plot <- p_ind_depth + p_loc_depth + p_ind_callrate
    depth_plot_path <- file.path(output_dir, paste0(run_name, "_depth_distribution.png"))
    # Increased width to accommodate the third plot
    ggsave(depth_plot_path, device = ragg::agg_png, depth_plot, width = 15, height = 5, dpi = 300)
    log_message("   - Depth and coverage distribution plot saved to: %s", depth_plot_path)

    # --- Create Data Output for Final Summary Table ---
    # This creates the exact file that your final summary code looks for!
    depth_summary_df <- data.frame(
         initial_n_loci = initial_locus_count,
        mean_depth_per_ind = mean(ind_depth, na.rm = TRUE),
        median_depth_per_ind = median(ind_depth, na.rm = TRUE),
        sd_depth_per_ind = sd(ind_depth, na.rm = TRUE)
    )
    depth_summary_path <- file.path(output_dir, paste0(run_name, "_depth_summary.csv"))
    write.csv(depth_summary_df, depth_summary_path, row.names = FALSE)
    log_message("   - Depth summary metrics saved to: %s", depth_summary_path)

} else {
    log_message("   - WARNING: Could not extract 'DP' field from VCF. Skipping depth analysis.", level = "WARN")
}

# --- Phase 2: Assigning and Merging Populations (Corrected Version) ---
log_message("--- Phase 2: Assigning and Merging Populations ---")

# First, check if the individual names appear to be numeric IDs
# This regex checks if all IDs are composed only of digits.
needs_conversion <- all(grepl("^[0-9]+$", indNames(genlight_obj)))

if (needs_conversion) {
    log_message("  - Numeric IDs detected. Converting to alphanumeric IDs using the master map...")
    master_id_map <- read.csv(MASTER_MAP_PATH)
    
    current_ids <- as.numeric(indNames(genlight_obj))
    new_names <- master_id_map$genotype[match(current_ids, master_id_map$targetid)]
    
    if(any(is.na(new_names))) {
        log_message("  - WARNING: %d individuals were not found in the master ID map and will be removed.", sum(is.na(new_names)), level="WARN")
        genlight_obj <- gl.keep.ind(genlight_obj, ind.list = indNames(genlight_obj)[!is.na(new_names)], verbose=0)
        new_names <- new_names[!is.na(new_names)]
    }
    
    indNames(genlight_obj) <- new_names
    log_message("  - Renaming complete.")
} else {
    log_message("  - Alphanumeric IDs detected. Proceeding without conversion.")
}

# Now, assign populations using the (now consistent) alphanumeric IDs for ALL runs
pop_map <- read.table(COMBINED_SAMPLES_POP_MAP_PATH, header = FALSE, sep = "\t", col.names = c("id", "pop"))

# Match individual names from the genlight object to the pop map's 'id' column,
# then assign the corresponding 'pop' value directly to the genlight object's population slot.
pop(genlight_obj) <- pop_map$pop[match(indNames(genlight_obj), pop_map$id)]

log_message("  - Initial populations assigned.")

# Check for any individuals that still couldn't be mapped
if (any(is.na(pop(genlight_obj)))) {
  unmapped_ids <- indNames(genlight_obj)[is.na(pop(genlight_obj))]
    log_message("  - WARNING: %d individuals could not be assigned a population and have NA.", length(unmapped_ids), level="WARN")
    log_message("  - DEBUG: Unmapped individual ID(s): %s", paste(unmapped_ids, collapse=", "), level="WARN")
}

#    --- Phase 3: Replicate Management and QC ---
log_message("--- Phase 3: Replicate Analysis and Data Filtering ---")

# --- Step 3.1: Identify Replicate Groups ---
log_message("  - Identifying technical replicate groups...")
all_inds <- indNames(genlight_obj)
replicate_inds <- grep("_rep[0-9]*$", all_inds, value = TRUE)
tr_list <- list() # Initialize an empty list

if (length(replicate_inds) > 0) {
    base_names <- unique(sub("_rep[0-9]*$", "", replicate_inds))
    tr_list <- lapply(base_names, function(base) {
        group_pattern <- paste0("^", base, "(_rep[0-9]*)?$")
        grep(group_pattern, all_inds, value = TRUE)
    })
    tr_list <- tr_list[sapply(tr_list, length) > 1]
}

# --- Step 3.2: Create Main Dataset by Removing Extra Replicates ---
if (length(tr_list) > 0) {
    log_message("  - Found %d replicate groups. Removing extras for main analysis.", length(tr_list))
    set.seed(123)
    ids_to_keep <- sapply(tr_list, function(group) sample(group, 1))
    ids_to_remove <- setdiff(unlist(tr_list), ids_to_keep)
    gl_no_reps <- gl.drop.ind(genlight_obj, ind.list = ids_to_remove)
    log_message("  - Removed %d technical replicates.", length(ids_to_remove))
} else {
    log_message("  - No replicate groups found. Proceeding with all individuals.")
    gl_no_reps <- genlight_obj
}

# --- Step 3.3: Perform QC on Replicates AFTER Filtering (Parallel Path) ---
if (length(tr_list) > 0) {
    log_message("  - Starting parallel QC analysis for replicates...")
    
    # Isolate only the replicate individuals from the original dataset
    gl_all_replicates <- gl.keep.ind(genlight_obj, ind.list = unlist(tr_list))
    
    # Apply the EXACT SAME filtering chain to this replicate-only dataset
    log_message("      - Filtering replicate-only dataset for fair comparison...")
    gl_replicates_filtered <- gl.filter.callrate(gl_all_replicates, method="loc", threshold=0.45, verbose=0)
    gl_replicates_filtered <- gl.filter.callrate(gl_replicates_filtered, method="ind", threshold=0.4, verbose=0)
    gl_replicates_filtered <- gl.filter.callrate(gl_replicates_filtered, method="loc", threshold=0.75, verbose=0)
    gl_replicates_filtered <- gl.filter.monomorphs(gl_replicates_filtered, verbose=0)
    log_message("      - Replicates filtered: %d inds, %d loc remain for QC.", nInd(gl_replicates_filtered), nLoc(gl_replicates_filtered))

    # Calculate Mismatch Rate
    replicate_mismatch_rates <- sapply(tr_list, function(group) {
        # Check which members of the group survived the filtering
        group_present <- group[group %in% indNames(gl_replicates_filtered)]
        
        if (length(group_present) > 1) {
            gl_sub <- gl.keep.ind(gl_replicates_filtered, ind.list = group_present)
            
            # 1. Get the number of loci in this specific comparison
            num_loci <- nLoc(gl_sub)
            if (num_loci > 0) {
                # 2. Calculate the COUNT of differing loci using method = "absolute"
                dist_mat_absolute <- as.matrix(gl.dist.ind(gl_sub, method = "absolute"))
                
                # 3. Calculate the average count of mismatches between pairs
                mean_mismatch_count <- mean(dist_mat_absolute[lower.tri(dist_mat_absolute)])
                
                # 4. Convert to a proportion (mismatch rate)
                return(mean_mismatch_count / num_loci)
            }
        }
        return(NA) # Return NA if not enough individuals to compare
    })

    replicate_qc_df <- data.frame(
        group = sapply(tr_list, paste, collapse = ", "),
        mean_mismatch_rate = replicate_mismatch_rates,
        mean_mismatch_percent = replicate_mismatch_rates * 100
    )
    
    rep_qc_path <- file.path(output_dir, paste0(run_name, "_replicate_QC_analysis.csv"))
    write.csv(replicate_qc_df, rep_qc_path, row.names = FALSE)
    log_message("  - Replicate Mismatch Rate QC report saved to: %s", rep_qc_path)
}

# ---  DIAGNOSTIC ---
log_message("--- DIAGNOSTIC: Populations present BEFORE merging ---")
log_message(paste(sort(unique(pop(gl_no_reps))), collapse = " | "))
# -----------------------------

#    --- Phase 3.4: Main Dataset Filtering ---
log_message("--- Applying main filtering to the replicate-cleaned dataset ---")
log_message("  - Before main filtering: %d ind, %d loc.", nInd(gl_no_reps), nLoc(gl_no_reps))

# Merge and drop populations from the main dataset ('gl_no_reps')
gl_merged <- gl.merge.pop(gl_no_reps, old = c("Muckle_Roe_Burn_of_Cooses", "Muckle_Roe_Burn_of_Scar_Jataing", "Muckle_Roe_Cross_Burn", "Muckle_Roe_Hours_Mill"), new = "Muckle_Roe", verbose=0)
gl_merged <- gl.merge.pop(gl_merged, old = c("Ronas_Hill_Swabie_Water", "Ronas_Hill_Twa_Roes"), new = "Ronas_Hill", verbose=0)
gl_merged <-gl.merge.pop(gl_merged,old=c("Alps_Grunwald_Karlalm_Donnersbach", "Alps_Stoderzinken", "Alps_Zirbitzkogel_Judenberg"), new="NE_Alps")
gl_merged <-gl.merge.pop(gl_merged,old=c("Alps_Gerlostal_Roller", "Alps_Wiederberger_Horn"), new="W_Alps")
gl_merged <-gl.merge.pop(gl_merged,old=c("Alps_Grunwald_Karlalm","Alps_Stoderzinken", "Alps_Zirbitzkogel"), new="NE_Alps")
gl_full2<-gl.merge.pop(gl_merged,old=c("Primorye_0_Sikhote_Alin","Primorye_A_Sikhote_Alin"), new="Primorye_Sikhote_Alin")
#gl_full <- gl.drop.pop(gl_merged, pop.list = "Malham_tarn", verbose=0)
log_message("  - Population merging complete.")

# Apply call rate and monomorph filters to the main dataset
gl_full2 <- gl.filter.callrate(gl_full, method="loc", threshold=0.45, verbose=0)
gl_full2 <- gl.filter.callrate(gl_full2, method="ind", threshold=0.4, verbose=0)
gl_full2 <- gl.filter.callrate(gl_full2, method="loc", threshold=0.75, verbose=0)
gl_full2 <- gl.filter.hwe(gl_full2, subset = "each", verbose = 0) 
gl_full2 <-gl.filter.monomorphs(gl_full2, verbose = 0)

# This log message and all of Phase 4 will now use the fully filtered 'gl_full2' object
log_message("  - After filtering: %d ind, %d loc remain for SNMF.", nInd(gl_full2), nLoc(gl_full2))

log_message("--- DIAGNOSTIC: Verifying population merges before SNMF ---")
final_pops <- unique(pop(gl_full2))

# List of populations that should have been merged and should no longer exist
pops_to_verify <- c(
    "Muckle_Roe_Burn_of_Cooses", "Muckle_Roe_Burn_of_Scar_Jataing", 
    "Muckle_Roe_Cross_Burn", "Muckle_Roe_Hours_Mill",
    "Ronas_Hill_Swabie_Water", "Ronas_Hill_Twa_Roes"
)

# Check if any of these old populations remain
remnant_pops <- intersect(final_pops, pops_to_verify)

if (length(remnant_pops) > 0) {
    log_message("  - VERIFICATION FAILED. The following populations were NOT merged correctly:", level = "WARN")
    log_message("    -> %s", paste(remnant_pops, collapse = ", "), level = "WARN")
} else {
    log_message("  - VERIFICATION SUCCESS: All specified sub-populations were merged successfully.")
}

log_message("  - Final populations for SNMF analysis: %s", paste(sort(final_pops), collapse = ", "))
# --- END DIAGNOSTIC ---

# ============================================================================
# --- Phase 3.5: UK-Specific Regional Analysis (Het, AMOVA, Fst) ---
# ============================================================================
#        --- Step 3.5.1: Subsetting and Filtering UK Data ---

log_message("--- Phase 3.5: Checking for UK samples to run regional analysis on final filtered data ---")

uk_individual_names <- indNames(gl_full2)[grepl("^UK_", indNames(gl_full2))]

if (length(uk_individual_names) > 0) { 
    gl_uk_only <- gl.keep.ind(gl_full2, ind.list = uk_individual_names, verbose = 0)
    
    log_message("       - Initial UK population counts before n=1 filter:")
    log_message("%s", paste(capture.output(print(table(pop(gl_uk_only)))), collapse="\n"))
    
    pop_counts <- table(pop(gl_uk_only))
    pops_to_keep <- names(pop_counts[pop_counts > 1])
    gl_uk_analysis <- gl.keep.pop(gl_uk_only, pop.list = pops_to_keep, verbose = 0)

    if (length(pops_to_keep) < length(pop_counts)) {
        pops_removed <- names(pop_counts[pop_counts <= 1])
        log_message("       - WARNING: For UK analysis, temporarily removing %d population(s) with only one individual: %s", 
                      length(pops_removed), paste(pops_removed, collapse=", "), level="WARN")
    }

    if (nPop(gl_uk_analysis) > 1 && nLoc(gl_uk_analysis) > 1) { 
        log_message("       - Proceeding with regional analysis on %d UK populations with n > 1.", nPop(gl_uk_analysis))
        
        region_map_path <- "/ssd0/jbaker/stacks/UK_population_regions.csv"
        master_region_map <- read.csv(region_map_path)
        current_pops <- data.frame(population = popNames(gl_uk_analysis)) 
        regional_assignments <- dplyr::left_join(current_pops, master_region_map, by = "population")

        if(any(is.na(regional_assignments$region))) {
            missing_pops <- regional_assignments$population[is.na(regional_assignments$region)]
            stop(paste("FATAL ERROR: The following populations are in your data but not in your region map file:", paste(missing_pops, collapse=", ")))
        }
        log_message("       - Dynamically assigned %d UK populations to %d regions.", nrow(regional_assignments), length(unique(regional_assignments$region)))
        
        log_message("         - Preparing strata for UK regional AMOVA...")
        amova_strata_uk <- data.frame(
            ID = indNames(gl_uk_analysis),
            population = as.character(pop(gl_uk_analysis))
         ) %>%
            dplyr::left_join(regional_assignments, by = "population")
        amova_strata_uk <- amova_strata_uk[match(indNames(gl_uk_analysis), amova_strata_uk$ID), ]
        strata(gl_uk_analysis) <- data.frame(
            region = factor(amova_strata_uk$region),
            population = factor(amova_strata_uk$population)
        )
        
        log_message("         - Running AMOVA on UK populations grouped by region...")
        amova_results_uk <- tryCatch({
            poppr::poppr.amova(gl_uk_analysis, hier = ~region/population, nperm = 999)
        }, error = function(e) {
            log_message("           - ERROR: UK regional AMOVA failed. Error: %s", as.character(e), level = "ERROR")
            return(NULL)
        })

         if (!is.null(amova_results_uk)) {
            amova_summary_csv_path <- file.path(output_dir, paste0(run_name, "_UK_regional_AMOVA_summary.csv"))
            con <- file(amova_summary_csv_path, "w")
            writeLines("AMOVA Results Table (Hierarchy: ~region/population)", con)
            write.csv(amova_results_uk$results, con, row.names = TRUE)
            writeLines("", con) 
            writeLines("Phi Statistics (Percent Variation)", con)
            write.csv(amova_results_uk$statphi, con, row.names = TRUE)
            writeLines("", con)
            writeLines("Significance of Phi-Statistics (P-values from permutations)", con)
            write.csv(amova_results_uk$signif, con, row.names = TRUE)
            close(con)
            log_message("           - UK regional AMOVA comprehensive summary saved to: %s", amova_summary_csv_path)
        }
        
        log_message("         - Calculating pairwise Fst within each defined region...")
        unique_regions <- unique(regional_assignments$region)

        for (current_region in unique_regions) {
            pops_in_region <- regional_assignments$population[regional_assignments$region == current_region]
            if (length(pops_in_region) > 1) {
                log_message("                 - Processing region: %s (%d pops)", current_region, length(pops_in_region))
                gl_subset <- gl.keep.pop(gl_uk_analysis, pop.list = pops_in_region, verbose = 0)
                if (nLoc(gl_subset) > 1) {
                    fst_matrix_region <- gl.fst.pop(gl_subset, nboots = 100, percent = 95, nclusters = 8, verbose=0)
                    safe_region_name <- gsub("[^a-zA-Z0-9_]", "_", current_region)
                    fst_output_path <- file.path(output_dir, paste0(run_name, "_Fst_within_", safe_region_name, ".csv"))
                    write.csv(fst_matrix_region$Fsts, fst_output_path, row.names = TRUE)
                    log_message("                       -> Saved Fst matrix to: %s", fst_output_path)
                } else {
                    log_message("                       -> SKIPPING Fst for region %s, not enough loci.", current_region)
                }
             } else {
                log_message("                 - Skipping region: %s (Only one population)", current_region)
            }
        }
        log_message("         - Completed within-region Fst calculations.")
        
        log_message("--- Phase 3.5.5: Preparing to export UK genetic distance matrix ---")

        tryCatch({
            target_pops <- COREHUNTER_UK_POPS[COREHUNTER_UK_POPS != ""]
            if (USE_ALL_POPS_FOR_COREHUNTER) {
                log_message("           - Using ALL available UK populations for the distance matrix.")
                gl_for_dist <- gl_uk_analysis 
            } else {
                log_message("           - Filtering to %d target populations for the distance matrix.", length(target_pops))
                gl_for_dist <- gl.keep.pop(gl_uk_analysis, pop.list = target_pops, verbose = 0)
            }
            
            if (nInd(gl_for_dist) < 2) {
                log_message("           - WARNING: Fewer than two individuals available. Skipping distance matrix calculation.", level="WARN")
            } else {
                log_message("           - Calculating pairwise Euclidean distance for %d individuals...", nInd(gl_for_dist))
                diss_mat <- as.matrix(gl.dist.ind(gl_for_dist, method = "euclidean", verbose = 0))
                dist_matrix_path <- file.path(output_dir, paste0(run_name, "_UK_distance_matrix.csv"))
                write.csv(diss_mat, file = dist_matrix_path, row.names = TRUE)
                log_message("           - Distance matrix successfully saved to: %s", dist_matrix_path)
            }
         }, error = function(e) {
            log_message("--- Distance matrix export FAILED: %s", as.character(e), level = "ERROR")
        })
    } # This brace properly closes the inner if-block.
} else {
    log_message("--- Phase 3.5: Skipped. No UK individuals found in the final filtered dataset. ---")
}

# ============================================================================
# --- Phase 4: SNMF Analysis and Plotting 
# ============================================================================

log_message("--- Phase 4: SNMF Analysis and Plotting ---")

#        --- Step 4.1: Pre-SNMF Data Cleaning ---

log_message("       - Checking for individuals with 100%% missing data before SNMF...")
call_rates <- rowMeans(!is.na(as.matrix(gl_full2)))
inds_to_remove <- indNames(gl_full2)[call_rates == 0]

gl_snmf <- gl_full2 # Start with the full object for this run
if (length(inds_to_remove) > 0) {
    log_message("       - WARNING: Found %d individual(s) with 0%% call rate. Removing for SNMF analysis.", length(inds_to_remove))
    log_message("       - Removing: %s", paste(inds_to_remove, collapse = ", "), level = "WARN")
    gl_snmf <- gl.drop.ind(gl_full2, ind.list = inds_to_remove, verbose = 0)
} else {
    log_message("       - No individuals with 0%% call rate found.")
}


log_message("       - Proceeding to SNMF with a sanitized object of %d individuals and %d loci.", nInd(gl_snmf), nLoc(gl_snmf))
#        --- Step 4.2: Run or Load SNMF Project ---
geno_path <- file.path(output_dir, paste0(run_name, ".geno"))
snmf_project_path <- file.path(output_dir, paste0(run_name, ".snmfProject"))

# This logic checks if an SNMF project already exists. If so, it loads it.
# If not, it runs the analysis from scratch.
if (file.exists(snmf_project_path)) {
    log_message("       - Found existing SNMF project. Loading from: %s", snmf_project_path)
    struc <- LEA::load.snmfProject(snmf_project_path)
} else {
    log_message("       - No SNMF project found. Running new analysis (can be slow)...")
    set.seed(12345)
    gl2geno(gl_snmf, outfile = run_name, outpath = output_dir)
    struc <- snmf(geno_path, K = 1:30, repetitions = 10, iterations = 10000, entropy = TRUE, ploidy = 2, project = "new", alpha = 50, CPU = 8)
}

# --- CRITICAL STEP: Always save the list of individuals after loading or running ---
# This ensures the standalone plotting script always has the file it needs.
log_message("       - Saving list of individuals for external plotting script.")
snmf_inds_path <- file.path(output_dir, paste0(run_name, "_snmf_ind_names.txt"))
write.table(indNames(gl_snmf), file = snmf_inds_path, row.names = FALSE, col.names = FALSE, quote = FALSE)


# log_message("       - Forcing new SNMF analysis (can be slow)...")
#     set.seed(12345)
#     # CRITICAL FIX: Ensure the cleaned gl_snmf object is used here.
#     log_message("       - Saving list of individuals for external plotting script.")
#     snmf_inds_path <- file.path(output_dir, paste0(run_name, "_snmf_ind_names.txt"))
#     write.table(indNames(gl_snmf), file = snmf_inds_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
#     gl2geno(gl_snmf, outfile = run_name, outpath = output_dir)
#     struc <- snmf(geno_path, K = 1:30, repetitions = 10, iterations = 10000, entropy = TRUE, ploidy = 2, project = "new", alpha = 50, CPU = 4)

#        --- Step 4.3: Determine Best K ---
log_message("       - Calculating cross-entropy and determining best K...")
all_ce_matrix <- sapply(2:30, function(k) cross.entropy(struc, K = k))
colnames(all_ce_matrix) <- 2:30
ce_summary <- as.data.frame(all_ce_matrix) %>%
    tidyr::pivot_longer(everything(), names_to = "K", values_to = "CrossEntropy") %>%
    mutate(K = as.numeric(K)) %>%
    group_by(K) %>%
    summarise(mean_ce = mean(CrossEntropy, na.rm = TRUE), .groups = "drop")

ce_summary_valid <- ce_summary %>%
    filter(!is.na(mean_ce) & !is.nan(mean_ce) & is.finite(mean_ce))

best_K <- NA 
if (nrow(ce_summary_valid) > 0) {
    best_K <- ce_summary_valid$K[which.min(ce_summary_valid$mean_ce)]
}

if (length(best_K) == 0 || is.na(best_K)) {
    log_message("       - FATAL FOR THIS RUN: Could not determine a valid best K from SNMF results.", level="ERROR")
    next 
}
log_message("       - SNMF analysis complete. Best K identified as %d.", best_K)
k_results_df <- data.frame(
    run = run_name,
    K = best_K,
    cross_entropy = ce_summary_valid$mean_ce[ce_summary_valid$K == best_K],
    n_loci_final = nLoc(gl_snmf) # Use the object that went into SNMF
)
kres_file_path <- file.path(output_dir, paste0(run_name, "_best_K_results.csv"))
write.csv(k_results_df, kres_file_path, row.names = FALSE)
log_message("         - Best K results saved for summary table: %s", kres_file_path)

#        --- Step 4.4: Generate SNMF Admixture Plots ---
ind_metrics_final <- data.frame(id = indNames(gl_snmf), pop = as.character(pop(gl_snmf)), stringsAsFactors = FALSE)
pop_coords_master <- read.csv(POP_COORDS_PATH)

q_matrix_check <- as.data.frame(Q(struc, K = best_K, run = which.min(cross.entropy(struc, K = best_K))))

# This check will now only fail if the user forgets to delete a stale project file.
if(nrow(ind_metrics_final) != nrow(q_matrix_check)){
    stop(sprintf("FATAL: Mismatch detected (%d vs %d). Delete the project file and retry: %s",
                 nrow(ind_metrics_final), nrow(q_matrix_check), snmf_project_path))
}

k_values_to_plot <- na.omit(unique(c(best_K - 1, best_K, best_K + 1)))
k_values_to_plot <- k_values_to_plot[k_values_to_plot >= 2]
all_final_plots <- list()


for (k in k_values_to_plot) {
    plot_title <- paste("K =", k, if (k == best_K) "(Best)" else "")
    q_matrix_k <- as.data.frame(Q(struc, K = k, run = which.min(cross.entropy(struc, K = k))))
    
    # This will now always have the same number of rows
    admixture_full_k <- dplyr::bind_cols(ind_metrics_final, q_matrix_k)
    
    
    # --- A. Generate the MAIN map ---
    if (!is_uk_run) { # Global run
         main_boundary <- c(xmin = -170, xmax = 170, ymin = 30, ymax = 72) 
        main_map <- generate_snmf_map(admixture_full_k, pop_coords_master, main_boundary, pie_size_val = 4, plot_title)
    } else { # UK-only run
        main_boundary <- c(xmin = -12, xmax = 2, ymin = 50, ymax = 62)
        main_map <- generate_snmf_map(admixture_full_k, pop_coords_master, main_boundary, pie_size_val = 2, plot_title)
    }
    
    if (is.null(main_map)) next
    
    # --- B. For combined runs, generate and combine the UK INSET ---
    if (!is_uk_run) {
        uk_boundary <- c(xmin = -12, xmax = 2, ymin = 50, ymax = 62)
        uk_pops_coords <- pop_coords_master %>% filter(Longitude >= uk_boundary["xmin"] & Longitude <= uk_boundary["xmax"], Latitude >= uk_boundary["ymin"] & Latitude <= uk_boundary["ymax"])
        
        # <<< Filter the already-combined data frame >>>
        admixture_uk_k <- admixture_full_k %>% filter(pop %in% uk_pops_coords$Site)
        
        if(nrow(admixture_uk_k) > 0) {
            uk_inset_map <- generate_snmf_map(admixture_uk_k, uk_pops_coords, uk_boundary, pie_size_val = .75, "")
            
            if (!is.null(uk_inset_map)) {
                uk_inset_map <- uk_inset_map +
                    theme_void() +
                    theme(legend.position = "none",
                          plot.title = element_blank(),
                          panel.background = element_rect(fill = "lightblue", color = "black", linewidth=0.5),
                          plot.margin = unit(c(0, 0, 0, 0), "lines"))
                
                # Combine main map and inset side-by-side
                    combined_plot_for_k <- main_map + uk_inset_map + plot_layout(widths = c(3, 1))
                all_final_plots[[as.character(k)]] <- combined_plot_for_k
            } else {
                all_final_plots[[as.character(k)]] <- main_map
            }
        } else {
            all_final_plots[[as.character(k)]] <- main_map
        }
    } else {
        all_final_plots[[as.character(k)]] <- main_map
    }
}

# --- C. Assemble Final Plot Panel ---
if (length(all_final_plots) > 0) {
    log_message("    - Assembling final composite plot...")
    entropy_plot <- create_entropy_plot(ce_summary, best_K, run_name)

    # Check if this is a UK-only run to adjust the layout for better legibility
    if (is_uk_run) {
        log_message("    - Using horizontal layout for UK-only maps.")
        
        # Arrange the K maps side-by-side in a single row
        map_panel <- wrap_plots(all_final_plots, nrow = 1)
        
        # Then, stack the row of maps on top of the entropy plot
        combined_output <- wrap_plots(map_panel, entropy_plot, ncol = 1, heights = c(1, 0.4))
        
        # Adjust image dimensions for this wider layout. 
        # Width is proportional to the number of maps, height is relatively fixed.
        final_width <- 5.5 * length(all_final_plots)
        final_height <- 8 
        
    } else {
        log_message("    - Using vertical layout for combined/global maps.")
        
        # Original logic: Create a single list containing all plots
        plots_to_assemble <- c(all_final_plots, list(entropy_plot))
        
        # Define relative heights for the single vertical column
        relative_heights <- c(rep(1, length(all_final_plots)), 0.5)
        
        # Stack everything into one column
        combined_output <- wrap_plots(plots_to_assemble, ncol = 1, heights = relative_heights)
        
        # Use original image dimensions
        final_width <- 12
        final_height <- 3 * sum(relative_heights)
    }

    output_png_path <- file.path(output_dir, paste0(run_name, "_SNMF_Composite_Mapmixture.png"))
    
    # Save the plot using the dynamically determined dimensions
    ggsave(
        output_png_path, 
        device = ragg::agg_png, 
        combined_output, 
        width = final_width, 
        height = final_height, 
        dpi = 300, 
        limitsize = FALSE
    )
    log_message("--- SNMF composite plot saved to: %s ---", output_png_path)

} else {
    log_message("WARNING: Skipping SNMF plot generation due to lack of valid maps.", level = "WARN")
}

#    --- Phase 4.5: AMOVA on SNMF Clusters ---

log_message("--- Phase 4.5: Running AMOVA on SNMF Clusters ---")
if (exists("struc") && exists("best_K") && best_K > 1 && !is.na(best_K)) {
    q_matrix <- as.data.frame(Q(struc, K = best_K, run = which.min(cross.entropy(struc, K = best_K))))
    colnames(q_matrix) <- paste0("K", 1:best_K)
    
    # CRITICAL FIX: Use the synchronized ind_metrics_final and gl_snmf objects.
    ind_data <- dplyr::bind_cols(ind_metrics_final, q_matrix)
    pop_to_cluster_map <- ind_data %>%
        dplyr::select(-id) %>%
        dplyr::group_by(pop) %>%
        dplyr::summarise(across(everything(), mean), .groups = "drop") %>%
        tidyr::pivot_longer(!pop, names_to = "snmf_cluster", values_to = "mean_q") %>%
        dplyr::group_by(pop) %>%
        dplyr::slice_max(order_by = mean_q, n = 1, with_ties = FALSE) %>%
        dplyr::select(pop, snmf_cluster)
        
    log_message("       - Assigned populations to predominant SNMF clusters.")
    
    amova_strata_snmf <- data.frame(ID = indNames(gl_snmf), Population = as.character(pop(gl_snmf))) %>%
        dplyr::left_join(pop_to_cluster_map, by = c("Population" = "pop"))
    amova_strata_snmf <- amova_strata_snmf[match(indNames(gl_snmf), amova_strata_snmf$ID), ]
    
    strata(gl_snmf) <- data.frame(
        snmf_cluster = factor(amova_strata_snmf$snmf_cluster),
        Population = factor(amova_strata_snmf$Population)
    )
    
    amova_results_snmf <- tryCatch({
        poppr::poppr.amova(gl_snmf, hier = ~snmf_cluster/Population, nperm = 999)
    }, error = function(e) {
        log_message("       - ERROR: SNMF-cluster AMOVA failed: %s", as.character(e), level = "ERROR")
        return(NULL)
    })

## save amovas
if (!is.null(amova_results_snmf)) {
    # Define a single, comprehensive CSV output path
    amova_summary_csv_path <- file.path(output_dir, paste0(run_name, "_SNMF_K", best_K, "_cluster_AMOVA_summary.csv"))

    # Use a file connection to write different sections
    con <- file(amova_summary_csv_path, "w")

    # Section 1: Main Results Table
    writeLines("AMOVA Results Table", con)
    # The poppr object stores results in a data.frame called 'results'
    write.csv(amova_results_snmf$results, con, row.names = TRUE)
    
    # Add a blank line for spacing
    writeLines("", con) 

    # Section 2: Phi Statistics (Variation components)
    writeLines("Phi Statistics (Percent Variation)", con)
    # This is stored in the 'statphi' table
    write.csv(amova_results_snmf$statphi, con, row.names = TRUE)

    # Add a blank line for spacing
    writeLines("", con)

    # Section 3: Significance Tests (from permutations)
    writeLines("Significance of Phi-Statistics (P-values from permutations)", con)
    # This is stored in the 'signif' table
    write.csv(amova_results_snmf$signif, con, row.names = TRUE)

    close(con) # Close the file connection
    
    log_message("         - SNMF cluster AMOVA comprehensive summary saved to: %s", amova_summary_csv_path)

    } else {
        log_message("    - SNMF results not available or Best K <= 1. Skipping AMOVA on clusters.")
    }
}

#    --- Phase 5: Geographic and Network Diversity Analysis ---
log_message("--- Phase 5: Calculating and Plotting Geographic Diversity Metrics ---")

# --- Step 1: Filter data for this Phase ---
log_message("       - Checking all populations for n=1 individuals...")
gl_for_het <- gl_full2
pop_counts_full <- table(pop(gl_for_het))
log_message("%s", paste(capture.output(print(pop_counts_full)), collapse="\n"))

if(any(pop_counts_full <= 1)) {
    pops_to_keep_full <- names(pop_counts_full[pop_counts_full > 1])
    pops_removed_full <- names(pop_counts_full[pop_counts_full <= 1])
    log_message("       - WARNING: Temporarily removing %d population(s) with n=1 for diversity calculations: %s",
                length(pops_removed_full), paste(pops_removed_full, collapse=", "), level="WARN")
    gl_for_het <- gl.keep.pop(gl_for_het, pop.list = pops_to_keep_full, verbose=0)
}
log_message("       - Proceeding with diversity calculations on %d populations (n>1).", nPop(gl_for_het))

# --- Step 2: Calculate Diversity Metrics ---

# A. Observed Heterozygosity (Ho) using the cleaned gl_for_het
log_message("       - Calculating Observed Heterozygosity (Ho)...")
het_stats <- gl.report.heterozygosity(gl_for_het, verbose = 0, nboots = 0, ncpus = 5)
pop_het <- het_stats %>%
    tibble::rownames_to_column("Site") %>%
    dplyr::select(Site, Ho)
het_file_path <- file.path(output_dir, paste0(run_name, "_heterozygosity_results.csv"))
write.csv(het_stats, het_file_path, row.names = TRUE) # Save the full stats table
log_message("         - Heterozygosity results saved for summary table: %s", het_file_path)

# B. Mean Pairwise Fst using the cleaned gl_for_het
log_message("       - Calculating Mean Pairwise Fst...")
fst_results_list <- gl.fst.pop(gl_for_het, nboots = 100, percent = 95, nclusters = 1, verbose = 0)
fst_matrix <- fst_results_list$Fsts
diag(fst_matrix) <- NA 
mean_fst <- data.frame(
    Site = rownames(fst_matrix),
    mean_Fst = rowMeans(fst_matrix, na.rm = TRUE)
)

# C. Private Alleles using the cleaned gl_for_het :

### METHOD FOR CALCULATING PRIVATE ALLELES EACH RUN 

log_message("         - Calculating Private Alleles...")
# The gl.report.pa function from dartR calculates private alleles per population
pa_results <- gl.report.pa(gl_for_het, verbose = 0)

# Save the results to the CSV file that the final summary code looks for
pa_file_path <- file.path(output_dir, paste0(run_name, "_private_alleles.csv"))
write.csv(pa_results, pa_file_path, row.names = FALSE)
log_message("         - Private allele counts saved for summary table: %s", pa_file_path)

# #### METHOD FOR SAVING AND LOADING PRIVATE ALLELE DATA TO DEBUG PLOTTING: 
# # C. Private Alleles using the cleaned gl_for_het
# pa_file_path <- file.path(output_dir, paste0(run_name, "_private_alleles.csv"))

# # Check if the private allele results file already exists to avoid re-calculation
# if (file.exists(pa_file_path)) {
#     log_message("         - Found existing private allele file. Loading from: %s", pa_file_path)
#     pa_results <- read.csv(pa_file_path)
# } else {
#     log_message("         - Calculating Private Alleles (this may take a while)...")
#     # The gl.report.pa function from dartR calculates private alleles per population 
#     pa_results <- gl.report.pa(gl_for_het, verbose = 0) [cite: 92]

#     # Save the results to the CSV file for future runs and the summary code 
#     write.csv(pa_results, pa_file_path, row.names = FALSE) [cite: 92]
#     log_message("         - Private allele counts calculated and saved to: %s", pa_file_path)
# }
# ##########

# # --- Step 3: Consolidate all metrics with coordinates ---
# log_message("       - Consolidating metrics into a data frame for plotting...")
# final_pops_for_plotting <- popNames(gl_for_het)
# plot_data <- pop_coords_master %>%
#     dplyr::filter(Site %in% final_pops_for_plotting) %>%
#     dplyr::left_join(pop_het, by = "Site") %>%
#     dplyr::left_join(mean_fst, by = "Site")


# log_message("DIAGNOSTIC: Number of populations with n > 1 available for diversity stats: %d", nPop(gl_for_het))

# # --- Step 4: Create and Assemble the Plots (Final Definitive Version) ---
# log_message("         - Generating final geographic and network diversity plots...")

# # Initialize the list that will hold all our final plots
# plot_list <- list()
# world_plot <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# # Define the master map theme for all plots
# master_map_theme <- theme_minimal() +
#     theme(
#         panel.background = element_rect(fill = "lightblue", color = NA),
#         panel.grid.major = element_blank(), # Removes dotted lines
#         panel.grid.minor = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank()
#     )

# # --- A. Shared Private Allele Network Plot (Geographic Version) ---
# shared_pa_df <- calculate_shared_private_alleles(gl_for_het, min_n = 5)
# if (!is.null(shared_pa_df) && nrow(shared_pa_df) > 0) {
#     # (Robustness checks for coordinates remain the same)
#     shared_pa_df_clean <- shared_pa_df %>% dplyr::filter(!is.na(pop1) & !is.na(pop2))
#     node_names <- unique(c(shared_pa_df_clean$pop1, shared_pa_df_clean$pop2))
#     node_coords <- pop_coords_master %>% dplyr::filter(Site %in% node_names)
#     missing_pops <- setdiff(node_names, node_coords$Site)
#     if(length(missing_pops) > 0){ log_message("           - WARNING: Pops missing from master_popmap.csv for network plot: %s", paste(missing_pops, collapse=", "), level="WARN") }
#     node_coords_clean <- node_coords %>% dplyr::filter(!is.na(Longitude) & !is.na(Latitude))
#     nodes_with_valid_coords <- node_coords_clean$Site
#     shared_pa_df_filtered <- shared_pa_df_clean %>% dplyr::filter(pop1 %in% nodes_with_valid_coords & pop2 %in% nodes_with_valid_coords)

#     if (nrow(shared_pa_df_filtered) > 0) {
#         edges_for_plot <- shared_pa_df_filtered %>%
#             left_join(node_coords_clean, by = c("pop1" = "Site")) %>% rename(x = Longitude, y = Latitude) %>%
#             left_join(node_coords_clean, by = c("pop2" = "Site")) %>% rename(xend = Longitude, yend = Latitude)

#         p_shared_pa <- ggplot(data = world_plot) +
#             geom_sf(fill = "gray80", color = "white") +
#             geom_segment(data = edges_for_plot, aes(x = x, y = y, xend = xend, yend = yend, linewidth = shared_pa_count), color = "gray40", alpha=0.6) +
#             geom_point(data = node_coords_clean, aes(x = Longitude, y = Latitude), size = 3, color = "skyblue") +
#             scale_linewidth_continuous(range = c(0.2, 1.5), name = "Shared Alleles") +
#             labs(title = "Shared Private Alleles (Pops n≥5)") +
#             master_map_theme # Apply the consistent theme
#         plot_list[['shared_pa']] <- p_shared_pa
#     }
# }

# # --- B. Geographic Diversity Plots (Ho and Fst) ---
# if (any(!is.na(plot_data$Ho))) {
#     p_ho <- ggplot(data = world_plot) +
#         geom_sf(fill = "gray80", color = "white") +
#         geom_point(data = plot_data, aes(x = Longitude, y = Latitude, size = Ho), color = "darkred", alpha = 0.7) +
#         scale_size_continuous(range = c(1, 8), name = "Obs. Het (Ho)") +
#         labs(title = "Observed Heterozygosity") +
#         master_map_theme
#     plot_list[['ho']] <- p_ho
# }
# if (any(!is.na(plot_data$mean_Fst))) {
#     p_fst <- ggplot(data = world_plot) +
#         geom_sf(fill = "gray80", color = "white") +
#         geom_point(data = plot_data, aes(x = Longitude, y = Latitude, size = mean_Fst), color = "darkblue", alpha = 0.7) +
#         scale_size_continuous(range = c(1, 8), name = "Mean Fst") +
#         labs(title = "Mean Pairwise Fst (Isolation)") +
#         master_map_theme
#     plot_list[['fst']] <- p_fst
# }

# # --- C. Add UK Insets and Final Projections to ALL plots ---
# if (length(plot_list) > 0) {
#     if(is_uk_run) {
#         uk_boundary <- c(xmin = -12, xmax = 2, ymin = 50, ymax = 62)
#         plot_list <- lapply(plot_list, function(p) p + coord_sf(xlim = uk_boundary[c("xmin", "xmax")], ylim = uk_boundary[c("ymin", "ymax")], expand = FALSE))
#     } else {
#         global_boundary <- c(xmin = -170, xmax = 170, ymin = 30, ymax = 72)
        
#         # Create the UK inset map once
#         uk_boundary <- c(xmin = -12, xmax = 2, ymin = 50, ymax = 62)
#         uk_inset_base <- ggplot(data = world_plot) +
#             geom_sf(fill = "gray80", color = "white") +
#             coord_sf(xlim = uk_boundary[c("xmin", "xmax")], ylim = uk_boundary[c("ymin", "ymax")], expand = FALSE) +
#             theme_void() +
#             theme(panel.background = element_rect(fill = "lightblue", color = "black"))

#         # Add the inset to each main plot before setting final coordinates
#         plot_list <- lapply(plot_list, function(p) {
#             p + annotation_custom(
#                     grob = ggplotGrob(uk_inset_base),
#                     xmin = 100, xmax = 170, ymin = 55, ymax = 72
#                 )
#         })

#         # Apply the final projection to all plots
#         plot_list <- lapply(plot_list, function(p) p + coord_sf(xlim = global_boundary[c("xmin", "xmax")], ylim = global_boundary[c("ymin", "ymax")], expand = FALSE))
#     }
# }

# # --- D. Assemble the final combined plot panel ---
# if (length(plot_list) > 0) {
#     diversity_plot_panel <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
#     diversity_plot_path <- file.path(output_dir, paste0(run_name, "_Diversity_Maps.png"))
#     ggsave(diversity_plot_path, device = ragg::agg_png, diversity_plot_panel, width = 7 * length(plot_list), height = 6, dpi = 300, limitsize=FALSE)
#     log_message("--- Geographic and network diversity plots saved to: %s ---", diversity_plot_path)
# } else {
#     log_message("--- WARNING: No diversity plots were generated as no valid data was found. ---", level = "WARN")
# }

}, error = function(e) {
    # This block runs if ANY error occurs for the current VCF file
    log_message("============ RUN FAILED ==========", level = "ERROR")
    log_message("An error occurred while processing run: '%s'", run_name, level = "ERROR")
    log_message("Error message: %s", as.character(e), level = "ERROR")
    # The loop will now automatically continue to the next VCF file
}) # This closes the tryCatch block
}


# ==============================================================================
# === FINAL SUMMARY TABLE CONSTRUCTION
# ==============================================================================

cat("\n=== Building and updating final summary table for all runs ===\n")

final_summary_path <- file.path(BASE_DIR, "all_runs_summary.csv")
existing_runs <- c()
old_summary_df <- data.frame()

if (file.exists(final_summary_path)) {
    cat("--- Found existing summary file. Reading for update ---\n")
    old_summary_df <- readr::read_csv(final_summary_path, show_col_types = FALSE)
    if ("run" %in% names(old_summary_df)) {
        existing_runs <- old_summary_df$run
    }
}

# Find all potential run directories by looking for depth summary files.
all_summary_files <- list.files(BASE_DIR, pattern = "_depth_summary.csv$", recursive = TRUE, full.names = TRUE)

# --- FIX 1: Apply the same filters from the control panel to the file list ---
if (DATASET_FILTER != "") all_summary_files <- all_summary_files[grepl(DATASET_FILTER, all_summary_files)]
if (TRUNCATION_FILTER != "") all_summary_files <- all_summary_files[grepl(TRUNCATION_FILTER, all_summary_files)]
if (PARAMS_FILTER != "") all_summary_files <- all_summary_files[grepl(PARAMS_FILTER, all_summary_files)]
all_summary_files <- all_summary_files[!grepl("/F1/", all_summary_files)]

run_names_found <- sub("_depth_summary.csv$", "", basename(all_summary_files))
new_runs_mask <- !(run_names_found %in% existing_runs)
summary_files_to_process <- all_summary_files[new_runs_mask]

if (length(summary_files_to_process) > 0) {
    cat(sprintf("--- Found %d new run(s) to add to the summary ---\n", length(summary_files_to_process)))

    new_summary_list <- lapply(summary_files_to_process, function(depth_file) {
        run_name <- sub("_depth_summary.csv$", "", basename(depth_file))
        dir <- dirname(depth_file)
        
        # Helper function to safely read a single value from a CSV
        safe_read_val <- function(file_path, col_name, row_num = 1) {
            if (!file.exists(file_path)) return(NA)
            tryCatch({
                df <- readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE)
                if(col_name %in% names(df)){ return(df[[col_name]][row_num]) } else { return(NA) }
            }, error = function(e) NA)
        }
        
        # Helper function to safely calculate the mean of a column
        safe_calc_mean <- function(file_path, col_name) {
            if (!file.exists(file_path)) return(NA)
            tryCatch({
                df <- readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE)
                if(col_name %in% names(df)){ return(mean(df[[col_name]], na.rm = TRUE)) } else { return(NA) }
            }, error = function(e) NA)
        }
        
        # --- FIX 2: Point to the files that are NOW being created ---
        
        # Read Depth (this was already correct)
        depth_summary <- readr::read_csv(depth_file, show_col_types = FALSE, progress = FALSE)
        mean_depth <- depth_summary$mean_depth_per_ind[1]
        median_depth <- depth_summary$median_depth_per_ind[1]
        sd_depth <- depth_summary$sd_depth_per_ind[1]

        # Read initial locus count from depth summary
        initial_loci <- safe_read_val(depth_file, "initial_n_loci")

        # Read Best K, Cross-Entropy, and Final Locus Count from the new file
        kres_file <- file.path(dir, paste0(run_name, "_best_K_results.csv"))
        best_K <- safe_read_val(kres_file, "K")
        ce <- safe_read_val(kres_file, "cross_entropy")
        n_loci_final <- safe_read_val(kres_file, "n_loci_final")

        # Read Heterozygosity from the new file
        het_file <- file.path(dir, paste0(run_name, "_heterozygosity_results.csv"))
        mean_Ho <- safe_calc_mean(het_file, "Ho")
        mean_He <- safe_calc_mean(het_file, "He") # Hs is Expected Het.

        # Read Replicate Mismatch Rate (this was already correct)
        rep_file <- file.path(dir, paste0(run_name, "_replicate_QC_analysis.csv"))
        mean_rep_dist <- safe_calc_mean(rep_file, "mean_mismatch_rate")

        # NOTE: Private allele calculation is not in the script. Returning NA.
        pa_file <- file.path(dir, paste0(run_name, "_private_alleles.csv"))
                priv_total <- NA
                if (file.exists(pa_file)) {
                    tryCatch({
                        pa_df <- readr::read_csv(pa_file, show_col_types = FALSE, progress = FALSE)
                        # The output from gl.report.pa has a column named 'pa' for the count
                        if("pa" %in% names(pa_df)) {
                            priv_total <- sum(pa_df$pa, na.rm = TRUE)
                        }
                    }, error = function(e) { priv_total <- NA })
                }
                
        data.frame(
            run = run_name,
            n_loci_initial = initial_loci,
            n_loci_final = n_loci_final,
            mean_depth = mean_depth,
            median_depth = median_depth,
            sd_depth = sd_depth,
            mean_Ho = mean_Ho,
            mean_He = mean_He,
            best_K = best_K,
            cross_ent = ce,
            rep_dist = mean_rep_dist,
            priv_alleles = priv_total,
            stringsAsFactors = FALSE
        )
    })
    
    new_summary_df <- dplyr::bind_rows(new_summary_list)
    final_summary_df <- dplyr::bind_rows(old_summary_df, new_summary_df)
    
    readr::write_csv(final_summary_df, final_summary_path)
    cat(sprintf("\n=== Summary updated. Total runs summarized: %d ===\n", nrow(final_summary_df)))

} else {
    cat("\n=== No new runs found. Summary file is already up-to-date ===\n")
}