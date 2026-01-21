# Loading packages 

library(dplyr)
library(readr)

set.seed(42)

# -------------------
# Settings
# -------------------

# Study window
xmin <- -10; xmax <- 45
ymin <-  30; ymax <- 60

# Slices
slice_baseline <- "slice_80000_70000"   # baseline slice
slice_bins_ref <- slice_baseline        # define longitude bin edges from baseline

# Targets (post-baseline)
slice_targets  <- c("slice_70000_60000", "slice_60000_50000", "slice_50000_40000")

# Rarefaction reps
nsim <- 1000

# Longitude bin width (degrees)
bin_width_deg <- 2

# Output folder
out_dir <- "3_output/road_analysis"

# Input (single combined CSV)
timeslices_file <- "3_output/road_analysis/time_slices/localities_all_slices.csv"


# -----------------------------
# Load all time-slice localities from single CSV
# -----------------------------

if (!file.exists(timeslices_file)) {
  stop("Cannot find localities_all_slices.csv at: ", timeslices_file)
}

all_loc <- readr::read_csv(timeslices_file, show_col_types = FALSE)

required_cols <- c("slice", "Locality", "X", "Y")
missing_cols <- setdiff(required_cols, names(all_loc))
if (length(missing_cols) > 0) {
  stop("Missing required columns in localities_all_slices.csv: ",
       paste(missing_cols, collapse = ", "))
}

# Build `timeslices` list (names are slice IDs like "slice_80000_70000")
timeslices <- split(all_loc, all_loc$slice)

# Ensure required slices exist
required_slices <- unique(c(slice_bins_ref, slice_baseline, slice_targets))
missing_slices <- setdiff(required_slices, names(timeslices))
if (length(missing_slices) > 0) {
  stop("Missing required slices in localities_all_slices.csv: ",
       paste(missing_slices, collapse = ", "),
       "\nAvailable slices: ", paste(sort(names(timeslices)), collapse = ", "))
}



# --------------------------------------------
# Helper: locality-level slice extractor
# ---------------------------------------------
get_slice_localities <- function(sname) {
  df <- timeslices[[sname]]
  if (is.null(df)) stop("Slice not found in timeslices: ", sname)
  
  df %>%
    distinct(Locality, .keep_all = TRUE) %>%
    transmute(
      Locality = Locality,
      lon = as.numeric(X),
      lat = as.numeric(Y),
      Country = if ("Country" %in% names(df)) Country else NA_character_
    ) %>%
    filter(
      !is.na(lon) & !is.na(lat),
      lon >= xmin & lon <= xmax,
      lat >= ymin & lat <= ymax
    )
}

# -----------------------------
# Load slice data
# -----------------------------
ref_df  <- get_slice_localities(slice_bins_ref)
base_df <- get_slice_localities(slice_baseline)

target_df_list <- lapply(slice_targets, get_slice_localities)
names(target_df_list) <- slice_targets

# Baseline N
N_base <- nrow(base_df)
if (N_base < 3) stop("Baseline slice has too few localities after clipping (N = ", N_base, ").")


# Diagnostic
all_lons_clipped <- unlist(lapply(c(slice_bins_ref, slice_baseline, slice_targets), function(s) {
  get_slice_localities(s)$lon
}))


# ==================================================
# (A) Longitude-bin coverage rarefaction
# ===================================================

# Define longitude bin edges using the REF slice (here: baseline 80–70 ka)
lon_min <- floor(min(ref_df$lon, na.rm = TRUE) / bin_width_deg) * bin_width_deg
lon_max <- ceiling(max(ref_df$lon, na.rm = TRUE) / bin_width_deg) * bin_width_deg
bin_edges <- seq(lon_min, lon_max, by = bin_width_deg)

assign_bins <- function(df, edges) {
  df %>%
    mutate(
      lon_bin = cut(lon, breaks = edges, include.lowest = TRUE, right = FALSE)
    ) %>%
    filter(!is.na(lon_bin))
}

ref_df  <- assign_bins(ref_df,  bin_edges)
base_df <- assign_bins(base_df, bin_edges)
target_df_list <- lapply(target_df_list, assign_bins, edges = bin_edges)

# Baseline occupied bins
base_bins <- sort(unique(base_df$lon_bin))
K_base <- length(base_bins)
all_bin_levels <- levels(base_df$lon_bin)


occupied_bins_count <- function(df) length(unique(df$lon_bin))

new_bins_vs_base <- function(df, base_bins) {
  b <- unique(df$lon_bin)
  sum(!(b %in% base_bins))
}

tv_distance <- function(dfA, dfB, all_bins_levels) {
  a <- table(factor(dfA$lon_bin, levels = all_bins_levels))
  b <- table(factor(dfB$lon_bin, levels = all_bins_levels))
  pa <- as.numeric(a) / sum(a)
  pb <- as.numeric(b) / sum(b)
  0.5 * sum(abs(pa - pb))
}

bin_sim_results <- list()
bin_ci_results  <- list()

for (sname in slice_targets) {
  
  dfT <- target_df_list[[sname]]
  N_T <- nrow(dfT)
  
  if (N_T < N_base) {
    warning("Skipping ", sname, ": target N (", N_T, ") < baseline N (", N_base, ") after clipping.")
    next
  }
  
  # Observed (full, not downsampled) for reference
  obs_bins <- occupied_bins_count(dfT)
  obs_new  <- new_bins_vs_base(dfT, base_bins)
  obs_tv   <- tv_distance(dfT, base_df, all_bin_levels)
  
  sim_occ_bins <- numeric(nsim)
  sim_new_bins <- numeric(nsim)
  sim_tv       <- numeric(nsim)
  
  nbins <- length(all_bin_levels)
  sim_bin_counts <- matrix(0L, nrow = nsim, ncol = nbins)
  colnames(sim_bin_counts) <- all_bin_levels
  
  for (i in seq_len(nsim)) {
    idx <- sample.int(N_T, size = N_base, replace = FALSE)
    samp <- dfT[idx, , drop = FALSE]
    
    sim_occ_bins[i] <- occupied_bins_count(samp)
    sim_new_bins[i] <- new_bins_vs_base(samp, base_bins)
    sim_tv[i]       <- tv_distance(samp, base_df, all_bin_levels)
    
    sim_bin_counts[i, ] <- as.integer(table(factor(samp$lon_bin, levels = all_bin_levels)))
  }
  
  # One-sided Monte Carlo p-values
  p_occ_more <- (1 + sum(sim_occ_bins <= K_base)) / (nsim + 1)
  p_new_more <- (1 + sum(sim_new_bins <= 0))     / (nsim + 1)
  p_tv_more  <- (1 + sum(sim_tv <= 0))           / (nsim + 1)
  
  bin_sim_results[[sname]] <- tibble(
    slice_target = sname,
    N_target_full = N_T,
    N_rarefied_to = N_base,
    
    bin_width_deg = bin_width_deg,
    baseline_bins_occupied = K_base,
    
    observed_full_bins_occupied = obs_bins,
    observed_full_new_bins_vs_baseline = obs_new,
    observed_full_tv_to_baseline = obs_tv,
    
    rarefied_bins_median = median(sim_occ_bins),
    rarefied_bins_q025   = as.numeric(quantile(sim_occ_bins, 0.025, names = FALSE)),
    rarefied_bins_q975   = as.numeric(quantile(sim_occ_bins, 0.975, names = FALSE)),
    
    rarefied_newbins_median = median(sim_new_bins),
    rarefied_newbins_q025   = as.numeric(quantile(sim_new_bins, 0.025, names = FALSE)),
    rarefied_newbins_q975   = as.numeric(quantile(sim_new_bins, 0.975, names = FALSE)),
    
    rarefied_tv_median = median(sim_tv),
    rarefied_tv_q025   = as.numeric(quantile(sim_tv, 0.025, names = FALSE)),
    rarefied_tv_q975   = as.numeric(quantile(sim_tv, 0.975, names = FALSE)),
    
    p_occ_more_than_baseline = p_occ_more,
    p_newbins_more_than_0    = p_new_more,
    p_tv_greater_than_0      = p_tv_more
  )
  
  # Per-bin CI table
  bin_ci_results[[sname]] <- tibble(
    slice_target = sname,
    lon_bin = all_bin_levels,
    base_count = as.integer(table(factor(base_df$lon_bin, levels = all_bin_levels))),
    target_full_count = as.integer(table(factor(dfT$lon_bin, levels = all_bin_levels))),
    rarefied_count_median = apply(sim_bin_counts, 2, median),
    rarefied_count_q025   = apply(sim_bin_counts, 2, quantile, probs = 0.025, names = FALSE),
    rarefied_count_q975   = apply(sim_bin_counts, 2, quantile, probs = 0.975, names = FALSE)
  )
}

bin_summary_df <- bind_rows(bin_sim_results)
write_csv(bin_summary_df, file.path(out_dir, "longitude_bin_rarefaction_summary_vs_80-70.csv"))
print(bin_summary_df)

bin_ci_df <- bind_rows(bin_ci_results)
write_csv(bin_ci_df, file.path(out_dir, "longitude_bin_rarefaction_binCIs_vs_80-70.csv"))


# =========================================================
# (B) Single-statistic rarefaction: longitude range
#     full range + central 95% range
# ========================================================

lon_range_full <- function(lon) max(lon) - min(lon)

lon_range_95 <- function(lon) {
  qs <- quantile(lon, probs = c(0.025, 0.975), names = FALSE)
  qs[2] - qs[1]
}

# Baseline values
base_range_full <- lon_range_full(base_df$lon)
base_range_95   <- lon_range_95(base_df$lon)

range_results <- list()

for (sname in slice_targets) {
  
  dfT_raw <- get_slice_localities(sname)
  N_T <- nrow(dfT_raw)
  
  if (N_T < N_base) {
    warning("Skipping ", sname, " for range test: target N (", N_T, ") < baseline N (", N_base, ").")
    next
  }
  
  obs_full_range <- lon_range_full(dfT_raw$lon)
  obs_95_range   <- lon_range_95(dfT_raw$lon)
  
  sim_full <- numeric(nsim)
  sim_95   <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    idx <- sample.int(N_T, size = N_base, replace = FALSE)
    lon_s <- dfT_raw$lon[idx]
    sim_full[i] <- lon_range_full(lon_s)
    sim_95[i]   <- lon_range_95(lon_s)
  }
  
  p_full_more <- (1 + sum(sim_full <= base_range_full)) / (nsim + 1)
  p_95_more   <- (1 + sum(sim_95   <= base_range_95))   / (nsim + 1)
  
  range_results[[sname]] <- tibble(
    slice_target = sname,
    N_rarefied_to = N_base,
    
    baseline_range_full = base_range_full,
    baseline_range_95   = base_range_95,
    
    observed_full_range_target = obs_full_range,
    observed_95_range_target   = obs_95_range,
    
    rarefied_full_median = median(sim_full),
    rarefied_full_q025   = as.numeric(quantile(sim_full, 0.025, names = FALSE)),
    rarefied_full_q975   = as.numeric(quantile(sim_full, 0.975, names = FALSE)),
    
    rarefied_95_median = median(sim_95),
    rarefied_95_q025   = as.numeric(quantile(sim_95, 0.025, names = FALSE)),
    rarefied_95_q975   = as.numeric(quantile(sim_95, 0.975, names = FALSE)),
    
    p_full_more_than_baseline = p_full_more,
    p_95_more_than_baseline   = p_95_more
  )
}

range_summary_df <- bind_rows(range_results)
write_csv(range_summary_df, file.path(out_dir, "longitude_range_rarefaction_vs_80-70.csv"))
print(range_summary_df)


#============================================================
#           Make mini-table that sums up results
#===========================================================

stopifnot(exists("range_summary_df"), exists("bin_summary_df"))

baseline_slice <- slice_baseline

mini_tbl <- range_summary_df %>%
  left_join(
    bin_summary_df %>%
      select(slice_target, new_bins_median = rarefied_newbins_median),
    by = "slice_target"
  ) %>%
  bind_rows(
    tibble::tibble(
      slice_target = baseline_slice,
      N_rarefied_to = unique(range_summary_df$N_rarefied_to)[1],
      baseline_range_full = unique(range_summary_df$baseline_range_full)[1],
      baseline_range_95   = unique(range_summary_df$baseline_range_95)[1],
      rarefied_full_median = unique(range_summary_df$baseline_range_full)[1],
      p_full_more_than_baseline = NA_real_,
      rarefied_95_median = unique(range_summary_df$baseline_range_95)[1],
      p_95_more_than_baseline = NA_real_,
      new_bins_median = 0
    )
  ) %>%
  mutate(
    Slice_ka = dplyr::case_when(
      slice_target == "slice_80000_70000" ~ "80–70",
      slice_target == "slice_70000_60000" ~ "70–60",
      slice_target == "slice_60000_50000" ~ "60–50",
      slice_target == "slice_50000_40000" ~ "50–40",
      TRUE ~ slice_target
    )
  ) %>%
  select(
    Slice_ka,
    Rarefied_N = N_rarefied_to,
    Full_longitude_range_median_deg = rarefied_full_median,
    p_full_gt_80_70 = p_full_more_than_baseline,
    Central95_range_median_deg = rarefied_95_median,
    p_95_gt_80_70 = p_95_more_than_baseline,
    New_longitude_bins_median = new_bins_median
  ) %>%
  mutate(
    Full_longitude_range_median_deg = round(Full_longitude_range_median_deg, 1),
    Central95_range_median_deg      = round(Central95_range_median_deg, 1),
    p_full_gt_80_70                 = ifelse(is.na(p_full_gt_80_70), NA, round(p_full_gt_80_70, 3)),
    p_95_gt_80_70                   = ifelse(is.na(p_95_gt_80_70), NA, round(p_95_gt_80_70, 3))
  ) %>%
  arrange(factor(Slice_ka, levels = c("80–70", "70–60", "60–50", "50–40")))

print(mini_tbl)
write_csv(mini_tbl, file.path(out_dir, "mini_table_vs_80-70.csv"))
