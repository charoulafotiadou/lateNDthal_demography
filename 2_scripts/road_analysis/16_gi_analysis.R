
# Loading Packages

library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)
library(units)
library(sfdep)
library(spdep)
library(rnaturalearth)

set.seed(42)  # reproducibility for all permutation tests


# Helper + I/O

read_txt_safe <- function(path) {
  # robust CSV2 reader (strings as chars, not factors)
  read.csv2(path, header = TRUE, colClasses = "character", stringsAsFactors = FALSE)
}
dir.create("3_output/road_analysis", showWarnings = FALSE, recursive = TRUE)  # ensure output folder exists


# Load data

data_1 <- read_txt_safe("1_data/road_analysis/askROAD/europe_130-60_middle_pal.txt")                      
data_2 <- read_txt_safe("1_data/road_analysis/askROAD/europe_asia_africa_130-40_humanremains_neanderthalensis.txt")

data_3_files <- list(
  "1_data/road_analysis/askROAD/data3_chatelperronian.txt",
  "1_data/road_analysis/askROAD/data3_denticulate_mousterian.txt",
  "1_data/road_analysis/askROAD/data3_late_mousterian.txt",
  "1_data/road_analysis/askROAD/data3_micoquian.txt",
  "1_data/road_analysis/askROAD/data3_mousterian_acheulean_tradition.txt",
  "1_data/road_analysis/askROAD/data3_mousterian_eurasia.txt"
)
data_3 <- dplyr::bind_rows(lapply(data_3_files, read_txt_safe))  # stack the “data_3_” files
assemblages <- dplyr::bind_rows(data_1, data_2, data_3)          # combine all sources


# Clean + filter

assemblages <- assemblages %>%
  mutate(
    Age_MIN = as.numeric(sub("#.*", "", Age_MIN)),
    Age_MAX = as.numeric(sub("#.*", "", Age_MAX)),
    age_diff = Age_MAX - Age_MIN,
    X = as.numeric(X),
    Y = as.numeric(Y)
  ) %>%
  rename(age_min = Age_MIN, age_max = Age_MAX) %>%
  filter(!is.na(X) & !is.na(Y))

# Age-span filter: exclude too-broad/too-narrow intervals
min_assemblage_age_diff <- 2000
max_assemblage_age_diff <- 30000
assemblages <- assemblages %>%
  filter(age_diff >= min_assemblage_age_diff & age_diff <= max_assemblage_age_diff)

# Remove duplicates: same locality + same age span
assemblages <- assemblages %>%
  distinct(Locality, age_min, age_max, .keep_all = TRUE)


# Time slices (10 ka) from 130–40 ka

form_age_max <- 130000   # oldest bound (years BP)
form_age_min <-  40000   # youngest bound (years BP)
form_timeslice_length <- 10000  # slice width (years)
slice_edges <- seq(form_age_max, form_age_min, by = -form_timeslice_length)  # descending edges

slice_names <- vapply(  # e.g., "slice_130000_120000", ..., "slice_50000_40000"
  seq_len(length(slice_edges) - 1),
  function(i) sprintf("slice_%d_%d", slice_edges[i], slice_edges[i+1]),
  FUN.VALUE = character(1)
)

# Assign assemblages to slices if their age interval overlaps the slice
timeslices <- setNames(vector("list", length(slice_names)), slice_names)
for (i in seq_along(slice_names)) {
  smax <- slice_edges[i]
  smin <- slice_edges[i+1]
  timeslices[[i]] <- assemblages %>%
    filter(age_min < smax & age_max > smin) %>%
    distinct(Locality, .keep_all = TRUE) %>%
    mutate(slice = slice_names[i])
}
names(timeslices) <- slice_names

# sf points per slice (WGS84 lon/lat)
timeslices_sf <- lapply(timeslices, function(df) st_as_sf(df, coords = c("X","Y"), crs = 4326))


# Study window + hex grid
win_wgs  <- st_as_sfc(st_bbox(c(xmin = -10, ymin = 30, xmax = 45, ymax = 60), crs = 4326))  # lon/lat box
win_laea <- st_transform(win_wgs, 3035) # project to 3035

# Tunable parameters for hex and neighbors
hex_size_km      <- 50     # hex “diameter” ~ 50 km
neighbor_band_km <- 150    # distance band for neighbors (Gi*)
gi_perm_nsim     <- 999    # permutations per slice for p-values (sfdep)
fdr_alpha        <- 0.10   # per-slice FDR threshold for hotspots

# Build hex grid + neighbor weights
cellsize_m <- hex_size_km * 1000
hex <- st_make_grid(win_laea, cellsize = cellsize_m, square = FALSE) |>
  st_as_sf() |>
  mutate(id = dplyr::row_number())
hex_cent <- st_centroid(st_geometry(hex))
nb <- st_dist_band(hex_cent, upper = neighbor_band_km * 1000)
lw <- spdep::nb2listw(nb, style = "W")


# Count points into hexes
count_slice <- function(pts_sf) {
  pts_laea <- st_transform(pts_sf, 3035)
  st_join(hex, pts_laea, join = st_contains, left = TRUE) |>
    st_drop_geometry() |>
    count(id, name = "n_sites") |>
    right_join(tibble::tibble(id = hex$id), by = "id") |>
    mutate(n_sites = tidyr::replace_na(n_sites, 0L)) |>
    arrange(id)
}

# Same logic but for an arbitrary subset of points (used in permutation tests)
make_counts_from_pts <- function(pts_subset){
  pts_laea <- st_transform(pts_subset, 3035)
  st_join(hex, pts_laea, join = st_contains, left = TRUE) |>
    st_drop_geometry() |>
    count(id, name = "n_sites") |>
    right_join(tibble::tibble(id = hex$id), by = "id") |>
    mutate(n_sites = tidyr::replace_na(n_sites, 0L)) |>
    arrange(id)
}

counts_by_slice <- lapply(timeslices_sf, count_slice)
names(counts_by_slice) <- names(timeslices_sf)


# Per-slice Gi*: p-values via sfdep; z-scores via spdep
gi_results_list <- lapply(names(counts_by_slice), function(sname){
  df <- counts_by_slice[[sname]] |> arrange(id)
  # Permutation-based p-values (right-tailed; sfdep computes its own weights)
  gi <- sfdep::local_gstar_perm(x = df$n_sites, nb = nb, wt = NULL,
                                nsim = gi_perm_nsim, alternative = "greater")
  # FDR hotspot flag
  df |>
    mutate(
      p_perm  = gi$p_value,
      q_fdr   = p.adjust(p_perm, method = "BH"),
      hotspot = !is.na(q_fdr) & q_fdr <= fdr_alpha,
      slice   = sname
    )
})
gi_results <- bind_rows(gi_results_list) |>
  left_join(hex, by = "id") |>
  st_as_sf(crs = 3035)

# Deterministic Gi* z-scores (spdep::localG) for intensity summaries and ΔGi*
if (!"gi_star" %in% names(gi_results)) {
  message("gi_star not found — computing z-scores with spdep::localG …")
  z_by_slice <- lapply(names(counts_by_slice), function(sname){
    df <- counts_by_slice[[sname]] |> arrange(id)
    gz <- spdep::localG(df$n_sites, listw = lw, zero.policy = TRUE)
    tibble::tibble(id = df$id, slice = sname, gi_star = as.numeric(gz))
  }) |> dplyr::bind_rows()
  gi_results <- gi_results |>
    mutate(slice = as.character(slice)) |>
    left_join(z_by_slice |> mutate(slice = as.character(slice)),
              by = c("id","slice"))
}

# ============================================================
# (1a) Slice-wise summary (global; focus on change through time)
# ============================================================
cell_area_km2 <- set_units(st_area(hex)[1], km^2) |> as.numeric()  # per-hex area (km2)

slice_summary <- gi_results |>
  st_drop_geometry() |>
  group_by(slice) |>
  summarize(
    n_cells         = dplyr::n(), # number of hexes in window
    hotspots_total  = sum(hotspot, na.rm=TRUE), # how many hotspots this slice
    sites_total     = sum(n_sites, na.rm=TRUE), # total sites counted into hexes
    .groups = "drop"
  ) |>
  mutate(
    area_hotspots_total_km2 = hotspots_total * cell_area_km2  # area proxy for hotspot footprint
  ) |>
  arrange(factor(slice, levels = names(timeslices_sf)))
write_csv(slice_summary, "3_output/road_analysis/gi_slice_summary_10ka_130-40.csv")

# (Optional) quick map panel to eyeball per-slice hotspots (binary: hotspot vs other)
label_from_slice <- function(s) paste0(gsub("_","–", gsub("^slice_", "", s)), " ka")
plot_list <- lapply(names(timeslices_sf), function(sname){
  df_hex <- gi_results |> filter(slice == sname)
  ggplot() +
    geom_sf(data = df_hex, aes(fill = ifelse(hotspot, "Hotspot", "Other")), color = NA) +
    scale_fill_manual(values = c("Hotspot" = "#c51b7d", "Other" = "#f0f0f0")) +
    geom_sf(data = st_transform(timeslices_sf[[sname]], 3035), size = 0.3, alpha = 0.6) +
    coord_sf(xlim = st_bbox(win_laea)[c("xmin","xmax")],
             ylim = st_bbox(win_laea)[c("ymin","ymax")]) +
    labs(title = label_from_slice(sname), fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", panel.grid = element_blank())
})
panel <- patchwork::wrap_plots(plot_list, ncol = 3)
ggsave("3_output/road_analysis/gi_hotspots_panel_10ka_130-40.png",
       panel, width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================
# (1b) Slice-wise global permutation test on hotspot counts
#       Null: with the same n_sites multiset, randomly assign
#       counts to hexes → is observed #hotspots unusually high?
# ============================================================
slice_perm_nsim <- 2000   # permutations per slice
zcut_global     <- 1.96   # hotspot rule for this test (z >= 1.96)

# Deterministic hotspot counter with fixed z cutoff
count_hotspots_zcut <- function(count_df, listw, zcut = zcut_global) {
  if (sum(count_df$n_sites, na.rm = TRUE) == 0) return(0L)  # no sites → 0 hotspots
  gz <- spdep::localG(count_df$n_sites, listw = listw, zero.policy = TRUE)
  gz <- as.numeric(gz)
  gz[is.na(gz)] <- -Inf
  sum(gz >= zcut, na.rm = TRUE)
}

# For each slice: (a) observed; (b) permute n_sites across hexes
slice_perm_df <- lapply(names(counts_by_slice), function(sname){
  df <- counts_by_slice[[sname]] |> arrange(id)
  obs_hot <- count_hotspots_zcut(df, lw, zcut = zcut_global) # observed
  perm_vals <- replicate(slice_perm_nsim, { # null by permuting counts across ids
    df_perm <- df
    df_perm$n_sites <- sample(df$n_sites, replace = FALSE)
    count_hotspots_zcut(df_perm, lw, zcut = zcut_global)
  })
  tibble::tibble(
    slice        = sname,
    obs_hotspots = obs_hot,
    mean_null    = mean(perm_vals),
    sd_null      = stats::sd(perm_vals),
    q05_null     = as.numeric(stats::quantile(perm_vals, 0.05, names = FALSE)),
    q95_null     = as.numeric(stats::quantile(perm_vals, 0.95, names = FALSE)),
    p_one_sided  = mean(perm_vals >= obs_hot) # Pr(null >= observed)
  )
}) |> bind_rows() |>
  arrange(factor(slice, levels = names(timeslices_sf)))
write_csv(slice_perm_df, "3_output/road_analysis/perm_slicewise_hotspot_counts.csv")

# ============================================================
# (1c) Consecutive (rolling) contrasts: each slice vs previous
#      → persist/new/lost/other + Jaccard + ΔGi*
# ============================================================
slice_order <- names(timeslices_sf)  # chronological order from 130–120 → 50–40

# Helper: build a one-row summary + a map-ready sf for a pair
compare_pair <- function(prev_s, curr_s) {
  prev_hex <- gi_results %>%
    st_drop_geometry() %>%
    filter(slice == prev_s) %>%
    select(id, gi_prev = gi_star, hot_prev = hotspot)
  curr_hex <- gi_results %>%
    st_drop_geometry() %>%
    filter(slice == curr_s) %>%
    select(id, gi_curr = gi_star, hot_curr = hotspot)
  
  j <- tibble(id = hex$id) %>% # ensure all hex ids appear (including zero-count)
    left_join(prev_hex, by = "id") %>%
    left_join(curr_hex, by = "id") %>%
    mutate(
      change_class = dplyr::case_when(
        hot_prev &  hot_curr ~ "persist_hot",
        hot_prev & !hot_curr ~ "lost_hot",
        !hot_prev & hot_curr ~ "new_hot",
        TRUE                 ~ "other"
      ),
      delta_gi = gi_curr - gi_prev
    )
  
  # one-row summary for this transition
  sum_row <- j %>%
    summarize(
      slice_prev         = prev_s,
      slice_curr         = curr_s,
      shared_hotspots    = sum(change_class == "persist_hot", na.rm = TRUE),
      new_hotspots       = sum(change_class == "new_hot",     na.rm = TRUE),
      lost_hotspots      = sum(change_class == "lost_hot",    na.rm = TRUE),
      prev_hotspots      = sum(isTRUE(hot_prev),              na.rm = TRUE),
      curr_hotspots      = sum(isTRUE(hot_curr),              na.rm = TRUE),
      jaccard_overlap    = shared_hotspots / (shared_hotspots + new_hotspots + lost_hotspots + 1e-9),
      mean_delta_gi      = mean(delta_gi, na.rm = TRUE),
      median_delta_gi    = median(delta_gi, na.rm = TRUE)
    )
  
  # map-ready sf for this transition
  j_sf <- j %>% left_join(hex, by = "id") %>% st_as_sf(crs = 3035)
  list(summary = sum_row, map_sf = j_sf)
}

# Build all consecutive pairs (slice i-1 vs slice i)
pairs_idx <- seq(2, length(slice_order))
pair_outputs <- lapply(pairs_idx, function(i) compare_pair(slice_order[i-1], slice_order[i]))

# Bind rolling summaries to a CSV
consec_summary <- bind_rows(lapply(pair_outputs, `[[`, "summary")) %>%
  mutate(
    pair_label = paste0(
      gsub("_","–", gsub("^slice_", "", slice_prev)), " ka → ",
      gsub("_","–", gsub("^slice_", "", slice_curr)), " ka"
    ),
    slice_prev = factor(slice_prev, levels = slice_order),
    slice_curr = factor(slice_curr, levels = slice_order)
  )
write_csv(consec_summary, "3_output/road_analysis/gi_consecutive_contrasts.csv")

# =================================
# (1d) Rolling change-class panel
# =================================

# Helper to create clean slice labels (in ka, not raw years)
label_from_slice <- function(s) {
  # pull the two numbers from e.g. "slice_130000_120000"
  nums <- as.numeric(unlist(regmatches(s, gregexpr("[0-9]+", s))))
  if (length(nums) == 2) {
    # integer division keeps it tidy: 130000 -> 130
    sprintf("%d–%d ka", nums[1] %/% 1000, nums[2] %/% 1000)
  } else {
    s
  }
}

# Plotting window
map_wgs  <- sf::st_as_sfc(sf::st_bbox(c(xmin = -10, ymin = 30, xmax = 40, ymax = 50), crs = 4326))
map_laea <- sf::st_transform(map_wgs, 3035)

# Basemap layers
world_outline <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  sf::st_transform(3035) |>
  suppressWarnings(sf::st_intersection(map_laea))
coastlines <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") |>
  sf::st_transform(3035) |>
  suppressWarnings(sf::st_intersection(map_laea))

# Fixed class levels & colors
change_levels <- c("persist_hot","new_hot","lost_hot")
change_colors <- c("persist_hot"="#6a51a3", "new_hot"="#238b45", "lost_hot"="#cb181d")

# Legend anchor using points
bb <- sf::st_bbox(map_laea)
cx <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
cy <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
legend_pts_df <- data.frame(
  x = c(cx - 15000, cx, cx + 15000),
  y = c(cy - 15000, cy, cy + 15000),
  change_class = factor(change_levels, levels = change_levels)
)

# Draw each rolling panel with consistent basemap + legend
plot_change_consec <- lapply(seq_along(pair_outputs), function(i){
  m <- pair_outputs[[i]]$map_sf
  m_plot <- m %>%
    dplyr::filter(change_class %in% change_levels) %>% # drop "other" to see basemap
    dplyr::mutate(change_class = factor(change_class, levels = change_levels))
  
  prev_s <- as.character(consec_summary$slice_prev[i])
  curr_s <- as.character(consec_summary$slice_curr[i])
  title_txt <- paste0(label_from_slice(prev_s), "  \u2192  ", label_from_slice(curr_s))
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = world_outline, color = "grey25", fill = "grey90", linewidth = 0.35) +
    ggplot2::geom_sf(data = coastlines,    color = "grey20", linewidth = 0.30) +
    ggplot2::geom_sf(data = m_plot, ggplot2::aes(fill = change_class), color = NA, alpha = 0.95) +
    ggplot2::geom_point(  # invisible anchor so all three legend keys always show
      data = legend_pts_df,
      mapping = ggplot2::aes(x = x, y = y, fill = change_class),
      inherit.aes = FALSE, shape = 22, size = 3, alpha = 0, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values       = change_colors,
      breaks       = change_levels,
      limits       = change_levels,
      drop         = FALSE,
      na.translate = FALSE,
      labels       = c("Persisting hotspot","New hotspot","Lost hotspot"),
      name         = NULL,
      guide        = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 4, shape = 22, colour = NA)  # uniform, solid keys
      )
    ) +
    ggplot2::coord_sf(
      xlim = sf::st_bbox(map_laea)[c("xmin","xmax")],
      ylim = sf::st_bbox(map_laea)[c("ymin","ymax")]
    ) +
    ggplot2::labs(
      title = title_txt,
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(size = 8)
    )
})
consec_panel <- patchwork::wrap_plots(plot_change_consec, ncol = 3)
ggplot2::ggsave("3_output/road_analysis/gi_change_vs_previous_panel.png",
                consec_panel, width = 13, height = 9, dpi = 300, bg = "white")

# ============================================================
# (2a) Baseline-anchored contrasts (global) vs 130–120 ka,
#      persist/new/lost + Jaccard + ΔGi* + matched panel
# ============================================================
baseline_slice <- "slice_130000_120000"  # reference slice

base_hex <- gi_results %>%
  st_drop_geometry() %>%
  filter(slice == baseline_slice) %>%
  select(id, gi_star_base = gi_star, hotspot_base = hotspot)

other_slices <- setdiff(names(timeslices_sf), baseline_slice)

# Summaries vs baseline
contrast_summary <- lapply(other_slices, function(sname){
  cur <- gi_results %>%
    st_drop_geometry() %>%
    filter(slice == sname) %>%
    select(id, gi_star, hotspot)
  
  j <- base_hex %>% left_join(cur, by = "id") %>%
    mutate(
      slice = sname,
      change_class = case_when(
        hotspot_base &  hotspot ~ "persist_hot",
        hotspot_base & !hotspot ~ "lost_hot",
        !hotspot_base & hotspot ~ "new_hot",
        TRUE                    ~ "other"
      ),
      delta_gi = gi_star - gi_star_base
    )
  
  sum_all <- j %>% summarize(
    shared_hotspots   = sum(change_class == "persist_hot", na.rm = TRUE),
    new_hotspots      = sum(change_class == "new_hot",     na.rm = TRUE),
    lost_hotspots     = sum(change_class == "lost_hot",    na.rm = TRUE),
    baseline_hotspots = sum(hotspot_base,                  na.rm = TRUE),
    current_hotspots  = sum(hotspot,                       na.rm = TRUE),
    jaccard_overlap   = shared_hotspots / (shared_hotspots + new_hotspots + lost_hotspots + 1e-9),
    mean_delta_gi     = mean(delta_gi, na.rm = TRUE),
    median_delta_gi   = median(delta_gi, na.rm = TRUE)
  )
  
  bind_cols(tibble(slice = sname), sum_all) %>%
    mutate(slice = factor(slice, levels = names(timeslices_sf)))
}) |> bind_rows() |>
  arrange(slice)

write_csv(contrast_summary, "3_output/road_analysis/gi_baseline_contrasts_vs_130-120ka.csv")

# Matched-style baseline-anchored panel
plot_change_baseline <- lapply(as.character(contrast_summary$slice), function(sname){
  cur <- gi_results %>%
    st_drop_geometry() %>%
    filter(slice == sname) %>%
    select(id, gi_star, hotspot)
  
  j <- base_hex %>% left_join(cur, by = "id") %>%
    mutate(
      change_class = dplyr::case_when(
        hotspot_base &  hotspot ~ "persist_hot",
        hotspot_base & !hotspot ~ "lost_hot",
        !hotspot_base & hotspot ~ "new_hot",
        TRUE                    ~ "other"
      )
    ) %>%
    left_join(hex, by = "id") %>%
    st_as_sf(crs = 3035)
  
  j_plot <- j %>% # drop "other" to see basemap
    dplyr::filter(change_class %in% change_levels) %>%
    dplyr::mutate(change_class = factor(change_class, levels = change_levels))
  
  title_txt <- paste0(label_from_slice(as.character(sname)), " vs ", label_from_slice(baseline_slice))
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = world_outline, color = "grey25", fill = "grey90", linewidth = 0.35) +
    ggplot2::geom_sf(data = coastlines,    color = "grey20", linewidth = 0.30) +
    ggplot2::geom_sf(data = j_plot, ggplot2::aes(fill = change_class), color = NA, alpha = 0.95) +
    ggplot2::geom_point(
      data = legend_pts_df,
      mapping = ggplot2::aes(x = x, y = y, fill = change_class),
      inherit.aes = FALSE, shape = 22, size = 3, alpha = 0, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values       = change_colors,
      breaks       = change_levels,
      limits       = change_levels,
      drop         = FALSE,
      na.translate = FALSE,
      labels       = c("Persisting hotspot","New hotspot","Lost hotspot"),
      name         = NULL,
      guide        = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 4, shape = 22, colour = NA)
      )
    ) +
    ggplot2::coord_sf(
      xlim = sf::st_bbox(map_laea)[c("xmin","xmax")],
      ylim = sf::st_bbox(map_laea)[c("ymin","ymax")]
    ) +
    ggplot2::labs(
      title = title_txt,
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(size = 8)
    )
})
baseline_panel <- patchwork::wrap_plots(plot_change_baseline, ncol = 3)
ggplot2::ggsave("3_output/road_analysis/gi_change_vs_baseline_panel.png",
                baseline_panel, width = 13, height = 9, dpi = 300, bg = "white")

# ============================================================
# (2b) Structured permutation test on baseline-anchored design
#      Test if 'new French hotspots vs baseline' > random
#     (France mask only used here)
# ============================================================
fr_poly   <- rnaturalearth::ne_countries(scale = "medium", country = "France", returnclass = "sf") |>
  st_transform(3035)
fr_mainland <- st_intersection(st_geometry(fr_poly), win_laea)  # clip to analysis window
fr_union    <- st_union(fr_mainland) # dissolved France polygon
in_fr_vec   <- st_intersects(hex, fr_union, sparse = FALSE)[,1] # hex-level France mask

nsim <- 2000 # permutations per slice
target_slices <- setdiff(names(timeslices_sf), baseline_slice)
results_perm <- vector("list", length(target_slices))

for (i in seq_along(target_slices)) {
  sname <- target_slices[i]
  
  # Combine points from baseline + current (projected to analysis CRS)
  pts_baseline <- st_transform(timeslices_sf[[baseline_slice]], 3035) %>% mutate(period = "Baseline")
  pts_current  <- st_transform(timeslices_sf[[sname]],          3035) %>% mutate(period = "Current")
  pts_all <- dplyr::bind_rows(pts_baseline, pts_current)
  
  # Statistic: # of new hotspots inside France (Current vs Baseline) using z ≥ 1.96
  stat_fun <- function(pts){
    pts_b <- pts[pts$period == "Baseline",]
    pts_c <- pts[pts$period == "Current",]
    cb <- make_counts_from_pts(pts_b)
    cc <- make_counts_from_pts(pts_c)
    gb <- spdep::localG(cb$n_sites, listw = lw, zero.policy = TRUE)
    gc <- spdep::localG(cc$n_sites, listw = lw, zero.policy = TRUE)
    hot_b <- (as.numeric(gb) >= 1.96)
    hot_c <- (as.numeric(gc) >= 1.96)
    sum((!hot_b & hot_c) & in_fr_vec, na.rm = TRUE) # “new” hotspots in France only
  }
  
  obs_stat <- stat_fun(pts_all) # observed new-FR hotspots
  
  perm_stats <- replicate(nsim, { # null: shuffle Baseline/Current labels across fixed locations
    pts_all$period <- sample(pts_all$period, replace = FALSE)
    stat_fun(pts_all)
  })
  
  results_perm[[i]] <- tibble::tibble(
    slice       = sname,
    obs         = obs_stat,
    mean_null   = mean(perm_stats),
    p_one_sided = mean(perm_stats >= obs_stat)
  )
}

perm_results_df <- dplyr::bind_rows(results_perm)
write_csv(perm_results_df, "3_output/road_analysis/perm_baseline_vs_slice_tests.csv")
