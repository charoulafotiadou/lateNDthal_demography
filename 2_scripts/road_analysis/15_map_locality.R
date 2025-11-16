# Loading Packages
library(dplyr)
library(sf)
library(ggplot2)
library(ggnewscale)
library(patchwork)
library(geodata)

# --- Load helper ---
read_txt_safe <- function(path) {
  read.csv2(path, header = TRUE, colClasses = "character", stringsAsFactors = FALSE)
}

# --- Load files ---
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

data_3 <- dplyr::bind_rows(lapply(data_3_files, read_txt_safe))
assemblages <- dplyr::bind_rows(data_1, data_2, data_3)

# --- Clean and rename ---
assemblages <- assemblages %>%
  mutate(
    Age_MIN = as.numeric(sub("#.*", "", Age_MIN)),
    Age_MAX = as.numeric(sub("#.*", "", Age_MAX)),
    age_diff = Age_MAX - Age_MIN
  ) %>%
  rename(age_min = Age_MIN, age_max = Age_MAX)

# --- Filter by age span ---
min_assemblage_age_diff <- 2000
max_assemblage_age_diff <- 30000

assemblages <- assemblages %>%
  filter(age_diff >= min_assemblage_age_diff & age_diff <= max_assemblage_age_diff)

# ---- Remove duplicates: same locality + same age span ----
assemblages <- assemblages %>%
  distinct(Locality, age_min, age_max, .keep_all = TRUE)

# --- Time slice parameters ---
form_age_max <- 130000
form_age_min <- 40000
form_timeslice_length <- 10000
num_timeslices <- ceiling((form_age_max - form_age_min) / form_timeslice_length)

# --- Build time slices ---
timeslices <- list()
for (i in 0:(num_timeslices - 1)) {
  slice_max <- form_age_max - form_timeslice_length * i
  slice_min <- max(form_age_max - form_timeslice_length * (i + 1), form_age_min)
  
  ts_df <- assemblages %>%
    filter(age_min < slice_max & age_max > slice_min) %>%
    distinct(Locality, .keep_all = TRUE)

  
  slice_name <- paste0("slice_", slice_max, "_", slice_min)
  timeslices[[slice_name]] <- ts_df
}



############## MAP #######################


# --- Convert all time slices to sf objects ---
# First ensure X and Y are numeric
timeslices <- lapply(timeslices, function(df) {
  df$X <- as.numeric(df$X)
  df$Y <- as.numeric(df$Y)
  df
})

# Convert each time slice to sf
timeslices_sf <- lapply(timeslices, function(df) {
  st_as_sf(df, coords = c("X", "Y"), crs = 4326)
})

# --- DEM Preparation ---
dem_global <- elevation_global(res = "2.5", path = "1_data/road_analysis/DEM_data")
dem_df <- as.data.frame(dem_global, xy = TRUE)
colnames(dem_df) <- c("Longitude", "Latitude", "Elevation")
dem_df <- dem_df %>% mutate(Elevation = ifelse(Elevation <= 0, NA, Elevation))  # Mask water

# --- Plot Settings ---
min_elevation <- 0
max_elevation <- 6000
min_density <- 0
max_density <- 0.01
xlim_bounds <- c(-10, 45)
ylim_bounds <- c(30, 60)

# --- Plot Function ---
create_density_map <- function(df, title) {
  ggplot() +
    geom_tile(data = dem_df, aes(x = Longitude, y = Latitude, fill = Elevation)) +
    scale_fill_gradient(
      low = "lightgrey", high = "black",
      na.value = "transparent", name = "Elevation",
      limits = c(min_elevation, max_elevation),
      oob = scales::squish,
      breaks = c(min_elevation, max_elevation),
      labels = c("Low", "High"),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    ggnewscale::new_scale_fill() +
    stat_density_2d(
      data = as.data.frame(st_coordinates(df)),
      aes(X, Y, fill = after_stat(level)),
      geom = "polygon", alpha = 0.5, h = c(5, 5), bins = 10
    ) +
    scale_fill_viridis_c(
      option = "magma", name = "Density",
      limits = c(min_density, max_density),
      breaks = c(min_density, max_density),
      labels = c("Low", "High"),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    geom_sf(data = df, shape = 21, fill = "white", color = "black",
            size = 2, stroke = 0.5, alpha = 1.0) +
    coord_sf(xlim = xlim_bounds, ylim = ylim_bounds) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.title = element_text(hjust = 0.5)
    ) +
    labs(title = title, x = NULL, y = NULL)
}

# --- Build Maps from sf Time Slices ---
# Define the order and labels
slice_names <- c(
  "slice_130000_120000", "slice_120000_110000", "slice_110000_1e+05",
  "slice_1e+05_90000",  "slice_90000_80000",   "slice_80000_70000",
  "slice_70000_60000",   "slice_60000_50000",   "slice_50000_40000"
)

map_titles <- c(
  "130 – 120 ka", "120 – 110 ka", "110 – 100 ka",
  "100 –  90 ka", "90 –  80 ka",  "80 –  70 ka",
  "70 –  60 ka",  "60 –  50 ka",  "50 –  40 ka"
)

# Generate all maps in a loop
maps <- mapply(
  function(name, title) {
    create_density_map(timeslices_sf[[name]], title)
  },
  name = slice_names,
  title = map_titles,
  SIMPLIFY = FALSE
)

# --- Combine in 3x3 Panel with Shared Legend ---
panel <- wrap_plots(maps, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

# --- Display Panel ---
panel

# --- Save Panel to File ---
ggsave(
  filename = "3_output/road_analysis/density_panel_locality.png",
  plot = panel,
  width = 12,
  height = 10,
  dpi = 300,
  units = "in",
  bg = "white"
)

