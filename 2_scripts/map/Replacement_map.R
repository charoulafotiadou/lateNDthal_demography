# Load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(readr)

# Load data
Neanderthals_database <- read_csv("lateNDthal_demography/1_data/map/Neanderthals_database.csv")

# Get world map
sf::sf_use_s2(FALSE)
worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
eurasia_africa <- st_crop(worldmap, xmin = -10, xmax = 90, ymin = 35, ymax = 72)

# Function to remove transparency from hex colors
opaque_hex <- function(hex_color, bg = c(255,255,255)) {
  hex_color <- gsub("#","", hex_color)
  if(nchar(hex_color) == 8) {
    r <- strtoi(substr(hex_color,1,2),16L)
    g <- strtoi(substr(hex_color,3,4),16L)
    b <- strtoi(substr(hex_color,5,6),16L)
    a <- strtoi(substr(hex_color,7,8),16L)/255
    r_opaque <- round(r*a + bg[1]*(1-a))
    g_opaque <- round(g*a + bg[2]*(1-a))
    b_opaque <- round(b*a + bg[3]*(1-a))
    sprintf("#%02X%02X%02X", r_opaque, g_opaque, b_opaque)
  } else if(nchar(hex_color) == 6) {
    paste0("#", hex_color)
  } else NA
}

# Remove transparency
Neanderthals_database <- Neanderthals_database %>%
  mutate(colour_no_alpha = sapply(colour, opaque_hex))

# Prepare separate offsets for each plot
Neanderthals_database <- Neanderthals_database %>%
  group_by(long, lat) %>%
  mutate(
    n = n(),
    # Zoom1: small spacing
    grid_x_zoom1 = ifelse(
      n == 1,
      long,
      long + (row_number() - (n + 1)/2) * 0.05
    ),
    grid_y_zoom1 = lat,
    # Zoom2: larger spacing
    grid_x_zoom2 = ifelse(
      n == 1,
      long,
      long + (row_number() - (n + 1)/2) * 0.25
    ),
    grid_y_zoom2 = lat,
    # Full map: largest spacing
    grid_x_full = ifelse(
      n == 1,
      long,
      long + (row_number() - (n + 1)/2) * 1.3
    ),
    grid_y_full = lat
  ) %>%
  ungroup()

# ----------------------------
# Zoomed maps
# ----------------------------
zoom1 <- ggplot() +
  geom_sf(data = worldmap) +
  geom_text(
    data = Neanderthals_database,
    aes(x = long, y = lat, label = arch_site),
    size = 5, vjust = -1
  ) +
  geom_point(
    data = Neanderthals_database,
    aes(x = grid_x_zoom1, y = grid_y_zoom1, fill = colour_no_alpha),
    shape = 21, color = "black", size = 5
  ) +
  scale_fill_identity() +
  coord_sf(xlim = c(3, 7), ylim = c(49, 51)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

zoom2 <- ggplot() +
  geom_sf(data = worldmap) +
  geom_text(
    data = Neanderthals_database,
    aes(x = long, y = lat, label = arch_site),
    size = 5, vjust = -1
  ) +
  geom_point(
    data = Neanderthals_database,
    aes(x = grid_x_zoom2, y = grid_y_zoom2, fill = colour_no_alpha),
    shape = 21, color = "black", size = 5
  ) +
  scale_fill_identity() +
  coord_sf(xlim = c(75, 90), ylim = c(48, 55)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# ----------------------------
# Full faceted map
# ----------------------------
full_map <- ggplot() +
  geom_sf(data = eurasia_africa) +
  geom_text(
    data = Neanderthals_database,
    aes(x = long, y = lat, label = arch_site),
    size = 4, vjust = +1.2
  ) +
  geom_point(
    data = Neanderthals_database,
    aes(x = grid_x_full, y = grid_y_full, fill = colour_no_alpha),
    shape = 21, color = "black", size = 4
  ) +
  scale_fill_identity() +
  facet_wrap(~ period) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 14, family = "Arial"),
    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, family = "Arial"),
    legend.position = "none"
  )

# ----------------------------
# Display plots
# ----------------------------
zoom1
zoom2
full_map
