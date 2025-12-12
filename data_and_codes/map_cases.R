
library(sf)
library(ggplot2)
library(dplyr)
library(tigris)
library(patchwork)
library(readr)
library(tidyr)

# Load Maryland county shapefiles
md_counties <- counties(state = "MD", cb = TRUE, resolution = "20m") %>%
  st_transform(crs = 4326)

# Clean county names
md_counties <- md_counties %>%
  mutate(NAME = case_when(
    GEOID == "24510" ~ "Baltimore City",
    GEOID == "24005" ~ "Baltimore",
    TRUE ~ gsub(" County", "", NAME)
  ))

# Load model-predicted and observed data
simulated_data <- read_csv("results/tables/MD_Lyme_Simulated_Cases.csv") %>%
  mutate(County = case_when(
    grepl("Baltimore city", County, ignore.case = TRUE) ~ "Baltimore City",
    grepl("Baltimore", County, ignore.case = TRUE) & !grepl("city", County, ignore.case = TRUE) ~ "Baltimore",
    TRUE ~ gsub(" County", "", County)
  ))

# Merge data with shapefile
md_counties <- md_counties %>%
  left_join(simulated_data, by = c("NAME" = "County")) %>%
  mutate(
    Cases_Current = replace_na(Cases_Current, 0),
    Cases_Cum = replace_na(Cases_Cum, 0)
  )

# Get centroids for labeling
county_centroids <- st_centroid(md_counties) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(NAME = md_counties$NAME)

# Set common color scale
common_min <- min(md_counties$Cases_Current, md_counties$Cases_Cum, na.rm = TRUE)
common_max <- max(md_counties$Cases_Current, md_counties$Cases_Cum, na.rm = TRUE)

# Panel (a): Model-predicted
map_model <- ggplot(md_counties) +
  geom_sf(aes(fill = Cases_Current), color = "white", size = 0.1) +
  geom_text(
    data = county_centroids,
    aes(x = X, y = Y, label = NAME),
    size = 5.5, fontface = "bold", color = "black"
  ) +
  scale_fill_gradient(
    low = "lightgreen",
    high = "red",
    limits = c(common_min, common_max),
    name = "Cumulative cases"
  ) +
  labs(title = "(a) Model-predicted") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    legend.title = element_text(face = "bold", size = 12),
    legend.key.height = unit(1.5, "cm"),
    legend.position = "right"
  )

# Panel (b): Observed
map_observed <- ggplot(md_counties) +
  geom_sf(aes(fill = Cases_Cum), color = "white", size = 0.1) +
  geom_text(
    data = county_centroids,
    aes(x = X, y = Y, label = NAME),
    size = 5.5, fontface = "bold", color = "black"
  ) +
  scale_fill_gradient(
    low = "lightgreen",
    high = "red",
    limits = c(common_min, common_max),
    name = "Cumulative cases"
  ) +
  labs(title = "(b) Observed") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    legend.title = element_text(face = "bold", size = 12),
    legend.key.height = unit(1.5, "cm"),
    legend.position = "right"
  )

# Print individual panels
print(map_model)
print(map_observed)

# Combine panels
combined_plot <- map_model | map_observed

# Print combined plot
print(combined_plot)

# Save plot (optional)
# ggsave("results/figures/MD_Lyme_Maps_Combined.png", 
#        combined_plot, width = 12, height = 8, dpi = 300)

