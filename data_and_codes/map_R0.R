

library(sf)
library(ggplot2)
library(dplyr)
library(tigris)
library(patchwork)
library(readr)

# Temperature-dependent functions
F_T <- function(T) {
  ifelse(T > 6 & T < 28, -24.58678 * T^2 + 835.9505 * T - 4105.579, 0)
}

s_E_T <- function(T) {
  T_min <- 6; T_max <- 28; s_min <- 0.1; s_max <- 0.8
  ifelse(T >= T_min & T <= T_max,
         s_min + (s_max - s_min) * (T - T_min) / (T_max - T_min), 0)
}

theta_L <- function(T) {
  ifelse(T > 10.8 & T < 30.2,
         (-0.0105 * T^2 + 0.4316 * T - 3.424) * (-1) *
           ((0.03116 - 0.007615 * DL + 0.00004469 * DL^2) /
              (1 - 0.1374 * DL + 0.004788 * DL^2)), 0)
}

F_LM <- function(T) { a * (M^b) * theta_L(T) }
F_NM <- F_LM
F_NH <- function(T) { a * (H^b) * theta_L(T) }

mu_L <- function(T) { ifelse(T > 10 & T < 30, (1 / (8.560 + 20.654 * (1 + (T / 17.759)^6.827)^-1)), 0) }
mu_N <- mu_L

b_M <- function(T) { ifelse(T > 10 & T < 25, (0.3 * exp(-((T - 18.5) / (25 - 15))^2)), 0) }

K_M <- function(T) { ifelse(T > 7 & T < 30, 150 * (1 - ((T - 18.5) / 23)^2), 0) }

mu <- function(T) { ifelse(T > 5 & T < 30, (0.012 + 0.4 * ((T - 18.5) / 25)^2), 0) }

# Parameters
a <- 0.0013; b <- 0.515; M <- 20000; H <- 600; DL <- 12

# County data
county_data <- data.frame(
  County = c("Allegany", "Anne Arundel", "Baltimore City", "Baltimore", "Calvert",
             "Caroline", "Carroll", "Cecil", "Charles", "Dorchester",
             "Frederick", "Garrett", "Harford", "Howard", "Kent",
             "Montgomery", "Prince George's", "Queen Anne's", "St. Mary's",
             "Somerset", "Talbot", "Washington", "Wicomico", "Worcester"),
  Population = c(68106, 588261, 585708, 854535, 92783,
                 33293, 172891, 103725, 166617, 32531,
                 271717, 28806, 260924, 332317, 19198,
                 1062061, 967201, 47798, 113777,
                 24620, 37526, 154705, 103588, 52460),
  Land_Area = c(1099, 1074, 209, 1549, 552,
                829, 1160, 896, 1186, 1400,
                1710, 1676, 1132, 648, 717,
                1285, 1250, 1097, 1093,
                832, 695, 1190, 974, 1212),
  Mean_Temp = c(15.3, 17.2, 16.9, 16.7, 18.3,
                17.7, 16.6, 17.3, 18.4, 17.9,
                16.1, 14.4, 16.8, 17.1, 17.5,
                17.0, 17.2, 17.6, 18.6,
                18.2, 17.8, 15.8, 18.0, 18.1),
  beta_NM = c(0.767924, 0.7353, 0.607473, 0.77486, 0.47175,
              0.457793, 0.7444, 0.6057433, 0.4567952, 0.42515,
              0.783, 0.7612, 0.727243, 0.7957, 0.45656,
              0.7867, 0.600576773, 0.486565, 0.46001,
              0.427665, 0.4757906, 0.687738, 0.4852, 0.42475167)
)

# R0 calculation function
calculate_R0 <- function(T, county_name) {
  tryCatch({
    county_info <- county_data %>% filter(County == county_name)
    
    L <- 1e5 * county_info$Land_Area
    Mt <- 100 * county_info$Land_Area
    Nt <- 1e2 * county_info$Land_Area
    Nh <- 6e3 * county_info$Population
    beta_NM <- county_info$beta_NM
    
    FLM <- F_LM(T); FNM <- F_NM(T); FNH <- F_NH(T)
    muN <- mu_N(T); d <- mu(T); KM <- K_M(T) * county_info$Land_Area; bM <- b_M(T)
    beta_L <- 1; d_M <- 0.015; m_L <- 27.8; m_N <- 4.5; h_N <- 2; s <- 1
    
    numerator <- sqrt(
      FLM *
        ((Mt * m_N + Nt) * (Nh * h_N + Nt) * muN +
           Nh * Mt * FNM * h_N +
           (Nh * FNH * m_N + Nt * FNM) * Mt + Nh * Nt * FNH) *
        beta_NM * Mt * FNM *
        (Mt * d_M + d * KM) * KM *
        (Mt * m_L + L) * beta_L * L * (Nh * h_N + Nt)
    )
    
    denominator <- ((m_N * muN + FNM) * Nt +
                      ((m_N * muN + FNM) * h_N + FNH * m_N) * Nh) * Mt +
      ((h_N * muN + FNH) * Nh + muN * Nt) * Nt
    denominator <- denominator * (Mt * d_M + d * KM) * (Mt * m_L + L)
    
    R0 <- s * numerator / denominator
    return(ifelse(is.finite(R0) && R0 >= 0, R0, 0))
  }, error = function(e) {
    message("Error calculating R0 for ", county_name, ": ", e$message)
    return(0)
  })
}

# Load and prepare shapefile
md_counties <- counties(state = "MD", cb = TRUE, resolution = "20m") %>%
  st_transform(crs = 4326) %>%
  mutate(NAME = case_when(
    GEOID == "24510" ~ "Baltimore City",
    GEOID == "24005" ~ "Baltimore",
    TRUE ~ gsub(" County", "", NAME)
  ))

# Calculate R0 values
simulated_data <- county_data %>%
  rowwise() %>%
  mutate(
    R0_Current = calculate_R0(Mean_Temp, County),
    R0_Increased_25 = calculate_R0(Mean_Temp + 2.5, County),
    R0_Increased_45 = calculate_R0(Mean_Temp + 4.5, County)
  ) %>%
  ungroup()

# Merge with spatial data
md_counties <- md_counties %>%
  left_join(simulated_data %>% select(County, R0_Current, R0_Increased_25, R0_Increased_45),
            by = c("NAME" = "County")) %>%
  mutate(across(starts_with("R0_"), ~replace_na(.x, 0)))

# Centroids for labeling
county_centroids <- st_centroid(md_counties) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(NAME = md_counties$NAME)

# Color palette
turbo_palette <- c("#30123B", "#4662D7", "#36AAF9", "#1AE4B6", 
                   "#72FE5E", "#FABA39", "red", "#7A0403")

# Plot function
make_r0_map <- function(data, fill_var, title_label) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), color = "white", size = 0.1) +
    geom_text(data = county_centroids, aes(x = X, y = Y, label = NAME),
              size = 3.75, fontface = "bold", color = "black") +
    scale_fill_gradientn(
      colors = turbo_palette,
      limits = c(0.5, 2.0),
      name = expression(R[0]),
      breaks = seq(0.5, 2.0, by = 0.5),
      na.value = "grey50"
    ) +
    labs(title = title_label) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0),
      legend.title = element_text(face = "bold", size = 20),
      legend.text = element_text(size = 20),
      legend.position = "right",
      legend.key.height = unit(0.1, "npc")
    )
}

# Create maps
map_r0_current <- make_r0_map(md_counties, "R0_Current", bquote((a) ~ R[0] ~ "Current"))
map_r0_25 <- make_r0_map(md_counties, "R0_Increased_25", bquote((b) ~ R[0] ~ "+2.5°C"))
map_r0_45 <- make_r0_map(md_counties, "R0_Increased_45", bquote((c) ~ R[0] ~ "+4.5°C"))

# Combine plots
combined_plot <- map_r0_current | map_r0_25 | map_r0_45

# Display
print(combined_plot)

# Save (optional)
# ggsave("results/figures/Maryland_Lyme_R0.png", 
#        combined_plot, width = 15, height = 5, dpi = 300)


