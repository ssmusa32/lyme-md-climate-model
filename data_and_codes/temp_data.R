
# Load necessary libraries
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(patchwork)

# File paths
temp_file <- "data/processed/MD_Statewide_Daily_Temp_2001_2022.csv"
cases_file <- "data/processed/MD_and_Counties_Cases.csv"
output_folder <- "results/figures"

# Load temperature data
temperature_data <- read_csv(temp_file)
temperature_data$date <- as.Date(temperature_data$date, format = "%m/%d/%Y")

# Create temperature plot
p_temp <- ggplot(temperature_data, aes(x = date, y = temperature_avg)) +
  geom_ribbon(aes(ymin = temp_min, ymax = temp_max), fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", linewidth = 0.8) +
  labs(title = "(a) ",
       x = NULL, y = "Temperature (°C)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10)
  )

# Load Lyme disease cases data
statewide_cases <- read_csv(cases_file)

# Prepare cases data
statewide_long <- statewide_cases %>%
  select(year, md_cases) %>%
  mutate(year = as.numeric(year))

# Smooth the cases data
n_points <- 10 * (max(statewide_long$year) - min(statewide_long$year) + 1)
spline_data <- data.frame(
  year = spline(statewide_long$year, statewide_long$md_cases, n = n_points)$x,
  md_cases = spline(statewide_long$year, statewide_long$md_cases, n = n_points)$y
)

# Create Lyme disease plot
p_lyme <- ggplot() +
  geom_path(data = spline_data, aes(x = year, y = md_cases), 
            lineend = "round", color = "red", linewidth = 1.2) +
  geom_point(data = statewide_long, aes(x = year, y = md_cases), 
             color = "red", size = 1.2) +
  labs(title = "(b) ",
       x = "Date (years)", y = "Reported Lyme disease cases") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10)
  )

# Combine plots
combined_plot <- p_temp / p_lyme +
  plot_layout(ncol = 1, heights = c(1, 1))

# Print and save
print(combined_plot)

static_plot_path <- file.path(output_folder, "MD_Temp_Lyme_Combined.png")
ggsave(filename = static_plot_path, plot = combined_plot, 
       width = 16, height = 12, dpi = 300)

# Interactive version
p_interactive <- ggplotly(combined_plot)
print(p_interactive)

# Save interactive (optional)


