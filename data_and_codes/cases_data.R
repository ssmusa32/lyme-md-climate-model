

# Load necessary libraries
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(patchwork)
library(tidyr)

# File paths
data_file <- "data/processed/MD_LD_cases_pop_temp.csv"
output_folder <- "results/figures"

# Load the Lyme disease cases data
cases_data <- read_csv(data_file)

# Convert the dataset from wide format to long format
cases_long <- cases_data %>%
  select(County, starts_with("Cases")) %>%
  pivot_longer(cols = starts_with("Cases"), names_to = "year", values_to = "cases") %>%
  mutate(year = as.numeric(sub("Cases", "", year))) %>%
  mutate(County = gsub(" County", "", County)) %>%
  mutate(County = factor(County))

# Function to create a yearly cases plot for a single county
create_cases_plot <- function(data, county_name, panel_letter, show_x_text = FALSE) {
  # Interpolate data for smoother line
  n_points <- 10 * (max(data$year) - min(data$year) + 1)
  spline_data <- data.frame(
    year = spline(data$year, data$cases, n = n_points)$x,
    cases = spline(data$year, data$cases, n = n_points)$y
  )
  
  p <- ggplot() +
    geom_path(data = spline_data, aes(x = year, y = cases), 
              lineend = "round", color = "red", linewidth = 1.2) +
    geom_point(data = data, aes(x = year, y = cases), 
               color = "red", size = 1.2) +
    labs(title = paste0("(", panel_letter, ") ", county_name),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(face = "bold", size = 16),
      axis.text.x = if (show_x_text) element_text(face = "bold", size = 16) else element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  return(p)
}

# Get unique counties
counties <- unique(cases_long$County)

# Create plots for all 24 counties
plots <- lapply(seq_along(counties), function(i) {
  county <- counties[i]
  county_data <- cases_long %>% filter(County == !!county)
  panel_letter <- letters[i]
  show_x_text <- i %in% 21:24
  create_cases_plot(county_data, county, panel_letter, show_x_text = show_x_text)
})

# Combine plots into a 6x4 grid
combined_plot <- wrap_plots(plots, ncol = 4, nrow = 6) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Date (years)",
    theme = theme(
      plot.caption = element_text(size = 18, hjust = 0.5, face = "bold", margin = margin(b = 10)),
      plot.margin = margin(t = 10, r = 10, b = 20, l = 40)
    )
  ) +
  geom_text(
    data = data.frame(x = -Inf, y = Inf),
    aes(x = x, y = y, label = "Reported Cases"),
    angle = 90, hjust = 1, vjust = -1, size = 18/2.845,
    inherit.aes = FALSE
  )

# Print the static combined plot
print(combined_plot)

# Save the static combined plot
static_plot_path <- file.path(output_folder, "Maryland_Yearly_Lyme_Cases_All_Counties.png")
ggsave(filename = static_plot_path, plot = combined_plot, 
       width = 16, height = 24, dpi = 300)

# Create interactive version
p_interactive <- ggplotly(combined_plot)
print(p_interactive)

# Save interactive version (optional)

