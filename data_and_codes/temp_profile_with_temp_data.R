



library(ggplot2)
library(dplyr)
library(readr)
library(lubridate)
library(patchwork)

# File path
temp_file <- "data/processed/MD_Statewide_Daily_Temp_2001_2022.csv"

# Load and clean temperature data
temperature_data <- read_csv(temp_file)

colnames(temperature_data) <- tolower(colnames(temperature_data))

temperature_data <- temperature_data %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         year = as.numeric(year),
         temperature_avg = as.numeric(temperature_avg),
         temp_max = as.numeric(temp_max),
         temp_min = as.numeric(temp_min)) %>%
  filter(!is.na(date) & !is.na(temperature_avg) & !is.na(temp_max) & !is.na(temp_min) & !is.na(year))

# Filter 2022 data
temp_filtered_data <- temperature_data %>%
  filter(year == 2022) %>%
  mutate(month = month(date))

# Month labels
month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# Parameter functions
F_T <- function(T) { ifelse(T > 6 & T < 28, -24.58678 * T^2 + 835.9505 * T - 4105.579, 0) }
b_T <- function(T) { ifelse(T > 6 & T < 28, 0.5 * F_T(T), 0) }
F_L <- function(T) { ifelse(T > 10.8 & T < 30.2, 0.1 * (T - 10), 0) }
F_H <- function(T) { ifelse(T > 10.8 & T < 30.2, 0.1 * (T - 10), 0) }
F_A <- function(T) { ifelse(T > 0 & T < 20.2, 0.08 * (T - 5), 0) }
mu_L <- function(T) { ifelse(T > 10 & T < 35, 1 / (8.56 + 20.65 * (1 + (T / 19.75)^6.82)^-1), 0) }
mu_A <- function(T) { ifelse(T > 10 & T < 35, -log(-0.000828 * T^2 + 0.0367 * T + 0.522), 0) }
b_M <- function(T) { ifelse(T > 5 & T < 25, 0.3 * exp(-((T - 20) / (25 - 10))^2), 0) }
K_T <- function(T) { ifelse(T > 7 & T < 30, 150 * (1 - ((T - 18.5) / 23)^2), 0) }
d_T <- function(T) { ifelse(T > 5 & T < 30, 0.012 + 0.4 * ((T - 18.5) / 25)^2, 0) }

# Functions and labels
functions_list <- list(b_T, F_L, F_H, F_A, mu_L, mu_A, b_M, K_T, d_T)
titles_list <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)")

y_labels <- list(
  expression(bold(b[T](T[A]))), 
  expression(bold(F[AD](T[A]))),
  expression(bold(F[LM](T[A]) ~ "or" ~ F[NM](T[A]))), 
  expression(bold(F[NH](T[A]))),
  expression(bold(mu[L](T[A]) ~ "or" ~ mu[N](T[A]))), 
  expression(bold(mu[A](T[A]))), 
  expression(bold(b[M](T[A]))), 
  expression(bold(K[M](T[A]))), 
  expression(bold(mu[M](T[A])))
)

# Create plots
plot_list <- list()

for (i in seq_along(functions_list)) {
  func <- functions_list[[i]]
  
  temp_filtered_data <- temp_filtered_data %>%
    mutate(main_curve = func(temperature_avg))
  
  monthly_avg_data <- temp_filtered_data %>%
    group_by(month) %>%
    summarise(main_avg = mean(main_curve, na.rm = TRUE)) %>%
    ungroup()
  
  p <- ggplot(monthly_avg_data, aes(x = month)) + 
    geom_smooth(aes(y = main_avg), color = "blue", method = "loess", linewidth = 1.5, se = FALSE) +  
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    annotate("text", x = 1, y = Inf, label = titles_list[i], 
             hjust = 0, vjust = 1.5, size = 5, fontface = "bold") +
    labs(x = "Months", y = y_labels[[i]]) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
  
  plot_list[[i]] <- p
}

# Combine plots
final_plot <- (plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) /
  (plot_list[[4]] + plot_list[[5]] + plot_list[[6]]) /
  (plot_list[[7]] + plot_list[[8]] + plot_list[[9]])

# Print and save
print(final_plot)

ggsave(
  filename = "results/figures/Combined_Temperature_Plots.png",
  plot = final_plot, width = 12, height = 12, dpi = 300
)
