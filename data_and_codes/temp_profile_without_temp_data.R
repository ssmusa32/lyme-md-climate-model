

library(ggplot2)
library(patchwork)

# Temperature-dependent functions
F_T <- function(T) {
  ifelse(T > 6 & T < 28, -24.58678 * T^2 + 835.9505 * T - 4105.579, 0)
}

s_E_T <- function(T) {
  T_min <- 6; T_max <- 28; s_min <- 0.1; s_max <- 0.8
  ifelse(T >= T_min & T <= T_max, s_min + (s_max - s_min) * (T - T_min) / (T_max - T_min), 0)
}
b_T <- function(T) { F_T(T) * s_E_T(T) }

a <- 0.0013; b <- 0.515 
M <- 1777.67
H <- 244
DL <- 12

theta_L <- function(T) {
  if (T > 10.8 & T < 30.2) {
    (-0.0105 * T^2 + 0.4316 * T - 3.424) * (-1) * 
      ((0.03116 - 0.007615 * DL + 0.00004469 * DL^2) / (1 - 0.1374 * DL + 0.004788 * DL^2))
  } else { 0 }
}

F_LM <- function(T) { a * (M^b) * theta_L(T) }
F_NH <- function(T) { a * (H^b) * theta_L(T) }

a_A <- 0.086; b_A <- 0.515; D <- 20

theta_A <- function(T) { 
  ifelse(T > 0 & T < 20.2, (-0.0095 * T^2 + 0.219 * T + 0.05)/2.332, 0) 
}

F_A <- function(T) { a_A * (D^b_A) * theta_A(T) }

mu_L <- function(T) {
  ifelse(T > 10 & T < 30, (1 / (8.560 + 20.654 * (1 + (T / 19.759)^6.827)^-1)), 0)
}

mu_N <- mu_L

mu_A <- function(T) {
  ifelse(T > 10 & T < 30, (-log(-0.000828 * T^2 + 0.0367 * T + 0.522)), 0)
}

b_M <- function(T) {
  ifelse(T > 10 & T < 25, (0.3 * exp(-((T - 20) / (25 - 10))^2)), 0)
}

K_M <- function(T) { ifelse(T > 7 & T < 30, 150 * (1 - ((T - 18.5) / 23)^2), 0) }

mu <- function(T) { ifelse(T > 5 & T < 30, (0.012 + 0.4 * ((T - 18.5) / 25)^2), 0) }

# Generate parameter data
data_list <- list(
  data.frame(Temperature = seq(6, 28, by = 0.1), Value = sapply(seq(6, 28, by = 0.1), b_T), Label = "(a) "),
  data.frame(Temperature = seq(0, 20.2, by = 0.1), Value = sapply(seq(0, 20.2, by = 0.1), F_A), Label = "(b) "),
  data.frame(Temperature = seq(10.8, 30.2, by = 0.1), Value = sapply(seq(10.8, 30.2, by = 0.1), F_LM), Label = "(c) "),
  data.frame(Temperature = seq(10.8, 30.2, by = 0.1), Value = sapply(seq(10.8, 30.2, by = 0.1), F_NH), Label = "(d) "),
  data.frame(Temperature = seq(10.1, 29.9, by = 0.1), Value = sapply(seq(10.1, 29.9, by = 0.1), mu_L), Label = "(e) "),
  data.frame(Temperature = seq(10.1, 29.9, by = 0.1), Value = sapply(seq(10.1, 29.9, by = 0.1), mu_A), Label = "(f) "),
  data.frame(Temperature = seq(10, 25, by = 0.1), Value = sapply(seq(10, 25, by = 0.1), b_M), Label = "(g) "),
  data.frame(Temperature = seq(7, 30, by = 0.1), Value = sapply(seq(7, 30, by = 0.1), K_M), Label = "(h) "),
  data.frame(Temperature = seq(5.1, 29.9, by = 0.1), Value = sapply(seq(5.1, 29.9, by = 0.1), mu), Label = "(i) ")
)

# Color palette
colors <- c("#0073C2FF", "#E69F00", "#009E73", "#CC79A7", "#F0E442",
            "#56B4E9", "#D55E00", "#6600CC", "#FF3300")

# Create individual plots
plot_list <- lapply(seq_along(data_list), function(i) {
  ggplot(data_list[[i]], aes(x = Temperature, y = Value)) +
    geom_line(color = colors[i], size = 1.8) +
    ggtitle(data_list[[i]]$Label[1]) +
    labs(
      x = ifelse(i %in% 7:9, "Temperature (°C)", ""),
      y = switch(i,
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
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
})

# Combine plots
combined_plot <- (plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) / 
  (plot_list[[4]] + plot_list[[5]] + plot_list[[6]]) / 
  (plot_list[[7]] + plot_list[[8]] + plot_list[[9]])

# Display plot
print(combined_plot)

# Save plot (optional)


