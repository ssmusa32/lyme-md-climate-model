

library(ggplot2)
library(dplyr)
library(tidyr)
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

a <- 0.0013; b <- 0.515; M <- 2454612; H <- 6612021; DL <- 12

theta_L <- function(T) {
  ifelse(T > 10.8 & T < 30.2,
         (-0.0105 * T^2 + 0.4316 * T - 3.424) * (-1) *
           ((0.03116 - 0.007615 * DL + 0.00004469 * DL^2) /
              (1 - 0.1374 * DL + 0.004788 * DL^2)), 0)
}

F_L <- function(T) { a * (M^b) * theta_L(T) }
F_N <- F_L
F_H <- function(T) { a * (H^b) * theta_L(T) }

a_A <- 0.086; b_A <- 0.515; D <- 19533

theta_A <- function(T) { ifelse(T > 0 & T < 20.2, -0.0095 * T^2 + 0.19 * T + 0.05, 0) }
F_A <- function(T) { a_A * (D^b_A) * theta_A(T) }

mu_L <- function(T) {
  ifelse(T > 10 & T < 30, (1 / (8.560 + 20.654 * (1 + (T / 19.759)^6.827)^-1)), 0)
}
mu_N <- mu_L

mu_A <- function(T) {
  ifelse(T > 10 & T < 30, (-log(-0.000828 * T^2 + 0.0367 * T + 0.522)), 0)
}

b_M <- function(T) {
  ifelse(T > 10 & T < 25, (0.3 * exp(-((T - 18.5) / (25 - 15))^2)), 0)
}

K <- function(T) { ifelse(T > 7 & T < 30, 150 * 25314 * (1 - ((T - 18.5) / 23)^2), 0) }

mu_M <- function(T) { ifelse(T > 5 & T < 30, (0.012 + 0.4 * ((T - 18.5) / 25)^2), 0) }

# R0 calculation with controls
R0_formula_temp <- function(C_p, epsilon_p, epsilon_c, epsilon_r, r, T) {
  F_L_val <- F_L(T)
  F_N_val <- F_N(T)
  F_H_val <- F_H(T)
  mu_N_val <- mu_N(T)
  mu_M_val <- mu_M(T)
  K_val <- K(T)
  b_M_val <- b_M(T)
  
  beta_L_val <- 1 * (1 - epsilon_c)
  beta_NM_val <- 1.75 * (1 - r * epsilon_r)
  F_H_control <- F_H_val * (1 - C_p * epsilon_p)
  
  m_L <- 27.8
  m_N <- 2.5
  m_A <- 239
  delta_M <- 0.015
  Pi_H <- 242
  mu_H <- 0.0000366
  h_N <- 2
  L_star <- 25983.98 * 25314
  S_N_star <- 240.87 * 25314
  
  N_M_star <- (b_M_val - mu_M_val) * K_val / delta_M
  N_H_star <- Pi_H / mu_H
  
  if (N_M_star <= 0 || !is.finite(N_M_star)) {
    return(NA)
  }
  
  N_N_star <- S_N_star
  
  A_1 <- beta_NM_val * F_N_val * (h_N * N_H_star + N_N_star) * N_M_star
  A_2 <- (N_M_star * m_N + N_N_star) * (N_H_star * h_N + N_N_star) * mu_N_val
  A_3 <- (F_N_val * h_N + F_H_control * m_N) * N_M_star * N_H_star
  A_4 <- (F_N_val + F_H_control) * N_N_star * N_M_star
  B_1 <- ((m_N * mu_N_val + F_N_val) * N_N_star + N_H_star * ((m_N * mu_N_val + F_N_val) * h_N + F_H_control * m_N)) * N_M_star
  B_2 <- ((h_N * mu_N_val + F_H_control) * N_H_star + mu_N_val * N_N_star) * N_N_star
  
  term_1 <- beta_L_val * L_star * F_L_val / ((m_L * N_M_star + L_star) * (mu_M_val + delta_M * N_M_star / K_val))
  term_2 <- (A_1 * (A_2 + A_3 + A_4)) / ((B_1 + B_2)^2)
  
  R0 <- sqrt(term_1 * term_2)
  
  return(ifelse(is.finite(R0) && R0 >= 0, R0, NA))
}

# Parameter ranges
r_values <- seq(0, 1, length.out = 50)
epsilon_r_values <- seq(0, 1, length.out = 50)

# Generate grids
grid_data_a <- expand.grid(r = r_values, epsilon_r = epsilon_r_values)
grid_data_b <- expand.grid(r = r_values, epsilon_r = epsilon_r_values)

# Calculate R0 values
grid_data_a$R0 <- mapply(
  function(r, epsilon_r) R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = epsilon_r, r = r, T = 18),
  r = grid_data_a$r, epsilon_r = grid_data_a$epsilon_r
)
grid_data_b$R0 <- mapply(
  function(r, epsilon_r) R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = epsilon_r, r = r, T = 20.5),
  r = grid_data_b$r, epsilon_r = grid_data_b$epsilon_r
)

# Calculate R_c at intersection
R_c_18 <- R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = 0.75, r = 0.82, T = 18)
R_c_20.5 <- R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = 0.75, r = 0.82, T = 20.5)

# Color palette
turbo_palette <- c("#30123B", "#4662D7", "#36AAF9", "#1AE4B6", "#72FE5E",
                   "#C8EF34", "#FABA39", "#F66B19", "#CB2A04", "#7A0403")

# Plot function
create_plot <- function(data, temp, R_c, label) {
  ggplot(data, aes(x = r, y = epsilon_r, fill = R0)) +
    geom_raster() +
    geom_segment(aes(x = 0.82, xend = 0.82, y = 0, yend = 0.75), 
                 color = "black", linetype = "dashed", linewidth = 1.5) +
    geom_segment(aes(x = 0, xend = 0.82, y = 0.75, yend = 0.75), 
                 color = "black", linetype = "dashed", linewidth = 1.5) +
    geom_point(aes(x = 0.82, y = 0.75), color = "black", size = 10, shape = 19) +
    geom_text(aes(x = 0.65, y = 0.77, label = sprintf("bold(R[c]^{(1)})~'='~%.2f", R_c)), 
              parse = TRUE, color = "black", size = 12, hjust = 0, vjust = 0, fontface = "bold") +
    geom_text(aes(x = 0.5, y = 1.075, label = sprintf("bold(T[A])==%.1f*degree*C", temp)), 
              parse = TRUE, color = "black", fontface = "bold", size = 14) +
    geom_text(aes(x = 0.015, y = 1.075, label = label), 
              color = "black", fontface = "bold", size = 14) +
    scale_fill_gradientn(colors = turbo_palette, name = expression(bold(R[c]^(1))), 
                         limits = c(0, 1.8)) +
    labs(x = expression(bold("Rodent bait coverage"~(r))), 
         y = if(label == "(a)") expression(bold("Bait efficacy"~(epsilon[r]))) else "") +
    theme_minimal() +
    theme(
      text = element_text(face = "bold", size = 35),
      axis.title = element_text(face = "bold", size = 35),
      axis.text = element_text(face = "bold", size = 35),
      legend.title = element_text(face = "bold", size = 35),
      legend.position = "right"
    )
}

# Create plots
plot_a <- create_plot(grid_data_a, 18, R_c_18, "(a)")
plot_b <- create_plot(grid_data_b, 20.5, R_c_20.5, "(b)")

# Combine plots
combined_plot <- plot_a + plot_b +
  plot_layout(ncol = 2, guides = "collect")

# Display
print(combined_plot)

# Print baseline values
R0_baseline_18 <- R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = 0, r = 0, T = 18)
cat(sprintf("Baseline R_c at 18°C: %.4f\n", R0_baseline_18))

R0_baseline_20.5 <- R0_formula_temp(C_p = 0, epsilon_p = 0, epsilon_c = 0, epsilon_r = 0, r = 0, T = 20.5)
cat(sprintf("Baseline R_c at 20.5°C: %.4f\n", R0_baseline_20.5))

# Print R_c values
cat(sprintf("R_c at (r=0.82, ε_r=0.75) for 18°C: %.2f\n", R_c_18))
cat(sprintf("R_c at (r=0.82, ε_r=0.75) for 20.5°C: %.2f\n", R_c_20.5))


