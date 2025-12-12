

library(deSolve)
library(readr)
library(ggplot2)
library(dplyr)

# Define temperature-dependent parameter functions
F_T <- function(T) {
  ifelse(T > 6 & T < 28, -24.58678 * T^2 + 835.9505 * T - 4105.579, 0)
}

s_E_T <- function(T) {
  T_min <- 6; T_max <- 28; s_min <- 0.1; s_max = 0.8
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
theta_A <- function(T) { ifelse(T > 0 & T < 25, -0.009 * T^2 + 0.19 * T + 0.05, 0) }
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

# Set temperature
T_A <- 18

# File paths
input_lyme_data <- "data/processed/MD_and_Counties_Cases.csv"
output_plot_lyme <- "results/figures/MD_fits_TC.png"
output_1_lyme <- "results/tables/MD_sim_TC.csv"
output_2_lyme <- "results/tables/MD_para_TC.csv"

# Load data
lyme_data <- read_csv(input_lyme_data)

# Demographic parameters
Nh <- 6e6
mu_h <- 0.000036
Pi_h <- Nh * mu_h
Nh_2020 <- 6e6

Nt <- 49092247.3
Mt <- 2454612
Deer <- 19533

# Initial states
H_S <- Nh * 0.7
H_E <- H_S * 0.25
H_I <- 608
H_H <- H_I * 0.1
H_R <- (H_I + H_H) * 0.05
L <- 1.665888479 * 10^7
N_S <- 0.7 * Nt
N_I <- N_S * 0.43
A <- 1.488740226 * 10^6
M_S <- 0.7 * Mt
M_I <- M_S * 0.1
C1 <- 600
t1 <- 0

# Temperature-dependent parameters
b_T_val <- b_T(T_A)
F_L_val <- F_L(T_A)
mu_L_val <- mu_L(T_A)
mu_N_val <- mu_N(T_A)
mu_A_val <- mu_A(T_A)
beta_L <- 1
b_M_val <- b_M(T_A)
F_N_val <- F_N(T_A)
m_L <- 27.8
m_N <- 2.5
d_d <- 0.02
d_M <- 0.015
K_val <- K(T_A)
m_A <- 239
F_A_val <- F_A(T_A)
F_H_val <- F_H(T_A)
h_N <- 2

vtime <- seq(1, nrow(lyme_data), 1)

# Model equations
lyme_model <- function(t, x, parameters) {
  with(as.list(c(x, parameters)), {
    dH_S <- Pi_h - (beta_H * F_H_val * N_I / (h_N * Nh + Nt)) * H_S + psi_H * H_R - mu_h * H_S
    dH_E <- (beta_H * F_H_val * N_I / (h_N * Nh + Nt)) * H_S - sigma_H * H_E - mu_h * H_E
    dH_I <- sigma_H * H_E - (tau_H + mu_h) * H_I
    dH_H <- (1 - theta_H) * tau_H * H_I - (gamma_H + mu_h) * H_H
    dH_R <- theta_H * tau_H * H_I + gamma_H * H_H - (psi_H + mu_h) * H_R
    dL <- b_T_val * F_A_val * Deer * A / (m_A * Deer + A) - F_L_val * Mt * L / (m_L * Mt + L) - mu_L_val * L
    dN_S <- F_L_val * (M_S + (1 - beta_L) * M_I) * L / (m_L * Mt + L) - mu_N_val * N_S
    dN_I <- beta_L * F_L_val * M_I * L / (m_L * Mt + L) - ((F_N_val * Mt * (N_I / (m_N * Mt + Nt))) + (F_H_val * Nh * N_I / (h_N * Nh + Nt))) - mu_N_val * N_I
    dA <- ((F_N_val * Mt * N_S / (m_N * Mt + Nt)) + (F_H_val * Nh * N_S / (h_N * Nh + Nt))) + ((F_N_val * Mt * N_I / (m_N * Mt + Nt)) + (F_H_val * N_I / (h_N * Nh + Nt))) - F_A_val * Deer * A / (m_A * Deer + A) - mu_A_val * A
    dM_S <- b_M_val * Mt - M_S * (d_d + d_M * Mt / K_val) - beta_M * F_N_val * M_S * N_I / (m_N * Mt + Nt)
    dM_I <- beta_M * F_N_val * M_S * N_I / (m_N * Mt + Nt) - M_I * (d_d + d_M * Mt / K_val)
    dC1 <- ((sigma_H * H_E) + ((1 - theta_H) * tau_H * H_I))^(rhom)
    dt1 <- 1 / 365
    
    list(c(dH_S, dH_E, dH_I, dH_H, dH_R, dL, dN_S, dN_I, dA, dM_S, dM_I, dC1, dt1))
  })
}

# Optimization function
optimize_model <- function(params) {
  parameters <- c(
    Pi_h = Pi_h, mu_h = mu_h, b_T_val = b_T_val, F_A_val = F_A_val, F_L_val = F_L_val, F_N_val = F_N_val,
    beta_H = params[1], sigma_H = params[2], psi_H = params[3],
    tau_H = params[4], gamma_H = params[5], theta_H = params[6], beta_M = params[7], rhom = params[8],
    beta_L = beta_L, mu_L_val = mu_L_val, mu_N_val = mu_N_val, mu_A_val = mu_A_val,
    d_d = d_d, d_M = d_M, K_val = K_val, m_L = m_L, m_N = m_N, m_A = m_A, F_H_val = F_H_val, h_N = h_N
  )
  
  inits <- c(H_S = H_S, H_E = H_E, H_I = H_I, H_H = H_H, H_R = H_R, L = L, 
             N_S = N_S, N_I = N_I, A = A, M_S = M_S, M_I = M_I, C1 = C1, t1 = t1)
  
  lyme_sim <- as.data.frame(lsoda(inits, vtime, lyme_model, parameters))
  observed <- lyme_data$md_ccases
  predicted <- lyme_sim$C1
  P_chi_sq <- sum((observed - predicted)^2 / predicted)
  return(P_chi_sq)
}

# Optimize parameters
init_params <- c(0.4, 0.15, 0.0055, 0.02, 0.4, 0.375, 0.75, 0.65)
opt_result <- optim(par = init_params, fn = optimize_model, method = "L-BFGS-B",
                    lower = c(0.03, 0.01, 0.0001, 0.001, 0.03, 0.005, 0.6, 0.1), 
                    upper = c(1.2, 0.2, 0.01, 0.03, 0.5, 0.7, 2.5, 1.0))

# Extract optimized parameters
optimized_params <- opt_result$par
best_beta_H <- optimized_params[1]
best_sigma_H <- optimized_params[2]
best_psi_H <- optimized_params[3]
best_tau_H <- optimized_params[4]
best_gamma_H <- optimized_params[5]
best_theta_H <- optimized_params[6]
best_beta_M <- optimized_params[7]
best_rhom <- optimized_params[8]

# Final parameters
parameters <- c(
  Pi_h = Pi_h, mu_h = mu_h, b_T_val = b_T_val, F_A_val = F_A_val, F_L_val = F_L_val, F_N_val = F_N_val,
  beta_H = best_beta_H, sigma_H = best_sigma_H, psi_H = best_psi_H,
  tau_H = best_tau_H, gamma_H = best_gamma_H, theta_H = best_theta_H,
  beta_M = best_beta_M, rhom = best_rhom, beta_L = beta_L, mu_L_val = mu_L_val, mu_N_val = mu_N_val, mu_A_val = mu_A_val,
  d_d = d_d, d_M = d_M, K_val = K_val, m_L = m_L, m_N = m_N, m_A = m_A, F_H_val = F_H_val, h_N = h_N
)

# Run simulation
inits <- c(H_S = H_S, H_E = H_E, H_I = H_I, H_H = H_H, H_R = H_R, L = L, 
           N_S = N_S, N_I = N_I, A = A, M_S = M_S, M_I = M_I, C1 = C1, t1 = t1)
lyme_sim <- as.data.frame(lsoda(inits, vtime, lyme_model, parameters))

# Save simulation results
write.csv(lyme_sim, output_1_lyme, row.names = FALSE)

# Calculate metrics
observed <- lyme_data$md_ccases
y_bfit <- lyme_sim$C1

MAPE <- mean(abs((observed - y_bfit) / observed)) * 100
SStot <- sum((observed - mean(observed))^2)
SSres <- sum((observed - y_bfit)^2)
R2 <- 1 - (SSres / SStot)

t_test_result <- t.test(observed, y_bfit, paired = TRUE)
residuals <- observed - y_bfit
mean_diff <- mean(residuals)
sd_diff <- sd(residuals)
n <- length(residuals)
z <- 1.96
CI_lower <- mean_diff - z * (sd_diff / sqrt(n))
CI_upper <- mean_diff + z * (sd_diff / sqrt(n))

disease_burden_2020 <- lyme_sim$H_I[20]

# Plot
narrowing_factor <- 0.04
widening_factor <- 0.052
ci_width <- (narrowing_factor * y_bfit) + (widening_factor * max(y_bfit))

p1 <- ggplot(lyme_data, aes(x = year, y = md_ccases)) +
  geom_ribbon(aes(ymin = y_bfit - ci_width, ymax = y_bfit + ci_width), 
              fill = "gray70", alpha = 0.4) +
  geom_point(aes(color = "Observed data"), size = 3.8, shape = 19) +
  geom_line(aes(y = y_bfit, color = "Model prediction"), linewidth = 2) +
  labs(
    title = " ",
    x = "Time (years)",
    y = "Cumulative cases (Maryland, 18°C)"
  ) +
  scale_color_manual(values = c("Observed data" = "red", "Model prediction" = "black")) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5)
  )

print(p1)

# Save plot
ggsave(output_plot_lyme, plot = p1, width = 8, height = 8, dpi = 300)

# Prepare output table
output_data <- data.frame(
  County = "Maryland",
  beta_H = round(best_beta_H, 4),
  sigma_H = round(best_sigma_H, 4),
  psi_H = round(best_psi_H, 4),
  tau_H = round(best_tau_H, 4),
  gamma_H = round(best_gamma_H, 4),
  theta_H = round(best_theta_H, 4),
  beta_M = round(best_beta_M, 4),
  rhom = round(best_rhom, 4),
  MAPE = round(MAPE, 2),
  R_squared = round(R2, 3),
  t_test = round(t_test_result$statistic, 4),
  p_value = round(t_test_result$p.value, 4),
  CI_lower = round(CI_lower, 2),
  CI_upper = round(CI_upper, 2),
  Estimated_Disease_Burden_2020 = round(disease_burden_2020, 2)
)

# Save results
write.csv(output_data, output_2_lyme, row.names = FALSE)

# Print parameters
cat("Estimated Parameters for Maryland (18°C):\n")
cat(sprintf("beta_H = %.4f\n", best_beta_H))
cat(sprintf("sigma_H = %.4f\n", best_sigma_H))
cat(sprintf("psi_H = %.4f\n", best_psi_H))
cat(sprintf("tau_H = %.4f\n", best_tau_H))
cat(sprintf("gamma_H = %.4f\n", best_gamma_H))
cat(sprintf("theta_H = %.4f\n", best_theta_H))
cat(sprintf("beta_M = %.4f\n", best_beta_M))
cat(sprintf("rhom = %.4f\n", best_rhom))
