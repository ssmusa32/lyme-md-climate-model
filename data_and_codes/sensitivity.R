

library(deSolve)
library(sensitivity)
library(ggplot2)
library(patchwork)

set.seed(123)

# R0 function
R0_lyme <- function(q_T, beta_NM, beta_NH, F_LM, F_NM, F_NH, h_N, mu_N, delta_M, mu_M, K_M, m_L, m_N, 
                    b_T, mu_L, mu_A, b_M, mu_H, Pi_H, sigma_H, g_H, tau_H, gamma_H, psi_H, m_A, F_AD) {
  
  R0_numerator <- sqrt(F_LM * ((m_N * mu_N + F_NM) * (h_N + 1) * mu_N +
                                 (F_NM * h_N + F_NH * m_N) + F_NH) *
                         beta_NM * F_NM * (delta_M + mu_M * K_M) * K_M * (m_L + 1) * q_T * (h_N + 1))
  
  R0_denominator <- ((m_N * mu_N + F_NM) + ((m_N * mu_N + F_NM) * h_N + F_NH * m_N)) + 
    ((h_N * mu_N + F_NH) + mu_N)
  
  R0_denominator <- R0_denominator * (delta_M + mu_M * K_M) * (m_L + 1)
  
  R0 <- R0_numerator / R0_denominator
  return(R0)
}

# Uncertainty analysis
generate_R0_samples <- function(n) {
  para_lyme <- data.frame(
    q_T = runif(n, 0.5, 2.0),
    beta_NM = runif(n, 0.5, 2.0),
    beta_NH = runif(n, 0.5, 2.0),
    F_LM = runif(n, 0.1, 0.5),
    F_NM = runif(n, 0.1, 0.5),
    F_NH = runif(n, 0.1, 0.5),
    h_N = runif(n, 20, 80),
    mu_N = runif(n, 0.02, 0.06),
    delta_M = runif(n, 0.01, 0.03),
    mu_M = runif(n, 0.001, 0.005),
    K_M = runif(n, 500000, 1500000),
    m_L = runif(n, 150, 250),
    m_N = runif(n, 1, 10),
    b_T = runif(n, 0.3, 0.5),
    mu_L = runif(n, 0.003, 0.008),
    mu_A = runif(n, 0.01, 0.03),
    b_M = runif(n, 0.01, 0.03),
    mu_H = runif(n, 0.00001, 0.0001),
    Pi_H = runif(n, 0.00001, 0.0001),
    sigma_H = runif(n, 0.01, 0.1),
    g_H = runif(n, 0.01, 0.3),
    tau_H = runif(n, 0.02, 0.2),
    gamma_H = runif(n, 0.02, 0.2),
    psi_H = runif(n, 0.0001, 0.01),
    m_A = runif(n, 200, 300),
    F_AD = runif(n, 0.1, 0.5)
  )
  
  R0_values <- apply(para_lyme, 1, function(row) {
    R0_lyme(row["q_T"], row["beta_NM"], row["beta_NH"], row["F_LM"], row["F_NM"], row["F_NH"],
            row["h_N"], row["mu_N"], row["delta_M"], row["mu_M"], row["K_M"], row["m_L"], row["m_N"], 
            row["b_T"], row["mu_L"], row["mu_A"], row["b_M"], row["mu_H"], row["Pi_H"], 
            row["sigma_H"], row["g_H"], row["tau_H"], row["gamma_H"], row["psi_H"], 
            row["m_A"], row["F_AD"])
  })
  
  return(R0_values)
}

sample_sizes <- seq(100, 1000, by = 100)
R0_all <- lapply(sample_sizes, generate_R0_samples)
names(R0_all) <- as.character(sample_sizes)

R0_df <- stack(R0_all)
colnames(R0_df) <- c("R0", "SampleSize")
R0_df$SampleSize <- factor(R0_df$SampleSize, levels = as.character(sample_sizes))

p_uncertainty <- ggplot(R0_df, aes(x = SampleSize, y = R0)) +
  geom_boxplot(outlier.colour = "red", color = "blue", fill = "white") +
  labs(
    x = "Runs",
    y = expression(bold("Basic reproduction number"~(R[0]))),
    title = expression(bold("(b)"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_uncertainty)

# PRCC analysis
n <- 1000
para_prcc <- data.frame(
  q_T = runif(n, 0.5, 2.0),
  beta_NM = runif(n, 0.5, 2.0),
  beta_NH = runif(n, 0.5, 2.0),
  F_LM = runif(n, 0.1, 0.5),
  F_NM = runif(n, 0.1, 0.5),
  F_NH = runif(n, 0.1, 0.5),
  h_N = runif(n, 20, 80),
  mu_N = runif(n, 0.02, 0.06),
  delta_M = runif(n, 0.01, 0.03),
  mu_M = runif(n, 0.001, 0.005),
  K_M = runif(n, 500000, 1500000),
  m_L = runif(n, 150, 250),
  m_N = runif(n, 1, 10),
  b_T = runif(n, 0.3, 0.5),
  mu_L = runif(n, 0.003, 0.008),
  mu_A = runif(n, 0.01, 0.03),
  b_M = runif(n, 0.01, 0.03),
  mu_H = runif(n, 0.00001, 0.0001),
  Pi_H = runif(n, 0.00001, 0.0001),
  sigma_H = runif(n, 0.01, 0.1),
  g_H = runif(n, 0.01, 0.3),
  tau_H = runif(n, 0.02, 0.2),
  gamma_H = runif(n, 0.02, 0.2),
  psi_H = runif(n, 0.0001, 0.01),
  m_A = runif(n, 200, 300),
  F_AD = runif(n, 0.1, 0.5)
)

R0_values_prcc <- apply(para_prcc, 1, function(row) {
  R0_lyme(row["q_T"], row["beta_NM"], row["beta_NH"], row["F_LM"], row["F_NM"], row["F_NH"],
          row["h_N"], row["mu_N"], row["delta_M"], row["mu_M"], row["K_M"], row["m_L"], row["m_N"], 
          row["b_T"], row["mu_L"], row["mu_A"], row["b_M"], row["mu_H"], row["Pi_H"], 
          row["sigma_H"], row["g_H"], row["tau_H"], row["gamma_H"], row["psi_H"], 
          row["m_A"], row["F_AD"])
})

prcc_result <- pcc(para_prcc, R0_values_prcc, nboot = 1000, rank = TRUE)

# Parameter labels
label_math <- c(
  q_T = "bold(q[T])", beta_NM = "bold(beta[NM])", beta_NH = "bold(beta[NH])",
  F_LM = "bold(F[LM])", F_NM = "bold(F[NM])", F_NH = "bold(F[NH])",
  h_N = "bold(h[N])", mu_N = "bold(mu[N])", delta_M = "bold(delta[M])", mu_M = "bold(mu[M])", 
  K_M = "bold(K[M])", m_L = "bold(m[L])", m_N = "bold(m[N])", b_T = "bold(b[T])",
  mu_L = "bold(mu[L])", mu_A = "bold(mu[A])", b_M = "bold(b[M])",
  mu_H = "bold(mu[H])", Pi_H = "bold(Pi[H])", sigma_H = "bold(sigma[H])",
  g_H = "bold(g[H])", tau_H = "bold(tau[H])", gamma_H = "bold(gamma[H])", 
  psi_H = "bold(psi[H])", m_A = "bold(m[A])", F_AD = "bold(F[AD])"
)

prcc_df <- data.frame(
  Parameter = names(label_math),
  PRCC = prcc_result$PRCC[, 1],
  Label = label_math
)
prcc_df <- prcc_df[order(abs(prcc_df$PRCC), decreasing = TRUE), ]

# PRCC plot
p_prcc <- ggplot(prcc_df, aes(x = PRCC, y = reorder(Label, PRCC))) +
  geom_col(fill = "darkblue", width = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = expression(bold("Partial Rank Correlation Coefficient (PRCC)")),
    y = "Parameter",
    title = expression(bold("(a)"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  scale_y_discrete(labels = function(x) parse(text = x))

print(p_prcc)

# Combine plots
combined_plot <- p_prcc + p_uncertainty
print(combined_plot)

