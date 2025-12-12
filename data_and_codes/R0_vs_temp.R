

library(ggplot2)
library(cowplot)
library(patchwork)

# Temperature-dependent functions
F_T <- function(T) {
  ifelse(T > 6 & T < 28, -24.58678 * T^2 + 835.9505 * T - 4105.579, 0)
}

s_E_T <- function(T) {
  T_min <- 6; T_max <- 28; s_min <- 0.1; s_max <- 0.8
  ifelse(T >= T_min & T <= T_max,
         s_min + (s_max - s_min) * (T - T_min) / (T_max - T_min), 0)
}

b_T <- function(T) { F_T(T) * s_E_T(T) }

a <- 0.0013; b <- 0.515; M <- 20000; H <- 600; DL <- 12

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

mu_A <- function(T) { ifelse(T > 10 & T < 30, (-log(-0.000828 * T^2 + 0.0367 * T + 0.522)), 0) }

b_M <- function(T) { ifelse(T > 10 & T < 25, (0.3 * exp(-((T - 18.5) / (25 - 15))^2)), 0) }

K_M <- function(T) { ifelse(T > 7 & T < 30, 150 * (1 - ((T - 18.5) / 23)^2), 0) }

mu <- function(T) { ifelse(T > 5 & T < 30, (0.012 + 0.4 * ((T - 18.5) / 25)^2), 0) }

# R0 calculation
calculate_R0 <- function(T) {
  L <- 1e5; Mt <- 100; Nt <- 1e2; Nh <- 6e3
  FLM <- F_LM(T); FNM <- F_NM(T); FNH <- F_NH(T)
  muN <- mu_N(T); d <- mu(T); KM <- K_M(T); bM <- b_M(T)
  beta_L <- 1; beta_M <- 0.75; d_M <- 0.015; m_L <- 27.8; m_N <- 4.5; h_N <- 2; s <- 1
  
  numerator <- sqrt(
    FLM *
      ((Mt * m_N + Nt) * (Nh * h_N + Nt) * muN +
         Nh * Mt * FNM * h_N +
         (Nh * FNH * m_N + Nt * FNM) * Mt + Nh * Nt * FNH) *
      beta_M * Mt * FNM *
      (Mt * d_M + d * KM) * KM *
      (Mt * m_L + L) * beta_L * L * (Nh * h_N + Nt)
  )
  
  denominator <- ((m_N * muN + FNM) * Nt +
                    ((m_N * muN + FNM) * h_N + FNH * m_N) * Nh) * Mt +
    ((h_N * muN + FNH) * Nh + muN * Nt) * Nt
  denominator <- denominator * (Mt * d_M + d * KM) * (Mt * m_L + L)
  
  R0 <- s * numerator / denominator
  return(ifelse(is.finite(R0) && R0 >= 0, R0, 0))
}

# Gamma calculation
calculate_Gamma <- function(T) {
  m_A <- 20
  b <- b_T(T)
  FAD <- theta_L(T)
  muL <- mu_L(T)
  Gamma <- (b * FAD) / (m_A * muL)
  return(ifelse(is.finite(Gamma) && Gamma >= 0, Gamma, NA))
}

# r0M calculation
calculate_r0M <- function(T) {
  bm <- b_M(T)
  mum <- mu(T)
  val <- ifelse(mum > 0, bm / mum, NA)
  return(val)
}

# Generate data
temp_values <- seq(0, 35, by = 0.5)

r0_data <- data.frame(
  Temperature = temp_values,
  R0 = sapply(temp_values, calculate_R0)
)
r0_data <- r0_data[r0_data$Temperature >= 5 & r0_data$Temperature <= 35 & r0_data$R0 <= 2, ]

gamma_data <- data.frame(
  Temperature = temp_values,
  Gamma = sapply(temp_values, calculate_Gamma)
)
gamma_data <- subset(gamma_data, Temperature >= 0 & Temperature <= 35 & !is.na(Gamma))

r0M_data <- data.frame(
  Temperature = temp_values,
  r0M = sapply(temp_values, calculate_r0M)
)
r0M_data <- subset(r0M_data, Temperature >= 5 & Temperature <= 35 & !is.na(r0M))

# Plotting function
create_plot <- function(data, y_var, color, y_lab, y_breaks = NULL, x_start = 5) {
  y_col <- data[[y_var]]
  
  p <- ggplot(data, aes_string(x = "Temperature", y = y_var)) +
    annotate("segment", x = x_start, xend = 35, y = 0, yend = 0, size = 1.4, color = "black") +
    annotate("segment", x = x_start, xend = x_start, y = 0, yend = max(y_col), size = 1.4, color = "black") +
    geom_ribbon(data = subset(data, y_col > 1),
                aes(ymin = 1, ymax = y_col), fill = "white", alpha = 0.15) +
    geom_line(color = color, linewidth = 1.8) +
    geom_point(color = color, fill = color, shape = 21, size = 3, stroke = 1.2) +
    labs(
      x = expression(bold("Mean monthly temperature") * " (°C)"),
      y = y_lab
    ) +
    scale_x_continuous(breaks = seq(x_start, 35, by = 5), limits = c(x_start, 35), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = y_breaks) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(face = "bold", size = 27),
      axis.text = element_text(face = "bold", size = 26),
      axis.ticks = element_line(color = "black", linewidth = 1.2),
      plot.margin = margin(10, 15, 10, 15)
    )
  return(p)
}

# Create individual plots
gamma_plot <- create_plot(gamma_data, "Gamma", "darkgreen",
                          expression(bold("Ticks reproduction number ("~Gamma(T[A])~")")),
                          NULL, 0)

r0M_plot <- create_plot(r0M_data, "r0M", "darkblue",
                        expression(bold("Mice reproduction number ("~r[0*M](T[A])~")")),
                        NULL, 5)

r0_plot <- create_plot(r0_data, "R0", "black",
                       expression(bold("Basic reproduction number ("~R[0](T[A])~")")),
                       seq(0, 2.25, by = 0.5), 5)

# Add labels
add_label <- function(plot, label_text) {
  ggdraw() +
    draw_plot(plot) +
    draw_label(label_text,
               x = 0.05, y = 0.95, hjust = -1.5, vjust = -0.5,
               fontface = "bold", size = 20)
}

gamma_plot_labeled <- add_label(gamma_plot, "(a) Ticks")
r0M_plot_labeled <- add_label(r0M_plot, "(b) Mice")
r0_plot_labeled <- add_label(r0_plot, "(c) Humans")

# Combine plots
combined_panel <- plot_grid(
  gamma_plot_labeled,
  r0M_plot_labeled,
  r0_plot_labeled,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

# Print
print(combined_panel)

# Save (optional)
# ggsave("results/figures/Combined_R0_Gamma_r0M_Temp_Panel.png", 
#        combined_panel, width = 18, height = 6.5, dpi = 300)

