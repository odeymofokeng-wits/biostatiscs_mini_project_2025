#------------------------------------------------------------
# 0) Install & Load Packages
#------------------------------------------------------------
# Installing required packages
if (!requireNamespace("copula", quietly = TRUE)) install.packages("copula")
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("compound.Cox", quietly = TRUE)) install.packages("compound.Cox")
if (!requireNamespace("HAC", quietly = TRUE)) install.packages("HAC")
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot") # For bootstrap CIs

library(copula)       # For generating copula-based data
library(survival)     # For Kaplan-Meier and Cox models
library(ggplot2)      # For enhanced plotting
library(compound.Cox) # For copula-graphic estimator
library(HAC)          # For hierarchical Archimedean copula
library(boot)         # For bootstrap resampling

#------------------------------------------------------------
# 1) Simulate Dependent Censoring Data
#------------------------------------------------------------
set.seed(255438)  # For reproducibility
n <- 1500  # 1500 patients per group

# Control group: Clayton copula
cop_ctl <- claytonCopula(0.5)  # alpha = 0.5, tau = 0.2
U_ctl <- rCopula(n, cop_ctl)
T_ctl <- qexp(U_ctl[,1], rate = 0.1)  # Control event times
C_ctl <- qexp(U_ctl[,2], rate = 0.1)  # Control censoring times
X_ctl <- pmin(T_ctl, C_ctl)
delta_ctl <- as.integer(T_ctl <= C_ctl)

# Treatment group: TRUE HAC (nested Clayton) in 3D with a latent S (e.g., side effects)
# Structure: inner node couples (T, C) with tau_inner; outer node couples ((T,C), S) with tau_outer
# Indices: 1 = T, 2 = C, 3 = S
tau_inner <- 0.5                     # strong T–C dependence
theta_inner <- 2 * tau_inner / (1 - tau_inner)   # = 2
tau_outer <- 1/3                     # moderate link of S with (T,C)
theta_outer <- 2 * tau_outer / (1 - tau_outer)   # = 1
hac_trt <- onacopula("Clayton", C(theta_outer, 3, C(theta_inner, c(1, 2))))
U_trt3 <- rnacopula(n, hac_trt)
T_trt <- qexp(U_trt3[, 1], rate = 0.08)  # Treatment event times
C_trt <- qexp(U_trt3[, 2], rate = 0.10)  # Treatment censoring times
S_trt_bin <- qbinom(U_trt3[, 3], size = 1, prob = 0.30)  # Bernoulli(0.3) for...fun
X_trt <- pmin(T_trt, C_trt)
delta_trt <- as.integer(T_trt <= C_trt)

# Combining into data frame
df <- data.frame(
  X = c(X_ctl, X_trt),
  delta = c(delta_ctl, delta_trt),
  Z = rep(0:1, each = n)
)


#------------------------------------------------------------
# 2) Standard Kaplan-Meier and Cox Analysis
#------------------------------------------------------------
fit_km <- survfit(Surv(X, delta) ~ Z, data = df)
fit_cox <- coxph(Surv(X, delta) ~ Z, data = df)
cat("Cox Model Summary:\n")
print(summary(fit_cox))
fit_cox_strat <- coxph(Surv(X, delta) ~ strata(Z), data = df)
cat("Stratified Cox Model Summary:\n")
print(summary(fit_cox_strat))

# Plot KM curves (Figure 1)
true_times <- seq(0, 20, by = 0.1)
S0_true <- pexp(true_times, rate = 0.1, lower.tail = FALSE)
S1_true <- pexp(true_times, rate = 0.08, lower.tail = FALSE)
km_summary <- summary(fit_km, times = true_times, extend = TRUE)
S0_km <- km_summary$surv[km_summary$strata == "Z=0"]
S1_km <- km_summary$surv[km_summary$strata == "Z=1"]
plot_df <- data.frame(
  Time = rep(true_times, 4),
  Survival = c(S0_true, S1_true, S0_km, S1_km),
  Group = rep(c("Control (True)", "Treatment (True)", "Control (KM)", "Treatment (KM)"), each = length(true_times)),
  Method = rep(c("True", "True", "KM", "KM"), each = length(true_times))
)
ggplot(plot_df, aes(x = Time, y = Survival, color = Group, linetype = Method)) +
  geom_step(linewidth = 1) +
  labs(title = "Figure 1: Kaplan-Meier vs. True Survival Curves", x = "Time", y = "Survival Probability") +
  scale_color_manual(values = c("blue", "red", "blue", "red")) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid")) +
  theme_light() + theme(panel.background = element_rect(fill = "white")) + theme(legend.position = "top")
ggsave("Figure1_KM_vs_True.png")

#------------------------------------------------------------
# 3) Copula-Graphic Estimator with compound.Cox ( I used AI for debugging)
#------------------------------------------------------------
cg_estimator <- function(t.vec, d.vec, tau = NULL, alpha = NULL, S.plot = FALSE) {
  if (!exists("CG.Clayton")) stop("CG.Clayton() not found: please load compound.Cox.")
  if (is.null(alpha)) {
    if (is.null(tau)) stop("Either 'tau' or 'alpha' must be supplied.")
    if (tau < 0 || tau >= 1) stop("'tau' should be in [0,1).")
    alpha <- if (tau == 0) 0 else 2 * tau / (1 - tau)
  }
  out <- CG.Clayton(t.vec, d.vec, alpha = alpha, S.plot = S.plot)
  return(list(time = out$time, surv = out$surv))
}

tau_ctl <- 0.2
tau_trt <- 0.5
alpha_ctl <- 2 * tau_ctl / (1 - tau_ctl)
alpha_trt <- 2 * tau_trt / (1 - tau_trt)
cg_ctl_comp <- CG.Clayton(df$X[df$Z == 0], df$delta[df$Z == 0], alpha = alpha_ctl, S.plot = FALSE)
cg_trt_comp <- CG.Clayton(df$X[df$Z == 1], df$delta[df$Z == 1], alpha = alpha_trt, S.plot = FALSE)

# Plot copula-graphic vs. KM vs. true for treatment group (Figure 2)
S1_cg <- approx(cg_trt_comp$time, cg_trt_comp$surv, xout = true_times, method = "constant", rule = 2)$y
plot_df_fig2 <- data.frame(
  Time = rep(true_times, 3),
  Survival = c(S1_true, S1_km, S1_cg),
  Method = rep(c("True", "KM", paste0("Copula-Graphic (tau = ", tau_trt, ")")), each = length(true_times))
)
ggplot(plot_df_fig2, aes(x = Time, y = Survival, color = Method)) +
  geom_step(linewidth = 1) +
  labs(title = "Figure 2: Treatment Group - KM vs. Copula-Graphic vs. True Survival", x = "Time", y = "Survival Probability") +
  scale_color_manual(values = c("red", "purple", "black")) +
  theme_light() + theme(panel.background = element_rect(fill = "white")) + theme(legend.position = "top")
ggsave("Figure2_Treatment_CG_vs_KM.png")

# CG-adjusted Cox model
df_cg <- df
df_cg$X[df$Z == 0] <- approx(cg_ctl_comp$time, cg_ctl_comp$surv, xout = df$X[df$Z == 0], method = "constant", rule = 2)$y
df_cg$X[df$Z == 1] <- approx(cg_trt_comp$time, cg_trt_comp$surv, xout = df$X[df$Z == 1], method = "constant", rule = 2)$y
fit_cox_cg <- coxph(Surv(X, delta) ~ Z, data = df_cg)
cat("CG-Adjusted Cox Model Summary:\n")
print(summary(fit_cox_cg))

#------------------------------------------------------------
# 4) Sensitivity Analysis with Bootstrap CIs
#------------------------------------------------------------
tau_values <- c(0, 0.2, 0.5, 0.8)
t0 <- 10
S0_true_t0 <- pexp(t0, rate = 0.1, lower.tail = FALSE)  # ~0.368
S1_true_t0 <- pexp(t0, rate = 0.08, lower.tail = FALSE) # ~0.527

# Bootstrap function for bias
boot_bias <- function(data, indices, tau, z) {
  d <- data[indices, ]
  cg <- cg_estimator(d$X[d$Z == z], d$delta[d$Z == z], tau = tau)
  # Removing duplicates and ensure unique times
  cg_time_unique <- cg$time[!duplicated(cg$time)]
  cg_surv_unique <- cg$surv[!duplicated(cg$time)]
  S_cg <- approx(cg_time_unique, cg_surv_unique, xout = t0, method = "constant", rule = 2)$y
  return(S_cg - (if (z == 0) S0_true_t0 else S1_true_t0))
}

#  results
sensitivity_results <- data.frame(
  Tau = numeric(),
  Tau_display = numeric(),
  Method = character(),
  S0_Bias = numeric(),
  S1_Bias = numeric(),
  S0_MSE = numeric(),
  S1_MSE = numeric(),
  S0_IMSE = numeric(),
  S1_IMSE = numeric(),
  S0_Bias_LCI = numeric(),
  S0_Bias_UCI = numeric(),
  S1_Bias_LCI = numeric(),
  S1_Bias_UCI = numeric(),
  stringsAsFactors = FALSE
)

set.seed(123) 
for (tau in tau_values) {
  # Bootstrap for control
  boot_ctl <- boot(data = df[df$Z == 0, ], statistic = boot_bias, R = 100, tau = tau, z = 0)
  ci_ctl <- boot.ci(boot_ctl, type = "perc")$percent[4:5]
  # Bootstrap for treatment
  boot_trt <- boot(data = df[df$Z == 1, ], statistic = boot_bias, R = 100, tau = tau, z = 1)
  ci_trt <- boot.ci(boot_trt, type = "perc")$percent[4:5]
  
  cg_ctl_sens <- cg_estimator(df$X[df$Z == 0], df$delta[df$Z == 0], tau = tau)
  cg_trt_sens <- cg_estimator(df$X[df$Z == 1], df$delta[df$Z == 1], tau = tau)
  # Removing duplicates for main estimates
  cg_ctl_time_unique <- cg_ctl_sens$time[!duplicated(cg_ctl_sens$time)]
  cg_ctl_surv_unique <- cg_ctl_sens$surv[!duplicated(cg_ctl_sens$time)]
  cg_trt_time_unique <- cg_trt_sens$time[!duplicated(cg_trt_sens$time)]
  cg_trt_surv_unique <- cg_trt_sens$surv[!duplicated(cg_trt_sens$time)]
  S0_cg_sens <- approx(cg_ctl_time_unique, cg_ctl_surv_unique, xout = t0, method = "constant", rule = 2)$y
  S1_cg_sens <- approx(cg_trt_time_unique, cg_trt_surv_unique, xout = t0, method = "constant", rule = 2)$y
  
  S0_cg_curve <- approx(cg_ctl_time_unique, cg_ctl_surv_unique, xout = true_times, method = "constant", rule = 2)$y
  S1_cg_curve <- approx(cg_trt_time_unique, cg_trt_surv_unique, xout = true_times, method = "constant", rule = 2)$y
  S0_imse <- mean((S0_cg_curve - S0_true)^2, na.rm = TRUE)
  S1_imse <- mean((S1_cg_curve - S1_true)^2, na.rm = TRUE)
  
  sensitivity_results <- rbind(sensitivity_results, data.frame(
    Tau = tau,
    Tau_display = tau,
    Method = "Copula-Graphic",
    S0_Bias = S0_cg_sens - S0_true_t0,
    S1_Bias = S1_cg_sens - S1_true_t0,
    S0_MSE = (S0_cg_sens - S0_true_t0)^2,
    S1_MSE = (S1_cg_sens - S1_true_t0)^2,
    S0_IMSE = S0_imse,
    S1_IMSE = S1_imse,
    S0_Bias_LCI = ci_ctl[1],
    S0_Bias_UCI = ci_ctl[2],
    S1_Bias_LCI = ci_trt[1],
    S1_Bias_UCI = ci_trt[2],
    stringsAsFactors = FALSE
  ))
}

# Add KM results
km_summary_t0 <- summary(fit_km, times = t0, extend = TRUE)
S0_km_t0 <- km_summary_t0$surv[km_summary_t0$strata == "Z=0"]
S1_km_t0 <- km_summary_t0$surv[km_summary_t0$strata == "Z=1"]
S0_km_curve <- km_summary$surv[km_summary$strata == "Z=0"]
S1_km_curve <- km_summary$surv[km_summary$strata == "Z=1"]
S0_km_imse <- mean((S0_km_curve - S0_true)^2, na.rm = TRUE)
S1_km_imse <- mean((S1_km_curve - S1_true)^2, na.rm = TRUE)

# Bootstrap for KM
boot_km_ctl <- boot(data = df[df$Z == 0, ], statistic = boot_bias, R = 100, tau = 0, z = 0)
ci_km_ctl <- boot.ci(boot_km_ctl, type = "perc")$percent[4:5]
boot_km_trt <- boot(data = df[df$Z == 1, ], statistic = boot_bias, R = 100, tau = 0, z = 1)
ci_km_trt <- boot.ci(boot_km_trt, type = "perc")$percent[4:5]

sensitivity_results <- rbind(
  data.frame(
    Tau = NA,
    Tau_display = -0.02,
    Method = "KM",
    S0_Bias = S0_km_t0 - S0_true_t0,
    S1_Bias = S1_km_t0 - S1_true_t0,
    S0_MSE = (S0_km_t0 - S0_true_t0)^2,
    S1_MSE = (S1_km_t0 - S1_true_t0)^2,
    S0_IMSE = S0_km_imse,
    S1_IMSE = S1_km_imse,
    S0_Bias_LCI = ci_km_ctl[1],
    S0_Bias_UCI = ci_km_ctl[2],
    S1_Bias_LCI = ci_km_trt[1],
    S1_Bias_UCI = ci_km_trt[2],
    stringsAsFactors = FALSE
  ),
  sensitivity_results
)

cat("Sensitivity Analysis Results at t = 10:\n")
print(sensitivity_results)

#------------------------------------------------------------
# 5) HAC Analysis Comparison
#------------------------------------------------------------
fit_km_hac <- fit_km
S0_km_hac <- S0_km_t0
S1_km_hac <- S1_km_t0
cg_ctl_hac <- cg_estimator(df$X[df$Z == 0], df$delta[df$Z == 0], tau = tau_ctl)
cg_trt_hac <- cg_estimator(df$X[df$Z == 1], df$delta[df$Z == 1], tau = tau_trt)
S0_cg_hac <- approx(cg_ctl_hac$time, cg_ctl_hac$surv, xout = t0, method = "constant", rule = 2)$y
S1_cg_hac <- approx(cg_trt_hac$time, cg_trt_hac$surv, xout = t0, method = "constant", rule = 2)$y
results_hac <- data.frame(
  Method = c("KM (HAC)", "Copula-Graphic (HAC)"),
  Control_Bias = c(S0_km_hac - S0_true_t0, S0_cg_hac - S0_true_t0),
  Treatment_Bias = c(S1_km_hac - S1_true_t0, S1_cg_hac - S1_true_t0),
  Control_MSE = c((S0_km_hac - S0_true_t0)^2, (S0_cg_hac - S0_true_t0)^2),
  Treatment_MSE = c((S1_km_hac - S1_true_t0)^2, (S1_cg_hac - S1_true_t0)^2),
  stringsAsFactors = FALSE
)
cat("Bias and MSE with HAC at t = 10:\n")
print(results_hac)

#------------------------------------------------------------
# 6) Boxplot of Bias (Figure 4)
#------------------------------------------------------------
bias_df <- rbind(
  data.frame(Group = "Control", Bias = sensitivity_results$S0_Bias, Method = sensitivity_results$Method, Tau = sensitivity_results$Tau),
  data.frame(Group = "Treatment", Bias = sensitivity_results$S1_Bias, Method = sensitivity_results$Method, Tau = sensitivity_results$Tau),
  data.frame(Group = "Treatment", Bias = results_hac$Treatment_Bias, Method = results_hac$Method, Tau = NA)
)
ggplot(bias_df, aes(x = Method, y = Bias, fill = Group)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 4: Bias in Survival Estimates at t=10", x = "Method", y = "Bias (S(10) - True S(10))") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  theme_light() + theme(panel.background = element_rect(fill = "white")) + theme(legend.position = "top")
ggsave("Figure4_Bias_Boxplot.png")