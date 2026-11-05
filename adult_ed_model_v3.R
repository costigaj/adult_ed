# =============================================================================
# Adult Ed Enrollment Flow Model — v3
#
# No weekly history is available, so rates are set from program knowledge.
#
#   1. Rates are set to a "low staff outreach" baseline.
#   2. Completion is scaled by completion_scale. By default it is fixed
#      (use_calibration = FALSE) so the model predicts the year-end total.
#      Optionally it can be solved so the low-staff run hits a known total.
#   3. Staff outreach works through two levers:
#        a. ns_start_pct: the % of non-starters who come to class within
#           their 90-day window (set directly as a percentage)
#        b. staff_mult: a multiplier on the other contact-driven rates
#           - return of dropped students (rho_dp, rho_dr)
#           - infrequent -> regular re-engagement (rho)
#      We don't yet know how big the staff effect is, so results are shown
#      across a range of scenarios.
#
# Time unit: weeks. Week 0 = start_date.
# =============================================================================

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- 0. Inputs you should set ---------------------------------------------------
typical_total  <- 24   # long-run average year-end total

# Calibration switch:
#   FALSE = keep rates fixed and let the model predict the year-end total
#   TRUE  = solve completion_scale so the low-staff run hits no_staff_total
use_calibration     <- FALSE
fixed_completion    <- 0.614   # completion_scale used when use_calibration = FALSE
no_staff_total      <- 26      # only used when use_calibration = TRUE

# Staff interaction scenarios
#   ns_start_pct = % of non-starters who come to class within their 90-day window
#   staff_mult   = boost to dropped-student returns and infrequent -> regular
staff_scenarios <- tibble::tribble(
  ~scenario,          ~staff_mult, ~ns_start_pct,
  "Low staff",         1.0,         0.01,
  "Moderate staff",    1.5,         0.10,
  "Good staff",        2.0,         0.20
)
ns_pct_low_staff <- staff_scenarios$ns_start_pct[1]

# ---- 1. Calendar --------------------------------------------------------------
start_date <- as.Date("2026-11-05")
end_date   <- as.Date("2027-06-30")
n_weeks    <- as.numeric(end_date - start_date) / 7

wk <- function(d) as.numeric(as.Date(d) - start_date) / 7
t_winter_slowdown <- wk("2026-11-30")
t_spring_pickup   <- wk("2027-01-25")

# ---- 2. Model -----------------------------------------------------------------
adult_ed_model_v3 <- function(t, y, p) {
  with(as.list(c(y, p)), {

    beta <- if (t < t_winter_slowdown) beta_fall
            else if (t < t_spring_pickup) beta_winter
            else beta_spring

    k <- 3 / weeks_to_12

    # Non-starter activation: convert "% who start within the window" into a
    # weekly rate. With expiry rate 1/window, share starting = g / (g + 1/window).
    g_ns  <- ns_start_pct / (1 - ns_start_pct) / ns_window

    # Other staff-affected (contact-driven) rates
    r_dp  <- rho_dp   * staff_mult
    r_dr  <- rho_dr   * staff_mult
    r_ir  <- rho      * staff_mult

    # Calibrated completion rates
    th_R <- theta_R * completion_scale
    th_I <- theta_I * completion_scale

    ns_activate <- g_ns * NS
    ns_expire   <- NS / ns_window
    back_pre    <- r_dp * Dp
    back_post   <- r_dr * Dr

    dS  <- beta - (alpha_start + alpha_noshow) * S
    dNS <- alpha_noshow * S - ns_activate - ns_expire
    dX  <- ns_expire

    dP1 <- alpha_start * S + ns_activate + back_pre - (k + delta_p) * P1
    dP2 <- k * P1 - (k + delta_p) * P2
    dP3 <- k * P2 - (k + delta_p) * P3

    dR  <- k * P3 + r_ir * I + back_post - (omega + delta_r + th_R) * R
    dI  <- omega * R - (r_ir + delta_i + th_I) * I

    dDp <- delta_p * (P1 + P2 + P3) - back_pre
    dDr <- delta_r * R + delta_i * I - back_post

    dA  <- th_R * R + th_I * I

    dRecruits   <- beta
    dNS_started <- ns_activate
    dReturned   <- back_pre + back_post

    list(c(dS, dNS, dX, dP1, dP2, dP3, dR, dI, dDp, dDr, dA,
           dRecruits, dNS_started, dReturned))
  })
}

# ---- 3. Initial conditions (Nov 5) ----------------------------------------------
state <- c(
  S  = 20, NS = 50, X = 0,
  P1 = 26, P2 = 26, P3 = 25,     # 77 pre-12
  R  = 38, I  = 63,
  Dp = 0,  Dr = 0,               # <-- current dropped-but-returnable counts
  A  = 3,
  Recruits = 0, NS_started = 0, Returned = 0
)
people <- c("S", "NS", "X", "P1", "P2", "P3", "R", "I", "Dp", "Dr", "A")

# ---- 4. Parameters: "no staff outreach" baseline (weekly rates) ----------------
params <- c(
  beta_fall    = 3.5,
  beta_winter  = 1.0,
  beta_spring  = 4.0,
  alpha_start  = 0.15,
  alpha_noshow = 0.08,
  ns_start_pct = ns_pct_low_staff,  # % of non-starters who start in window
  ns_window    = 90 / 7,
  weeks_to_12  = 3,
  delta_p      = 0.15,
  omega        = 0.12,
  rho          = 0.06,     # infrequent -> regular, unprompted
  delta_r      = 0.01,
  delta_i      = 0.05,
  rho_dp       = 0.01,     # dropped (pre-12) returning, unprompted
  rho_dr       = 0.015,    # dropped (post-12) returning, unprompted
  theta_R      = 0.012,
  theta_I      = 0.002,
  completion_scale = 1,    # solved below
  staff_mult       = 1     # 1 = no staff outreach
)

# ---- 5. Helpers -------------------------------------------------------------------
run_model <- function(p = params, y0 = state, times = seq(0, n_weeks, by = 1)) {
  out <- as.data.frame(ode(y = y0, times = times, func = adult_ed_model_v3,
                           parms = p, method = "lsoda", hmax = 0.25))
  out$date <- start_date + out$time * 7
  last <- tail(out, 1)
  stopifnot(abs(sum(last[people]) - sum(y0[people]) - last$Recruits) < 1e-4)
  out
}

year_end_A <- function(p) tail(run_model(p, times = c(0, n_weeks)), 1)$A

set_p <- function(p, ...) { ch <- c(...); p[names(ch)] <- ch; p }

# ---- 6. Calibrate: no-staff run must hit last year's total ---------------------
calibrate <- function(target, p = params) {
  uniroot(function(cs) year_end_A(set_p(p, completion_scale = cs, staff_mult = 1,
                                        ns_start_pct = ns_pct_low_staff)) - target,
          interval = c(0.05, 5))$root
}

if (use_calibration) {
  cs <- calibrate(no_staff_total)
  cat(sprintf("Calibrated completion_scale = %.3f (low-staff run hits %d)\n\n",
              cs, no_staff_total))
} else {
  cs <- fixed_completion
  cat(sprintf("Fixed completion_scale = %.3f (no calibration)\n\n", cs))
}
params["completion_scale"] <- cs

# ---- 7. Staff-effect scenarios ----------------------------------------------------
summarise_run <- function(out, label) {
  last <- tail(out, 1)
  tibble(
    scenario          = label,
    outcomes_total    = last$A,
    new_outcomes      = last$A - state[["A"]],
    ns_who_started    = last$NS_started,
    ns_expired        = last$X,
    returned_events   = last$Returned,
    dropped_now       = last$Dp + last$Dr,
    active_june30     = last$P1 + last$P2 + last$P3 + last$R + last$I
  ) |>
    mutate(across(where(is.numeric), \(x) round(x, 1)))
}

scenarios <- staff_scenarios |>
  mutate(scenario = sprintf("%s (%d%% of NS start)", scenario, round(ns_start_pct * 100)))

run_scenario <- function(i) {
  run_model(set_p(params, staff_mult   = scenarios$staff_mult[i],
                          ns_start_pct = scenarios$ns_start_pct[i]))
}

results <- bind_rows(lapply(seq_len(nrow(scenarios)), \(i)
  summarise_run(run_scenario(i), scenarios$scenario[i])))

cat("--- Year-End Projection (June 30, 2027) ---\n")
print(results, width = Inf)

# ---- 8. Sensitivity to the calibration anchor (calibration mode only) ----------
if (use_calibration) {
  anchor_check <- expand.grid(anchor = no_staff_total + c(-2, 0, 2), staff_mult = c(1, 2)) |>
    rowwise() |>
    mutate(outcomes = round(year_end_A(set_p(params,
                                              completion_scale = calibrate(anchor),
                                              staff_mult = staff_mult)), 1)) |>
    ungroup()
  cat("\n--- Sensitivity to last year's total ---\n")
  print(anchor_check)
}

# ---- 9. Plot: cumulative outcomes by scenario -------------------------------
bind_rows(lapply(seq_len(nrow(scenarios)), \(i)
  run_scenario(i) |> transmute(date, A, scenario = scenarios$scenario[i]))) |>
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) |>
  ggplot(aes(date, A, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = typical_total, linetype = "dashed", colour = "grey50") +
  annotate("text", x = start_date + 7, y = typical_total + 0.8,
           label = "Average year", hjust = 0, colour = "grey40", size = 3.5) +
  labs(title = "Projected cumulative outcomes through June 30",
       subtitle = if (use_calibration) "Calibrated so the low-staff run matches last year"
                  else "Fixed rates (no calibration)",
       x = NULL, y = "Outcomes achieved", colour = NULL) +
  theme_minimal()
