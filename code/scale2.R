### script to run two experiments on scale1 to show the possibility of dispersal being harmful
output_dir <- here("output", "figures", "note3", "scale2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# m: proportion of individuals moving to the other patch after local growth
exchange <- function(growth, m) {
  stopifnot(length(growth) == 2L, 0 <= m, m <= 0.5)
  c(A = (1 - m) * growth[1L] + m * growth[2L],    # number of individuals migrate from A to B is in proportion of population in A
    B = m * growth[1L] + (1 - m) * growth[2L])    # vice versa
}

##################### Experiment 1: a prolonged deterioration creates a sink in B #####################
# Beverton-Holt, no Allee effect
logistic <- function(N, r, K, allee_threshold) {
  r * N / (1 + (r - 1) * N / K)
}

# strong Allee effect. match the Beverton-Holt slope at K, while fixing the Allee threshold at 10
cubic_allee <- function(N, r, K, allee_threshold) {
  eta <- (1 - 1 / r) / (K / allee_threshold - 1)
  N + eta * N * (1 - N / K) * (N / allee_threshold - 1)
}

run_press <- function(m, growth_fun, press_steps = 20L, after_steps = 60L,   # patch B remains unsuitable during the first 20 steps
                      K = 100, r = 1.25, press_factor = 0.2,
                      allee_threshold = 10) {
  stopifnot(press_steps >= 1L, after_steps >= 1L, press_factor >= 0, press_factor < 1, r > 1, 0 < allee_threshold, allee_threshold < K)
  n_steps <- press_steps + after_steps
  N <- matrix(NA_real_, nrow = n_steps + 1L, ncol = 2L, dimnames = list(NULL, c("A", "B")))
  N[1L, ] <- c(K, K)
  for (t in seq_len(n_steps)) {
    growth <- growth_fun(N[t, ], r = r, K = K, allee_threshold = allee_threshold)
    if (t <= press_steps) {
      growth[2L] <- press_factor * growth[2L]      # the disturbed population is multiplied by press_factor in each disturbance year
    }
    N[t + 1L, ] <- exchange(growth, m)
  }
  data.frame(time = 0:n_steps, A = N[, 1L], B = N[, 2L],
             total = rowSums(N), m = m)
}

press_steps <- 20L
m_press <- c(0, 0.08, 0.35)
press_cols <- c("#9ECAE1", "#4292C6", "#08519C")
growth_models <- list(no_allee = logistic, strong_allee = cubic_allee)
K <- 100
r <- 1.25
allee_threshold <- 10    # where per capita growth is 0. growth rate lower than that is negative.
### for total population 
# N_grid <- seq(0, 1.2 * K, length.out = 500L)
# rate_ylim <- range(c(logistic(N_grid, r, K, allee_threshold) - N_grid, cubic_allee(N_grid, r, K, allee_threshold) - N_grid)) 
### for per capita
N_grid <- seq(1e-6, 1.2 * K, length.out = 500L)
rate_ylim <- range(c((logistic(N_grid, r, K, allee_threshold) - N_grid) / N_grid, (cubic_allee(N_grid, r, K, allee_threshold) - N_grid) / N_grid))

plot_press <- function(press_runs, growth_fun, file, plot_title) {
  pdf(file.path(output_dir, file), width = 10, height = 7.5)
  on.exit(dev.off())
  par(mfrow = c(2, 2), mar = c(5.5, 5.5, 3.5, 1.5), oma = c(0, 0, 2, 0), cex.lab = 1.5)
  
  for (variable in c("A", "B", "total")) {
    ymax <- if (variable == "total") 210 else 110
    plot(0, 0, type = "n",
         xlim = range(press_runs[[1L]]$time),
         ylim = c(0, ymax),
         xlab = "time step", ylab = "population size",
         main = switch(variable,
                       A = "productive patch A",
                       B = "deteriorated patch B",
                       total = "metapopulation"))
    abline(v = press_steps, lty = "dotdash", col = "black", lwd = 0.7)
    
    for (j in seq_along(press_runs)) {
      x <- press_runs[[j]]
      lines(x$time, x[[variable]], col = press_cols[j], lwd = 2)
    }
    
    if (variable == "total") {
      abline(h = 200, lty = "longdash", col = "black", lwd = 0.7)
    }
  }
  
  ### for total population 
  # growth_rate <- growth_fun(N_grid, r, K, allee_threshold) - N_grid
  # plot(N_grid, growth_rate, type = "l", col = "black", lwd = 2,
  #      ylim = rate_ylim, xlab = "Population size N",
  #      ylab = "proliferation rate (N[t+1] - N[t])", main = "Local growth without press")
  ### for per capita
  growth_rate <- (growth_fun(N_grid, r, K, allee_threshold) - N_grid) / N_grid
  plot(N_grid, growth_rate, type = "l", col = "black", lwd = 2,
       ylim = rate_ylim, xlab = "population size",
       ylab = "per-capita growth rate", main = "local growth without press")
  abline(h = 0, lty = 3, col = "grey55")
  abline(v = K, lty = 3, col = "grey55")
  if (identical(growth_fun, cubic_allee)) {
    abline(v = allee_threshold, lty = 2, col = "grey55")
  }
    legend("topright", legend = names(press_runs), col = press_cols, lty = 1, lwd = 2, bty = "n", cex = 1.3, title = "dispersal rate per time step")
  mtext(plot_title, outer = TRUE, side = 3, line = 0.5)
}

for (model in names(growth_models)) {
  press_runs <- lapply(m_press, run_press,
                       growth_fun = growth_models[[model]],
                       press_steps = press_steps, K = K, r = r,
                       allee_threshold = allee_threshold)
  names(press_runs) <- paste0("m = ", m_press)
  plot_press(press_runs, growth_models[[model]], paste0("press_", model, ".pdf"),
             if (model == "no_allee") "Without Allee effect"
             else "With strong Allee effect")
  cat("\n", model, "\n", sep = "")
  for (x in press_runs) {
    end <- x[x$time == press_steps, ]
    recovered <- x$time[x$time > press_steps & x$total >= 190]  # the time after disturbance for the total population to return to 95% capacity
    cat("m =", x$m[1L], ": end of press =", round(end$total, 2),
        "; steps to 95% =",
        if (length(recovered)) recovered[1L] - press_steps else NA_integer_,
        "\n")
  }
}



##################### Experiment 2: migration synchronizes initially offset cycles #####################
# use Ricker dynamics to generate inherent oscillation
ricker <- function(N, r, K) {
  N * exp(r * (1 - N / K))
}

run_cycles <- function(m, n_steps = 20L, r = 2.4, K = 100,   # 2<r<2.52 generates sustained two-cycle oscillations
                       initial = c(A = 50, B = 180)) {
  stopifnot(n_steps >= 1L, length(initial) == 2L, all(initial > 0))
  N <- matrix(NA_real_, nrow = n_steps + 1L, ncol = 2L, dimnames = list(NULL, c("A", "B")))
  N[1L, ] <- initial
  for (t in seq_len(n_steps)) {
    growth <- ricker(N[t, ], r = r, K = K)
    N[t + 1L, ] <- exchange(growth, m)
  }
  data.frame(time = 0:n_steps, A = N[, 1L], B = N[, 2L], total = rowSums(N), m = m)
}

regional_threshold <- 100      # illustrative risk threshold
patch_cols <- c(A = "#8E7CC3", B = "#C99A6B")

plot_settings <- list(
  list(m_cycles = c(0.1, 0.2), n_steps = 20L, file = "cycles_synchrony_20.pdf"),          # to show the effect of low and high migration rate on population synchronization
  list(m_cycles = c(0.105, 0.106), n_steps = 300L, file = "cycles_synchrony_3000.pdf"))   # to show the critical point of population synchronization

for (settings in plot_settings) {
  cycles_runs <- lapply(settings$m_cycles, run_cycles, n_steps = settings$n_steps)
  names(cycles_runs) <- paste0("m = ", settings$m_cycles)
  pdf(file.path(output_dir, settings$file), width = 10, height = 7.5)
  par(mfrow = c(2, 2), mar = c(5.5, 5.5, 3.5, 1.5), oma = c(0, 0, 1.6, 0), cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.2)
  
  for (x in cycles_runs) {
    plot(x$time, x$A, type = "l", lwd = 2, col = patch_cols["A"], ylim = c(0, 190),
         xlab = "time step", ylab = "local population size", main = paste0("patches; m = ", x$m[1L]))
    lines(x$time, x$B, col = patch_cols["B"], lwd = 2)
    if (x$m[1L] == settings$m_cycles[1L]) {
      legend("topright", legend = c("patch A", "patch B"), col = patch_cols, lty = 1, lwd = 2, bty = "n", cex = 1.15)
    }
  }
  
  for (x in cycles_runs) {
    plot(x$time, x$total, type = "l", lwd = 2, col = "black", ylim = c(0, 350),
         xlab = "time step", ylab = "metapopulation size", main = paste0("metapopulation; m = ", x$m[1L]))
    abline(h = regional_threshold, col = "#b2182b", lty = 2, lwd = 1.5)
    if (x$m[1L] == settings$m_cycles[1L]) {
      legend("topright", legend = c("metapopulation size", "illustrative risk threshold"),
             col = c("black", "#b2182b"), lty = c(1, 2), lwd = 1.5,  bty = "n", cex = 1.05)
    }
  }
  dev.off()
}

message("Figures saved to ", normalizePath(output_dir))

