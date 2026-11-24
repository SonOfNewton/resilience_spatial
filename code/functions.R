##########################  functions for note 2 plots ########################## 
################# for figure 2 #################
plot_pcf <- function( result, collision = NULL, curve_colour = "black") {
  ok <- is.finite(result$trans)
  upper_limit <- max(2, quantile(result$trans[ok], 0.99, na.rm = TRUE))
  plot(result$r, result$trans, type = "n", xlim = c(0, 20), ylim = c(0, upper_limit), xlab = "Distance r", ylab = expression(g(r)), main = "")
  
  # plot collision distance
  if (!is.null(collision)) {
    rect(xleft = 0, ybottom = 0, xright = collision, ytop = upper_limit, col = adjustcolor("grey70", alpha.f = 0.3), border = NA)
  }
  # CSR baseline
  abline(h = 1, col = "black", lty = 2, lwd = 2)
  
  # empirical PCF
  lines( result$r, result$trans, col = curve_colour, lwd = 2)
  
  if (!is.null(collision)) {
    abline(v = collision, col = "#2166AC", lty = 3, lwd = 2)
  }
  legend("topright", inset = c(0, 0),
  legend = c("Empirical PCF", "CSR: g(r) = 1"), col = c(curve_colour, "black"), lty = c(1, 2), lwd = 2, bty = "n"
  )
}

################# for figure 3 #################
plot_population <- function(marks, title_text) {
  plot( c(0, side), c(0, side), type = "n", asp = 1, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  box()
  title(main = title_text, line = 1)
  points( xy[marks == "survive", , drop = FALSE], pch = 16, col = "grey45", cex = 0.8)
  points( xy[marks == "removed", , drop = FALSE], pch = 16, col = "orange", cex = 1)
}

################# for figure 4, 5, 6 #################
standardize_field <- function(r) {
  mu <- terra::global(r, "mean", na.rm = TRUE)[1, 1]
  sd_r <- terra::global(r, "sd", na.rm = TRUE)[1, 1]
  (r - mu) / sd_r
}

# break the Eurasia map down into several regions for stratified random sampling
region_of <- function(lon, lat) {
  region <- rep("Other", length(lon))
  region[
    lon >= -12 & lon <= 40 & lat >= 35 & lat <= 72
  ] <- "Europe"
  region[
    lon >= 100 & lon <= 140 & lat >= 20 & lat <= 55
  ] <- "EastAsia"
  region[
    lon >= 40 & lat >= 60
  ] <- "NorthernRussia"
  south_asia <- (
    lon >= 68 & lon <= 90 & lat >= 5 & lat <= 30
  )
  southeast_asia <- (
    lon >= 95 & lon <= 141 & lat >= -12 & lat <= 22
  )
  region[south_asia | southeast_asia] <- "SouthSoutheastAsia"

  region
}

# lat long -> X Y Z (in km)
to_sphere <- function(lon, lat) {
  a <- lon * pi / 180
  b <- lat * pi / 180
  data.frame(
    X = 6371 * cos(b) * cos(a),
    Y = 6371 * cos(b) * sin(a),
    Z = 6371 * sin(b)
  )
}

scores <- function(y, pred) {
  c(RMSE = sqrt(mean((y - pred)^2)), MAE = mean(abs(y - pred)))
}

# main function for plots
draw_map <- function(r, colors, unit, filename = NULL, limits = NULL, title = "", overlay = c("plain", "species"), circle_radii = NULL) {
  overlay <- match.arg(overlay)    # detailed plot with the five examples or plain plot without that emphasis
  if (!is.null(filename)) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    png(filename, width = 2800, height = 1500, res = 220)
    on.exit(dev.off(), add = TRUE)
  }
  terra::plot(r, col = colors, range = limits, main = title, xlab = "", ylab = "", mar = c(3, 3, 2.5, 6), 
              cex.main = 1.25, pax = list(cex.axis = 1.25), plg = list(title = unit, cex = 1.25, title.cex = 1.25))
  
  border_col <- "grey35"
  terra::lines(eurasia, col = border_col, lwd = 0.4)
  
  if (overlay == "plain") {
    # minimal notation of samples
    points(dat$lon, dat$lat, pch = 21, bg = "white", col = "black", cex = 0.55)
    legend_labels <- "100 observed sites"
    legend_colors <- "white"
  } else {
    # detailed notation of samples: colored according to species and annotate the examples
    point_colors <- unname(species_colors[as.character(dat$S)])
    idx <- match(labels_xy$population, dat$population)
    stopifnot(!anyNA(point_colors), !anyNA(idx))
    focus <- dat[idx, ]
    points(dat$lon, dat$lat, pch = 21, bg = "white", col = "white", cex = 1.65)
    points(dat$lon, dat$lat, pch = 21, bg = point_colors, col = "grey15", cex = 1.2, lwd = 0.65)
    end_lon <- focus$lon - 0.07 * (focus$lon - labels_xy$label_lon)
    end_lat <- focus$lat - 0.07 * (focus$lat - labels_xy$label_lat)
    arrows(labels_xy$label_lon, labels_xy$label_lat, end_lon, end_lat, length = 0.12, angle = 25, col = "white", lwd = 4)
    arrows(labels_xy$label_lon, labels_xy$label_lat, end_lon, end_lat,length = 0.12, angle = 25, col = "black", lwd = 1.4)
    points(focus$lon, focus$lat, pch = 21, bg = "white", col = "white", cex = 2.1)
    points(focus$lon, focus$lat, pch = 21, bg = point_colors[idx], col = "black", cex = 1.65, lwd = 1)
    half_w <- strwidth(labels_xy$population, cex = 1.5, font = 2) * 0.85
    half_h <- strheight("A", cex = 1.5, font = 2) * 0.85
    rect(
      labels_xy$label_lon - half_w, labels_xy$label_lat - half_h,
      labels_xy$label_lon + half_w, labels_xy$label_lat + half_h,
      col = "white", border = NA
    )
    text(labels_xy$label_lon, labels_xy$label_lat, labels = labels_xy$population, cex = 1.5, font = 2)
    legend_labels <- c(
      "species 1 (5 kg)",
      "species 2 (10 kg)",
      "species 3 (20 kg)"
    )
    legend_colors <- unname(species_colors)
  }
  
  if (is.null(circle_radii)) {
    legend("bottomleft", inset = 0.1, legend = legend_labels, pch = 21, pt.bg = legend_colors, 
           col = "black", pt.cex = 1.4, cex = 1.05, bg = "white", box.col = "grey65")
  } else {
    stopifnot(overlay == "plain", length(circle_radii) == 2)
    draw_range_circles(circle_radii[1], circle_radii[2])
  }
}

# INLA helper function: create empty prediction row
prediction_stack <- function(A, T_values, tag) {
  n <- nrow(A)
  stopifnot(length(T_values) == n)
  INLA::inla.stack(
    data = list(R = rep(NA_real_, n)),
    A = list(A, 1),
    effects = list(w_index, data.frame(Intercept = rep(0, n), T = T_values, S = rep(0, n))),
    tag = tag
  )
}

# INLA helper function: recover the posterior of the linear predictor
stack_summary <- function(tag, expected_rows) {
  rows <- INLA::inla.stack.index(stack_all, tag = tag)$data
  ans <- fit_inla$summary.linear.predictor[rows, , drop = FALSE]
  stopifnot(nrow(ans) == expected_rows)
  ans
}

# INLA helper function: link values to map cell ids
to_map <- function(x, label) {
  stopifnot(length(x) == n_map)
  out <- T_map
  v <- rep(NA_real_, terra::ncell(out))
  v[map_sites$cell] <- x
  terra::values(out) <- v
  names(out) <- label
  out
}

# INLA helper function: summarise posterior density
summarise_marginal <- function(m) {
  mu <- INLA::inla.emarginal(function(x) x, m)
  second <- INLA::inla.emarginal(function(x) x^2, m)
  c(
    mean = mu,
    sd = sqrt(max(0, second - mu^2)),
    q025 = INLA::inla.qmarginal(0.025, m),
    median = INLA::inla.qmarginal(0.5, m),
    q975 = INLA::inla.qmarginal(0.975, m)
  )
}

# tool for figure 5: draw circles to help visualize posterior range
draw_range_circles <- function(r1, r2) {
  centres <- data.frame(name = c("C1", "C2"), lon = c(93, 46), lat = c(21, 50))
  rad <- c(r1, r2)
  circle_cols <- c("#B8860B", "#1B9E77")
  bearings <- matrix(c(0, 90, 180, 270), nrow = 2, byrow = TRUE)
  
  label_box <- function(x, y, label, colour) {
    w <- strwidth(label, cex = 1.1, font = 2) / 2
    h <- strheight(label, cex = 1.1, font = 2) / 2
    rect(x - w - 1, y - h - 1, x + w + 1, y + h + 1, col = "white", border = "grey30")
    text(x, y, label, col = colour, cex = 1.1, font = 2)
  }
  
  for (i in seq_len(nrow(centres))) {
    origin <- c(centres$lon[i], centres$lat[i])
    for (j in seq_along(rad)) {
      colour <- circle_cols[j]
      circle <- geosphere::destPoint(origin, b = seq(0, 360, by = 0.5), d = rad[j] * 1000)
      # separate lines at the antimeridian
      breaks <- c(0L, which(abs(diff(circle[, 1])) > 180), nrow(circle))
      for (k in seq_len(length(breaks) - 1L)) {
        ids <- seq.int(breaks[k] + 1L, breaks[k + 1L])
        if (length(ids) < 2L) next
        lines(circle[ids, 1], circle[ids, 2], col = "white", lwd = 4)
        lines(circle[ids, 1], circle[ids, 2], col = colour, lwd = 2.2)
      }
      bearing <- bearings[i, j]
      edge <- geosphere::destPoint(origin, b = bearing, d = rad[j] * 1000)
      lines(c(origin[1], edge[1]), c(origin[2], edge[2]), col = "white", lwd = 4, lty = 2)
      lines(c(origin[1], edge[1]), c(origin[2], edge[2]), col = colour, lwd = 2.2, lty = 2)
      
      # label radii at C1 only
      if (i == 1L) {
        middle <- geosphere::destPoint(origin, b = bearing, d = rad[j] * 500)
        label_box(
          middle[1] + (if (j == 1L) 8 else 0),
          middle[2] + (if (j == 2L) -4 else 0),
          sprintf("r%d = %.0f km", j, rad[j]), colour
        )
      }
    }
    points(origin[1], origin[2], pch = 21, bg = "white", col = "black", cex = 1.8, lwd = 1.2)
    text(origin[1] + 2, origin[2] + 2, centres$name[i], cex = 1.2, font = 2)
  }
  legend("bottomleft", inset = c(0.06, 0.04),
    legend = c(
      sprintf("r1 = %.0f km (corr. ~ 0.14)", r1),
      sprintf("r2 = %.0f km (corr. = 0.50)", r2),
      "100 observed sites"
    ),
    col = c(circle_cols, "grey20"), lty = c(1, 1, NA), lwd = c(3, 3, NA),
    pch = c(NA, NA, 1), pt.cex = 1.3, cex = 1.05, bg = "white", box.col = "grey50"
  )
  terra::sbar(
    d = 2000, xy = c(140, -5), type = "bar", divs = 2,
    labels = c(0, 1000, 2000), below = "km", lonlat = TRUE, bg = "white", cex = 1.15
  )
}





##########################  function for note 3 experiments ########################## 
############# function for scale1_nonstructure_dd_spatial ############# 
empty_xy <- function() matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("x", "y")))

# periodic boundary
torus_d2 <- function(points, focal, L) {
  dx <- abs(points[, 1L] - focal[1L])
  dy <- abs(points[, 2L] - focal[2L])
  pmin(dx, L - dx)^2 + pmin(dy, L - dy)^2
}

# accelerate neighbour searches using spatial buckets (only search within the surrounding 3x3 buckets)
make_grid <- function(radius, L) {
  nc <- max(1L, as.integer(floor(L / radius)))
  near <- vector("list", nc^2)
  for (iy in 0:(nc - 1L)) for (ix in 0:(nc - 1L)) {
    xs <- unique((ix + (-1:1)) %% nc)
    ys <- unique((iy + (-1:1)) %% nc)
    near[[1L + ix + nc * iy]] <- as.integer(as.vector(outer(xs, ys, function(x, y) 1L + x + nc * y)))
  }
  list(nc = nc, width = L / nc, near = near)
}

cell_ids <- function(xy, grid) {
  ix <- pmin(grid$nc - 1L, floor(xy[, 1L] / grid$width))
  iy <- pmin(grid$nc - 1L, floor(xy[, 2L] / grid$width))
  as.integer(1L + ix + grid$nc * iy)
}

make_buckets <- function(ids, grid) {
  unname(split(seq_along(ids), factor(ids, levels = seq_len(grid$nc^2))))
}

# put initial individuals at random locations
initial_population <- function(z) {
  xy <- matrix(runif(2L * z$initial_n, 0, z$L), ncol = 2L)
  colnames(xy) <- c("x", "y")
  xy
}

# competition kernel
make_kernel <- function(z) {
  list(
    cutoff2 = (z$cutoff_sd * z$Sw)^2,
    twice_variance = 2 * z$Sw^2,
    normalizer = if (z$kernel_normalized)
      2 * pi * z$Sw^2 * (1 - exp(-z$cutoff_sd^2 / 2)) else 1
  )
}

# overall competition strength
kernel_sum <- function(points, focal, L, kernel) {
  if (!nrow(points)) return(0)
  d2 <- torus_d2(points, focal, L)
  d2 <- d2[d2 <= kernel$cutoff2]
  sum(exp(-d2 / kernel$twice_variance)) / kernel$normalizer
}

# adult survival probability: "seeds" do not contribute, and focal adults do not compete with themselves
adult_survival_prob <- function(adults, z, grid, kernel) {
  n <- nrow(adults)
  if (!n) return(numeric(0))
  ids <- cell_ids(adults, grid)
  buckets <- make_buckets(ids, grid)
  crowding <- numeric(n)
  for (i in seq_len(n)) {
    neighbours <- unlist(buckets[grid$near[[ids[i]]]], use.names = FALSE)
    neighbours <- neighbours[neighbours != i]
    if (!length(neighbours)) next
    crowding[i] <- kernel_sum(adults[neighbours, , drop = FALSE], adults[i, ], z$L, kernel)
  }
  z$d * exp(-z$d_prime * crowding)           # nonlinear dd
  #z$d * pmax(0, 1 - z$d_prime * crowding)    # linear dd
}

# establishment of "seeds"
establish_seeds <- function(survivors, seeds, z, grid, kernel) {
  ns <- nrow(survivors)
  nb <- nrow(seeds)
  if (!nb) return(list(adults = survivors, mean_prob = NA_real_))
  seeds <- seeds[sample.int(nb), , drop = FALSE]    # globally shuffle all seeds, avoiding systematic priority by parent
  seed_cells <- cell_ids(seeds, grid)
  buckets <- make_buckets(cell_ids(survivors, grid), grid)
  occupied <- matrix(NA_real_, nrow = ns + nb, ncol = 2L)
  if (ns) occupied[seq_len(ns), ] <- survivors
  n_live <- ns
  sum_prob <- 0
  
  for (j in seq_len(nb)) {
    cell <- seed_cells[j]
    neighbours <- unlist(buckets[grid$near[[cell]]], use.names = FALSE)
    crowding <- kernel_sum(occupied[neighbours, , drop = FALSE], seeds[j, ], z$L, kernel)   # only compete with "adults" and survived "seeds" before it
    prob <- z$g * exp(-z$d_double_prime * crowding)    # nonlinear dd
    #prob <- z$g * pmax(0, 1 - z$d_double_prime * crowding)   # linear dd
    sum_prob <- sum_prob + prob
    if (runif(1L) < prob) {
      n_live <- n_live + 1L
      occupied[n_live, ] <- seeds[j, ]
      buckets[[cell]] <- c(buckets[[cell]], n_live)
    }
  }
  answer <- occupied[seq_len(n_live), , drop = FALSE]
  colnames(answer) <- c("x", "y")
  list(adults = answer, mean_prob = sum_prob / nb)
}

# run one time step
step_population <- function(adults, z, grid, kernel, keep_details = FALSE) {
  n <- nrow(adults)
  # step 1: each adult produces independently Poisson(b) seeds
  offspring <- rpois(n, lambda = z$b)
  parents <- rep(seq_len(n), times = offspring)
  n_seeds <- length(parents)
  seeds <- empty_xy()
  if (n_seeds) {
    displacement <- matrix(rnorm(2L * n_seeds, sd = z$Sm), ncol = 2L)
    seeds <- (adults[parents, , drop = FALSE] + displacement) %% z$L
  }
  # step 2: adult survival
  adult_prob <- adult_survival_prob(adults, z, grid, kernel)
  survivors <- adults[runif(n) < adult_prob, , drop = FALSE]
  # step 3: seed survival (sequential)
  recruitment <- establish_seeds(survivors, seeds, z, grid, kernel)
  next_adults <- recruitment$adults
  n_recruits <- nrow(next_adults) - nrow(survivors)
  answer <- list(
    adults = next_adults,
    diagnostics = data.frame(
      n_before = n, seeds = n_seeds,
      adult_deaths = n - nrow(survivors), recruits = n_recruits,
      n_after = nrow(next_adults),
      mean_adult_survival_prob = if (n) mean(adult_prob) else NA_real_,
      mean_seed_survival_prob = recruitment$mean_prob
    )
  )
  if (keep_details) {
    answer$snapshot <- list(
      adults_before = adults,
      seeds = seeds,
      adult_survivors = survivors,
      recruits = next_adults[nrow(survivors) + seq_len(n_recruits), , drop = FALSE],
      adults_after = next_adults
    )
  }
  answer
}

# run for designated time steps
simulate <- function(adults, steps, z, grid, kernel, keep_first = FALSE) {
  counts <- integer(steps + 1L)
  counts[1L] <- nrow(adults)
  diagnostic_rows <- vector("list", steps)
  first_snapshot <- NULL
  for (t in seq_len(steps)) {
    ans <- step_population(adults, z, grid, kernel, keep_details = keep_first && t == 1L)
    if (keep_first && t == 1L) first_snapshot <- ans$snapshot
    adults <- ans$adults
    counts[t + 1L] <- nrow(adults)
    diagnostic_rows[[t]] <- cbind(step = t, ans$diagnostics)
  }
  list(adults = adults, counts = counts,
       diagnostics = do.call(rbind, diagnostic_rows),
       first_snapshot = first_snapshot)
}

# run full experiment: before and after disturbance
run_experiment <- function(z, n_cores = 1L) {
  stopifnot(z$b >= 0, z$d >= 0, z$d <= 1, z$d_prime >= 0,
            z$d_double_prime >= 0, z$Sm > 0, z$Sw > 0, z$L > 0,
            z$p >= 0, z$p <= 1,
            z$cutoff_sd > 0, z$cutoff_sd * z$Sw <= z$L / 2,
            z$initial_n >= 1, z$n_rep >= 1,
            z$pre_steps >= 1, z$post_steps >= 1)
  stopifnot(length(n_cores) == 1L, is.numeric(n_cores), is.finite(n_cores), n_cores >= 1, n_cores == floor(n_cores))
  n_cores <- as.integer(min(n_cores, z$n_rep))
  grid <- make_grid(z$cutoff_sd * z$Sw, z$L)
  kernel <- make_kernel(z)
  run_one <- function(r) {
    set.seed(z$seed + 3L * r)
    before <- simulate(initial_population(z), z$pre_steps, z, grid, kernel)
    pre <- before$counts
    diagnostics <- vector("list", 3L)
    diagnostics[[1L]] <- cbind(replicate = r, phase = "pre", before$diagnostics)
    shared_state <- before$adults
    n <- nrow(shared_state)
    removed <- as.integer(floor(z$p * n + 0.5))
    # LC removes the nearest m adults to the centre of the study area
    clear_ids <- head(order(torus_d2(shared_state, c(z$L / 2, z$L / 2), z$L)), removed)
    keep_lc <- rep(TRUE, n)
    keep_lc[clear_ids] <- FALSE
    lc_start <- shared_state[keep_lc, , drop = FALSE]
    set.seed(z$seed + 3L * r + 1L)
    thin_ids <- if (removed) sample.int(n, removed) else integer(0)
    keep_dt <- rep(TRUE, n)
    keep_dt[thin_ids] <- FALSE
    dt_start <- shared_state[keep_dt, , drop = FALSE]
    stopifnot(nrow(dt_start) == nrow(lc_start))
    dt_run <- simulate(dt_start, z$post_steps, z, grid, kernel, keep_first = TRUE)
    set.seed(z$seed + 3L * r + 2L)
    lc_run <- simulate(lc_start, z$post_steps, z, grid, kernel, keep_first = TRUE)
    DT <- dt_run$counts
    LC <- lc_run$counts
    snapshots <- list(DT = dt_run$first_snapshot, LC = lc_run$first_snapshot)
    stopifnot(
      nrow(snapshots$DT$adults_before) == DT[1L],
      nrow(snapshots$LC$adults_before) == LC[1L],
      nrow(snapshots$DT$adults_after) == DT[2L],
      nrow(snapshots$LC$adults_after) == LC[2L]
    )
    diagnostics[[2L]] <- cbind(replicate = r, phase = "DT", dt_run$diagnostics)
    diagnostics[[3L]] <- cbind(replicate = r, phase = "LC", lc_run$diagnostics)
    disturbance_row <- data.frame(
      replicate = r, before = n, removed = removed,
      after_DT = nrow(dt_start), after_LC = nrow(lc_start))
    message(sprintf("Repeat %d/%d: N%d=%d; removed=%d; DT%d=%d; LC%d=%d",
                    r, z$n_rep, z$pre_steps, n, removed,
                    z$pre_steps + z$post_steps, tail(dt_run$counts, 1L),
                    z$pre_steps + z$post_steps, tail(lc_run$counts, 1L)))
    list(pre = pre, DT = DT, LC = LC, diagnostics = diagnostics,
         disturbance = disturbance_row, snapshots = snapshots)
  }
  if (n_cores == 1L) {
    runs <- lapply(seq_len(z$n_rep), run_one)
  } else {
    cl <- parallel::makePSOCKcluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("empty_xy", "torus_d2", "cell_ids", "make_buckets", "initial_population",
                                  "kernel_sum", "adult_survival_prob", "establish_seeds", "step_population", "simulate"), envir = environment())
    parallel::clusterCall(cl, function(kind) do.call(RNGkind, as.list(kind)), RNGkind())
    runs <- parallel::parLapply(cl, seq_len(z$n_rep), run_one)
  }
  pre <- do.call(cbind, lapply(runs, `[[`, "pre"))
  DT <- do.call(cbind, lapply(runs, `[[`, "DT"))
  LC <- do.call(cbind, lapply(runs, `[[`, "LC"))
  diagnostics <- do.call(c, lapply(runs, `[[`, "diagnostics"))
  disturbance_rows <- lapply(runs, `[[`, "disturbance")
  snapshots <- lapply(runs, `[[`, "snapshots")
  diagnostic_table <- do.call(rbind, diagnostics)
  disturbance_table <- do.call(rbind, disturbance_rows)
  mean_trajectory <- rbind(
    data.frame(phase = "pre", time = 0:z$pre_steps, time_since_disturbance = NA_integer_, mean_N = rowMeans(pre)),
    data.frame(phase = "DT", time = z$pre_steps + 0:z$post_steps, time_since_disturbance = 0:z$post_steps, mean_N = rowMeans(DT)),
    data.frame(phase = "LC", time = z$pre_steps + 0:z$post_steps, time_since_disturbance = 0:z$post_steps, mean_N = rowMeans(LC))
  )
  list(pre = pre, DT = DT, LC = LC, mean_trajectory = mean_trajectory,
       diagnostics = diagnostic_table, disturbances = disturbance_table,
       snapshots = snapshots, params = z)
}

# spatial scatter plots (optional)
plot_spatial_snapshot <- function(result, rep_id = 1L, stage = c("seeds", "next")) {
  stage <- match.arg(stage)
  stopifnot(length(rep_id) == 1L, is.finite(rep_id),    # only plot for the first replicate, plot immediately after disturbance
            rep_id == as.integer(rep_id), rep_id >= 1L,
            rep_id <= length(result$snapshots))
  z <- result$params
  op <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(op) }, add = TRUE)
  # a map for LC, a map for DT, a shared legend
  layout(matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE),
         heights = c(1, 0.14))
  par(mar = c(1, 1, 3.8, 1), oma = c(0, 0, 1.5, 0))
  green <- adjustcolor("#008000", alpha.f = 0.45)
  red <- "#E31A1C"
  for (type in c("DT", "LC")) {
    snap <- result$snapshots[[rep_id]][[type]]
    plot(NA_real_, NA_real_, type = "n", xlim = c(0, z$L),
         ylim = c(0, z$L), asp = 1, xaxs = "i", yaxs = "i",
         xlab = "", ylab = "", axes = FALSE)
    box()
    if (stage == "seeds") {
      # plot all seeds in green, then overwrite successful seed positions in red
      points(snap$seeds, pch = 16, col = green, cex = 0.45)
      points(snap$adults_before, pch = 16, col = "black", cex = 0.65)
      points(snap$recruits, pch = 16, col = red, cex = 0.65)
      # title(main = type, cex.main = 1.5)
      mtext(sprintf("adults = %d; seeds = %d; established = %d",
                    nrow(snap$adults_before), nrow(snap$seeds),
                    nrow(snap$recruits)), side = 3, line = 0.6, cex = 1.5)
    } else {
      points(snap$adult_survivors, pch = 16, col = "black", cex = 0.65)
      points(snap$recruits, pch = 16, col = red, cex = 0.65)
      #title((main = type), cex.main = 1.5)
      mtext(sprintf("total = %d; old survivors = %d; new adults = %d",
                    nrow(snap$adults_after), nrow(snap$adult_survivors),
                    nrow(snap$recruits)), side = 3, line = 0.6, cex = 0.8)
    }
  }
  par(mar = rep(0, 4))
  plot.new()
  if (stage == "seeds") {
    legend("center", c("\"adults\" immediately after disturbance", "\"seeds\"",
                       "successfully established \"seeds\""),
           col = c("black", "#008000", red), pch = 16, horiz = TRUE,
           bty = "n", cex = 1.5)
  } else {
    legend("center", c("surviving old adults", "new adults from seeds"),
           col = c("black", red), pch = 16, horiz = TRUE, bty = "n", cex = 1.5)
  }
  invisible(result$snapshots[[rep_id]])
}

save_spatial_snapshots <- function(result, scenario_name) {
  out_dir <- result$params$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rep_id <- 1L
  # only plot the first repetition
  for (stage in c("seeds", "next")) {
    pdf(file.path(out_dir, sprintf("%s_spatial_%s_rep%02d.pdf", scenario_name, stage, rep_id)),
        width = 12, height = 7
    )
    tryCatch(
      plot_spatial_snapshot(result, rep_id = rep_id, stage = stage),
      finally = dev.off()
    )
  }
  invisible(out_dir)
}

plot_trajectories <- function(result, font_scale = 1.5) {
  z <- result$params
  pre_mean <- rowMeans(result$pre)
  dt_mean <- rowMeans(result$DT)
  lc_mean <- rowMeans(result$LC)
  n_rep <- ncol(result$pre)
  stopifnot(n_rep >= 2L)
  pre_se <- apply(result$pre, 1L, sd) / sqrt(n_rep)
  dt_se <- apply(result$DT, 1L, sd) / sqrt(n_rep)
  lc_se <- apply(result$LC, 1L, sd) / sqrt(n_rep)
  # collapse blank spaces
  all_mean <- c(pre_mean, dt_mean, lc_mean)
  all_se <- c(pre_se, dt_se, lc_se)
  y_range <- range(c(pmax(0, all_mean - 1.96 * all_se), all_mean + 1.96 * all_se), finite = TRUE)
  y_pad <- max(diff(y_range) * 0.08, 1)
  y_limits <- c(
    max(0, y_range[1L] - y_pad),
    y_range[2L] + y_pad
  )
  op <- par(
    mar = c(5.5, 6, 1.5, 1),
    mgp = c(3.5, 1, 0),
    cex.axis = font_scale,
    cex.lab = font_scale,
    cex.main = font_scale
  )
  on.exit(par(op), add = TRUE)
  plot(
    0:z$pre_steps, pre_mean,
    type = "n", col = "#666666", lwd = 2,
    xlim = c(0, z$pre_steps + z$post_steps),
    ylim = y_limits, yaxs = "i",
    xlab = "time step",
    ylab = "mean population size",
    main = ""
  )
  draw_band <- function(time, avg, se, col) {
    polygon(
      c(time, rev(time)),
      c(pmax(0, avg - 1.96 * se), rev(avg + 1.96 * se)),
      col = adjustcolor(col, alpha.f = 0.14), border = NA
    )
  }
  draw_band(0:z$pre_steps, pre_mean, pre_se, "#666666")
  draw_band(z$pre_steps + 0:z$post_steps, dt_mean, dt_se, "#1b9e77")
  draw_band(z$pre_steps + 0:z$post_steps, lc_mean, lc_se, "orange")
  lines(0:z$pre_steps, pre_mean, col = "#666666", lwd = 2)
  abline(v = z$pre_steps, lty = 3, col = "grey50")
  post_x <- c(z$pre_steps, z$pre_steps + 0:z$post_steps)
  lines(post_x, c(tail(pre_mean, 1L), dt_mean), col = "#1b9e77", lwd = 2)
  lines(post_x, c(tail(pre_mean, 1L), lc_mean), col = "orange", lwd = 2)
  legend("bottomright",c("before disturbance","dispersed thinning (DT)","localised clearance (LC)"),
    col = c("#666666", "#1b9e77", "orange"),lty = 1, lwd = 2, bty = "n",cex = font_scale)
}




############# additional function for scale1_structure_dd_spatial ############# 
st_initial_population <- function(z) {
  na <- as.integer(round(z$initial_n * z$adult_fraction))
  xy <- matrix(runif(2L * z$initial_n, 0, z$L), ncol = 2L)
  colnames(xy) <- c("x", "y")
  list(J = xy[seq_len(z$initial_n - na), , drop = FALSE],
       A = xy[z$initial_n - na + seq_len(na), , drop = FALSE])
}

st_counts <- function(pop) {
  c(J = nrow(pop$J), A = nrow(pop$A), N = nrow(pop$J) + nrow(pop$A))
}

# competition_mode: 1 = old adults only; 2 = all old individuals; exclude the focal individual
st_crowding <- function(pop, z, grid, kernel) {
  nj <- nrow(pop$J)
  xy <- rbind(pop$J, pop$A)
  n <- nrow(xy)
  if (!n) return(numeric(0))
  ids <- cell_ids(xy, grid)
  buckets <- make_buckets(ids, grid)
  crowding <- numeric(n)
  for (i in seq_len(n)) {
    neighbours <- unlist(buckets[grid$near[[ids[i]]]], use.names = FALSE)
    if (z$competition_mode == 1L) neighbours <- neighbours[neighbours > nj]
    neighbours <- neighbours[neighbours != i]
    crowding[i] <- kernel_sum(xy[neighbours, , drop = FALSE], xy[i, ], z$L, kernel)
  }
  crowding
}

# one synchronous time step, matching the two-stage MPM
st_step <- function(pop, z, grid, kernel) {
  nj <- nrow(pop$J)
  na <- nrow(pop$A)
  # evaluate competition using the population at the start of the time step
  crowding <- st_crowding(pop, z, grid, kernel)
  crowding_J <- crowding[seq_len(nj)]
  crowding_A <- crowding[nj + seq_len(na)]
  # potential dd on all four vital rates 
  sigma1 <- z$sigma1_max * exp(-z$k1 * crowding_J)
  sigma2 <- z$sigma2_max * exp(-z$k2 * crowding_A)
  gamma <- z$gamma_max * exp(-z$k3 * crowding_J)
  phi <- z$phi_max * exp(-z$k4 * crowding_A)
  # step 1: adults reproduce 
  offspring <- rpois(na, phi)
  parents <- rep(seq_len(na), times = offspring)
  seeds <- empty_xy()
  if (length(parents)) {
    displacement <- matrix(rnorm(2L * length(parents), sd = z$Sm), ncol = 2L)
    seeds <- (pop$A[parents, , drop = FALSE] + displacement) %% z$L
  }
  # step 2:  survival of juveniles, then surviving juveniles may mature
  alive <- runif(nj) < sigma1
  mature <- alive & (runif(nj) < gamma)
  # step 3: survival of adults
  adults_alive <- runif(na) < sigma2
  list(J = rbind(pop$J[alive & !mature, , drop = FALSE], seeds),
       A = rbind(pop$A[adults_alive, , drop = FALSE], pop$J[mature, , drop = FALSE]))
}

st_simulate <- function(pop, steps, z, grid, kernel) {
  counts <- matrix(0L, nrow = steps + 1L, ncol = 3L, dimnames = list(NULL, c("J", "A", "N")))
  counts[1L, ] <- st_counts(pop)
  for (t in seq_len(steps)) {
    if (!sum(st_counts(pop)[1:2])) break
    pop <- st_step(pop, z, grid, kernel)
    counts[t + 1L, ] <- st_counts(pop)
    if (counts[t + 1L, "N"] > z$max_population) stop("Population exceeded max_population.")
  }
  list(pop = pop, counts = counts)
}

ss_cases <- function(z) {
  if (isTRUE(z$stage_specific_lc)) {
    cases <- expand.grid(space = c("LC_A", "LC_J"), stage = "mixed", factor = z$abundance_factors, stringsAsFactors = FALSE)
    cases$key <- paste(cases$stage, cases$space, cases$factor, sep = "_")
    return(cases)
  }
  cases <- expand.grid(space = c("DT", "LC"), stage = c("J", "A"), factor = z$abundance_factors, stringsAsFactors = FALSE)
  cases$key <- paste(cases$stage, cases$space, cases$factor, sep = "_")
  cases
}

ss_uniform_xy <- function(n, L) {
  xy <- matrix(runif(2L * n, 0, L), ncol = 2L)
  colnames(xy) <- c("x", "y")
  xy
}

# relocate the points inside the central circle, do not change total population size
ss_relocate <- function(xy, z) {
  radius2 <- z$cleared_area * z$L^2 / pi
  inside <- which(torus_d2(xy, c(z$L / 2, z$L / 2), z$L) < radius2)
  while (length(inside)) {
    xy[inside, ] <- ss_uniform_xy(length(inside), z$L)
    inside <- inside[torus_d2(xy[inside, , drop = FALSE],
                              c(z$L / 2, z$L / 2), z$L) < radius2]
  }
  xy
}

ss_initial_states <- function(K_empirical, z, pop = NULL) {
  if (isTRUE(z$stage_specific_lc)) {
    stopifnot(!is.null(pop), z$abundance_factors == 1)
    cases <- ss_cases(z)
    states <- setNames(vector("list", nrow(cases)), cases$key)
    for (i in seq_len(nrow(cases))) {
      states[[i]] <- pop
      target <- if (cases$space[i] == "LC_A") "A" else "J"
      states[[i]][[target]] <- ss_relocate(pop[[target]], z)
    }
    return(states)
  }
  n_full <- as.integer(floor(K_empirical + 0.5))
  dt <- ss_uniform_xy(n_full, z$L)
  lc <- ss_relocate(dt, z)
  indices <- lapply(z$abundance_factors, function(f) {
    n <- as.integer(floor(f * K_empirical + 0.5))
    if (n == n_full) seq_len(n_full) else sort(sample.int(n_full, n))
  })
  cases <- ss_cases(z)
  states <- setNames(vector("list", nrow(cases)), cases$key)
  for (i in seq_len(nrow(cases))) {
    take <- indices[[match(cases$factor[i], z$abundance_factors)]]
    xy <- if (cases$space[i] == "DT") dt[take, , drop = FALSE] else lc[take, , drop = FALSE]
    states[[i]] <- list(J = empty_xy(), A = empty_xy())
    states[[i]][[cases$stage[i]]] <- xy
  }
  states
}

ss_run_one <- function(job) {
  z <- job$z
  r <- job$rep
  seed <- z$seed + 1000L * r
  grid <- make_grid(z$cutoff_sd * z$Sw, z$L)
  kernel <- make_kernel(z)
  set.seed(seed)
  before <- st_simulate(st_initial_population(z), z$pre_steps, z, grid, kernel)
  times <- seq.int(z$equilibrium_window[1L], z$equilibrium_window[2L])
  baseline <- colMeans(before$counts[times + 1L, , drop = FALSE])
  set.seed(seed + 1L)
  starts <- ss_initial_states(baseline["N"], z, pop = before$pop)
  cases <- ss_cases(z)
  post <- setNames(vector("list", nrow(cases)), cases$key)
  for (i in seq_len(nrow(cases))) {
    pair <- if (isTRUE(z$stage_specific_lc)) 1L else
      match(cases$stage[i], c("J", "A")) + 2L * (match(cases$factor[i], z$abundance_factors) - 1L)
    set.seed(seed + 10L + pair)
    post[[i]] <- st_simulate(starts[[i]], z$post_steps, z, grid, kernel)$counts
  }
  centre <- floor(mean(z$equilibrium_window))
  early <- mean(before$counts[(times[times < centre]) + 1L, "N"])
  late <- mean(before$counts[(times[times > centre]) + 1L, "N"])
  log <- data.frame(rep = r, seed = seed, K_empirical = unname(baseline["N"]),
                    J_empirical = unname(baseline["J"]), A_empirical = unname(baseline["A"]),
                    N_at_pre_end = tail(before$counts[, "N"], 1L),
                    early_window_mean = early, late_window_mean = late,
                    fractional_window_change = if (early > 0) (late - early) / early else NA_real_)
  list(pre = before$counts, post = post, log = log, starts = if (r == 1L) starts else NULL)
}

ss_summary_path <- function(paths) {
  avg <- Reduce("+", paths) / length(paths)
  sd <- sqrt(Reduce("+", lapply(paths, function(x) (x - avg)^2)) / (length(paths) - 1L))
  list(mean = avg, se = sd / sqrt(length(paths)))
}

ss_run_experiment <- function(z, n_cores = 1L) {
  rate_names <- c("sigma1_max", "sigma2_max", "gamma_max", "phi_max", "k1", "k2", "k3", "k4")
  stopifnot(all(vapply(z[rate_names], function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }, logical(1))))
  stopifnot(length(n_cores) == 1L, is.finite(n_cores), n_cores >= 1,
            n_cores == floor(n_cores), z$n_rep >= 2, z$n_rep == floor(z$n_rep),
            z$sigma1_max >= 0, z$sigma1_max <= 1, z$gamma_max >= 0, z$gamma_max <= 1,
            z$sigma2_max >= 0, z$sigma2_max <= 1, z$phi_max >= 0,
            z$k1 >= 0, z$k2 >= 0, z$k3 >= 0, z$k4 >= 0,
            length(z$competition_mode) == 1L, is.numeric(z$competition_mode),
            is.finite(z$competition_mode), z$competition_mode %in% c(1L, 2L),
            z$L > 0, z$Sm > 0, z$Sw > 0, z$kernel_normalized,
            z$cutoff_sd > 0, z$cutoff_sd * z$Sw <= z$L / 2,
            z$cleared_area > 0, z$cleared_area <= pi / 4,
            length(z$abundance_factors) == 1L, is.numeric(z$abundance_factors),
            is.finite(z$abundance_factors), z$abundance_factors > 0, z$abundance_factors <= 1,
            length(z$equilibrium_window) == 2L, z$equilibrium_window[1L] >= 0,
            diff(z$equilibrium_window) >= 4, z$equilibrium_window[2L] <= z$pre_steps,
            z$initial_n >= 1, z$initial_n == floor(z$initial_n),
            z$adult_fraction >= 0, z$adult_fraction <= 1,
            z$post_steps >= 1, z$max_population > z$initial_n,
            z$seed >= 0, z$seed + 1000L * z$n_rep + 20L < .Machine$integer.max)
  jobs <- lapply(seq_len(z$n_rep), function(r) list(z = z, rep = r))
  n_cores <- as.integer(min(n_cores, z$n_rep))
  message("Running ", z$n_rep, " replicates x ", nrow(ss_cases(z)), " initial states on ", n_cores, " workers.")
  if (n_cores == 1L) {
    runs <- lapply(jobs, ss_run_one)
  } else {
    cl <- parallel::makePSOCKcluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl))
    env <- environment(ss_run_experiment)
    helpers <- c("empty_xy", "torus_d2", "make_grid", "cell_ids", "make_buckets", "make_kernel", "kernel_sum",
                 "st_initial_population", "st_counts", "st_crowding", "st_step", "st_simulate")
    parallel::clusterExport(cl, c(helpers, ls(env, pattern = "^ss_")), envir = env)
    parallel::clusterCall(cl, function(kind) do.call(RNGkind, as.list(kind)), RNGkind())
    runs <- parallel::parLapply(cl, jobs, ss_run_one)
  }
  cases <- ss_cases(z)
  pre <- ss_summary_path(lapply(runs, `[[`, "pre"))
  post <- setNames(lapply(cases$key, function(key) {
    ss_summary_path(lapply(runs, function(x) x$post[[key]]))
  }), cases$key)
  calibration <- do.call(rbind, lapply(runs, `[[`, "log"))
  endpoint <- cases
  for (i in seq_len(nrow(cases))) {
    paths <- lapply(runs, function(x) x$post[[cases$key[i]]])
    endpoint$initial_N[i] <- post[[i]]$mean[1L, "N"]
    endpoint$mean_final_N[i] <- tail(post[[i]]$mean[, "N"], 1L)
    endpoint$mean_last_20_N[i] <- mean(tail(post[[i]]$mean[, "N"], min(20L, z$post_steps)))
    endpoint$extinct_replicates[i] <- sum(vapply(paths, function(x) tail(x[, "N"], 1L) == 0, logical(1)))
  }
  calibration_summary <- data.frame(
    n_rep = z$n_rep, mean_K_empirical = mean(calibration$K_empirical),
    sd_K_empirical = sd(calibration$K_empirical),
    mean_N_pre_end = mean(calibration$N_at_pre_end),
    mean_early_window = mean(calibration$early_window_mean),
    mean_late_window = mean(calibration$late_window_mean),
    extinct_at_pre_end = sum(calibration$N_at_pre_end == 0))
  list(params = z, n_cores = n_cores, cases = cases, runs = runs,
       pre = pre, post = post, calibration = calibration,
       calibration_summary = calibration_summary, endpoint_summary = endpoint)
}

ss_band <- function(time, avg, se, col) {
  polygon(c(time, rev(time)), c(pmax(0, avg - 1.96 * se), rev(avg + 1.96 * se)), col = adjustcolor(col, alpha.f = 0.14), border = NA)
}

ss_plot_burnin <- function(result) {
  z <- result$params
  time <- 0:z$pre_steps
  avg <- result$pre$mean
  se <- result$pre$se
  op <- par(mar = c(4.5, 4.7, 3.2, 1))
  on.exit(par(op))
  ylim <- c(0, 1.08 * max(avg + 1.96 * se))
  plot(time, avg[, "N"], type = "n", ylim = ylim, xlab = "burn-in time step",
       ylab = "mean population size", main = paste0("Empirical equilibrium calibration (", z$n_rep, " replicates)"))
  rect(z$equilibrium_window[1L], 0, z$equilibrium_window[2L], ylim[2L], col = "#EEF1F5", border = NA)
  ss_band(time, avg[, "N"], se[, "N"], "#333333")
  abline(h = mean(result$calibration$K_empirical), lty = 3, col = "#666666")
  matlines(time, avg[, c("N", "J", "A")], lty = c(1, 1, 1), lwd = c(2.4, 1.8, 1.8),
           col = c("#333333", "#B24058", "#366B93"))
  legend("bottomright", c("total", "juveniles", "adults", paste0("mean over steps ", paste(z$equilibrium_window, collapse = "-"))),
         col = c("#333333", "#B24058", "#366B93", "#666666"),
         lty = c(1, 1, 1, 3), lwd = c(2.4, 1.8, 1.8, 1), bty = "n", cex = 0.85)
}

ss_plot_trajectories <- function(result, component = "N", time_max = NULL) {
  z <- result$params
  f <- z$abundance_factors
  stopifnot(length(f) == 1L)
  component <- match.arg(component, c("N", "J", "A"))
  time <- 0:z$post_steps
  if (is.null(time_max)) time_max <- z$post_steps
  keys <- if (isTRUE(z$stage_specific_lc)) {
    paste(c("mixed_LC_A", "mixed_LC_J"), f, sep = "_")
  } else {
    paste(c("A_LC", "J_LC", "A_DT", "J_DT"), f, sep = "_")
  }
  stopifnot(all(keys %in% names(result$post)))
  paths <- result$post[keys]
  colours <- c("orange", "orange", "#1b9e77", "#1b9e77")
  line_types <- c(1, 2, 1, 2)
  labels <- c("LC + all adults", "LC + all juveniles",
              "DT + all adults", "DT + all juveniles")
  if (isTRUE(z$stage_specific_lc)) {
    colours <- c("#B24058", "#366B93")
    line_types <- c(1, 1)
    labels <- c("LC on adults only", "LC on juveniles only")
  }
  y_range <- range(c(
    unlist(lapply(paths, function(x) {
      c(pmax(0, x$mean[, component] - 1.96 * x$se[, component]), x$mean[, component] + 1.96 * x$se[, component])
    })),
    mean(result$calibration$K_empirical)
  ), finite = TRUE)
  y_pad <- max(diff(y_range) * 0.08, 1)
  y_limits <- c(max(0, y_range[1L] - y_pad), y_range[2L] + y_pad)
  op <- par(mfrow = c(1, 1), mar = c(4.5, 4.7, 3.5, 1), oma = c(0, 0, 0, 0))
  on.exit(par(op))
  plot(time, paths[[1L]]$mean[, component], type = "n",
       ylim = y_limits, yaxs = "i", xlim = c(0, time_max),
       xlab = "time step", ylab = "mean population size")
  abline(h = mean(result$calibration$K_empirical), lty = 3, col = "#666666")
  for (i in seq_along(paths)) {
    ss_band(time, paths[[i]]$mean[, component], paths[[i]]$se[, component], colours[i])
  }
  for (i in seq_along(paths)) {
    lines(time, paths[[i]]$mean[, component], col = colours[i], lty = line_types[i], lwd = 2)
  }
  legend("topright", labels, col = colours, lty = line_types, lwd = 2, bty = "n", cex = 0.9)
}

ss_plot_initial_states <- function(result) {
  z <- result$params
  if (isTRUE(z$stage_specific_lc)) {
    op <- par(mfrow = c(1, 2), mar = c(2.2, 2.2, 3.2, 0.5), oma = c(0, 0, 2, 0))
    on.exit(par(op))
    for (i in seq_len(nrow(result$cases))) {
      pop <- result$runs[[1L]]$starts[[result$cases$key[i]]]
      label <- if (result$cases$space[i] == "LC_A") "LC on adults only" else "LC on juveniles only"
      plot(NA_real_, NA_real_, type = "n", xlim = c(0, z$L), ylim = c(0, z$L), asp = 1,
           xlab = "x", ylab = "y", main = paste0(label, "\nJ=", nrow(pop$J), "; A=", nrow(pop$A)))
      points(pop$J, pch = 16, cex = 0.6, col = "#B24058")
      points(pop$A, pch = 16, cex = 0.6, col = "#366B93")
      symbols(z$L / 2, z$L / 2, circles = sqrt(z$cleared_area * z$L^2 / pi),
              inches = FALSE, add = TRUE, fg = "#666666", lty = 2)
      legend("topright", c("juveniles", "adults"), col = c("#B24058", "#366B93"), pch = 16, bty = "n", cex = 0.8)
    }
    mtext("Initial states (replicate 1)", outer = TRUE)
    return(invisible(NULL))
  }
  op <- par(mfrow = c(2, 4), mar = c(2.2, 2.2, 3.2, 0.5), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  for (stage in c("J", "A")) for (f in z$abundance_factors) for (space in c("DT", "LC")) {
    pop <- result$runs[[1L]]$starts[[paste(stage, space, f, sep = "_")]]
    xy <- pop[[stage]]
    plot(xy, xlim = c(0, z$L), ylim = c(0, z$L), asp = 1, pch = 16, cex = 0.6,
         col = if (stage == "J") "#B24058" else "#366B93", xlab = "x", ylab = "y",
         main = paste(stage, space, paste0(f, "K"), paste0("N=", nrow(xy))))
    if (space == "LC") symbols(z$L / 2, z$L / 2, circles = sqrt(z$cleared_area * z$L^2 / pi),
                               inches = FALSE, add = TRUE, fg = "#666666", lty = 2)
  }
  mtext("Initial states (replicate 1)", outer = TRUE)
}

ss_save_result <- function(result, figure_dir, scenario_name = "scenario") {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(figure_dir, paste0(scenario_name, "_total.pdf")), width = 8, height = 5.5)
  tryCatch(ss_plot_trajectories(result, component = "N"), finally = dev.off())
  pdf(file.path(figure_dir, paste0(scenario_name, "_burnin.pdf")), width = 8, height = 5.5)
  tryCatch(ss_plot_burnin(result), finally = dev.off())
  pdf(file.path(figure_dir, paste0(scenario_name, "_initial.pdf")), width = 9, height = 5.5)
  tryCatch(ss_plot_initial_states(result), finally = dev.off())
}


