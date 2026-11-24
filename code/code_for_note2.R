### code for generating all figures and tables in note 2
output_dir_figures <- here("output", "figures", "note2")
output_dir_tables <- here("output", "tables", "note2")
dir.create(output_dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_tables, recursive = TRUE, showWarnings = FALSE)

################# figure 1 #################
### panel a
# data from the first quadrat of PLANTPOPNET_Y0_V1.3_2020-11-18
pts_df <- data.frame(
  x_coord = c(30, 19, 0, 31, 35, 45, 14, 8, 18, 25, 28, 35, 2),
  y_coord = c(8, 40, 24, 32, 34, 46, 30, 32, 12, 12, 12, 13, 37)
)

window <- owin(xrange = c(0, 50), yrange = c(0, 50))    # 50x50 cm quadrat
pts_ppp <- as.ppp(pts_df, W = window)    # transform to point pattern
pts_voronoi <- spatstat.geom::dirichlet(pts_ppp)    # Voronoi

pdf(file.path(output_dir_figures, "scale1.pdf"), width = 5, height = 5)
plot(pts_voronoi, main = "", border = "grey", lwd = 1.5)
plot(pts_ppp, add = TRUE, pch = 16, col = "black", cex = 1.2)
dev.off()


### panel c
# data from PLANTPOPNET
sites <- read.csv(here("data", "Coordinates_Oct2020_site_level.csv"))
sites <- sites[is.finite(sites$longitude) & is.finite(sites$latitude),]

# temperature data from WorldClim 2.1
tavg <- geodata::worldclim_global(var = "tavg", res = 10, path = here("data", "worldclim"))
temperature <- mean(tavg)

pdf(file.path(output_dir_figures, "scale3.pdf"), width = 10, height = 5)
plot(temperature, col = hcl.colors(100, "RdYlBu", rev = TRUE))
# add national borders
maps::map("world", add = TRUE, col = "grey25", lwd = 0.4)
# add population sites
points(sites$longitude, sites$latitude, pch = 21, bg = "green", col = "black", cex = 0.9)
# add scale
terra::sbar(d = 5000, xy = "bottomleft", type = "bar", below = "km", lonlat = TRUE)
dev.off()


################################## scale 1 ##################################
################# figure 2 #################
set.seed(123)
W <- owin(xrange = c(0, 100), yrange = c(0, 100))

### simulate point process
# process 1: inhibition
collision_distance <- 6
X_regular <- rSSI(r = collision_distance, n = 120, win = W, giveup = 100000)
# process 2: Thomas cluster
X_cluster <- rThomas(
  kappa = 0.0037,  # intensity of cluster centre
  scale = 5,       # deffusion scale of each cluster
  mu = 4,          # intensity of each cluster
  win = W
)

### pcf
r_values <- seq(0, 20, by = 0.05)
pcf_regular <- pcf(X_regular, r = r_values, correction = "translate", divisor = "d", kernel = "epanechnikov", bw = 0.3)
pcf_cluster <- pcf(X_cluster, r = r_values, correction = "translate", divisor = "d", kernel = "epanechnikov", bw = 0.3)
pdf(file.path(output_dir_figures, "pcf.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# inhibition
plot(X_regular, chars = 16, cols = "#2166AC", cex = 0.7)
plot_pcf(pcf_regular, collision = collision_distance, curve_colour = "#2166AC")
# Thomas cluster 
plot(X_cluster, chars = 16, cols = "#B2182B", cex = 0.7)
plot_pcf(pcf_cluster, collision = NULL, curve_colour = "#B2182B")
par(mfrow = c(1, 1))
dev.off()

npoints(X_regular)
npoints(X_cluster)


################# figure 3 #################
set.seed(1)

side <- 100
n_individuals <- 300
min_distance <- 2  # collision

### simulate point process
xy <- matrix(numeric(0), ncol = 2)
colnames(xy) <- c("x", "y")
while (nrow(xy) < n_individuals) {
  candidate <- runif(2, min = 0, max = side)
  if (nrow(xy) == 0 || all(rowSums(sweep(xy, 2, candidate, "-")^2) >= min_distance^2)) {
    xy <- rbind(xy, candidate)
  }
}

### assign marks (deterministic for visualization)
# dispersed removal
grid <- expand.grid(
  x = (seq_len(10) - 0.5) * side / 10,
  y = (seq_len(10) - 0.5) * side / 10
)
removed_dispersed <- integer(0)
for (i in seq_len(nrow(grid))) {
  available <- setdiff(seq_len(n_individuals), removed_dispersed)
  d2 <- (xy[available, "x"] - grid$x[i])^2 + (xy[available, "y"] - grid$y[i])^2
  removed_dispersed <- c(removed_dispersed, available[which.min(d2)])
}
# clustered removal
patch_centres <- matrix(
  c(22, 25,    # specify two centres at random
    80, 57),
  ncol = 2, byrow = TRUE
)
removed_patches <- integer(0)
for (i in seq_len(nrow(patch_centres))) {
  available <- setdiff(seq_len(n_individuals), removed_patches)
  d2 <- (xy[available, "x"] - patch_centres[i, 1])^2 + (xy[available, "y"] - patch_centres[i, 2])^2
  removed_patches <- c(removed_patches, available[order(d2)[seq_len(50)]])
}

marks_dispersed <- rep("survive", n_individuals)
marks_dispersed[removed_dispersed] <- "removed"
marks_patches <- rep("survive", n_individuals)
marks_patches[removed_patches] <- "removed"
stopifnot(
  sum(marks_dispersed == "removed") == 100,
  sum(marks_patches == "removed") == 100
)

### plot marked point process
pdf(file.path(output_dir_figures, "markedpp.pdf"), width = 9, height = 5)
layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
par(mar = c(2, 2, 4, 1)) 
plot_population(marks_dispersed, "dispersed thinning")
plot_population(marks_patches, "local clearance")
legend("topright", legend = c("survived", "removed"), col = c("grey45", "orange"), 
       pch = 16,  horiz = TRUE,  bty = "n", xpd = TRUE, inset = c(0, -0.18))
dev.off()



################################## scale 3 ##################################
################# figure 4 & table 2 #################
set.seed(123)
##### data part #####
### import data
# map of Eurasia
area <- terra::ext(-12, 180, -12, 82)
eurasia <- rnaturalearth::ne_countries(scale = 110, continent = c("Europe", "Asia"), returnclass = "sv")
eurasia <- terra::crop(eurasia, area)

# monthly temperature and wind speed data from WorldClim
tavg <- geodata::worldclim_global(var = "tavg", res = 10, path = here("data", "worldclim"))
wind <- geodata::worldclim_global(var = "wind", res = 10, path = here("data", "worldclim"))
T_field <- terra::mask(terra::crop(mean(tavg), area), eurasia)
wind_field <- terra::mask(terra::crop(mean(wind), area), eurasia)
names(T_field) <- "T"
names(wind_field) <- "wind"

### simulate "historical disturbance D" using wind speed as template and add noise
variation_deg <- 25   # range of random patch, can be tuned to simulate different patterns of D
coarse <- terra::rast(
  nrows = ceiling(94 / variation_deg),
  ncols = ceiling(192 / variation_deg),
  xmin = -12, xmax = 180,
  ymin = -12, ymax = 82,
  crs = terra::crs(T_field)
)
# sample an iid standard normal distribution for each coarse grid cell and smooth between cells
terra::values(coarse) <- rnorm(terra::ncell(coarse))
smooth_noise <- terra::mask(terra::resample(coarse, T_field, method = "bilinear"), T_field)
# standardized wind speed + standardized spatial noise
D_field <- standardize_field(0.8 * standardize_field(wind_field) +0.7 * standardize_field(smooth_noise))
names(D_field) <- "D"
# D_field is visualized at the end of this section

### simulate samples across Eurasia 
# fix five sites as examples
fixed <- data.frame(
  population = LETTERS[1:5],
  lon = c(-0.13, -3.70, 117, 117.9, 113.92),
  lat = c(51.51, 40.40, 36.67, 33.6, -2.21),
  S = c(20, 10, 10, 5, 5)   # here S is species size (denote as G in the note)
)
fixed$region <- region_of(fixed$lon, fixed$lat) 
# randomly sample the other 95 sites
candidate <- as.data.frame(c(T_field, D_field), xy = TRUE, cells = TRUE, na.rm = TRUE)    # candidates of other samples
candidate$region <- region_of(candidate$x, candidate$y)
fixed_cells <- terra::cellFromXY(T_field, as.matrix(fixed[, c("lon", "lat")]))
candidate <- candidate[!candidate$cell %in% fixed_cells, ]
# assign heterogeneous density among regions
targets <- c(
  Europe = 60,
  EastAsia = 15,
  NorthernRussia = 5,
  SouthSoutheastAsia = 15,
  Other = 5
)
# other samples
n_extra <- as.integer(targets) - as.integer(
  table(factor(fixed$region, levels = names(targets)))
)
names(n_extra) <- names(targets)
stopifnot(all(n_extra >= 0), sum(n_extra) == 95)
extra <- do.call(rbind, lapply(names(targets), function(g) {
  possible <- which(candidate$region == g)
  chosen <- candidate[sample(possible, n_extra[[g]]), ]
  data.frame(lon = chosen$x, lat = chosen$y, region = g)
}))
extra$population <- sprintf("P%03d", 6:100)
extra$S <- sample(rep(c(5, 10, 20), length.out = 95))
# merge samples
cols <- c("population", "lon", "lat", "region", "S")
dat <- rbind(fixed[, cols], extra[, cols])

# values of T and D at sample locations
climate_at_sites <- terra::extract(
  c(T_field, D_field),
  as.matrix(dat[, c("lon", "lat")])
)
dat$T <- climate_at_sites$T
dat$D <- climate_at_sites$D
stopifnot(nrow(dat) == 100, !anyNA(dat$T), !anyNA(dat$D))

### simulate resilience at sampled locations
# prespecified model parameters
beta0 <- 40
beta_T <- 0.4
beta_S <- 0.25
beta_D <- 5
sigma_eps <- 1.2

# underlying mechanism of generating resilience: R_i=40+0.4*T_i+0.25*S_i+5*D_i+epsilon_i
dat$R <- with(
  dat,
  beta0 + beta_T * T + beta_S * S + beta_D * D +
    rnorm(nrow(dat), mean = 0, sd = sigma_eps)
)
dat <- cbind(dat, to_sphere(dat$lon, dat$lat))

# train test split 
train_idx <- sample.int(100, size = 70)    # 70 samples for training
dat$set <- "test"
dat$set[train_idx] <- "train"
train <- dat[dat$set == "train", ]
test  <- dat[dat$set == "test", ]
write.csv(dat, file.path(output_dir_tables, "toy_data.csv"), row.names = FALSE)

##### model part #####
### fit model
# M0
fit_lm <- lm(R ~ T + S, data = train) 

# M1
fit_spatial <- mgcv::gam(
  R ~ T + S + 
    # latent field implemented as a smoothing function using mgcv
    # here we choose this representation because it is simple and more comparable to the linear model M0
    # but it is non-Bayesian, so we will need to change to the standard INLA-SPDE approach for figure 5, 6
    s(X, Y, Z,    # 3D coordinates in km
    bs = "gp",    # Gaussian process smoother
    m = c(-3, 1200),    # parameters of the correlation function
    k = 35    # maximum complexity
  ),
  data = train,
  method = "REML",
  family = gaussian()
)

### prediction
test$pred_lm <- as.numeric(predict(fit_lm, newdata = test))
test$pred_spatial <- as.numeric(predict(fit_spatial, newdata = test))

### evaluate on test set
linear_scores  <- scores(test$R, test$pred_lm)
spatial_scores <- scores(test$R, test$pred_spatial)
linear_coef  <- coef(fit_lm)
spatial_coef <- coef(fit_spatial)

comparison_table <- data.frame(
  model = c("linear", "spatial_GP"),
  RMSE = c(linear_scores["RMSE"], spatial_scores["RMSE"]),
  MAE  = c(linear_scores["MAE"],  spatial_scores["MAE"]),
  beta_T_true = beta_T,
  beta_S_true = beta_S,  
  beta_T_estimate = c(linear_coef["T"], spatial_coef["T"]),
  beta_S_estimate = c(linear_coef["S"], spatial_coef["S"])
)

comparison_table$beta_T_error_pct <- 100 * abs(
  comparison_table$beta_T_estimate - comparison_table$beta_T_true
) / abs(comparison_table$beta_T_true)

comparison_table$beta_S_error_pct <- 100 * abs(
  comparison_table$beta_S_estimate - comparison_table$beta_S_true
) / abs(comparison_table$beta_S_true)
comparison_table[-1] <- lapply(comparison_table[-1], round, digits = 3)
print(comparison_table, row.names = FALSE)
# save table
write.csv(comparison_table, file.path(output_dir_tables, "model_comparison.csv"), row.names = FALSE)

### detailed map (figure 4)
# three species
species_colors <- c(
  "5"  = "green",  
  "10" = "blue",  
  "20" = "pink"   
)
point_colors <- unname(species_colors[as.character(dat$S)])

# annotate the examples
labels_xy <- data.frame(
  population = LETTERS[1:5],
  label_lon = c(-6, -6, 137, 141, 133),
  label_lat = c(68, 23, 50, 25, 6)
)
idx <- match(labels_xy$population, dat$population)
stopifnot(!anyNA(idx))
focus <- dat[idx, ]

draw_map(T_field, hcl.colors(100, "RdYlBu", rev = TRUE), "\u00b0C", file.path(output_dir_figures, "temperature_100_sites.png"), overlay = "species")

### visualize latent field
draw_map(D_field, hcl.colors(100, "RdYlBu", rev = TRUE), "R units", file.path(output_dir_figures, "D_field.png"), overlay = "plain")


################# figure 5 & 6 #################
### data preparation
# use the same toy example data as in figure 4
# create mesh
earth_radius_km <- 6371
loc_obs <- as.matrix(dat[, c("X", "Y", "Z")]) / earth_radius_km    # unit-sphere coordinates for INLA
mesh_globe <- 25L
mesh_cutoff <- 0.02
mesh <- fmesher::fm_rcdt_2d_inla(
  loc = loc_obs, globe = mesh_globe, cutoff = mesh_cutoff,
  crs = fmesher::fm_crs("sphere")
)

# aggregate WorldClim 10 arc-min cells to obtain a 0.5-degree map grid 
map_factor <- 3L
T_map <- terra::aggregate(T_field, fact = map_factor, fun = "mean", na.rm = TRUE)
T_map <- terra::mask(T_map, eurasia)
names(T_map) <- "T"
map_sites <- as.data.frame(T_map, xy = TRUE, cells = TRUE, na.rm = TRUE)
n_map <- nrow(map_sites)
T_ref <- mean(dat$T)             # mean of sampled T, used as reference in L

##### INLA-SPDE model #####
### parameters for prior
# range of latent field
rho0_km <- 1000
p_rho <- 0.5     # P(rho < 1000 km) = 0.5
# marginal standard deviation of latent field
sigmaW0 <- 10
p_sigmaW <- 0.05    # P(sigmaW > 10) = 0.05
# standard deviation of i.i.d. error
sigmaE0 <- 5
p_sigmaE <- 0.05    # P(sigmaE > 5) = 0.05
# select a member of the Matern family
alpha_spde <- 2L    # smoothness parameter: nu = alpha - dimension/2 = 1 on S^2
# prior of beta_0: N(0,50^2)
beta0_mean <- 0
beta0_sd <- 50
# prior of beta_T: N(0,1^2)
betaT_mean <- 0
betaT_sd <- 1
# prior of beta_S: N(0,1^2)
betaS_mean <- 0
betaS_sd <- 1

### model construction
spde <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  alpha = alpha_spde,
  constr = FALSE,              
  prior.range = c(rho0_km / earth_radius_km, p_rho),
  prior.sigma = c(sigmaW0, p_sigmaW)
)
w_index <- INLA::inla.spde.make.index("w", n.spde = spde$n.spde)
A_obs <- INLA::inla.spde.make.A(mesh, loc = loc_obs)

loc_field <- as.matrix(to_sphere(map_sites$x, map_sites$y)) / earth_radius_km
A_field <- INLA::inla.spde.make.A(mesh, loc = loc_field)

stack_obs <- INLA::inla.stack(
  data = list(R = dat$R),
  A = list(A_obs, 1),
  effects = list(
    w_index,
    data.frame(Intercept = rep(1, nrow(dat)), T = dat$T, S = dat$S)   # use all 100 samples to fit the model
  ),
  tag = "observations"
)
stack_W <- prediction_stack(A_field, rep(0, nrow(A_field)), "W")
stack_L <- prediction_stack(A_field, map_sites$T - T_ref, "L")
stack_all <- INLA::inla.stack(stack_obs, stack_W, stack_L)

### fit model
fit_inla <- INLA::inla(
  R ~ 0 + Intercept + T + S + f(w, model = spde),
  family = "gaussian",
  data = INLA::inla.stack.data(stack_all),
  control.predictor = list(A = INLA::inla.stack.A(stack_all), compute = TRUE),
  control.fixed = list(
    mean = list(Intercept = beta0_mean, T = betaT_mean, S = betaS_mean, default = 0),
    prec = list(Intercept = 1 / beta0_sd^2, T = 1 / betaT_sd^2, S = 1 / betaS_sd^2, default = 1)
  ),
  control.family = list(
    hyper = list(prec = list(
      prior = "pc.prec", param = c(sigmaE0, p_sigmaE)
    ))
  )
)

### posteriors and predictions
W_summary <- stack_summary("W", n_map)
L_summary <- stack_summary("L", n_map)
map_rows <- seq_len(n_map)
# raw maps of figure 5, 6
W_mean_map <- to_map(W_summary[map_rows, "mean"], "Posterior mean of W(s)")
W_sd_map <- to_map(W_summary[map_rows, "sd"], "Posterior SD of W(s)")
L_mean_map <- to_map(L_summary[, "mean"], "Posterior mean of L(s)")
L_sd_map <- to_map(L_summary[, "sd"], "Posterior SD of L(s)")

spde_result <- INLA::inla.spde.result(fit_inla, "w", spde)
rho_posterior <- INLA::inla.tmarginal(
  function(x) earth_radius_km * x,
  spde_result$marginals.range.nominal[[1]]
)
sigmaW_posterior <- INLA::inla.tmarginal(
  sqrt, spde_result$marginals.variance.nominal[[1]]
)
sigmaE_posterior <- INLA::inla.tmarginal(
  function(precision) 1 / sqrt(precision),
  fit_inla$marginals.hyperpar[["Precision for the Gaussian observations"]]
)

posterior_marginals <- list(
  beta_0 = fit_inla$marginals.fixed[["Intercept"]],
  beta_T = fit_inla$marginals.fixed[["T"]],
  beta_S = fit_inla$marginals.fixed[["S"]],
  rho_km = rho_posterior,
  sigma_W = sigmaW_posterior,
  sigma_epsilon = sigmaE_posterior
)
post <- t(vapply(posterior_marginals, summarise_marginal, numeric(5)))
write.csv(post, file.path(output_dir_tables, "posteriors.csv"))

### figures (posterior mean and sd of W and L)
# minor adjustments are made to unify the two figures
# this causes slight numerical differences compared to the note
r1_km <- post["rho_km", "mean"]
matern_cor <- function(distance_km, rho_km) {
  z <- sqrt(8) * distance_km / rho_km
  z * besselK(z, nu = 1)
}
r2_km <- uniroot(
  function(d) matern_cor(d, r1_km) - 0.5,
  interval = c(1e-6, r1_km)
)$root

mean_lim <- max(abs(c(
  W_summary[map_rows, "mean"], L_summary[, "mean"]
)), na.rm = TRUE)
if (!is.finite(mean_lim) || mean_lim < .Machine$double.eps) mean_lim <- 1
mean_limits <- c(-mean_lim, mean_lim)

draw_map(    # figure 6 a
  L_mean_map,
  colors = hcl.colors(100, "RdBu", rev = TRUE),
  unit = "R units",
  filename = file.path(output_dir_figures, "posterior_mean_L.png"),
  limits = mean_limits
)

draw_map(    # figure 6 b
  L_sd_map,
  colors = hcl.colors(100, "YlOrRd"),
  unit = "R units",
  filename = file.path(output_dir_figures, "posterior_sd_L.png")
)

draw_map(
  W_mean_map,
  colors = hcl.colors(100, "RdBu", rev = TRUE),
  unit = "R units",
  filename = file.path(output_dir_figures, "posterior_mean_W.png"),
  limits = mean_limits
)

draw_map(
  W_sd_map,
  colors = hcl.colors(100, "YlOrRd"),
  unit = "R units",
  filename = file.path(output_dir_figures, "posterior_sd_W.png")
)

draw_map(    # figure 5
  W_mean_map,
  colors = hcl.colors(100, "RdBu", rev = TRUE),
  unit = "R units",
  filename = file.path(output_dir_figures, "posterior_mean_W_with_circles.png"),
  limits = mean_limits,
  circle_radii = c(r1_km, r2_km)
)

