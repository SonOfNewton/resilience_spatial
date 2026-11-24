### script to run experiments on scale1 nonstructured density dependent spatial models
output_dir <- here("output", "figures", "note3", "scale1_nonstructure_dd_spatial")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# parallel processing
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- 1L
n_cores <- max(1L, detected_cores - 1L)

# runs experiment on 12 sets of parameter settings
scenarios <- list(
# Slow life history, long dispersal. Expect LC to recover faster than DT.
scenario_1 = list(
  b = 5,                      # birth rate per individual
  d = 0.95,                    # "adult" survival rate (note that it's not death rate)
  d_prime = 0.001,             # competition coefficient for "adults"
  d_double_prime = 0.3,        # competition coefficient for "seeds"
  g = 0.2,                     # baseline survival rate for "seeds"
  Sm = 0.5,                    # spatial parameter for dispersal
  Sw = 0.01,                   # spatial parameter for competition
  L = 1,                       # side length of simulation region
  initial_n = 400L,             # initial number of individuals
  pre_steps = 20L,             # run steps before disturbance
  post_steps = 20L,            # run steps after disturbance
  n_rep = 50L,                 # number of repeated experiments
  p = 0.2,                     # disturbance intensity
  seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir),   # cutoff the normal distribution at +- 3 sd

# Slow life history, short dispersal. Expect DT to recover faster than LC.
scenario_2 = list(
  b = 5, d = 0.95, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,    # tests show that the stable equilibrium is around 500-600 individuals
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

# Fast life history, long dispersal. Expect LC and DT to have similar recovery.
scenario_3 = list(
  b = 25, d = 0.2, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.5, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

# Fast life history, short dispersal. Expect DT to recover faster than LC, but the difference to be smaller than scenario 2.
scenario_4 = list(
  b = 25, d = 0.2, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

# scenario 5-8 mimics scenario 1-4, but alters the population to be seed limited (1-4 are establishment limited). Expect DT to almost always recover faster than LC.
scenario_5 = list(
  b = 0.2, d = 0.95, d_prime = 0.001, d_double_prime = 0, g = 1,    # compare with scenario 1: b = 10, d_double_prime = 0.3, g = 0.1
  Sm = 0.5, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

scenario_6 = list(
  b = 0.2, d = 0.95, d_prime = 0.001, d_double_prime = 0, g = 1,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3,output_dir = output_dir
),

scenario_7 = list(   # scenario 7 and 8 are quite unstable, likely due to reduced effect of the negative feedback of dd on death rate because of high baseline death rate 
  b = 0.9, d = 0.2, d_prime = 0.01, d_double_prime = 0, g = 1,
  Sm = 0.5, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

scenario_8 = list(
  b = 0.9, d = 0.2, d_prime = 0.01, d_double_prime = 0, g = 1,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 20L, n_rep = 50L,
  p = 0.2, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

# scenario 9-12 mimics scenario 1-4, but alters the disturbance intensity.
scenario_9 = list(
  b = 5, d = 0.95, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,   
  Sm = 0.5, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 30L, n_rep = 50L,
  p = 0.7, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir   # compare with scenario 1: p = 0.2
),

scenario_10 = list(
  b = 5, d = 0.95, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 30L, n_rep = 50L,
  p = 0.7, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

scenario_11 = list(
  b = 25, d = 0.2, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.5, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 30L, n_rep = 50L,
  p = 0.7, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
),

scenario_12 = list(
  b = 25, d = 0.2, d_prime = 0.001, d_double_prime = 0.3, g = 0.2,
  Sm = 0.05, Sw = 0.01, L = 1, initial_n = 400L,
  pre_steps = 20L, post_steps = 30L, n_rep = 50L,
  p = 0.7, seed = 1, kernel_normalized = TRUE, cutoff_sd = 3, output_dir = output_dir
)
)

all_results <- list()

# run
for (scenario_name in names(scenarios)) {
  message("\nRunning ", scenario_name)
  params <- scenarios[[scenario_name]]
  results <- run_experiment(params, n_cores=n_cores)
  all_results[[scenario_name]] <- results
  # save the mean trajectory
  figure_file <- file.path(params$output_dir, paste0(scenario_name, "_mean_trajectory.pdf"))
  pdf(figure_file, width = 10, height = 6.3)
  tryCatch(plot_trajectories(results), finally = dev.off())
  
  # (optional) visualize the situation right after disturbance
  save_spatial_snapshots(results, scenario_name = scenario_name)
}
message("Figures saved to ", output_dir)
