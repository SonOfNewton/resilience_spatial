### script to run experiments on scale1 structured density dependent spatial models
output_dir <- here("output", "figures", "note3", "scale1_structure_dd_spatial")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# parallel processing
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- 1L
n_cores <- max(1L, detected_cores - 1L)

# runs experiment on 12 sets of parameter settings
scenarios <- list(
  # scenario 1 to 4: standard parameters, dd on one parameter at a time
  scenario_1 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,    # MPM parameters
    k1 = 0.01, k2 = 0, k3 = 0, k4 = 0,                                     # density dependence of each MPM parameter
    competition_mode = 2L,                                                 # values: 1 means only consider adults as competition source; 2 means consider adults and juveniles 
    Sm = 0.5, Sw = 0.01,                                                   # spatial parameters
    L = 1, initial_n = 100L, adult_fraction = 0.2,                         # set initial state
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                     # burn-in period
    post_steps = 30L, n_rep = 50L,                                         # after disturbance
    abundance_factors = 1, cleared_area = 0.5,                             # start with how many individuals; intensity of LC
    stage_specific_lc = FALSE,                                             # FALSE: only adults / only juveniles; TRUE: keep J/A proportion, do LC on only J / only A
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_2 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0.01, k3 = 0, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.5, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_3 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0, k3 = 0.01, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.5, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_4 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0, k3 = 0, k4 = 0.01,                                 
    competition_mode = 2L,        
    Sm = 0.5, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  # dd on all four parameters
  scenario_5 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0.01, k2 = 0.01, k3 = 0.01, k4 = 0.01,                                 
    competition_mode = 2L,        
    Sm = 0.5, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  # compare with scenario 1, changed Sm from 0.5 to 0.1
  scenario_6 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0.01, k2 = 0, k3 = 0, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.1, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  # compare with scenario 1, changed Sm from 0.5 to 0.1
  scenario_7 = list(
    sigma1_max = 0.9, sigma2_max = 0.9, gamma_max = 0.95, phi_max = 1,  
    k1 = 0.01, k2 = 0, k3 = 0, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.1, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  # compare with scenario 1, changed competition_mode from 2 to 1
  scenario_8 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0.01, k2 = 0, k3 = 0, k4 = 0,                                 
    competition_mode = 1L,        
    Sm = 0.5, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = FALSE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  # scenarios 9 to 12 mimics scenarios 1 to 4 respectively, changing stage_specific_lc from FALSE to TRUE
  scenario_9 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0.01, k2 = 0, k3 = 0, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.05, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = TRUE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_10 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0.01, k3 = 0, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.05, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = TRUE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_11 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0, k3 = 0.01, k4 = 0,                                 
    competition_mode = 2L,        
    Sm = 0.05, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = TRUE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  ),
  
  scenario_12 = list(
    sigma1_max = 0.7, sigma2_max = 0.6, gamma_max = 0.8, phi_max = 1.2,  
    k1 = 0, k2 = 0, k3 = 0, k4 = 0.01,                                 
    competition_mode = 2L,        
    Sm = 0.05, Sw = 0.01,  
    L = 1, initial_n = 100L, adult_fraction = 0.2,                      
    pre_steps = 50L, equilibrium_window = c(40L, 50L),                
    post_steps = 30L, n_rep = 50L,                                     
    abundance_factors = 1, cleared_area = 0.5, stage_specific_lc = TRUE,                         
    kernel_normalized = TRUE, cutoff_sd = 3,  
    seed = 1, max_population = 20000L
  )
)

results <- list()

# run
for (scenario_name in names(scenarios)) {
  z <- scenarios[[scenario_name]]
  message("Running ", scenario_name)
  results[[scenario_name]] <- ss_run_experiment(z, n_cores = n_cores)
  ss_save_result(results[[scenario_name]], output_dir, scenario_name)
}
message("Figures saved to ", output_dir)
