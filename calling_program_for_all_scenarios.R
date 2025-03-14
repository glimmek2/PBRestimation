library(clustermq)
source("simulation_functions_myseed.R")
source("analysis_functions_sep_robustvar_fix.R")
source("oneSimRun.R")

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Define the wrapper function for multiple scenarios
run_multiple_simulations <- function(scenarios, num_cores, num_sims_per_core, alpha) {
  # Create a folder for the simulations if it doesn't exist
  if (!dir.exists("SIMRES")) {
    dir.create("SIMRES")
  }
  
  for (scenario in scenarios) {
    scenario_C <- scenario$scenario_C
    scenario_T <- scenario$scenario_T
    scenario_name <- scenario$scenario_name
    
    scenario <- list(scenario_C = scenario_C, scenario_T = scenario_T, scenario_name = scenario_name)
    
    beginsim <- Sys.time()
    sim_results <- run_simulation_test(scenario=scenario, fun_simArm=simulate1trialarm, getPBIRcurve=getPBIRcurve, robustvar=robustvar,
                                       num_cores=num_cores, num_sims_per_core=num_sims_per_core, alpha=alpha)
    endsim <- Sys.time()
    cat("Simulation duration for scenario", scenario_name, ":", endsim - beginsim, "\n")
    
    # Create subfolder for the current scenario
    folder_name <- file.path(getwd(), "SIMRES", paste0("scenario_", scenario$scenario_name))
    if (!dir.exists(folder_name)) {
      dir.create(folder_name)
    }
    
    # Generate filenames for the CSV files
    file1 <- paste0(folder_name, "/scensetting.csv")
    file2 <- paste0(folder_name, "/simoutput.csv")
    file3 <- paste0(folder_name, "/sim_summaries.csv")
    
    # Write data to CSV files
    write.csv(scenario, file1, row.names = FALSE)
    write.csv(sim_results$allruns, file2, row.names = FALSE)
    write.csv(sim_results$summary_sims, file3, row.names = FALSE)
  }
}

# Define your scenarios
# base scenarios (without sample size)
base_scenarios <- list(
  list(scenario_name = "N1",
       lambda01 = 1, lambda02 = 0.6, lambda12raw = 0.5, slope12 = 0.4, lambda_cens = 0.3),
  list(scenario_name = "N2",
       lambda01 = 1, lambda02 = 0.6, lambda12raw = 0.5, slope12 = 0.0, lambda_cens = 0.3),
  list(scenario_name = "N3",
       lambda01 = 1, lambda02 = 0.6, lambda12raw = 0.5, slope12 = 0.8, lambda_cens = 0.3),
  list(scenario_name = "N4",
       lambda01 = 0.5, lambda02 = 0.6, lambda12raw = 0.8, slope12 = 0.4, lambda_cens = 0.3),
  list(scenario_name = "N5",
       lambda01 = 1.5, lambda02 = 0.6, lambda12raw = 0.8, slope12 = 0.4, lambda_cens = 0.3),
  list(scenario_name = "N6",
       lambda01 = 1.5, lambda02 = 0.6, lambda12raw = 0.2, slope12 = 0.1, lambda_cens = 0.3),
  list(scenario_name = "A1",
       scenario_C = list(lambda01 = 1, lambda02 = 0.6, lambda12raw = 0.5, slope12 = 0.4, lambda_cens = 0.3),
       scenario_T = list(lambda01 = 1.2, lambda02 = 0.5, lambda12raw = 0.3, slope12 = 0.2, lambda_cens = 0.4)),
  list(scenario_name = "A2",
       scenario_C = list(lambda01 = 1, lambda02 = 0.6, lambda12raw = 0.5, slope12 = 0.0, lambda_cens = 0.3),
       scenario_T = list(lambda01 = 1.5, lambda02 = 0.5, lambda12raw = 0.2, slope12 = 0.1, lambda_cens = 0.3)),
  list(scenario_name = "A3",
       scenario_C = list(lambda01 = 1.4, lambda02 = 0.6, lambda12raw = 0.7, slope12 = 0.2, lambda_cens = 0.3),
       scenario_T = list(lambda01 = 0.9, lambda02 = 0.6, lambda12raw = 0.15, slope12 = 0.1, lambda_cens = 0.3)),
  list(scenario_name = "A4",
       scenario_C = list(lambda01 = 0.8, lambda02 = 0.6, lambda12raw = 0.4, slope12 = 0.3, lambda_cens = 0.3),
       scenario_T = list(lambda01 = 1.5, lambda02 = 0.6, lambda12raw = 0.7, slope12 = 0.2, lambda_cens = 0.3))
)

# Define sample sizes
sample_sizes <- c(100, 200, 300)

# Function to add sample size to base scenario
generate_scenario <- function(base_scenario, sample_size) {
  if (startsWith(base_scenario$scenario_name, "N")) {
    scenario_C <- data.frame(npatients = sample_size, trtname = "C", 
                             lambda01 = base_scenario$lambda01, lambda02 = base_scenario$lambda02, 
                             lambda12raw = base_scenario$lambda12raw, slope12 = base_scenario$slope12, 
                             lambda_cens = base_scenario$lambda_cens)
    scenario_T <- data.frame(npatients = sample_size, trtname = "T", 
                             lambda01 = base_scenario$lambda01, lambda02 = base_scenario$lambda02, 
                             lambda12raw = base_scenario$lambda12raw, slope12 = base_scenario$slope12, 
                             lambda_cens = base_scenario$lambda_cens)
    scenario_name <- paste0(base_scenario$scenario_name, "_", sample_size)
  } else {
    scenario_C <- data.frame(npatients = sample_size, trtname = "C", 
                             lambda01 = base_scenario$scenario_C$lambda01, lambda02 = base_scenario$scenario_C$lambda02, 
                             lambda12raw = base_scenario$scenario_C$lambda12raw, slope12 = base_scenario$scenario_C$slope12, 
                             lambda_cens = base_scenario$scenario_C$lambda_cens)
    scenario_T <- data.frame(npatients = sample_size, trtname = "T", 
                             lambda01 = base_scenario$scenario_T$lambda01, lambda02 = base_scenario$scenario_T$lambda02, 
                             lambda12raw = base_scenario$scenario_T$lambda12raw, slope12 = base_scenario$scenario_T$slope12, 
                             lambda_cens = base_scenario$scenario_T$lambda_cens)
    scenario_name <- paste0(base_scenario$scenario_name, "_", sample_size)
  }
  list(scenario_name = scenario_name, scenario_C = scenario_C, scenario_T = scenario_T)
}

# Generate all scenarios with different sample sizes
scenarios <- unlist(
  lapply(base_scenarios, function(base_scenario) {
    lapply(sample_sizes, function(sample_size) {
      generate_scenario(base_scenario, sample_size)
    })
  }),
  recursive = FALSE
)

# Run the wrapper function
run_multiple_simulations(scenarios = scenarios, num_cores = 10, num_sims_per_core = 1000, alpha = 0.025)