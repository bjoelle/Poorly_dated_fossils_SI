# simulating the main datasets can be done by adjusting parameters in main_sim.R then:
source("aux_functions.R")
source("sim.deposit.values.R")
source("main_sim.R")

save_folder = "path/to/folder/to/write/simulated/data"
run_simulation(save_folder)

# other additional datasets can be simulated with 
source("aux_functions.R")
source("sim.deposit.values.R")
source("extended_main_sim.R")

save_folder = "path/to/folder/to/write/simulated/data"
run_simulation_extended(save_folder)