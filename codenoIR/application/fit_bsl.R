#### File used to fit the baseline model involving 4C Score 

###############
# 1. Copying files in the correct directory
###############

file.copy("datas/sim_surv_monolix.txt", "intermediate_results/bsl_fit/sim_surv_monolix.txt")
file.copy("modfiles/mod_bsl.txt", "intermediate_results/bsl_fit/mod_bsl.txt")
file.copy("modfiles/mod_bsl.mlxtran", "intermediate_results/bsl_fit/mod_bsl.mlxtran")

###############
# 2. Run estimation 
###############

loadProject("intermediate_results/bsl_fit/mod_bsl.mlxtran")
runPopulationParameterEstimation()
runStandardErrorEstimation()

