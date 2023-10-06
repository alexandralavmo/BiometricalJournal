#### File used to perform the backward on stringent scenario 

# Only one biomarker selected

## Copying files in the correct directory 
file.copy("intermediate_results/uni_fits/lin/mod_7.txt", "intermediate_results/multi_fits/stringent/mod_7.txt", overwrite = T)
file.copy("intermediate_results/uni_fits/datas/data7.txt", "intermediate_results/multi_fits/stringent/data7.txt", overwrite = T)

## change data file name in the mlxtran file
mod = readLines("intermediate_results/uni_fits/lin/mod7.mlxtran")
mod[4] = "file = 'data7.txt'"
write.table(mod, file = "intermediate_results/multi_fits/stringent/mod7.mlxtran", quote = F, row.names = F, col.names = F)

## run estimation 
loadProject(paste0("intermediate_results/multi_fits/stringent/mod7.mlxtran"))
runPopulationParameterEstimation()
