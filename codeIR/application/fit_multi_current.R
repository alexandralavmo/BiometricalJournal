#### File used to perform the backward on current scenario 

###############
# 1. Loading and formating data to keep the three biomarkers involved 
###############

data = read.table("datas/alldata.txt", header = T)

load("intermediate_results/uni_fits/selec_lin.RData")
load("intermediate_results/uni_fits/selec_nlin.RData")

data3biom = data[data$ytype %in% c(selec_lin, selec_nlin, 60, 61), ]
data3biom$ytype[data3biom$ytype == 7] = 1
data3biom$ytype[data3biom$ytype == 24] = 2
data3biom$ytype[data3biom$ytype == 43] = 3
data3biom$ytype[data3biom$ytype == 60] = 4
data3biom$ytype[data3biom$ytype == 61] = 5
write.table(data3biom, "intermediate_results/multi_fits/current/back3/data.txt", row.names = F)


###############
# 2. Preparing the back3 directory and fit full model
###############

file.copy("modfiles/mod_123.txt", "intermediate_results/multi_fits/current/back3/mod_123.txt")
file.copy("modfiles/mod123.mlxtran", "intermediate_results/multi_fits/current/back3/mod.mlxtran")

loadProject("intermediate_results/multi_fits/current/back3/mod.mlxtran")
runPopulationParameterEstimation()
runStandardErrorEstimation()

###############
# 3. Prepare all possible models in each directory (e.g. /back2 will contain all possible models with 2 biomarkers
###############

sapply(1:2, create_combtxt_curr)

###############
# 3. Backward selection procedure
###############

# Step with 2 biomarkers and next

# create data, model txt and mlxtran files given the results of the backward with 3 biomarkers 
iter_back_curr(2)
# fit the model with 2 biomarkers
fit_all_models_curr(2)

# create data, model txt and mlxtran files given the results of the backward with 2 biomarkers 
iter_back_curr(1)
# fit the model with one biomarker
fit_all_models_curr(1)