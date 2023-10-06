# This code performs the backward strategy for the first scenario of correlation

###############
# 1. Prepare all possible models in each directory (e.g. /back3 will contain all possible models with 3 biomarkers
###############

# Fills all directories with all possible model files of the correct size
sapply(1:6, create_combtxt_scen1) # for back6 to back1 directories


###############
# 2. Preparing the back7 directory and fit full model to every simulated datasets
###############

#copying data files
file.copy(paste0("data/", list.files("data/")), "intermediate_results/backward/back7/")

# creating all model text files
sapply(1:p, function(j) adapt_txt_scen1(c(1, 2, 3, 4, 5, 6, 7), j))

# creating all the mlxtran files
sapply(1:p, setmodel_scen1)

# fit all models with 7 biomarkers 
sapply(1:p, fit_fullmodels_scen1)


###############
# 3. Backward selection procedure
###############

# Step with 6 biomarkers and next
iter_back_scen1(6) # create all mlxtran files from the previous step depending on which biomarker has been removed
fit_all_models_scen1(6) # fits all resulting models

iter_back_scen1(5)
fit_all_models_scen1(5)

iter_back_scen1(4)
fit_all_models_scen1(4)

iter_back_scen1(3)
fit_all_models_scen1(3)

iter_back_scen1(2)
fit_all_models_scen1(2)

iter_back_scen1(1)
fit_all_models_scen1(1)

iter_back_scen1(0)
