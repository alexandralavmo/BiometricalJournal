#### File used to derive rocAUC for all landmark and horizon times for the final joint model

###############
# 1. Generation landmark datasets
###############
# landmark data set (Monolix format)
v = expand.grid(lmk = c(3, 6, 9))
mapply(gen_landmark_curr, v$lmk)

# landmark survival data set (R format)
v = expand.grid(lmk = c(3, 6, 9))
mapply(gen_survdata_curr, v$lmk)


###############
# 2. Bayesian individual dynamic predictions : simulation of individual parameters for all datasets and all landmarks
###############

v = expand.grid(lmk = c(3, 6, 9))
mapply(gen_indv_params_curr, v$lmk)

###############
# 3. AUC derivation 
###############

# Computing ROC AUC for each landmark and each horizon time 
roc3 = sapply(4:30, function(j) getAUC(lmk = 3, hrz = j))
roc6 = sapply(7:30, function(j) getAUC(lmk = 6, hrz = j))
roc9 = sapply(10:30, function(j) getAUC(lmk = 9, hrz = j))

# Saving results
save(roc3, file = "intermediate_results/evalAUC/rocJM3.RData")
save(roc6, file = "intermediate_results/evalAUC/rocJM6.RData")
save(roc9, file = "intermediate_results/evalAUC/rocJM9.RData")


