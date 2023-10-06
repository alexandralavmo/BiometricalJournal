# This code is used to compute ROC AUC for both true and final joint models

###############
# 1. Generation of landmark datasets
###############

v = expand.grid(lmk = 9, j = 1:p)
mapply(gen_landmark, v$lmk, v$j)

v = expand.grid(lmk = 9, j = 1:p)
mapply(gen_survdata, v$lmk, v$j)

###############
# 2. Bayesian individual dynamic predictions : simulation of individual parameters for all datasets and all landmarks
###############

v = expand.grid(lmk = 9, j = 1:p, truemod = c(T, F))
mapply(gen_indv_params, v$lmk, v$j, v$truemod)

###############
# 3. AUC derivation 
###############

# To reproduce Figure 7, need to get AUC for landmark time = 9 and horizon time = 30
trueroc = sapply(1:p, function(j) getAUC(lmk = 9, hrz = 30, j = j, truemod = T))
finalroc = sapply(1:p, function(j) getAUC(lmk = 9, hrz = 30, j = j, truemod = F))

# Saving results
save(trueroc, file = "intermediate_results/evalAUC/trueroc.RData")
save(finalroc, file = "intermediate_results/evalAUC/finalroc.RData")
