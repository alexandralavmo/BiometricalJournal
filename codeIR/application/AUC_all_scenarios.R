#### File to compute ROC AUC for landmark = D9 and horizon = D30 for all scenarios of the sensibility analysis on thresholds used in univariate modeling

###############
# 1. generating landmark data sets
###############

gen_landmark_flex(9)
gen_landmark_strin(9)
gen_survdata_comp(9)

###############
# 2. Simulating individual parameters for flexible and stringent scenario
############### 

gen_indv_params_flexible(9)
gen_indv_params_stringent(9)

###############
# 3. Computing ROC AUC for landmark = D9 and horizon time = D30 with 95% confidence intervals 
###############

roc_flex = getAUC_flex(lmk = 9, hrz = 30)
roc_strin = getAUC_strin(lmk = 9, hrz = 30)
roc_cur = getAUC_cur(lmk = 9, hrz = 30)

# saving results
save(roc_flex, file = "intermediate_results/compAUC/roc_flex.RData")
save(roc_strin, file = "intermediate_results/compAUC/roc_strin.RData")
save(roc_cur, file = "intermediate_results/compAUC/roc_cur.RData")
