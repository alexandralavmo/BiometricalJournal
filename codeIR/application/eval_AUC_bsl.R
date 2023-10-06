#### File used to derive rocAUC for all landmark and horizon times for the baseline model

###############
# 1. Generation of directory and landmark datasets
###############

surv_table = read.table("datas/sim_survival.txt", h = T)

roc_bsl3 = sapply(4:30, function(j) getAUC_b(lmk = 3, hrz = j))
roc_bsl6 = sapply(7:30, function(j) getAUC_b(lmk = 6, hrz = j))
roc_bsl9 = sapply(10:30, function(j) getAUC_b(lmk = 9, hrz = j))

# save the results 
save(roc_bsl3, file = "intermediate_results/evalAUC/rocBS3.RData")
save(roc_bsl6, file = "intermediate_results/evalAUC/rocBS6.RData")
save(roc_bsl9, file = "intermediate_results/evalAUC/rocBS9.RData")