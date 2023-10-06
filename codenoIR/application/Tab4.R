#### File to derive ROC AUC and comparison tests for baseline and joint models
#### Do not reproduce exactly Table 4 because real data is not shared. 

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Warning: the intermediate results are not provided in this version of the code
# Please, refer to the README file for further details

# Executes the script that fits all the univariate joint models with the linear modeling (writes in subfolder 
# "intermediate_results/uni_fits/lin")
# Warning: time consuming (several days)... to be run on a computing cluster
source("fit.R")  

# Executes the script that fits the univariate joint models with the nonlinear modeling (writes in subfolder 
# "intermediate_results/uni_fits/nonlin")
# Warning: time consuming (several days)... to be run on a computing cluster
source("fit_nl.R") 

# Executes the script that performs the backward strategy under current scenario (uses subfolder "intermediate_results/uni_fits/" 
# and writes in subfolder "intermediate_results/multi_fits/current")
source("fit_multi_current.R")

# Executes the script that simulates individual parameters and derives ROC AUC for all landmark and horizon times under the current joint model 
# (uses subfolder "intermediate_results/multi_fits/current" and writes in subfolder intermediate_results/evalAUC")
source("eval_AUC_current.R")

# Executes the script that fits the baseline model (writes in subfolder "intermediate_results/bsl_fit")
source("fit_bsl.R") 

###############
# 2. Creating Table 4 
###############

# load survival data
surv_table = read.table("datas/sim_survival.txt", header = T)

# evaluating ROC AUC for each landmark times
roc_cur_l3 = getAUC_cur(lmk = 3, hrz = 30)
roc_cur_l6 = getAUC_cur(lmk = 6, hrz = 30)
roc_cur_l9 = getAUC_cur(lmk = 9, hrz = 30)
roc_bsl_l3 = getAUC_bsl(lmk = 3, hrz = 30)
roc_bsl_l6 = getAUC_bsl(lmk = 6, hrz = 30)
roc_bsl_l9 = getAUC_bsl(lmk = 9, hrz = 30)

# get ROC AUC value 
rcl3 = round(roc_cur_l3$AUC_2[2], 2)
rcl6 = round(roc_cur_l6$AUC_2[2], 2)
rcl9 = round(roc_cur_l9$AUC_2[2], 2)
rbl3 = round(roc_bsl_l3$AUC_2[2], 2)
rbl6 = round(roc_bsl_l6$AUC_2[2], 2)
rbl9 = round(roc_bsl_l9$AUC_2[2], 2)

# get 95% confidence interval 
CIc_l3 = paste0("[", round(confint(roc_cur_l3)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_cur_l3)$CI_AUC_2[2]/100, 2), "]")
CIc_l6 = paste0("[", round(confint(roc_cur_l6)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_cur_l6)$CI_AUC_2[2]/100, 2), "]")
CIc_l9 = paste0("[", round(confint(roc_cur_l9)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_cur_l9)$CI_AUC_2[2]/100, 2), "]")
CIb_l3 = paste0("[", round(confint(roc_bsl_l3)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_bsl_l3)$CI_AUC_2[2]/100, 2), "]")
CIb_l6 = paste0("[", round(confint(roc_bsl_l6)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_bsl_l6)$CI_AUC_2[2]/100, 2), "]")
CIb_l9 = paste0("[", round(confint(roc_bsl_l9)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_bsl_l9)$CI_AUC_2[2]/100, 2), "]")

# to print results in the table, paste ROC value + confidence interval 
roc_cl3 = paste0(rcl3, CIc_l3)
roc_cl6 = paste0(rcl6, CIc_l6)
roc_cl9 = paste0(rcl9, CIc_l9)
roc_bl3 = paste0(rbl3, CIb_l3)
roc_bl6 = paste0(rbl6, CIb_l6)
roc_bl9 = paste0(rbl9, CIb_l9)

# comparing ROc AUC for all landmark times
test_l3 = signif(compare(roc_cur_l3, roc_bsl_l3)$p_values_AUC_2[2], 2)
test_l6 = signif(compare(roc_cur_l6, roc_bsl_l6)$p_values_AUC_2[2], 2)
test_l9 = signif(compare(roc_cur_l9, roc_bsl_l9)$p_values_AUC_2[2], 2)

# formatting data 
label = c("AUC[95%CI] - baseline model", "AUC[95%CI] - joint model", "p-value")
tab4 = data.frame(HorizonTime30 = label, LandmarkDay3 = c(roc_bl3, roc_cl3, test_l3), LandmarkDay6 = c(roc_bl6, roc_cl6, test_l6), LandmarkDay9 = c(roc_bl9, roc_cl9, test_l9))

save(tab4, file = "../results/Tab4.RData")