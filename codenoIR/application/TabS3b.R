#### File to derive results of sensibility analysis on thresholds for univariate modeling
#### Do not reproduce exactly Table S3b because real data is not shared. 

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# set interresults = TRUE if you want to rely on intermediate results (beware that the intermediate results are not provided in this version)
# the version with intermediate results is available at https://github.com/alexandralavmo/BiometricalJournal/tree/main/codeIR)

interresults <- FALSE
if (!interresults) {
source("fit.R") # Executes this script that fits all the univariate joint models with the linear modeling (writes in subfolder "intermediate_results/uni_fits/lin")
		# Warning: time consuming... to be run on a computing cluster
source("fit_nl.R") # Executes the script that fits the univariate joint models with the nonlinear modeling (writes in subfolder "intermediate_results/uni_fits/nonlin")
		   # Warning: time consuming... to be run on a computing cluster
source("fit_multi_current.R") # Executes the script that performs the backward strategy under current scenario (uses subfolder "intermediate_results/uni_fits/" and 
			      # writes in subfolder "intermediate_results/multi_fits/current")
source("fit_multi_flex.R") # Executes the script that performs the backward strategy under the current scenario (uses subfolder "intermediate_results/uni_fits/" and writes 
			   # in subfolder "intermediate_results/multi_fits/current")
			   # Warning: time consuming... to be run on a computing cluster
source("eval_AUC_current.R")  # Executes the script that simulates individual parameters and derives ROC AUC for all landmark and horizon times under the current joint model 
			      # (uses subfolder "intermediate_results/multi_fits/current" and writes in subfolder "intermediate_results/evalAUC")
source("fit_multi_strin.R") # Executes the script that performs the backward strategy under the stringent scenario (uses subfolder "intermediate_results/uni_fits/" and writes 
			    # in subfolder "intermediate_results/multi_fits/stringent")
source("AUC_all_scenarios.R") # Executes the script that computes ROC AUC for all scenarios of the sensibility analysis (uses subfolder "intermediate_results/multi_fits/" and 
			      # writes in subfolder "intermediate_results/compAUC")
}
 

###############
# 2. Creating Table S3b 
###############

# Loading results 

load("intermediate_results/compAUC/roc_flex.RData")
load("intermediate_results/compAUC/roc_strin.RData")
load("intermediate_results/compAUC/roc_cur.RData")

rf = round(roc_flex$AUC_2[2], 2)
rc = round(roc_cur$AUC_2[2], 2)
rs = round(roc_strin$AUC_2[2], 2)

ci_f = paste0("[", round(confint(roc_flex)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_flex)$CI_AUC_2[2]/100, 2), "]")
ci_c = paste0("[", round(confint(roc_cur)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_cur)$CI_AUC_2[2]/100, 2), "]")
ci_s = paste0("[", round(confint(roc_strin)$CI_AUC_2[1]/100, 2), ", ", round(confint(roc_strin)$CI_AUC_2[2]/100, 2), "]")

roc_f = paste0(rf, ci_f)
roc_c = paste0(rc, ci_c)
roc_s = paste0(rs, ci_s)

# comparing ROc AUC with current scenario
test_cf = compare(roc_cur, roc_flex)
test_cs = compare(roc_cur, roc_strin)


# Load the results of univariate analysis
# Creating data using results of univariate analysis

load("intermediate_results/uni_fits/selec_lin.RData")
load("intermediate_results/uni_fits/selec_nlin.RData")
sl = length(selec_lin)
snl = length(selec_nlin)
load("intermediate_results/uni_fits/selec_lin_flex.RData")
load("intermediate_results/uni_fits/selec_nlin_flex.RData")
sl_f = length(selec_lin)
snl_f = length(selec_nlin)
load("intermediate_results/uni_fits/selec_lin_stringent.RData")
load("intermediate_results/uni_fits/selec_nlin_stringent.RData")
sl_s = length(selec_lin)
snl_s = length(selec_nlin)


# some manual steps (number of biomarkers included in the final model, and their names)

biom_lin = c("Number of linear biomarkers ", sl_s, sl, sl_f) # biomarker selected with linear model for stringent, current and flexible scenario
biom_nonlin = c("Number of nonlinear biomarkers", snl_s, snl, snl_f) # biomarker selected with nonlinear model for stringent, current and flexible scenario
biom_multi = c("Number of biomarkers selected for the multivariable analysis", 1, 3, 6)
biom_final = c("Number of biomarkers included in the final model", 1, 2, 2)
name_biom = c("Biomarkers included in the final model", "Biomarker 7", "Biomarkers 24 and 43", "Biomarkers 24 and 46")
roc = c("ROC AUC [95%] CI(lmk = D9, horizon = D30)", roc_s, roc_c, roc_f)
pval = c("p-value comparing ROC AUC with the current scenario", signif(test_cs$p_values_AUC_2[2], 2), "", signif(test_cf$p_values_AUC_2[2], 2))

table = data.frame(NULL = NA, Stringent_scenario = NA, Current_scenario = NA, Flexible_scenario = NA)
table = rbind(table, biom_lin, biom_nonlin, biom_multi, biom_final, name_biom, roc, pval)
tabS3b = table[-1, ]
save(tabS3b, file = "../results/TabS3b.RData")