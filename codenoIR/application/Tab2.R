#### Used to reproduce baseline model parameter estimates 
#### Do not reproduce excatly the Table 2 because real data is not shared. 

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

# Executes the script that fits the baseline model (writes in subfolder "intermediate_results/bsl_fit")
source("fit_bsl.R")  

###############
# 2. Creating Table 2 
###############

### Name of parameters ###

name_surv = c("h1", "beta1", "h2", "beta2")

### Load result file ###
data = read.table("intermediate_results/bsl_fit/mod_bsl/populationParameters.txt", header = T, sep = ',')

### Re-organized data ###
datanew = data[c(1, 4, 2, 6), ]

colnames(datanew) = c("Parameter", "Value", "SE", "RSE(%)", "p-value")
datanew[, 1] = name_surv

datanew$`p-value` = ifelse(datanew$`p-value`<10e-05, "<10e-05", datanew$`p-value`)

datanew$Value = signif(datanew$Value, digits = 2)
datanew$SE = signif(datanew$SE, digits = 2)
datanew$`RSE(%)` = signif(datanew$`RSE(%)`, digits = 2)
tab2 = datanew

save(tab2, file = "../results/Tab2.RData")