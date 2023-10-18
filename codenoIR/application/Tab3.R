####  File used to derive parameter estimates (current multivariate joint model)
#### Do not reproduce excatly the Table 3 because real data is not shared

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
}

###############
# 2. Creating Table 2 
###############

### Name of parameters ###

name_bm2 = c("mu0n", "mu1n", "omega0n", "omega1n", "sigmaan")
name_bm3 = c("mu0c", "mu1c", "omega0c", "omega1c", "sigmaac")
name_surv = c("h1", "alpha1p", "alpha1c", "beta1", "h2", "alpha2p", "alpha2c", "beta2")

### Load result file ###
data = read.table("intermediate_results/multi_fits/current/back2/mod/populationParameters.txt", header = T, sep = ',')

### Re-organized data ###
datanew = data[c(1:2, 15:16, 19, 3:4, 17:18, 20, 5, 7:8, 12, 6, 9:10, 14), ]

colnames(datanew) = c("Parameter", "Value", "SE", "RSE(%)", "p-value")
datanew[, 1] = c(name_bm2, name_bm3, name_surv)

datanew$Value = signif(datanew$Value, digits = 2)
datanew$SE = signif(datanew$SE, digits = 2)
datanew$`RSE(%)` = signif(datanew$`RSE(%)`, digits = 2)

### Compute p-value of link coefficients (not done in Monolix)
subdata = datanew[datanew$Parameter %in% c("alpha1p", "alpha1c", "alpha2p", "alpha2c"), ] 
subdata$`p-value` = signif(2*pnorm(abs(subdata$Value)/subdata$SE, lower.tail = F), 2)
datanew[datanew$Parameter %in% c("alpha1p", "alpha1c", "alpha2p", "alpha2c"), ] = subdata

tab3 = datanew

save(tab3, file = "../results/Tab3.RData")