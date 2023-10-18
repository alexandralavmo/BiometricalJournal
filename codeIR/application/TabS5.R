#### File to get estimates of full initial joint model (current scenario)
#### Do not reproduce exactly Table S5 because real data is not shared. 

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# set interresults = FALSE if you want to (re-) compute the intermediate results

interresults <- TRUE
if (!interresults) {
source("fit.R") # Executes this script that fits all the univariate joint models with the linear modeling (writes in subfolder "intermediate_results/uni_fits/lin")
		# Warning: time consuming... to be run on a computing cluster
source("fit_nl.R") # Executes the script that fits the univariate joint models with the nonlinear modeling (writes in subfolder "intermediate_results/uni_fits/nonlin")
		   # Warning: time consuming... to be run on a computing cluster
source("fit_multi_current.R") # Executes the script that performs the backward strategy under current scenario (uses subfolder "intermediate_results/uni_fits/" and 
			      # writes in subfolder "intermediate_results/multi_fits/current")
}

###############
# 2. Creating Table S5
###############

### Name of parameters ###

name_bm1 = c("mu0n", "mu1n", "omega0n", "omega1n", "sigmaan")
name_bm2 = c("mu0p", "mu1p", "omega0p", "omega1p", "sigmaap")
name_bm3 = c("mu0c", "mu1c", "omega0c", "omega1c", "sigmaac")
name_surv = c("h1", "alpha1n", "alpha1p", "alpha1c", "beta1", "h2", "alpha2n", "alpha2p", "alpha2c", "beta2")

### Load result file ###
data = read.table("intermediate_results/multi_fits/current/back3/mod/populationParameters.txt", header = T, sep = ',')

### Re-organized data ###
datanew = data[c(1:2, 19:20, 25, 3:4, 21:22, 26, 5:6, 23:24, 27, 7, 9:11, 16, 8, 12:14, 18), ]

colnames(datanew) = c("Parameter", "Value", "SE", "RSE(%)", "p-value")
datanew[, 1] = c(name_bm1, name_bm2, name_bm3, name_surv)

datanew$Value = signif(datanew$Value, digits = 2)
datanew$SE = signif(datanew$SE, digits = 2)
datanew$`RSE(%)` = signif(datanew$`RSE(%)`, digits = 2)

### Compute p-value of link coefficients (not done in Monolix)
subdata = datanew[datanew$Parameter %in% c("alpha1n", "alpha1p", "alpha1c", "alpha2n", "alpha2p", "alpha2c"), ] 
subdata$`p-value` = signif(2*pnorm(abs(subdata$Value)/subdata$SE, lower.tail = F), 2)
datanew[datanew$Parameter %in% c("alpha1n", "alpha1p", "alpha1c", "alpha2n", "alpha2p", "alpha2c"), ] = subdata

tabS5 = datanew

save(tabS5, file = "../results/TabS5.RData")