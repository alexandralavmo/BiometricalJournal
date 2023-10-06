# Required libraries	
library(mvtnorm)
library(ggplot2)
library(cowplot)
library(ComplexUpset)
library(survival)
library(stringr)
library(timeROC)
library(R6)
library(RJSONIO)
library(MlxConnectors)
# !!! Change directory to specify correct Monolix localization !!! #
initializeMlxConnectors(software = "monolix", mlxDirectory = "C:/ProgramData/Lixoft/MonolixSuite2018R2")

# General simulation parameters
N = 300 # Number of patients
p = 100 # Number of simulation replicates

# First biomarker (nonlinear mixed model, see manuscript for parameter interpretation)
params_nonlin = c(b01 = 4.6, b11 = -0.15, b21 = -0.16, a1 = log(5.3), omega_b01 = 2**2, omega_b11 = 0.1**2, omega_b21 = 0.07**2, omega_a1 = 0.8**2, sigma_b1 = 0.3)
params_nonlin = data.frame(as.list(params_nonlin))

# Other biomarkers (linear mixed model, 6 biomarkers)
b0 = c(7.4, 4.2, 5.2, 7, 5.2, 338)			# Population intercept
b1 = c(0.003, -0.16, 0.02, -0.05, -0.07, -5)		# Population slope
omega_b0 = c(0.04, 0.9, 1.7, 0.8, 0.7, 107)**2	# individual variabiliy on intercept (variances)
omega_b1 = c(0.005, 0.15, 0.08, 0.1, 0.09, 5)**2	# individual variabiliy on slope (variances)
sigma_a = c(0.05, 0.7, 0.7, 0, 0, 0)		# residual error (additive part)
sigma_b = c(0, 0, 0, 0.08, 0.1, 0.2)		# residual error (proportional part)

params_lin = data.frame(b0, b1, omega_b0, omega_b1, sigma_a, sigma_b)
params = list(nonlinear = params_nonlin, linear = params_lin)
list_nl = unlist(params_nonlin[1:4])
list_l = unlist(params_lin[, c("b0", "b1")])
names(list_l) = c("b02", "b03", "b04", "b05", "b06", "b07", "b12", "b13", "b14", "b15", "b16", "b17")
pop_means = c(list_nl, list_l)

# Covariance matrix of random effects for all biomarkers
list_nl_omega = unlist(params_nonlin[5:8])
list_l_omega = unlist(params_lin[, c("omega_b0", "omega_b1")])
names(list_l_omega) = c("omega_b02", "omega_b03", "omega_b04", "omega_b05", "omega_b06", "omega_b07", "omega_b12", "omega_b13", "omega_b14", "omega_b15", "omega_b16", "omega_b17")
v = c(list_nl_omega, list_l_omega)
mcov = diag(v)
rownames(mcov) = colnames(mcov) = names(v)
rho_int = 0
rho_slp = 0.8
mcov["omega_b02", "omega_b04"] = mcov["omega_b04", "omega_b02"] = rho_int*sqrt(params_lin[1, "omega_b0"])*sqrt(params_lin[3, "omega_b0"])
mcov["omega_b03", "omega_b05"] = mcov["omega_b05", "omega_b03"] = rho_int*sqrt(params_lin[2, "omega_b0"])*sqrt(params_lin[4, "omega_b0"])
mcov["omega_b12", "omega_b14"] = mcov["omega_b14", "omega_b12"] = rho_slp*sqrt(params_lin[1, "omega_b1"])*sqrt(params_lin[3, "omega_b1"])
mcov["omega_b13", "omega_b15"] = mcov["omega_b15", "omega_b13"] = rho_slp*sqrt(params_lin[2, "omega_b1"])*sqrt(params_lin[4, "omega_b1"])

# Sampling times for longitudinal markers
tt1 = tt3 = seq(0, 30, 2)	# One observation every two days until day 30 for biomarker 1 and 3
tt2 = seq(0, 30, 1.5)
tt4 = tt5 = tt6 = tt7 = seq(0, 30, 3)

# Survival parameters 
p1 = 0.05	# proportion of event 1 at infinite time for covariates equal to 0
g1 = 0.1	# rate of event 1 occurence
alpha = c(0.14, -11, 0.6, 0, 0, 0, 0) # Link coefficient for all biomarkers

# To test scalability 
n_scalable = c(100, 200, 300, 400, 800, 1600) # number of patients included in the joint model (3 biomarkers fixed)
k_scalable = 1:7                         # number of biomarkers involved in the joint model (300 patients fixed)
