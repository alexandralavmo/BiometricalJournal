# File to define parameters used in scenario of strong correlation

# all parameters are the same as previously
source("parameters.R")

# only the correlation on slope parameters change
rho_int = 0.95
rho_slp = 0.95
mcov["omega_b02", "omega_b04"] = mcov["omega_b04", "omega_b02"] = rho_int*sqrt(params_lin[1, "omega_b0"])*sqrt(params_lin[3, "omega_b0"])
mcov["omega_b03", "omega_b05"] = mcov["omega_b05", "omega_b03"] = rho_int*sqrt(params_lin[2, "omega_b0"])*sqrt(params_lin[4, "omega_b0"])
mcov["omega_b12", "omega_b14"] = mcov["omega_b14", "omega_b12"] = rho_slp*sqrt(params_lin[1, "omega_b1"])*sqrt(params_lin[3, "omega_b1"])
mcov["omega_b13", "omega_b15"] = mcov["omega_b15", "omega_b13"] = rho_slp*sqrt(params_lin[2, "omega_b1"])*sqrt(params_lin[4, "omega_b1"])
