# File to simulate data set under first scenario of correlation

###############
# 1. Generation of all datasets
###############
set.seed(1996)
sapply(1:p, sim)

###############
# 4. Simulation procedure to test scalability depending the number of patients (3 biomarkers fixed)
###############
set.seed(1996)
# Start the simulation with the model including 1600 patients 
# (multivariate gaussian distribution because of a possible correlations between random effects)
# each row represents a patient and each column an individual parameter
indv_params = as.data.frame(rmvnorm(1600, pop_means, mcov))

# full longitudinal trajectories
longi1 = sim_nonlin(1600, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), params_nonlin$sigma_b1, tt1) # exp(a0) because this parameter follow a lognormal distribution
longi2 = sim_lin(1600, indv_params$b02, indv_params$b12, params_lin$sigma_a[1], params_lin$sigma_b[1], tt2)
longi3 = sim_lin(1600, indv_params$b03, indv_params$b13, params_lin$sigma_a[2], params_lin$sigma_b[2], tt3)
	
# median of biomarker observations 
med1 = median(longi1$obs)
med2 = median(longi2$obs)
med3 = median(longi3$obs)

# Generation of all event types and times
survobs = mapply(gettimes, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), indv_params$b02, indv_params$b12, indv_params$b03, indv_params$b13)
tte.data = data.frame(id = 1:1600, time = survobs[1, ], obs = survobs[2, ])	# dataset of survival data
	
# Administrative censoring at D30
tte.data$obs[tte.data$time>30 | is.na(tte.data$time) == T] = 0	
tte.data$time[tte.data$time>30 | is.na(tte.data$time) == T] = 30

# Longitudinal observations after the event time are removed
longi1 = longi1[longi1$time<= sapply(longi1$id, function(x) tte.data$time[tte.data$id == x]), ]
longi2 = longi2[longi2$time<= sapply(longi2$id, function(x) tte.data$time[tte.data$id == x]), ]
longi3 = longi3[longi3$time<= sapply(longi3$id, function(x) tte.data$time[tte.data$id == x]), ] 

# Adding ytype column (before merging datasets)
longi1$ytype = 1
longi2$ytype = 2
longi3$ytype = 3

# creating time to event data in Monolix format
tte.data = do.call(rbind, mapply(makesurv, tte.data$id, tte.data$time, tte.data$obs, SIMPLIFY = F))
	
# merging all datasets and writing to output file
alldata = rbind(longi1, longi2, longi3, tte.data)
write.table(alldata, file = "data_scalability/data_n1600.txt", row.names = F)


###############
# 3. Generation of all datasets to test scalability 
###############

sapply(n_scalable, sim_n)
sapply(k_scalable, sim_k)