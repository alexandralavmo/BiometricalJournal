##############
# 1. Data simulation (first scenario of correlation)
##############

### a. Simulation fonctions for linear and nonlinear model

# sim_lin : returns a data frame with columns id, sampling time and observation.
# argument N represent the number of patients 
# Observations are derived from a linear mixed model with parameters a and b (vector of individual parameters for intercept and slope) 
# and with a residual error composed of an additive part (of variance sigma_a^2) and a proportional part (of variance sigma_b^2)
# argument times represent observations times for all patients
sim_lin = function(N, a, b, sigma_a, sigma_b, times){
	ids = as.vector(sapply(1:N, function(x) rep(x, length(times))))
	longi = data.frame(id = ids, time = times)
	longi$obs = a[longi$id] + b[longi$id]*longi$time		
	longi$obs = longi$obs + (sigma_a + longi$obs*sigma_b)*rnorm(nrow(longi), sd = 1) 	
	longi
}

# sim_nonlin : returns a data frame with columns id, sampling time and observation.
# argument N represent the number of patients 
# Observations are derived from a nonlinear mixed model 
# b0, b1, b2 and a are vectors of individual parameters (see paper for parameter explanation)
# and a residual proportional error of variance sigma^2 
# argument times represent observations times for all patients
sim_nonlin = function(N, b0, b1, b2, a, sigma, times){
	ids = as.vector(sapply(1:N, function(x) rep(x, length(times))))
	longi = data.frame(id = ids, time = times)
	longi$obs = b0[longi$id] + a[longi$id] *(exp(b1[longi$id]*longi$time)-exp(b2[longi$id]*longi$time)) 
	longi$obs = longi$obs*(1+rnorm(nrow(longi), sd = sigma))	# proportional error
	longi
}


### b. Simulation procedure

# function to save the data set corresponding to a simulation it 
sim = function(it) {

	# Simulation of individual parameters 
	# (multivariate gaussian distribution because of a possible correlations between random effects)
	# each row represents a patient and each column an individual parameter
	indv_params = as.data.frame(rmvnorm(N, pop_means, mcov))

	# full longitudinal trajectories
	longi1 = sim_nonlin(N, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), params_nonlin$sigma_b1, tt1) # exp(a0) because this parameter follow a lognormal distribution
	longi2 = sim_lin(N, indv_params$b02, indv_params$b12, params_lin$sigma_a[1], params_lin$sigma_b[1], tt2)
	longi3 = sim_lin(N, indv_params$b03, indv_params$b13, params_lin$sigma_a[2], params_lin$sigma_b[2], tt3)
	longi4 = sim_lin(N, indv_params$b04, indv_params$b14, params_lin$sigma_a[3], params_lin$sigma_b[3], tt4)
	longi5 = sim_lin(N, indv_params$b05, indv_params$b15, params_lin$sigma_a[4], params_lin$sigma_b[4], tt5)
	longi6 = sim_lin(N, indv_params$b06, indv_params$b16, params_lin$sigma_a[5], params_lin$sigma_b[5], tt6)
	longi7 = sim_lin(N, indv_params$b07, indv_params$b17, params_lin$sigma_a[6], params_lin$sigma_b[6], tt7)

	# median of biomarker observations 
  	med1 = median(longi1$obs)
  	med2 = median(longi2$obs)
  	med3 = median(longi3$obs)
	med4 = median(longi4$obs)
	med5 = median(longi5$obs)
	med6 = median(longi6$obs)
	med7 = median(longi7$obs)

	# gettimes: function for simulating event type (1 or 2) and event time from a specified joint model
	# vector of individual parameters for biomarkers 1 to 3 are passed as argument
	# (the others biomarkers are independant of the survival process)
	# The function proceeds by constructing both subdistribution functions
	gettimes = function(r1, r2, r3, r4, r5, r6, r7, r8){
		tt = seq(0, 1000, length.out = 10000) # sequence of times for numerical integration
		pas = tt[2]-tt[1]
		# h denotes the subsitribution hazard for event 1 at each time of the sequence tt
		h = (p1*g1*exp(-g1*tt)/(1-p1*(1-exp(-g1*tt))))*exp(alpha[1]*(r1+r4*(exp(r2*tt)-exp(r3*tt))-med1)+alpha[2]*(r5+r6*tt-med2)+alpha[3]*(r7+r8*tt-med3))
		F1 = 1-exp(-cumsum(h)*pas) # numerical integration for the derivation of the subdistribution of event 1
		F2 = (1-rev(F1)[1])*(1-exp(-tt/10)) # subdistribution function of event 2
		prob1 = rev(F1)[1] # marginal probability of event 1
		prob2 = 1-p1
		event = ifelse(runif(1)<prob1, 1, 2) # event type simulation based on the marginal probability
		if(event == 1) z = which(F1/prob1>runif(1))[1]*pas  # event time after conditionning on event type
		if(event == 2) z = which(F2/prob2>runif(1))[1]*pas 
		c(z, event)
	}

	# Generation of all event types and times
	survobs = mapply(gettimes, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), indv_params$b02, indv_params$b12, indv_params$b03, indv_params$b13)
	tte.data = data.frame(id = 1:N, time = survobs[1, ], obs = survobs[2, ])	# dataset of survival data
	
	# Administrative censoring at D30
	tte.data$obs[tte.data$time>30 | is.na(tte.data$time) == T] = 0	
	tte.data$time[tte.data$time>30 | is.na(tte.data$time) == T] = 30

	# Longitudinal observations after the event time are removed
	longi1 = longi1[longi1$time<= sapply(longi1$id, function(x) tte.data$time[tte.data$id == x]), ]
	longi2 = longi2[longi2$time<= sapply(longi2$id, function(x) tte.data$time[tte.data$id == x]), ]
	longi3 = longi3[longi3$time<= sapply(longi3$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi4 = longi4[longi4$time<= sapply(longi4$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi5 = longi5[longi5$time<= sapply(longi5$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi6 = longi6[longi6$time<= sapply(longi6$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi7 = longi7[longi7$time<= sapply(longi7$id, function(x) tte.data$time[tte.data$id == x]), ]

	# Adding ytype column (before merging datasets)
	longi1$ytype = 1
	longi2$ytype = 2
	longi3$ytype = 3
	longi4$ytype = 4
	longi5$ytype = 5
	longi6$ytype = 6
	longi7$ytype = 7

	# Monolix requirement for fitting competing events: the survival dataset must include
	# 2 lines for each patient and each event (i.e. 4 lines per patient)
	# a first line with time = 0 and obs = 0 for both events
	# a second line with t = 30 and obs = 0 (for censored observations or competing event occurred)
	# or with event type and obs = 1 if an event has been observed
	# event types are denoted 8 and 9 for the first and second event respectively
	makesurv = function(id, ev_time, obs){
		if(obs == 0) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, 30), obs = 0, ytype = c(8, 8, 9, 9))
		if(obs == 1) surv = data.frame(id = rep(id, 4), time = c(0, ev_time, 0, 30), obs = c(0, 1, 0, 0), ytype = c(8, 8, 9, 9))
		if(obs == 2) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, ev_time), obs = c(0, 0, 0, 1), ytype = c(8, 8, 9, 9))
		surv
	}
	tte.data = do.call(rbind, mapply(makesurv, tte.data$id, tte.data$time, tte.data$obs, SIMPLIFY = F))
	
	# merging all datasets and writing to output file
	alldata = rbind(longi1, longi2, longi3, longi4, longi5, longi6, longi7, tte.data)
	write.table(alldata, file = paste0("data/data", it, ".txt"), row.names = F)

	# keeping the median of each biomarker for estimation
	med = read.table("modfiles/med_biom_sim.txt", header = T)
	med[it, 1] = it; med[it, 2] = med1; med[it, 3] = med2; med[it, 4] = med3; med[it, 5] = med4; med[it, 6] = med5; med[it, 7] = med6; med[it, 8] = med7
	write.table(med, file = paste0("modfiles/med_biom_sim.txt"), row.names = F)
}


### c. Simulation functions to test scalability 

# function for simulating event type (1 or 2) and event time from a specified joint model
# vector of individual parameters for biomarkers 1 to 3 are passed as argument
# (the others biomarkers are independant of the survival process)
# The function proceeds by constructing both subdistribution functions
gettimes = function(r1, r2, r3, r4, r5, r6, r7, r8){
	tt = seq(0, 1000, length.out = 10000) # sequence of times for numerical integration
	pas = tt[2]-tt[1]
	# h denotes the subsitribution hazard for event 1 at each time of the sequence tt
	h = (p1*g1*exp(-g1*tt)/(1-p1*(1-exp(-g1*tt))))*exp(alpha[1]*(r1+r4*(exp(r2*tt)-exp(r3*tt))-med1)+alpha[2]*(r5+r6*tt-med2)+alpha[3]*(r7+r8*tt-med3))
	F1 = 1-exp(-cumsum(h)*pas) # numerical integration for the derivation of the subdistribution of event 1
	F2 = (1-rev(F1)[1])*(1-exp(-tt/10)) # subdistribution function of event 2
	prob1 = rev(F1)[1] # marginal probability of event 1
	prob2 = 1-p1
	event = ifelse(runif(1)<prob1, 1, 2) # event type simulation based on the marginal probability
	if(event == 1) z = which(F1/prob1>runif(1))[1]*pas  # event time after conditionning on event type
	if(event == 2) z = which(F2/prob2>runif(1))[1]*pas 
	c(z, event)
}

# function giving the survival observations in Monolix format for individual id, knowing his event time ev_time and its event indicator obs 
# Monolix requirement for fitting competing events: the survival dataset must include
# 2 lines for each patient and each event (i.e. 4 lines per patient)
# a first line with time = 0 and obs = 0 for both events
# a second line with t = 30 and obs = 0 (for censored observations or competing event occurred)
# or with event type and obs = 1 if an event has been observed
# event types are denoted 8 and 9 for the first and second event respectively
makesurv = function(id, ev_time, obs){
	if(obs == 0) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, 30), obs = 0, ytype = c(8, 8, 9, 9))
	if(obs == 1) surv = data.frame(id = rep(id, 4), time = c(0, ev_time, 0, 30), obs = c(0, 1, 0, 0), ytype = c(8, 8, 9, 9))
	if(obs == 2) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, ev_time), obs = c(0, 0, 0, 1), ytype = c(8, 8, 9, 9))
	surv
}

# function to keep the n first patients in initial data set with 1600 patients
# save data in "data_scalability" folder
sim_n = function(n) {
	data = read.table(file = "data_scalability/data_n1600.txt", header = T)
	seq = 1:n 
	subdata = data[data$id %in% seq, ]
	write.table(subdata, file = paste0("data_scalability/data_n", n, ".txt"), row.names = F)
}

# function to keep the k first biomarkers from the first simulated data set
# save results in "data_scalability" folder
sim_k = function(k) {
	data = read.table("data/data1.txt", header = T)
	seq = c(1:k, 8, 9)
	subdata = data[data$ytype %in% seq, ]
	write.table(subdata, file = paste0("data_scalability/data_bm", k, ".txt"), row.names = F)
}


##############
# 2. True model fit 
##############

# function to import data referring to simulation i, create corresponding model text and mlxtran files
create_model_true = function(i){
	# First filter data file in order to keep only three biomarkers
	dat = read.table(paste0("data/data", i, ".txt"), header = T, sep = ' ')
	dat = dat[dat$ytype %in% c(1, 2, 3, 8, 9), ]
	write.table(dat, file = paste0("intermediate_results/fit_truemodel/data", i, ".txt"), row.names = F)

	# Then, modify the txt model file (with the correct median values)
	med = read.table("modfiles/med_biom_sim.txt", header = T)
	obj = readLines("modfiles/mod_123.txt")
	obj[12] = paste0("haz1 = (p1*g1*exp(-g1*t)/(1-p1*(1-exp(-g1*t)))) * exp( alpha11*(m1-", med$med1[i], ") +alpha12*(m2-", med$med2[i], ") +alpha13*(m3-", med$med3[i], ") )")
	obj[13] = paste0("haz2 = (p2*g2*exp(-g2*t)/(1-p2*(1-exp(-g2*t)))) * exp( alpha21*(m1-", med$med1[i], ") +alpha22*(m2-", med$med2[i], ") +alpha23*(m3-", med$med3[i], ") )")
	writeLines(obj, paste0("intermediate_results/fit_truemodel/mod_123_", i, ".txt"))
	
	# Finally, modify the mlxtran parameter file 
	obj = readLines("modfiles/ref123.mlxtran")
	obj[4] = paste0("file = 'data", i, ".txt'")
	obj[42] = paste0("file = 'mod_123_", i, ".txt'")
	obj[95] = paste0("exportpath = 'mod", i, "'")
	write.table(obj, file = paste0("intermediate_results/fit_truemodel/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
}

# functions to fit the true model on simulation i 
# and save results in "intermediate_results/fit_truemodel" folder
estimate_true = function(i){
	loadProject(paste0("intermediate_results/fit_truemodel/mod", i, ".mlxtran"))
	runPopulationParameterEstimation()
	fit = data.frame(par = names(getEstimatedPopulationParameters()), val = getEstimatedPopulationParameters())
	write.table(fit, paste0("intermediate_results/fit_truemodel/pop", i, ".txt"), row.names = F)
}


##############
# 3. Backward strategy on first scenario of correlations
##############

# function to create model text files with a specified list of biomarkers passed as argument
# create_modtxt(i) creates a txt model including biomarkers i (i can be a vector)
create_modtxt = function(i){
  
  num = setdiff(c("1", "2", "3", "4", "5", "6", "7"), i)
  ref = readLines("modfiles/mod_1234567.txt")
  warn = 0
  l4 = unlist(strsplit(ref[4], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(num)){
    indb = c(grep(paste0("b0", num[j]), l4), grep(paste0("b1", num[j]), l4), grep(paste0("b2", num[j]), l4), grep(paste0("a", num[j], "_pop"), l4), grep(paste0("alpha1", num[j]), l4), grep(paste0("alpha2", num[j]), l4), grep(paste0("a", num[j], ", "), l4))
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l4)-1) warn = 1
  l4 = l4[-c(ind)]
  
  if (warn == 1) {der = paste0(str_sub(rev(l4)[2], 1, 7)) ; l4[c(length(l4)-1)] = der}
  ref[4] = paste(l4, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 8:14){
    for (j in num){
      l = unlist(strsplit(ref[ligne], split = " "))
      if (l[1] == paste0("m", j)) ind_ligne = c(ind_ligne, ligne)  
    }
  }
  
  ind = c()
  indb = c()
  l16 = unlist(strsplit(ref[16], split = " "))
  for (j in 1:length(num)){
    indb = grep(paste0("alpha1", num[j]), l16)
    ind = c(ind, indb)
  }
  l16 = l16[-c(ind)]
  ref[16] = paste(l16, collapse = " ")
  
  l17 = unlist(strsplit(ref[17], split = " "))
  l17 = l17[-c(ind)]
  ref[17] = paste(l17, collapse = " ")
  
  ind = c()
  indb = c()
  l25 = unlist(strsplit(ref[25], split = " "))
  for (j in 1:length(num)){
    indb = grep(num[j], l25)
    ind = c(ind, indb)
  }
  l25 = l25[-c(ind)]
  ref[25] = paste(l25, collapse = " ")
  
  ref = ref[-c(ind_ligne)]
  return(ref)
}


# function to create all possible model text files with k biomarkers
# and write them into the appropriate directory
create_combtxt_scen1 = function(k){
  create_combination = function(col) {
    txt = create_modtxt(as.character(col))
    path = paste0("intermediate_results/backward/back", k, "/mod_", paste(col, collapse = ''), ".txt")
    writeLines(txt, path)
  }
  apply(combn(1:7, k), 2, create_combination)
}

# Function to adapt the median values in the model files for each simulation 
# the arguments are num: the vector of biomamarkers involved and j: the simulation j
adapt_txt_scen1 = function(num, j){
  med = read.table("modfiles/med_biom_sim.txt", header = T)
  if (length(num) == 7) {ref = readLines("modfiles/mod_1234567.txt")}
  else {ref = readLines(paste0("intermediate_results/backward/back", length(num), "/mod_", paste0(num, collapse = ""), ".txt"))}
  l = unlist(strsplit(ref[length(num)+9], split = " "))
  l[which(l == "+alpha11*(m1-4.58)")] = paste0("+alpha11*(m1-", med$med1[med$sim == j], ")")
  l[which(l == "+alpha12*(m2-7.44)")] = paste0("+alpha12*(m2-", med$med2[med$sim == j], ")")
  l[which(l == "+alpha13*(m3-2.36)")] = paste0("+alpha13*(m3-", med$med3[med$sim == j], ")")
  l[which(l == "+alpha14*(m4-5.28)")] = paste0("+alpha14*(m4-", med$med4[med$sim == j], ")")
  l[which(l == "+alpha15*(m5-6.39)")] = paste0("+alpha15*(m5-", med$med5[med$sim == j], ")")
  l[which(l == "+alpha16*(m6-4.41)")] = paste0("+alpha16*(m6-", med$med6[med$sim == j], ")")
  l[which(l == "+alpha17*(m7-256.6)")] = paste0("+alpha17*(m7-", med$med7[med$sim == j], ")")
  l = paste(l, collapse = " ")
  m = unlist(strsplit(ref[length(num)+10], split = " "))
  m[which(m == "+alpha21*(m1-4.58)")] = paste0("+alpha21*(m1-", med$med1[med$sim == j], ")")
  m[which(m == "+alpha22*(m2-7.44)")] = paste0("+alpha22*(m2-", med$med2[med$sim == j], ")")
  m[which(m == "+alpha23*(m3-2.36)")] = paste0("+alpha23*(m3-", med$med3[med$sim == j], ")")
  m[which(m == "+alpha24*(m4-5.28)")] = paste0("+alpha24*(m4-", med$med4[med$sim == j], ")")
  m[which(m == "+alpha25*(m5-6.39)")] = paste0("+alpha25*(m5-", med$med5[med$sim == j], ")")
  m[which(m == "+alpha26*(m6-4.41)")] = paste0("+alpha26*(m6-", med$med6[med$sim == j], ")")
  m[which(m == "+alpha27*(m7-256.6)")] = paste0("+alpha27*(m7-", med$med7[med$sim == j], ")")
  m = paste(m, collapse = " ")
  ref[length(num)+9] = l
  ref[length(num)+10] = m
  path = paste0("intermediate_results/backward/back", length(num), "/mod_", paste(num, collapse = ''), "_", j, ".txt")
  writeLines(ref, path)
}

# function to adapt the mlxtran file of simulation i involving 7 biomarker (full initial joint model)
setmodel_scen1 = function(i) {
  obj = readLines("modfiles/ref1234567.mlxtran")
  obj[4] = paste0("file = 'data", i, ".txt'")
  obj[58] = paste0("file = 'mod_1234567_", i, ".txt'")
  obj[143] = paste0("exportpath = 'mod", i, "'")
  write.table(obj, file = paste0("intermediate_results/backward/back7/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
}

# function to fit the full initial joint model on simulation i 
fit_fullmodels_scen1 = function(i){
  loadProject(paste0("intermediate_results/backward/back7/mod", i, ".mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}

# function to create mlxtran files given a list of biomarkers and for a given dataset
# create_mlx_scen1(i, j) creates a mlxtran file removing biomarkers i and setting dataset j as associated data
create_mlx_scen1 = function(i, j){
  mod_ind = 1:7
  mod_ind = paste0(mod_ind[-c(as.integer(i))], collapse = "")
  mod_ind2 = seq(1, 7)[-c(as.integer(i))]
  adapt_txt_scen1(mod_ind2, j)
  
  ref7 = readLines("modfiles/ref1234567.mlxtran")
  ref7[4] = paste0("file = 'data", j, ".txt'") # Set correct data file
  ref7[58] = paste0("file = 'mod_", mod_ind, "_", j, ".txt'")# Set correct model file
  ref7[143] = paste0("exportpath = 'mod", j, "'")
  
  l11 = unlist(strsplit(ref7[11], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    ind = grep(i[j], l11)
    indb = c(indb, ind)
  }
  ind = c(indb, which(l11 == "continuous,")[1:length(i)])
  l11 = l11[-c(ind)]
  ref7[11] = paste(l11, collapse = " ") 
  
  l17 = unlist(strsplit(ref7[17], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = c(grep(paste0("b0", i[j]), l17), grep(paste0("b1", i[j]), l17), grep(paste0("b2", i[j]), l17), grep(paste0("a", i[j], "_pop"), l17), grep(paste0("omega_a", i[j]), l17), grep(paste0("alpha1", i[j]), l17), grep(paste0("alpha2", i[j]), l17))
    ind = c(ind, indb)
  }
  warn = 0
  if (rev(ind)[1] == length(l17)) warn = 1
  l17 = l17[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l17)[1], 1, 11), "}") ; l17[length(l17)] = der}
  ref7[17] = paste(l17, collapse = " ")
  
  l56 = unlist(strsplit(ref7[56], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l56)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l56)) warn = 1
  l56 = l56[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l56)[1], 1, 3), "}") ; l56[length(l56)] = der}
  ref7[56] = paste(l56, collapse = " ") 
  
  l70 = unlist(strsplit(ref7[70], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l70)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l70)) warn = 1
  l70 = l70[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l70)[1], 1, 2), "}") ; l70[length(l70)] = der}
  ref7[70] = paste(l70, collapse = " ") 
  
  l71 = unlist(strsplit(ref7[71], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l71)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l71)) warn = 1
  l71 = l71[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l71)[1], 1, 2), "}") ; l71[length(l71)] = der}
  ref7[71] = paste(l71, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 1:length(ref7)){
    l = unlist(strsplit(ref7[ligne], split = " "))
    for (j in 1:length(i)){
      if (length(grep(paste0("b0", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b1", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b2", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("alpha1", i[j], "_pop"), l)) != 0 | length(grep(paste0("alpha2", i[j], "_pop"), l)) != 0 | length(grep(paste0("b_", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("a_", i[j]), l, fixed = T)) != 0 | length(grep(paste0("omega_a", i[j]), l, fixed = T)) != 0 | length(grep(paste0("a", i[j], "_pop"), l, fixed = T)) != 0) ind_ligne = c(ind_ligne, ligne)
    }
  }
  ref7 = ref7[-c(ind_ligne)]
  return(ref7)
}

# function to create mlxtran files and data files for iteration k of the backward 
# load the results found at the previous iteration and create the new model and data files
iter_back_scen1 = function(k){
  for (i in 1:p){ 
    if( !file.exists(paste0("intermediate_results/backward/back", k+1, "/mod", i, "/populationParameters.txt"))) next
    tab = read.table(paste0("intermediate_results/backward/back", k+1, "/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    dat = read.table(paste0("intermediate_results/backward/back", k+1, "/data", i, ".txt"), header = T, sep = ' ')
    stab = tab[c(which(tab$parameter %in% c("alpha11_pop", "alpha12_pop", "alpha13_pop", "alpha14_pop", "alpha15_pop", "alpha16_pop", "alpha17_pop"))), ]
    max = stab$parameter[stab$rse_sa == max(stab$rse_sa[is.nan(stab$rse_sa) == F])&is.nan(stab$rse_sa) == F]
    prev = substr(stab$parameter, 7, 7)
    num = substr(max, 7, 7)
    cur = prev[-which(prev == num)]
    numb = setdiff(c("1", "2", "3", "4", "5", "6", "7"), cur)
    if (tab$rse_sa[tab$parameter == max]>50){
      mlx = create_mlx_scen1(numb, i)
      write.table(mlx, file = paste0("intermediate_results/backward/back", k, "/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
      dat = dat[-c(which(dat$ytype == as.integer(num))), ]
      write.table(dat, file = paste0("intermediate_results/backward/back", k, "/data", i, ".txt"), row.names = F)
      if (k == 0){
        print(paste0("No biomarker kept for simulation ", i))
        tab = data.frame(parameter = NULL, value = NULL, se_sa = NULL, rse_sa = NULL)
        write.table(tab, paste0("intermediate_results/backward/final/mod", i, ".txt"), row.names = F) 
      }
    }
    else{
      indnan = which(is.nan(stab$rse_sa) == T)
      if(length(indnan)>0){
        min = stab$parameter[stab$value == min(stab$value[indnan])]
        prev = substr(stab$parameter, 7, 7)
        num = substr(min, 7, 7)
        cur = prev[-which(prev == num)]
        numb = setdiff(c("1", "2", "3", "4", "5", "6", "7"), cur)
        mlx = create_mlx_scen1(numb, i)
        write.table(mlx, file = paste0("intermediate_results/backward/back", k, "/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
        dat = dat[-c(which(dat$ytype == as.integer(num))), ]
        write.table(dat, file = paste0("intermediate_results/backward/back", k, "/data", i, ".txt"), row.names = F)
      } 
      else{
        print(paste0("End for simulation ", i))
        write.table(tab, paste0("intermediate_results/backward/final/mod", i, ".txt"), row.names = F)
        write.table(dat, paste0("intermediate_results/backward/final/data", i, ".txt"), row.names = F)
	file.copy(paste0("intermediate_results/backward/back", k+1, "/mod", i, ".mlxtran"), paste0("intermediate_results/backward/final/mod", i, ".mlxtran"))
        file.copy(paste0("intermediate_results/backward/back", k+1, "/mod_", paste(prev, collapse = ''), "_", i, ".txt"), paste0("intermediate_results/backward/final/mod_", paste(prev, collapse = ''), "_", i, ".txt"))
        file.copy(paste0("intermediate_results/backward/back", k+1, "/mod", i), paste0("intermediate_results/backward/final"), recursive = T) 
      }
    }
  }
}

# function to fit the model involving k biomarkers for a simulation i 
fit_all_models_scen1 = function(k){
  fitmodel = function(i){
    if( !file.exists(paste0("intermediate_results/backward/back", k, "/mod", i, ".mlxtran"))) return(NA)
    loadProject(paste0("intermediate_results/backward/back", k, "/mod", i, ".mlxtran"))
    runPopulationParameterEstimation()
    runStandardErrorEstimation()
  }
  sapply(1:p, fitmodel)
}

# creating indicator vector of length 7 (referring to each biomarker) for simulation i
# 1 if the biomarker is selected, 0 otherwise 
getdata = function(i){
  tab = read.table(paste0("intermediate_results/backward/final/mod", i, ".txt"), header = T, sep = '')
  1*!is.na(match(vars, tab$parameter))
}

##############
# 4. Data simulation (strong scenario of correlation)
##############

# function to save the data set corresponding to a simulation it 
sim2 = function(it) {

	# Simulation of individual parameters 
	# (multivariate gaussian distribution because of a possible correlations between random effects)
	# each row represents a patient and each column an individual parameter
	indv_params = as.data.frame(rmvnorm(N, pop_means, mcov))

	# full longitudinal trajectories
	longi1 = sim_nonlin(N, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), params_nonlin$sigma_b1, tt1) # exp(a0) because this parameter follow a lognormal distribution
	longi2 = sim_lin(N, indv_params$b02, indv_params$b12, params_lin$sigma_a[1], params_lin$sigma_b[1], tt2)
	longi3 = sim_lin(N, indv_params$b03, indv_params$b13, params_lin$sigma_a[2], params_lin$sigma_b[2], tt3)
	longi4 = sim_lin(N, indv_params$b04, indv_params$b14, params_lin$sigma_a[3], params_lin$sigma_b[3], tt4)
	longi5 = sim_lin(N, indv_params$b05, indv_params$b15, params_lin$sigma_a[4], params_lin$sigma_b[4], tt5)
	longi6 = sim_lin(N, indv_params$b06, indv_params$b16, params_lin$sigma_a[5], params_lin$sigma_b[5], tt6)
	longi7 = sim_lin(N, indv_params$b07, indv_params$b17, params_lin$sigma_a[6], params_lin$sigma_b[6], tt7)

	# median of biomarker observations 
  med1 = median(longi1$obs)
  med2 = median(longi2$obs)
  med3 = median(longi3$obs)
	med4 = median(longi4$obs)
	med5 = median(longi5$obs)
	med6 = median(longi6$obs)
	med7 = median(longi7$obs)

	# gettimes: function for simulating event type (1 or 2) and event time from a specified joint model
	# vector of individual parameters for biomarkers 1 to 3 are passed as argument
	# (the others biomarkers are independant of the survival process)
	# The function proceeds by constructing both subdistribution functions
	gettimes = function(r1, r2, r3, r4, r5, r6, r7, r8){
		tt = seq(0, 1000, length.out = 10000) # sequence of times for numerical integration
		pas = tt[2]-tt[1]
		# h denotes the subsitribution hazard for event 1 at each time of the sequence tt
		h = (p1*g1*exp(-g1*tt)/(1-p1*(1-exp(-g1*tt))))*exp(alpha[1]*(r1+r4*(exp(r2*tt)-exp(r3*tt))-med1)+alpha[2]*(r5+r6*tt-med2)+alpha[3]*(r7+r8*tt-med3))
		F1 = 1-exp(-cumsum(h)*pas) # numerical integration for the derivation of the subdistribution of event 1
		F2 = (1-rev(F1)[1])*(1-exp(-tt/10)) # subdistribution function of event 2
		prob1 = rev(F1)[1] # marginal probability of event 1
		prob2 = 1-p1
		event = ifelse(runif(1)<prob1, 1, 2) # event type simulation based on the marginal probability
		if(event == 1) z = which(F1/prob1>runif(1))[1]*pas  # event time after conditionning on event type
		if(event == 2) z = which(F2/prob2>runif(1))[1]*pas 
		c(z, event)
	}

	# Generation of all event types and times
	survobs = mapply(gettimes, indv_params$b01, indv_params$b11, indv_params$b21, exp(indv_params$a1), indv_params$b02, indv_params$b12, indv_params$b03, indv_params$b13)
	tte.data = data.frame(id = 1:N, time = survobs[1, ], obs = survobs[2, ])	# dataset of survival data
	
	# Administrative censoring at D30
	tte.data$obs[tte.data$time>30 | is.na(tte.data$time) == T] = 0	
	tte.data$time[tte.data$time>30 | is.na(tte.data$time) == T] = 30

	# Longitudinal observations after the event time are removed
	longi1 = longi1[longi1$time<= sapply(longi1$id, function(x) tte.data$time[tte.data$id == x]), ]
	longi2 = longi2[longi2$time<= sapply(longi2$id, function(x) tte.data$time[tte.data$id == x]), ]
	longi3 = longi3[longi3$time<= sapply(longi3$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi4 = longi4[longi4$time<= sapply(longi4$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi5 = longi5[longi5$time<= sapply(longi5$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi6 = longi6[longi6$time<= sapply(longi6$id, function(x) tte.data$time[tte.data$id == x]), ] 
	longi7 = longi7[longi7$time<= sapply(longi7$id, function(x) tte.data$time[tte.data$id == x]), ]

	# Adding ytype column (before merging datasets)
	longi1$ytype = 1
	longi2$ytype = 2
	longi3$ytype = 3
	longi4$ytype = 4
	longi5$ytype = 5
	longi6$ytype = 6
	longi7$ytype = 7

	# Monolix requirement for fitting competing events: the survival dataset must include
	# 2 lines for each patient and each event (i.e. 4 lines per patient)
	# a first line with time = 0 and obs = 0 for both events
	# a second line with t = 30 and obs = 0 (for censored observations or competing event occurred)
	# or with event type and obs = 1 if an event has been observed
	# event types are denoted 8 and 9 for the first and second event respectively
	makesurv = function(id, ev_time, obs){
		if(obs == 0) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, 30), obs = 0, ytype = c(8, 8, 9, 9))
		if(obs == 1) surv = data.frame(id = rep(id, 4), time = c(0, ev_time, 0, 30), obs = c(0, 1, 0, 0), ytype = c(8, 8, 9, 9))
		if(obs == 2) surv = data.frame(id = rep(id, 4), time = c(0, 30, 0, ev_time), obs = c(0, 0, 0, 1), ytype = c(8, 8, 9, 9))
		surv
	}
	tte.data = do.call(rbind, mapply(makesurv, tte.data$id, tte.data$time, tte.data$obs, SIMPLIFY = F))
	
	# merging all datasets and writing to output file
	alldata = rbind(longi1, longi2, longi3, longi4, longi5, longi6, longi7, tte.data)
	write.table(alldata, file = paste0("data2/data", it, ".txt"), row.names = F)

	# keeping the median of each biomarker for estimation
	med = read.table("modfiles/med_biom_sim2.txt", header = T)
	med[it, 1] = it; med[it, 2] = med1; med[it, 3] = med2; med[it, 4] = med3; med[it, 5] = med4; med[it, 6] = med5; med[it, 7] = med6; med[it, 8] = med7
	write.table(med, file = paste0("modfiles/med_biom_sim2.txt"), row.names = F)
}


##############
# 5. Backward strategy on the strong scenario of correlations
##############

# function to create model text files with a specified list of biomarkers passed as argument
# create_modtxt_scen2(i) creates a txt model including biomarkers i (i can be a vector)
create_modtxt_scen2 = function(i){
  
  num = setdiff(c("1", "2", "3", "4", "5", "6", "7"), i)
  ref = readLines("modfiles/mod_1234567.txt")
  warn = 0
  l4 = unlist(strsplit(ref[4], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(num)){
    indb = c(grep(paste0("b0", num[j]), l4), grep(paste0("b1", num[j]), l4), grep(paste0("b2", num[j]), l4), grep(paste0("a", num[j], "_pop"), l4), grep(paste0("alpha1", num[j]), l4), grep(paste0("alpha2", num[j]), l4), grep(paste0("a", num[j], ", "), l4))
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l4)-1) warn = 1
  l4 = l4[-c(ind)]
  
  if (warn == 1) {der = paste0(str_sub(rev(l4)[2], 1, 7)) ; l4[c(length(l4)-1)] = der}
  ref[4] = paste(l4, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 8:14){
    for (j in num){
      l = unlist(strsplit(ref[ligne], split = " "))
      if (l[1] == paste0("m", j)) ind_ligne = c(ind_ligne, ligne)  
    }
  }
  
  ind = c()
  indb = c()
  l16 = unlist(strsplit(ref[16], split = " "))
  for (j in 1:length(num)){
    indb = grep(paste0("alpha1", num[j]), l16)
    ind = c(ind, indb)
  }
  l16 = l16[-c(ind)]
  ref[16] = paste(l16, collapse = " ")
  
  l17 = unlist(strsplit(ref[17], split = " "))
  l17 = l17[-c(ind)]
  ref[17] = paste(l17, collapse = " ")
  
  ind = c()
  indb = c()
  l25 = unlist(strsplit(ref[25], split = " "))
  for (j in 1:length(num)){
    indb = grep(num[j], l25)
    ind = c(ind, indb)
  }
  l25 = l25[-c(ind)]
  ref[25] = paste(l25, collapse = " ")
  
  ref = ref[-c(ind_ligne)]
  return(ref)
}

# function to create all possible model text files with k biomarkers
# and write them into the appropriate directory
create_combtxt_scen2 = function(k){
  create_combination = function(col) {
    txt = create_modtxt_scen2(as.character(col))
    path = paste0("intermediate_results/backward2/back", k, "/mod_", paste(col, collapse = ''), ".txt")
    writeLines(txt, path)
  }
  apply(combn(1:7, k), 2, create_combination)
}

# Function to adapt the median values in the model files for each simulation 
# the arguments are num: the biomamarkers involved and j: the simulation j
adapt_txt_scen2 = function(num, j){
  med = read.table("modfiles/med_biom_sim2.txt", header = T)
  if (length(num) == 7) {ref = readLines("modfiles/mod_1234567.txt")}
  else {ref = readLines(paste0("intermediate_results/backward2/back", length(num), "/mod_", paste0(num, collapse = ""), ".txt"))}
  l = unlist(strsplit(ref[length(num)+9], split = " "))
  l[which(l == "+alpha11*(m1-4.58)")] = paste0("+alpha11*(m1-", med$med1[med$sim == j], ")")
  l[which(l == "+alpha12*(m2-7.44)")] = paste0("+alpha12*(m2-", med$med2[med$sim == j], ")")
  l[which(l == "+alpha13*(m3-2.36)")] = paste0("+alpha13*(m3-", med$med3[med$sim == j], ")")
  l[which(l == "+alpha14*(m4-5.28)")] = paste0("+alpha14*(m4-", med$med4[med$sim == j], ")")
  l[which(l == "+alpha15*(m5-6.39)")] = paste0("+alpha15*(m5-", med$med5[med$sim == j], ")")
  l[which(l == "+alpha16*(m6-4.41)")] = paste0("+alpha16*(m6-", med$med6[med$sim == j], ")")
  l[which(l == "+alpha17*(m7-256.6)")] = paste0("+alpha17*(m7-", med$med7[med$sim == j], ")")
  l = paste(l, collapse = " ")
  m = unlist(strsplit(ref[length(num)+10], split = " "))
  m[which(m == "+alpha21*(m1-4.58)")] = paste0("+alpha21*(m1-", med$med1[med$sim == j], ")")
  m[which(m == "+alpha22*(m2-7.44)")] = paste0("+alpha22*(m2-", med$med2[med$sim == j], ")")
  m[which(m == "+alpha23*(m3-2.36)")] = paste0("+alpha23*(m3-", med$med3[med$sim == j], ")")
  m[which(m == "+alpha24*(m4-5.28)")] = paste0("+alpha24*(m4-", med$med4[med$sim == j], ")")
  m[which(m == "+alpha25*(m5-6.39)")] = paste0("+alpha25*(m5-", med$med5[med$sim == j], ")")
  m[which(m == "+alpha26*(m6-4.41)")] = paste0("+alpha26*(m6-", med$med6[med$sim == j], ")")
  m[which(m == "+alpha27*(m7-256.6)")] = paste0("+alpha27*(m7-", med$med7[med$sim == j], ")")
  m = paste(m, collapse = " ")
  ref[length(num)+9] = l
  ref[length(num)+10] = m
  path = paste0("intermediate_results/backward2/back", length(num), "/mod_", paste(num, collapse = ''), "_", j, ".txt")
  writeLines(ref, path)
}

# function to adapt the mlxtran file of simulation i involving 7 biomarker (full initial joint model)
setmodel_scen2 = function(i) {
  obj = readLines("modfiles/ref1234567.mlxtran")
  obj[4] = paste0("file = 'data", i, ".txt'")
  obj[58] = paste0("file = 'mod_1234567_", i, ".txt'")
  obj[143] = paste0("exportpath = 'mod", i, "'")
  write.table(obj, file = paste0("intermediate_results/backward2/back7/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
}

# function to fit the full initial joint model on simulation i
fit_fullmodels_scen2 = function(i){
  loadProject(paste0("intermediate_results/backward2/back7/mod", i, ".mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}

# function to create mlxtran files given a list of biomarkers and for a given dataset
# create_mlx_scen2(i) creates a mlxtran files removing biomarkers i and setting dataset j as associated data
create_mlx_scen2 = function(i, j){
  mod_ind = 1:7
  mod_ind = paste0(mod_ind[-c(as.integer(i))], collapse = "")
  mod_ind2 = seq(1, 7)[-c(as.integer(i))]
  adapt_txt_scen2(mod_ind2, j)
  
  ref7 = readLines("modfiles/ref1234567.mlxtran")
  ref7[4] = paste0("file = 'data", j, ".txt'") # Set correct data file
  ref7[58] = paste0("file = 'mod_", mod_ind, "_", j, ".txt'")# Set correct model file
  ref7[143] = paste0("exportpath = 'mod", j, "'")
  
  l11 = unlist(strsplit(ref7[11], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    ind = grep(i[j], l11)
    indb = c(indb, ind)
  }
  ind = c(indb, which(l11 == "continuous,")[1:length(i)])
  l11 = l11[-c(ind)]
  ref7[11] = paste(l11, collapse = " ") 
  
  l17 = unlist(strsplit(ref7[17], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = c(grep(paste0("b0", i[j]), l17), grep(paste0("b1", i[j]), l17), grep(paste0("b2", i[j]), l17), grep(paste0("a", i[j], "_pop"), l17), grep(paste0("omega_a", i[j]), l17), grep(paste0("alpha1", i[j]), l17), grep(paste0("alpha2", i[j]), l17))
    ind = c(ind, indb)
  }
  warn = 0
  if (rev(ind)[1] == length(l17)) warn = 1
  l17 = l17[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l17)[1], 1, 11), "}") ; l17[length(l17)] = der}
  ref7[17] = paste(l17, collapse = " ")
  
  l56 = unlist(strsplit(ref7[56], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l56)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l56)) warn = 1
  l56 = l56[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l56)[1], 1, 3), "}") ; l56[length(l56)] = der}
  ref7[56] = paste(l56, collapse = " ") 
  
  l70 = unlist(strsplit(ref7[70], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l70)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l70)) warn = 1
  l70 = l70[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l70)[1], 1, 2), "}") ; l70[length(l70)] = der}
  ref7[70] = paste(l70, collapse = " ") 
  
  l71 = unlist(strsplit(ref7[71], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l71)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l71)) warn = 1
  l71 = l71[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l71)[1], 1, 2), "}") ; l71[length(l71)] = der}
  ref7[71] = paste(l71, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 1:length(ref7)){
    l = unlist(strsplit(ref7[ligne], split = " "))
    for (j in 1:length(i)){
      if (length(grep(paste0("b0", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b1", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b2", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("alpha1", i[j], "_pop"), l)) != 0 | length(grep(paste0("alpha2", i[j], "_pop"), l)) != 0 | length(grep(paste0("b_", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("a_", i[j]), l, fixed = T)) != 0 | length(grep(paste0("omega_a", i[j]), l, fixed = T)) != 0 | length(grep(paste0("a", i[j], "_pop"), l, fixed = T)) != 0) ind_ligne = c(ind_ligne, ligne)
    }
  }
  ref7 = ref7[-c(ind_ligne)]
  return(ref7)
}



# function to create mlxtran files and data files for iteration k of the backward 
# load the results found at the previous iteration and create the new model and data files
iter_back_scen2 = function(k){
  for (i in 1:p){ 
    if( !file.exists(paste0("intermediate_results/backward2/back", k+1, "/mod", i, "/populationParameters.txt"))) next
    tab = read.table(paste0("intermediate_results/backward2/back", k+1, "/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    dat = read.table(paste0("intermediate_results/backward2/back", k+1, "/data", i, ".txt"), header = T, sep = ' ')
    stab = tab[c(which(tab$parameter %in% c("alpha11_pop", "alpha12_pop", "alpha13_pop", "alpha14_pop", "alpha15_pop", "alpha16_pop", "alpha17_pop"))), ]
    max = stab$parameter[stab$rse_sa == max(stab$rse_sa[is.nan(stab$rse_sa) == F])&is.nan(stab$rse_sa) == F]
    prev = substr(stab$parameter, 7, 7)
    num = substr(max, 7, 7)
    cur = prev[-which(prev == num)]
    numb = setdiff(c("1", "2", "3", "4", "5", "6", "7"), cur)
    if (tab$rse_sa[tab$parameter == max]>50){
      mlx = create_mlx_scen2(numb, i)
      write.table(mlx, file = paste0("intermediate_results/backward2/back", k, "/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
      dat = dat[-c(which(dat$ytype == as.integer(num))), ]
      write.table(dat, file = paste0("intermediate_results/backward2/back", k, "/data", i, ".txt"), row.names = F)
      if (k == 0){
        print(paste0("No biomarker kept for simulation ", i))
        tab = data.frame(parameter = NULL, value = NULL, se_sa = NULL, rse_sa = NULL)
        write.table(tab, paste0("intermediate_results/backward2/final/mod", i, ".txt"), row.names = F) 
      }
    }
    else{
      indnan = which(is.nan(stab$rse_sa) == T)
      if(length(indnan)>0){
        min = stab$parameter[stab$value == min(stab$value[indnan])]
        prev = substr(stab$parameter, 7, 7)
        num = substr(min, 7, 7)
        cur = prev[-which(prev == num)]
        numb = setdiff(c("1", "2", "3", "4", "5", "6", "7"), cur)
        mlx = create_mlx_scen2(numb, i)
        write.table(mlx, file = paste0("intermediate_results/backward2/back", k, "/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
        dat = dat[-c(which(dat$ytype == as.integer(num))), ]
        write.table(dat, file = paste0("intermediate_results/backward2/back", k, "/data", i, ".txt"), row.names = F)
      } 
      else{
        print(paste0("End for simulation ", i))
        write.table(tab, paste0("intermediate_results/backward2/final/mod", i, ".txt"), row.names = F)     
      }
    }
  }
}

# function to fit the model involving k biomarkers for a simulation i 
fit_all_models_scen2 = function(k){
  fitmodel = function(i){
    if( !file.exists(paste0("intermediate_results/backward2/back", k, "/mod", i, ".mlxtran"))) return(NA)
    loadProject(paste0("intermediate_results/backward2/back", k, "/mod", i, ".mlxtran"))
    runPopulationParameterEstimation()
    runStandardErrorEstimation()
  }
  sapply(1:p, fitmodel)
}

# creating indicator vector of length 7 (referring to each biomarker) for simulation i
# 1 if the biomarker is selected, 0 otherwise 
getdata_scen2 = function(i){
  tab = read.table(paste0("intermediate_results/backward2/final/mod", i, ".txt"), header = T, sep = '')
  1*!is.na(match(vars, tab$parameter))
}

##############
# 6. Evaluating ROC AUC under first and strong scenario of correlation
##############

# generates landmark dataset for dataset j and landmark time lmk
gen_landmark = function(lmk, j) {
  dat = read.table(paste0("data/data", j, ".txt"), header = T)
  
  # Only keep patients still at risk at t = lmk
  ind = sapply(1:N, function(i) setdiff(dat$time[dat$id == i & dat$ytype %in% 8:9], c(0, 30))<lmk)
  ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
  dat = dat[dat$id %in% ind, ]
  
  # Remove observations after landmark (bayesian prediction at t = lmk)
  ind = (dat$ytype %in% 1:7) & (dat$time>lmk)
  dat = dat[!ind, ]
  
  # Censoring all patients at t = lmk
  dat$obs[dat$ytype %in% 8:9] = 0
  dat$time[(dat$ytype %in% 8:9) & (dat$time>0)] = lmk
  
  # Export dataset into the appropriate directory
  write.table(dat, file = paste0("intermediate_results/evalAUC/data_lmk", j, ".txt"), row.names = F)
}

# generates survival data in R format for patients of simation j that are in the landmark dataset for time = lmk
gen_survdata = function(lmk, j) {
  dat = read.table(paste0("data/data", j, ".txt"), header = T)
  
  # Only keep patients still at risk at t = lmk
  ind = sapply(1:N, function(i) setdiff(dat$time[dat$id == i & dat$ytype %in% 8:9], c(0, 30))<lmk)
  ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
  dat = dat[dat$id %in% ind, ]
  
  # Translate survival data into R standard format
  dat = dat[dat$ytype %in% 8:9, ]	# only keep survival informations
  dat = dat[dat$time>0, ]	# removing Monolix-imposed starting observation time
  
  gen_Rformat = function(id){
    dd = dat[dat$id == id, ]
    if (all(dd$time == 30)) return(c(id = id, time = 30, obs = 0))
    if (dd$time[dd$ytype == 8] == 30) return(c(id = id, time = dd$time[dd$ytype == 9], obs = 2))
    if (dd$time[dd$ytype == 9] == 30) return(c(id = id, time = dd$time[dd$ytype == 8], obs = 1))
  }
  dat = data.frame(t(sapply(unique(dat$id), gen_Rformat)))
  
  # Export dataset into the appropriate directory
  write.table(dat, file = paste0("intermediate_results/evalAUC/data_surv", j, ".txt"), row.names = F)
}

# function for simulating individual parameters from a landmark time lmk, a simulation dataset j and a model
# truemod = T : the considered model is correct (biomarkers 1, 2, 3). Population parameters used are still estimated from the data
# truemod = F : the considered model is taken as the result of the backward selection procedure
gen_indv_params = function(lmk, j, truemod) {
  
  # First duplicate project files and directories
  if(truemod) {
    path2 = "intermediate_results/fit_truemodel"
    file.copy(paste0(path2, "/mod_123_", j, ".txt"), paste0("intermediate_results/evalAUC/mod_123_", j, ".txt"), overwrite = T)
  } else { 
    path2 = "intermediate_results/backward/final"
    file.copy(paste0(path2, "/", list.files("intermediate_results/backward/final", pattern = c(paste0("_", j, ".txt")))), "intermediate_results/evalAUC/", overwrite = T)
  }
  file.copy(paste0(path2, "/data", j, ".txt"), paste0("intermediate_results/evalAUC/data", j, ".txt"), overwrite = T)
  file.copy(paste0(path2, "/mod", j, ".mlxtran"), paste0("intermediate_results/evalAUC/mod", j, ".mlxtran"), overwrite = T)
  file.copy(paste0(path2, "/mod", j), "intermediate_results/evalAUC", recursive = T, overwrite = T)
  
  # load project and fix all population parameters
  loadProject(paste0("intermediate_results/evalAUC/mod", j, ".mlxtran"))
  setInitialEstimatesToLastEstimates()
  popparams = getPopulationParameterInformation()
  popparams$method <- "FIXED"
  setPopulationParameterInformation(popparams)
  
  # Replace data file
  BaseData = getData()
  setData(paste0("intermediate_results/evalAUC/data_lmk", j, ".txt"), BaseData$headerTypes, BaseData$observationTypes)
  
  # Simulation under the conditionnal distribution
  runPopulationParameterEstimation() # does nothing since all parameters are fixed (but still mandatory for further steps)
  setConditionalDistributionSamplingSettings(nbsimulatedparameters = 10)
  runConditionalDistributionSampling()
  
  # Export dataset into the appropriate directory
  write.table(getSimulatedIndividualParameters(), file = paste0("intermediate_results/evalAUC/data_params_", truemod, "_", j, ".txt"), row.names = F, sep = "\t")
}

# function for the computation of time-dependent AUC based on
# a landmark time lmk, a vector of horizon times hrz 
# a model (true model if truemod = T or result of backward selection if truemod = F)
# and a simulation replicate j
getAUC = function(lmk, hrz, j, truemod){
  # Reading data files (survival table, simulated individual parameters and median table)
  surv_table = read.table(paste0("intermediate_results/evalAUC/data_surv", j, ".txt"), h = T)
  sim = read.csv(paste0("intermediate_results/evalAUC/data_params_", truemod, "_", j, ".txt"), sep = "\t")
  med = read.table("modfiles/med_biom_sim.txt", header = T)
  
  calc_pred = function(x) {
    # Setting default values for parameters not estimated
    if (is.na(x["b01"])) x["b01"] = x["a1"] = x["b11"] = x["b21"] = x["alpha11"] = x["alpha21"] = 0
    if (is.na(x["b02"])) x["b02"] = x["b12"] = x["alpha12"] = x["alpha22"]  = 0
    if (is.na(x["b03"])) x["b03"] = x["b13"] = x["alpha13"] = x["alpha23"] = 0
    if (is.na(x["b04"])) x["b04"] = x["b14"] = x["alpha14"] = x["alpha24"] = 0
    if (is.na(x["b05"])) x["b05"] = x["b15"] = x["alpha15"] = x["alpha25"] = 0
    if (is.na(x["b06"])) x["b06"] = x["b16"] = x["alpha16"] = x["alpha26"] = 0
    if (is.na(x["b07"])) x["b07"] = x["b17"] = x["alpha17"] = x["alpha27"] = 0
    
    # Deriving individual biomarker trajectories and hazard functions for both risks 
    bm1 = function(tt) x["b01"]+x["a1"]*(exp(x["b11"]*tt)-exp(x["b21"]*tt))
    bm2 = function(tt) x["b02"]+x["b12"]*tt
    bm3 = function(tt) x["b03"]+x["b13"]*tt
    bm4 = function(tt) x["b04"]+x["b14"]*tt
    bm5 = function(tt) x["b05"]+x["b15"]*tt
    bm6 = function(tt) x["b06"]+x["b16"]*tt
    bm7 = function(tt) x["b07"]+x["b17"]*tt
    h1 = function(tt) {
      a = x["p1"]*x["g1"]*exp(-x["g1"]*tt)/(1-x["p1"]*(1-exp(-x["g1"]*tt))) 
      a = a*exp( x["alpha11"]*(bm1(tt)-med$med1[j]) )
      a = a*exp( x["alpha12"]*(bm2(tt)-med$med2[j]) )
      a = a*exp( x["alpha13"]*(bm3(tt)-med$med3[j]) )
      a = a*exp( x["alpha14"]*(bm4(tt)-med$med4[j]) )
      a = a*exp( x["alpha15"]*(bm5(tt)-med$med5[j]) )
      a = a*exp( x["alpha16"]*(bm6(tt)-med$med6[j]) )
      a = a*exp( x["alpha17"]*(bm7(tt)-med$med7[j]) )
      a
    }
    h2 = function(tt) {
      a = x["p2"]*x["g2"]*exp(-x["g2"]*tt)/(1-x["p2"]*(1-exp(-x["g2"]*tt))) 
      a = a*exp( x["alpha21"]*(bm1(tt)-med$med1[j]) )
      a = a*exp( x["alpha22"]*(bm2(tt)-med$med2[j]) )
      a = a*exp( x["alpha23"]*(bm3(tt)-med$med3[j]) )
      a = a*exp( x["alpha24"]*(bm4(tt)-med$med4[j]) )
      a = a*exp( x["alpha25"]*(bm5(tt)-med$med5[j]) )
      a = a*exp( x["alpha26"]*(bm6(tt)-med$med6[j]) )
      a = a*exp( x["alpha27"]*(bm7(tt)-med$med7[j]) )
      a
    }
    
    # Cumulative incidence functions
    F1 = function(tt) 1-exp(-(integrate(h1, 0, tt, stop.on.error = F)$value))
    F2 = function(tt) 1-exp(-(integrate(h2, 0, tt, stop.on.error = F)$value))
    F1 = Vectorize(F1)
    F2 = Vectorize(F2)
    
    # Compute survival probabilities for the desired horizon times
    if( (h1(1000) == Inf)){
      if (h1(hrz) == Inf){
        1-(1+F2(1000)-1-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk))
      }
      else if(h2(1000) == Inf){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))}
      else{1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk))}
        
    } else if((h2(1000) == Inf)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk)) 
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  # Making individual predictions for each patient and each required time
  preds = data.frame(id = sim$id, rep = sim$rep, pred = apply(sim, 1, calc_pred))
  preds = aggregate(.~id, preds, mean)	# Aggregating individual predictions across simulation replicates
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  # Returning ROC AUC for the required time horizons
  timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)$AUC_2[2]
}


##############
# 7. Scalability of the algorithm 
##############

# function to create mlxtran files involving N = i patients and 3 biomarkers 
create_model_scaN = function(i){
	# First copy data file in the new directory
	file.copy(paste0("data_scalability/data_n", i, ".txt"), paste0("intermediate_results/fit_scalability/data_n", i, ".txt"))
	file.copy("modfiles/mod_123.txt", "intermediate_results/fit_scalability/mod_123.txt")

	# Then, modify the mlxtran parameter file
	obj = readLines("modfiles/ref123.mlxtran")
	obj[4] = paste0("file = 'data_n", i, ".txt'")
	obj[95] = paste0("exportpath = 'mod", i, "'")
	write.table(obj, file = paste0("intermediate_results/fit_scalability/mod_n", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
}

# function to create mlxtran files involving N = 300 patients and K = i biomarkers
create_model_scaK = function(i){
	seq = paste0(1:i, collapse = "")
	# First copy data file in the new directory
	file.copy(paste0("data_scalability/data_bm", i, ".txt"), paste0("intermediate_results/fit_scalability/data_bm", i, ".txt"))
	file.copy(paste0("modfiles/mod_", seq, ".txt"), paste0("intermediate_results/fit_scalability/mod_", seq, ".txt"))

	# Then, modify the mlxtran parameter file
	obj = readLines(paste0("modfiles/ref", seq, ".mlxtran"))
	obj[4] = paste0("file = 'data_bm", i, ".txt'")
	write.table(obj, file = paste0("intermediate_results/fit_scalability/mod_k", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
}

# function to run the model involving N = i patients and K = 3 biomarkers
cputime_n = function(i){
	loadProject(paste0("intermediate_results/fit_scalability/mod_n", i, ".mlxtran"))
	runScenario()
}

# function to run the model involving N = 300 patients and K = k biomarkers
cputime_k = function(i){
	loadProject(paste0("intermediate_results/fit_scalability/mod_k", i, ".mlxtran"))
	runScenario()
}

# function to read the summary file of the model fit
# and returns the cpu time spent to estimate parameters and the Fisher Information Matrix
get_cpu = function(i){
  res = readLines(paste0("intermediate_results/fit_scalability/mod", i, "/summary.txt"))
  l1 = unlist(strsplit(res[47], split = " "))
  l2 = unlist(strsplit(res[88], split = " "))
  l3 = unlist(strsplit(res[87], split = " "))
  cpu_par = as.numeric(rev(l1)[1])
  cpu_fim = as.numeric(rev(l2)[1])
  if (is.na(cpu_fim == T)) cpu_fim = as.numeric(rev(l3)[1])
  c(cpu_par, cpu_fim)
}