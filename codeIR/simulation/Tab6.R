###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# set interresults = FALSE if you want to (re-) compute the intermediate results

interresults <- TRUE
if (!interresults) {
source("simulation.R") # Executes the script that simulates datasets under the first scenario of correlation (writes in subfolder "data")
		       # Warning: some differences in the 13rd digit after the coma point have been reported for some simulated values 
		       # depending on the operating system and may affect the results 
source("fit.R") # Executes the script that fits the true joint model on all the simulated data sets (uses data stored in subfolder "data" 
		# and writes in subfolder "intermediate_results/fit_truemodel")
		# Warning: very time consuming... to be run on a computing cluster
}

###############
# 2. Importation of results
###############

est = sapply(1:p, function(i) read.table(paste0("intermediate_results/fit_truemodel/pop", i, ".txt"), sep = ' ', header = T)[, 2])
rownames(est) = read.table("intermediate_results/fit_truemodel/pop1.txt", sep = ',', header = T)[, 1]
est = data.frame(t(est))

###############
# 3. Derivation of relative bias
###############

# for fixed effects
rb01 = (mean(est$b01)-params_nonlin$b01)/params_nonlin$b01
rb11 = (mean(est$b11)-params_nonlin$b11)/params_nonlin$b11
rb21 = (mean(est$b21)-params_nonlin$b21)/params_nonlin$b21
rba1 = (mean(est$a1)-exp(params_nonlin$a1))/exp(params_nonlin$a1)
rb02 = (mean(est$b02)-b0[1])/b0[1]
rb12 = (mean(est$b12)-b1[1])/b1[1]
rb03 = (mean(est$b03)-b0[2])/b0[2]
rb13 = (mean(est$b13)-b1[2])/b1[2]

# for random effects
rbomega01 = (mean(est$omega_b01**2)-params_nonlin$omega_b01)/params_nonlin$omega_b01
rbomega11 = (mean(est$omega_b11**2)-params_nonlin$omega_b11)/params_nonlin$omega_b11
rbomega21 = (mean(est$omega_b21**2)-params_nonlin$omega_b21)/params_nonlin$omega_b21
rbomegaa1 = (mean(est$omega_a1**2)-params_nonlin$omega_a1)/params_nonlin$omega_a1
rbomega02 = (mean(est$omega_b02**2)-omega_b0[1])/omega_b0[1]
rbomega12 = (mean(est$omega_b12**2)-omega_b1[1])/omega_b1[1]
rbomega03 = (mean(est$omega_b03**2)-omega_b0[2])/omega_b0[2]
rbomega13 = (mean(est$omega_b13**2)-omega_b1[2])/omega_b1[2]

# for error parameters
rbsigmab1 = (mean(est$b_1)-params_nonlin$sigma_b1)/params_nonlin$sigma_b1
rbsigmaa2 = (mean(est$a_2)-sigma_a[1])/sigma_a[1]
rbsigmaa3 = (mean(est$a_3)-sigma_a[2])/sigma_a[2]

# for survival parameters
rbp1 = (mean(est$p1)-p1)/p1
rbg1 = (mean(est$g1)-g1)/g1
rbalpha1 = (mean(est$alpha11)-alpha[1])/alpha[1]
rbalpha2 = (mean(est$alpha12)-alpha[2])/alpha[2]
rbalpha3 = (mean(est$alpha13)-alpha[3])/alpha[3]


###############
# 4. Derivation of relative root mean square errors 
###############

rrmseb01 = sqrt(mean((est$b01-params_nonlin$b01)^2))/params_nonlin$b01
rrmseb11 = sqrt(mean((est$b11-params_nonlin$b11)^2))/params_nonlin$b11
rrmseb21 = sqrt(mean((est$b21-params_nonlin$b21)^2))/params_nonlin$b21
rrmsea1 = sqrt(mean((est$a1-exp(params_nonlin$a1))^2))/exp(params_nonlin$a1)
rrmseb02 = sqrt(mean((est$b02-b0[1])^2))/b0[1]
rrmseb12 = sqrt(mean((est$b12-b1[1])^2))/b1[1]
rrmseb03 = sqrt(mean((est$b03-b0[2])^2))/b0[2]
rrmseb13 = sqrt(mean((est$b13-b1[2])^2))/b1[2]


rrmsesigmab1 = sqrt(mean((est$b_1-params_nonlin$sigma_b1)^2))/params_nonlin$sigma_b1
rrmsesigmaa2 = sqrt(mean((est$a_2-sigma_a[1])^2))/sigma_a[1]
rrmsesigmaa3 = sqrt(mean((est$a_3-sigma_a[2])^2))/sigma_a[2]

rrmseomega01 = sqrt(mean((est$omega_b01**2-params_nonlin$omega_b01)^2))/params_nonlin$omega_b01
rrmseomega11 = sqrt(mean(((est$omega_b11**2-params_nonlin$omega_b11)/params_nonlin$omega_b11)^2))
rrmseomega21 = sqrt(mean((est$omega_b21**2-params_nonlin$omega_b21)^2))/params_nonlin$omega_b21
rrmseomegaa1 = sqrt(mean((est$omega_a1-params_nonlin$omega_a1)^2))/params_nonlin$omega_a1
rrmseomega02 = sqrt(mean((est$omega_b02**2-omega_b0[1])^2))/omega_b0[1]
rrmseomega12 = sqrt(mean((est$omega_b12**2-omega_b1[1])^2))/omega_b1[1]
rrmseomega03 = sqrt(mean((est$omega_b03**2-omega_b0[2])^2))/omega_b0[2]
rrmseomega13 = sqrt(mean((est$omega_b13**2-omega_b1[2])^2))/omega_b1[2]


rrmsep1 = sqrt(mean((est$p1-p1)^2))/p1
rrmseg1 = sqrt(mean((est$g1-g1)^2))/g1
rrmsealpha1 = sqrt(mean((est$alpha11-alpha[1])^2))/alpha[1]
rrmsealpha2 = sqrt(mean((est$alpha12-alpha[2])^2))/alpha[2]
rrmsealpha3 = sqrt(mean((est$alpha13-alpha[3])^2))/alpha[3]


###############
# 5. Generate the Table 6 
###############

name_bm1 = c("mu01", "mu11", "mu21", "mua1", "omega01", "omega11", "omega21", "omegaa1", "sigmab1")
name_bm2 = c("mu02", "mu12", "omega02", "omega12", "sigmaa2")
name_bm3 = c("mu03", "mu13", "omega03", "omega13", "sigmaa3")
name_surv = c("p1", "g1", "alpha1", "alpha2", "alpha3")

true_bm1 = as.numeric(params_nonlin)
true_bm1[4] = exp(true_bm1[4])
true_bm1[5:8] = sqrt(true_bm1[5:8])
true_bm2 = c(b0[1], b1[1], sqrt(omega_b0[1]), sqrt(omega_b1[1]), sigma_a[1])
true_bm3 = c(b0[2], b1[2], sqrt(omega_b0[2]), sqrt(omega_b1[2]), sigma_a[2])
true_surv = c(p1, g1, alpha[1:3])

RBbm1 = signif(c(rb01, rb11, rb21, rba1, rbomega01, rbomega11, rbomega21, rbomegaa1, rbsigmab1)*100, digits = 2) 
RBbm2 = signif(c(rb02, rb12, rbomega02, rbomega12, rbsigmaa2)*100, digits = 2) 
RBbm3 = signif(c(rb03, rb13, rbomega03, rbomega13, rbsigmaa3)*100, digits = 2) 
RBsurv = signif(c(rbp1, rbg1, rbalpha1, rbalpha2, rbalpha3)*100, digits = 2)

RRMSEbm1 = signif(abs(c(rrmseb01, rrmseb11, rrmseb21, rrmsea1, rrmseomega01, rrmseomega11, rrmseomega21, rrmseomegaa1, rrmsesigmab1))*100, digits = 2)  
RRMSEbm2 = signif(abs(c(rrmseb02, rrmseb12, rrmseomega02, rrmseomega12, rrmsesigmaa2))*100, digits = 2) 
RRMSEbm3 = signif(abs(c(rrmseb03, rrmseb13, rrmseomega03, rrmseomega13, rrmsesigmaa3))*100, digits = 2) 
RRMSEsurv = signif(abs(c(rrmsep1, rrmseg1, rrmsealpha1, rrmsealpha2, rrmsealpha3))*100, digits = 2)

tablong1 = data.frame(Parameter = name_bm1, Truevalue = true_bm1, RelativeBias = RBbm1, RRMSE = RRMSEbm1)
tablong2 = data.frame(Parameter = name_bm2, Truevalue = true_bm2, RelativeBias = RBbm2, RRMSE = RRMSEbm2)
tablong3 = data.frame(Parameter = name_bm3, Truevalue = true_bm3, RelativeBias = RBbm3, RRMSE = RRMSEbm3)
tabsurv =  data.frame(Parameter = name_surv, Truevalue = true_surv, RelativeBias = RBsurv, RRMSE = RRMSEsurv)

tab6 = rbind(c("Longitudinal parameters for bm1", "", "", "", ""), tablong1,
             c("Longitudinal parameters for bm2", "", "", "", ""), tablong2,
             c("Longitudinal parameters for bm3", "", "", "", ""), tablong3,
             c("Survival parameters", "", "", "", ""), tabsurv)

save(tab6, file = "../results/Tab6.RData")