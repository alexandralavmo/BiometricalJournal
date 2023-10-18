###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# set interresults = TRUE if you want to rely on intermediate results (beware that the intermediate results are not provided in this version)
# the version with intermediate results is available at https://github.com/alexandralavmo/BiometricalJournal/tree/main/codeIR)

interresults <- FALSE
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
# 3. Derivation of relative errors
###############

est$e.b01 = (est$b01-params_nonlin$b01)/params_nonlin$b01
est$e.b11 = (est$b11-params_nonlin$b11)/params_nonlin$b11
est$e.b21 = (est$b21-params_nonlin$b21)/params_nonlin$b21
est$e.a1 = (est$a1-exp(params_nonlin$a1))/exp(params_nonlin$a1)
est$e.b02 = (est$b02-b0[1])/b0[1]
est$e.b12 = (est$b12-b1[1])/b1[1]
est$e.b03 = (est$b03-b0[2])/b0[2]
est$e.b13 = (est$b13-b1[2])/b1[2]

est$e.p1 = (est$p1-p1)/p1
est$e.g1 = (est$g1-g1)/g1

est$e.alpha11 = (est$alpha11-alpha[1])/alpha[1]
est$e.alpha12 = (est$alpha12-alpha[2])/alpha[2]
est$e.alpha13 = (est$alpha13-alpha[3])/alpha[3]

est$e.b1 = (est$b_1-params_nonlin$sigma_b1)/params_nonlin$sigma_b1 
est$e.a2 = (est$a_2-sigma_a[1])/sigma_a[1]
est$e.a3 = (est$a_3-sigma_a[2])/sigma_a[2]

est$e.omega_b01 = (est$omega_b01**2-params_nonlin$omega_b01)/params_nonlin$omega_b01
est$e.omega_b11 = (est$omega_b11**2-params_nonlin$omega_b11)/params_nonlin$omega_b11
est$e.omega_b21 = (est$omega_b21**2-params_nonlin$omega_b21)/params_nonlin$omega_b21
est$e.omega_a1 = (est$omega_a1**2-params_nonlin$omega_a1)/params_nonlin$omega_a1
est$e.omega_b02 = (est$omega_b02**2-omega_b0[1])/omega_b0[1]
est$e.omega_b12 = (est$omega_b12**2-omega_b1[1])/omega_b1[1]
est$e.omega_b03 = (est$omega_b03**2-omega_b0[2])/omega_b0[2]
est$e.omega_b13 = (est$omega_b13**2-omega_b1[2])/omega_b1[2]


###############
# 4. Violin plot (supplementary Figure S2)
###############

# Longitudinal parameters for biomarker 1 (nonlinear model)
liste1 = c(rep("μ01", nrow(est)), rep("μ11", nrow(est)), rep("μ21", nrow(est)), rep("μa1", nrow(est)), 
         rep("omega_01", nrow(est)), rep("omega_11", nrow(est)), 
         rep("omega_21", nrow(est)), rep("omega_a1", nrow(est)), 
         rep("sigma_b1", nrow(est)))
values1 = 100*unlist(est[, c("e.b01", "e.b11", "e.b21", "e.a1", "e.omega_b01", "e.omega_b11", "e.omega_b21", "e.omega_a1", "e.b1")])
tb1 = data.frame(param = liste1, val = values1)
tb1$param = factor(tb1$param, levels = c("μ01", "μ11", "μ21", "μa1", "omega_01", "omega_11", "omega_21", "omega_a1", "sigma_b1"))

# violin plot for biomarker 1 (nonlinear model)
vioplot1 = ggplot(data = tb1, aes(x = param, y = val)) +
  geom_violin(fill = "#99CCFF", width = 0.6, scale = "width") + geom_boxplot(width = 0.07, outlier.size = 0.5, size = 0.2) + 
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed', lwd = 0.2)+
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 1")+
  scale_x_discrete(labels = c("μ01" = expression(mu["01"]), "μ11" = expression(mu["11"]), "μ21" = expression(mu["21"]), 
                              "μa1" = expression(mu["a1"]), "omega_01" = expression(omega["01"]), "omega_11" = expression(omega["11"]), 
                              "omega_21" = expression(omega["21"]), "omega_a1" = expression(omega["a1"]), "sigma_b1" = expression(sigma["b1"])))+
  scale_y_continuous(limits = c(-100, 100))+
  theme(axis.text.x = element_text(hjust = 1, size = 8), axis.title = element_text(size = 10), axis.text.y = element_text(size = 8))

# Longitudinal parameters for biomarker 2 (linear model)
liste2 = c(rep("μ02", nrow(est)), rep("μ12", nrow(est)), 
         rep("omega_02", nrow(est)), rep("omega_12", nrow(est)), 
         rep("sigma_a2", nrow(est)))
values2 = 100*unlist(est[, c("e.b02", "e.b12", "e.omega_b02", "e.omega_b12", "e.a2")])
tb2 = data.frame(param = liste2, val = values2)
tb2$param = factor(tb2$param, levels = c("μ02", "μ12", "omega_02", "omega_12", "sigma_a2"))

# violin plot for biomarker 2 (linear model)
vioplot2 = ggplot(data = tb2, aes(x = param, y = val)) +
  geom_violin(fill = "#99CCFF", width = 0.6, scale = "width") + geom_boxplot(width = 0.07, outlier.size = 0.5, size = 0.2) + 
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed', lwd = 0.2)+
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 2")+
  scale_x_discrete(labels = c("μ02" = expression(mu["02"]), "μ12" = expression(mu["12"]), "omega_02" = expression(omega["02"]), 
                              "omega_12" = expression(omega["12"]), "sigma_a2" = expression(sigma["a2"])))+
  scale_y_continuous(limits = c(-100, 100))+
  theme(axis.text.x = element_text(hjust = 1, size = 8), axis.title = element_text(size = 10), axis.text.y = element_text(size = 8))

# Longitudinal parameters for biomarker 3 (linear model)
liste3 = c(rep("μ03", nrow(est)), rep("μ13", nrow(est)), 
         rep("omega_03", nrow(est)), rep("omega_13", nrow(est)), 
         rep("sigma_a3", nrow(est)))
values3 = 100*unlist(est[, c("e.b03", "e.b13", "e.omega_b03", "e.omega_b13", "e.a3")])
tb3 = data.frame(param = liste3, val = values3)
tb3$param = factor(tb3$param, levels = c("μ03", "μ13", "omega_03", "omega_13", "sigma_a3"))

# violin plot for biomarker 3 (linear model)
vioplot3 = ggplot(data = tb3, aes(x = param, y = val)) +
  geom_violin(fill = "#99CCFF", width = 0.6, scale = "width") + geom_boxplot(width = 0.07, outlier.size = 0.5, size = 0.2) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed', lwd = 0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 3")+
  scale_y_continuous(limits = c(-100, 100))+
  scale_x_discrete(labels = c("μ03" = expression(mu["03"]), "μ13" = expression(mu["13"]), "omega_03" = expression(omega["03"]), 
                              "omega_13" = expression(omega["13"]), "sigma_a3" = expression(sigma["a3"])))+
  theme(axis.text.x = element_text(hjust = 1, size = 8), axis.title = element_text(size = 10), axis.text.y = element_text(size = 8))

# Survival parameters
liste4 = c(rep("p1", nrow(est)), rep("g1", nrow(est)), 
         rep("alpha11", nrow(est)), rep("alpha12", nrow(est)), 
         rep("alpha13", nrow(est)))
values4 = 100*unlist(est[, c("e.p1", "e.g1", "e.alpha11", "e.alpha12", "e.alpha13")])
tb4 = data.frame(param = liste4, val = values4)
tb4$param = factor(tb4$param, levels = c("p1", "g1", "alpha11", "alpha12", "alpha13"))

# violin plot for survival parameters
vioplot4 = ggplot(data = tb4, aes(x = param, y = val)) +
  geom_violin(fill = "#99CCFF", width = 0.6, scale = "width") + geom_boxplot(width = 0.07, outlier.size = 0.5, size = 0.2) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed', lwd = 0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Survival parameters")+
  scale_y_continuous(limits = c(-100, 200))+
  scale_x_discrete(labels = c("p1" = expression(p["1"]), "g1" = expression(g["1"]), "alpha11" = expression(alpha["11"]), 
                              "alpha12" = expression(alpha["12"]), "alpha13" = expression(alpha["13"])))+
  theme(axis.text.x = element_text(hjust = 1, size = 8), axis.title = element_text(size = 10), axis.text.y = element_text(size = 8))

# arranging the plots and saving figure 
plot_grid(vioplot1, vioplot2, vioplot3, vioplot4, nrow = 4)

ggsave("FigS3.png", device = "png", path = "../results", width = 15, height = 20, units = "cm", dpi = 1000)