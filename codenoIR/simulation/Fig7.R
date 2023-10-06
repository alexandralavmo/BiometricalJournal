#### To reproduce the violin plots for the distribution of AUC under true and final models in Figure 7

###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Warning: the intermediate results are not provided in this version of the code
# Please, refer to the README file for further details

# Executes the script that simulates datasets under the first scenario of correlation (writes in subfolder "data") 
# Warning: some differences in the 13rd digit after the coma point have been reported for some simulated values 
# depending on the operated system and may affect the results 
source("simulation.R")

# Executes the script that fits the true joint model on all the simulated data sets (uses data stored in subfolder "data" and 
# writes in subfolder "intermediate_results/fit_truemodel")
# Warning: very time consuming (several days)... to be run on a computing cluster
source("fit.R") 

# Executes the script that performs the backward selection process on all simulated datasets under the first scenario of correlation 
# (uses data stored in folder "data" and writes in subfolder "intermediate_results/backward")
# Warning: very time consuming (several days)... to be run on a computing cluster
source("backward.R") 

# Executes the script that computes ROC AUC under the true and final models (uses subfolders "intermediate_results/fit_truemodel" 
# and "intermediate_results/backward", and writes in subfolder "intermediate_results/eval_AUC") 
source("eval_AUC.R")

###############
# 1. Importation of results
###############

load("intermediate_results/evalAUC/trueroc.RData")
load("intermediate_results/evalAUC/finalroc.RData")


###############
# 2. Data summarizing results
###############

box = data.frame(sim = rep(1:100, 2), roc = c(trueroc, finalroc), grp = c(rep("True model", 100), rep("Final model", 100)), paired = rep(1:100, 2))
box$grpp = ifelse(finalroc == trueroc, 1, 0)
box$grp = factor(box$grp, levels = c("True model", "Final model"))
box$grpp = factor(box$grpp, levels = c("1", "0"))

###############
# 3. t-test
###############
signif(t.test(trueroc, finalroc, paired = T)$p.value, 1)

###############
# 3. Creating and saving Figure 7
###############
ggplot(data = box, aes(x = grp, y = roc))+geom_violin(aes(fill = grp))+scale_x_discrete(name = " ")+
  geom_boxplot(width = 0.4)+
  scale_fill_manual(values = c("#99CCFF", "#999CFF"))+
  geom_line(aes(group = paired, col = grpp), show.legend = F, lwd = 0.4)+
  geom_point(aes(group = paired, col = grpp), show.legend = F, size = 0.4)+
  scale_color_manual(values = c("#CCCCCC", "red"))+
  scale_y_continuous(name = "ROC AUC")+theme_classic()+
  theme(legend.position = "none")+
  annotate(geom = 'text', x = 0.7, y = 1, label = "Paired t-test: p = 0.007", size = 3)


ggsave("Fig7.png", device = "png", path = "../results", width = 14.5, height = 8, units = "cm", dpi = 800)