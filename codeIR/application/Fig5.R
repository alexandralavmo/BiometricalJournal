#### File to save plot for time-dependant AUC (baseline and joint model) 
#### Do not reproduce exactly Figure 5 because real data is not shared. 

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
source("eval_AUC_current.R")  # Executes the script that simulates individual parameters and derives ROC AUC for all landmark and horizon times under the current joint model 
			      # (uses subfolder "intermediate_results/multi_fits/current" and writes in subfolder "intermediate_results/evalAUC")
source("fit_bsl.R") # Executes the script that fits the baseline model (writes in subfolder "intermediate_results/bsl_fit")
source("eval_AUC_bsl.R") # Executes the scripts that derives ROC AUC for all landmark and horizon times under the baseline model (uses subfolder "intermediate_results/bsl_fit" 
			 # and writes in subfolder "intermediate_results/evalAUC")
}

###############
# 2. Creating Figure 5
###############

# Loading results 
load("intermediate_results/evalAUC/rocBS3.RData")
load("intermediate_results/evalAUC/rocBS6.RData")
load("intermediate_results/evalAUC/rocBS9.RData")
load("intermediate_results/evalAUC/rocJM3.RData")
load("intermediate_results/evalAUC/rocJM6.RData")
load("intermediate_results/evalAUC/rocJM9.RData")

# Formatting data set to plot results 
dataAUC<-data.frame(temps = c(4:30, 7:30, 10:30, 4:30, 7:30, 10:30), 
                    AUC = c(roc_bsl3, roc_bsl6, roc_bsl9, roc3, roc6, roc9), 
                    grp = c(rep("score3", 27), rep("score6", 24), rep("score9", 21), rep("3", 27), rep("6", 24), rep("9", 21)))
dataAUC$grp<-factor(dataAUC$grp, levels = c("score3", "score6", "score9", "3", "6", "9"))
dataAUC$grplmk = ifelse(dataAUC$grp == "score3"|dataAUC$grp == '3', "3", ifelse(dataAUC$grp == "score6"|dataAUC$grp == '6', "6", "9"))
dataAUC$grplmk<-factor(dataAUC$grplmk, levels = c("3", "6", "9"))
dataAUC$grpline = ifelse(dataAUC$grp == "score3"|dataAUC$grp == 'score6'|dataAUC$grp == 'score9', "baseline score", "baseline score + longitudinal neutrophils and CRP")
dataAUC$grpline<-factor(dataAUC$grpline, levels = c("baseline score", "baseline score + longitudinal neutrophils and CRP"))
dataAUC$opacity = c(rep(1, 72), rep(0, 6), rep(1, 21), rep(0, 10), rep(1, 14), rep(0, 13), rep(1, 8))
dataAUC$opacity = factor(dataAUC$opacity, levels = c("1", "0"))

# Computing number of at-risk patients
data.surv = read.table("datas/sim_survival.txt", header =  T)
N3 = length(data.surv$id[data.surv$time >= 3])
N6 = length(data.surv$id[data.surv$time >= 6])
N9 = length(data.surv$id[data.surv$time >= 9])

# Plot results 
gpauc<-ggplot(data = dataAUC, aes(x = temps, y = AUC, group = grp, color = grplmk, fill = grplmk))+
  scale_color_manual(values = c("#99CCFF", "#FF9900", "red")) +
  scale_fill_manual(values = c("#99CCFF", "#FF9900", "red"))+
  geom_line(lwd = 0.5, aes(linetype = grpline))+ylim(0.2, 1)+ scale_linetype_manual(values = c("dashed", "solid"))+
  theme_classic()+
  labs(colour = "Landmark time", linetype = "")+
  geom_vline(xintercept = c(3, 6, 9), color = c("#99CCFF", "#FF9900", "red"), lwd = 0.3)+
  scale_x_continuous(name = "Time (days)", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.box = "vertical", 
        legend.text = element_text(size = 6), element_line(size = 0.5), axis.text = element_text(size = 9), 
        axis.title = element_text(size = 9), legend.key.size = unit(0.75, 'cm'), panel.grid.major.y = element_line())+
  annotate(geom = "text", x = 3, y = 1, label = paste0("N = ", N3), size = 1)+ 
  annotate(geom = "text", x = 6, y = 1, label = paste0("N = ", N6), size = 1)+ 
  annotate(geom = "text", x = 9, y = 1, label = paste0("N = ", N9), size = 1) 
gpauc


ggsave("Fig5.png", device = "png", path = "../results", width = 14.5, height = 8, units = "cm", dpi = 800)



