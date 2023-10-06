#### File to save ggplot objects for individual dynamic predictions for a given patient of the data set 
#### Do not reproduce exactly Figure 4 because real data is not shared. 

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Warning: the intermediate results are not provided in this version of the code
# Please, refer to the README file for further details

# Executes the script that fits all the univariate joint models with the linear modeling (writes in subfolder 
# "intermediate_results/uni_fits/lin")
# Warning: time consuming (several days)... to be run on a computing cluster
source("fit.R")  

# Executes the script that fits the univariate joint models with the nonlinear modeling (writes in subfolder 
# "intermediate_results/uni_fits/nonlin")
# Warning: time consuming (several days)... to be run on a computing cluster
source("fit_nl.R") 

# Executes the script that performs the backward strategy under current scenario (uses subfolder "intermediate_results/uni_fits/" 
# and writes in subfolder "intermediate_results/multi_fits/current")
source("fit_multi_current.R")

# Executes the script that simulates individual parameters and derives ROC AUC for all landmark and horizon times under the current joint model 
# (uses subfolder "intermediate_results/multi_fits/current" and writes in subfolder intermediate_results/evalAUC")
source("eval_AUC_current.R")
                         
###############
# 2. Creating Figure 4  
###############

# load the data sets
dataobs_all = read.table("datas/alldata.txt", header =  T)
data.surv = read.table("datas/sim_survival.txt", header =  T)

### Dynamic predictions with two random patients

# patA = "181"  is dead
# patB = "286"  is discharged 

# for patient A, creation of the ggplot objects (for each biomarker and survival probability, for each landmark)
plotA_n1 = landmark_neutro(3, 181) 
plotA_n2 = landmark_neutro(6, 181) 
plotA_n3 = landmark_neutro(9, 181) 
plotA_c1 = landmark_crp(3, 181) 
plotA_c2 = landmark_crp(6, 181) 
plotA_c3 = landmark_crp(9, 181) 

plotA_surv1 = landmark_survie(3, 181) 
plotA_surv2 = landmark_survie(6, 181) 
plotA_surv3 = landmark_survie(9, 181) 

# for patient B, creation of the ggplot objects (for each biomarker and survival probability)
plotB_n1 = landmark_neutro(3, 286) 
plotB_n2 = landmark_neutro(6, 286) 
plotB_n3 = landmark_neutro(9, 286) 
plotB_c1 = landmark_crp(3, 286) 
plotB_c2 = landmark_crp(6, 286) 
plotB_c3 = landmark_crp(9, 286) 
 
plotB_surv1 = landmark_survie(3, 286) 
plotB_surv2 = landmark_survie(6, 286) 
plotB_surv3 = landmark_survie(9, 286) 

# creating the legend of the figure
gleg = ggplot(data = NULL)+geom_rect(aes(xmin = 0, ymin = 0, xmax = 0.1, ymax = 0.05), col = 'grey', fill = 'grey')+xlim(0, 1)+ylim(0, 1.1)+
  theme_classic()+ annotate(geom = 'text', x = 0.3, y = 0.025, label = "95% prediction interval", size = 1.3)+
  geom_segment(aes(x = 0.05, xend = 0.05, y = 0.1, yend = 0.2), linetype = 'dashed', lwd = 0.2)+
  annotate(geom = 'text', x = 0.3, y = 0.15, label = "discharge", size = 1.3)+
  geom_segment(aes(x = 0.05, xend = 0.05, y = 0.25, yend = 0.35), lwd = 0.2)+
  annotate(geom = 'text', x = 0.3, y = 0.275, label = "death", size = 1.3)+
  geom_segment(aes(x = 0, xend = 0.1, y = 0.5, yend = 0.5), lwd = 0.2)+
  annotate(geom = 'text', x = 0.3, y = 0.5, label = "predicted survival", size = 1.3)+
  
  geom_segment(aes(x = 0, xend = 0.1, y = 0.7, yend = 0.7), lwd = 0.2, col = '#99FF66')+
  annotate(geom = 'text', x = 0.3, y = 0.7, label = "predicted neutrophil value", size = 1.3)+
  geom_segment(aes(x = 0, xend = 0.1, y = 0.6, yend = 0.6), lwd = 0.2, col = '#9933FF')+
  annotate(geom = 'text', x = 0.3, y = 0.6, label = "predicted CRP value", size = 1.3)+
  geom_point(aes(x = 0.05, y = 0.8), col = 'black', show.legend = F, lwd = 0.15)+
  annotate(geom = 'text', x = 0.3, y = 0.8, label = "observed marker value", size = 1.3)+
  geom_point(aes(x = 0.05, y = 0.9), col = 'blue', show.legend = F, lwd = 0.15)+
  annotate(geom = 'text', x = 0.3, y = 0.9, label = "future marker value", size = 1.3)+
  geom_segment(aes(x = 0.05, xend = 0.05, y = 0.95, yend = 1.05), lwd = 0.2, col = 'red')+
  annotate(geom = 'text', x = 0.3, y = 1, label = "landmark time", size = 1.3)+
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank()) 
gleg

# arranging the plots 
a = plot_grid(plotA_n1, plotA_n2, plotA_n3, plotA_c1, plotA_c2, plotA_c3, plotA_surv1, plotA_surv2, plotA_surv3, align = "h", nrow = 3, ncol = 3)
b = plot_grid(plotB_n1, plotB_n2, plotB_n3, plotB_c1, plotB_c2, plotB_c3, plotB_surv1, plotB_surv2, plotB_surv3, align = "h", nrow = 3, ncol = 3)
plot_grid(a, b, gleg, nrow = 2, ncol = 2, labels = c("A", "B", ""), rel_heights = c(0.75, 0.25))
# saving the figure
ggsave("Fig4.png", path = "../results", device = "png", width = 14.5, height = 14, units = "cm", dpi = 800)
