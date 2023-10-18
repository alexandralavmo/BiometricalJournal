#### File to derive goodness of fit plots (3 plots: observations vs predictions, iwres and parameter distributions)
#### Do not reproduce exactly Figures S1abc because real data is not shared. 

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# set interresults = TRUE if you want to rely on intermediate results (beware that the intermediate results are not provided in this version)
# the version with intermediate results is available at https://github.com/alexandralavmo/BiometricalJournal/tree/main/codeIR)

interresults <- FALSE
if (!interresults) {
source("fit.R") # Executes this script that fits all the univariate joint models with the linear modeling (writes in subfolder "intermediate_results/uni_fits/lin")
		# Warning: time consuming... to be run on a computing cluster
source("fit_nl.R") # Executes the script that fits the univariate joint models with the nonlinear modeling (writes in subfolder "intermediate_results/uni_fits/nonlin")
		   # Warning: time consuming... to be run on a computing cluster
source("fit_multi_current.R") # Executes the script that performs the backward strategy under current scenario (uses subfolder "intermediate_results/uni_fits/" and 
			      # writes in subfolder "intermediate_results/multi_fits/current")
}

###############
# 2. Creating FigureS1abc 
###############

#### Observations vs predictions plots #####

# linear model definition 
mod_l = function(b0, b1, t) b0+b1*t

# Sampling times for longitudinal markers
tt2 = seq(0, 30, 2) # One observation every two days until day 30 for biomarker 2 
tt3 = seq(0, 29) # One observation every day until day 29 for biomarker 3

# Load the data sets 
data = read.table("intermediate_results/multi_fits/current/back2/data.txt", sep = ' ', header = T)
par = read.table("intermediate_results/multi_fits/current/back2/mod/IndividualParameters/estimatedIndividualParameters.txt", sep = ',', header = T)
parsim = read.table("intermediate_results/multi_fits/current/back2/mod/IndividualParameters/simulatedIndividualParameters.txt", sep = ',', header = T)
pop = read.table("intermediate_results/multi_fits/current/back2/mod/populationParameters.txt", sep = ',', header = T)

# seperate accoring to biomarkers
data_bm2 = data[data$ytype == 2, ]
data_bm3 = data[data$ytype == 3, ]

# Compute the predictions for each individual i and biomarker 2
for (i in unique(data_bm2$id)){
  b02 = par$b02_SAEM[par$id == i]
  b12 = par$b12_SAEM[par$id == i]
  t2 = data_bm2$time[data_bm2$id == i]
  data_bm2$pred[data_bm2$id == i] = mod_l(b02, b12, t2)
}

# Compute the predictions for each individual i and biomarker 3
for (i in unique(data_bm3$id)){
  b03 = par$b03_SAEM[par$id == i]
  b13 = par$b13_SAEM[par$id == i]
  t3 = data_bm3$time[data_bm3$id == i]
  data_bm3$pred[data_bm3$id == i] = mod_l(b03, b13, t3)
}

# spline for the relation between observations of bm2 and predictions
spl_bm2 = smooth.spline(x = data_bm2$pred, y = data_bm2$obs, spar = 1)

gp_bm2 = ggplot(data = NULL, aes(x = data_bm2$pred, y = data_bm2$obs))+geom_point(lwd = 0.07)+theme_classic()+geom_abline(intercept = 0, slope = 1, col = 'red', lwd = 0.4)+
  scale_x_continuous(limits = c(min(data_bm2$obs), max(data_bm2$obs)), name = bquote("Predicted neutrophils (x"*10^9*"/L)"))+
  scale_y_continuous(limits = c(min(data_bm2$obs), max(data_bm2$obs)), name = bquote("Observed neutrophils (x"*10^9*"/L)"))+
  geom_line(data = NULL, aes(x = spl_bm2$x, y = spl_bm2$y), linetype = 'dashed', col = 'red', lwd = 0.4)+
  theme(element_line(size = 0.5), axis.text = element_text(size = 5), axis.title = element_text(size = 7))

# spline for the relation between observations of bm3 and predictions
spl_bm3 = smooth.spline(x = data_bm3$pred, y = data_bm3$obs, spar = 1)

gp_bm3 = ggplot(data = NULL, aes(x = data_bm3$pred, y = data_bm3$obs))+geom_point(lwd = 0.07)+theme_classic()+geom_abline(intercept = 0, slope = 1, col = 'red', lwd = 0.4)+
  scale_x_continuous(limits = c(min(data_bm3$obs), max(data_bm3$obs)), name = "Predicted CRP (log(mg/L))")+
  scale_y_continuous(limits = c(min(data_bm3$obs), max(data_bm3$obs)), name = "Observed CRP (log(mg/L)")+
  geom_line(data = NULL, aes(x = spl_bm3$x, y = spl_bm3$y), linetype = 'dashed', col = 'red', lwd = 0.4)+
  theme(element_line(size = 0.5), axis.text = element_text(size = 5), axis.title = element_text(size = 7))


plot_grid(gp_bm2, gp_bm3, nrow = 1, labels = c("A", "B"), label_size = 9)
ggsave("FigS1a.png", path = "../results", device = "png", width = 15, height = 5, units = "cm", dpi = 500)


#### Individual Weighted Residuals ####

# get the residual error parameters 
sigma_a2 = pop$value[pop$parameter == 'a_2']
sigma_a3 = pop$value[pop$parameter == 'a_3']
# compute the individual weighted residuals 
iwres_bm2 = data.frame(time = data_bm2$time, iwres = (data_bm2$obs-data_bm2$pred)/sigma_a2, pred = data_bm2$pred)
iwres_bm3 = data.frame(time = data_bm3$time, iwres = (data_bm3$obs-data_bm3$pred)/sigma_a3, pred = data_bm3$pred)

# computed empirical individual weighted residuals for each biomarker and times from 0 to 30
emp_bm2 = sapply(0:30, emp_res, 2)
emp_bm3 = sapply(0:30, emp_res, 3)

gp_2 = ggplot(data = NULL)+geom_point(aes(x = iwres_bm2$time, y = iwres_bm2$iwres), size = 0.3, col = '#666666')+
  theme_classic()+
  theme(element_line(size = 0.5), axis.text = element_text(size = 9), 
        axis.title = element_text(size = 9))+
  geom_hline(yintercept = 0)+geom_hline(yintercept = 1.96)+geom_hline(yintercept = -1.96)+
  geom_line(data = NULL, aes(x = tt2, y = emp_bm2[1, ][tt2+1]), linetype = 'dashed', col = 'red')+
  geom_line(data = NULL, aes(x = tt2, y = emp_bm2[2, ][tt2+1]), linetype = 'dashed', col = 'red')+
  geom_line(data = NULL, aes(x = tt2, y = emp_bm2[3, ][tt2+1]), linetype = 'dashed', col = 'red')+
  scale_y_continuous(name = 'IWRES', limits = c(-5, 5))+scale_x_continuous(name = 'Time (days)')+
  geom_segment(aes(x = 0, xend = 2, y = 5, yend = 5), lwd = 0.3)+
  annotate(geom = 'text', x = 9, y = 5, label = "Theorical percentiles", size = 1.3)+
  geom_segment(aes(x = 0, xend = 2, y = 4.5, yend = 4.5), lwd = 0.3, linetype = 'dashed', col = 'red')+
  annotate(geom = 'text', x = 9, y = 4.5, label = "Empirical percentiles", size = 1.3)

gp_3 = ggplot(data = NULL)+geom_point(aes(x = iwres_bm3$time, y = iwres_bm3$iwres), size = 0.3, col = '#666666')+
  theme_classic()+
  theme(element_line(size = 0.5), axis.text = element_text(size = 9), 
        axis.title = element_text(size = 9))+
  geom_hline(yintercept = 0)+geom_hline(yintercept = 1.96)+geom_hline(yintercept = -1.96)+
  geom_line(data = NULL, aes(x = tt3, y = emp_bm3[1, ][tt3+1]), linetype = 'dashed', col = 'red')+
  geom_line(data = NULL, aes(x = tt3, y = emp_bm3[2, ][tt3+1]), linetype = 'dashed', col = 'red')+
  geom_line(data = NULL, aes(x = tt3, y = emp_bm3[3, ][tt3+1]), linetype = 'dashed', col = 'red')+
  scale_y_continuous(name = 'IWRES', limits = c(-5, 5))+scale_x_continuous(name = 'Time (days)')+
  geom_segment(aes(x = 0, xend = 2, y = 5, yend = 5), lwd = 0.3)+
  annotate(geom = 'text', x = 9, y = 5, label = "Theorical percentiles", size = 1.3)+
  geom_segment(aes(x = 0, xend = 2, y = 4.5, yend = 4.5), lwd = 0.3, linetype = 'dashed', col = 'red')+
  annotate(geom = 'text', x = 9, y = 4.5, label = "Empirical percentiles", size = 1.3)


plot_grid(gp_2, gp_3, nrow = 1, labels = c("A", "B"), label_size = 9)
ggsave("FigS1b.png", path = "../results", device = "png", width = 15, height = 7, units = "cm", dpi = 500)


#### Distribution of individual parameters ####

# creating parameter distribution plots
dist("b02");dist("b12");dist("b03");dist("b13")

# arranging plots 
plot2 = plot_grid(gp_b02, gp_b12, nrow = 1, ncol = 2, byrow = T)
plot3 = plot_grid(gp_b03, gp_b13, nrow = 1, ncol = 2, byrow = T)

plot_grid(plot2, plot3, nrow = 2, labels = c("A", "B"), label_size = 9)
ggsave("FigS1c.png", path = "../results", device = "png", width = 15, height = 7, units = "cm", dpi = 500)
