#### This file is used to derived the flowchart in Figure 2 
#### Do not reproduce excatly the flowcharts because real data is not shared. 

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
}

###############
# 2. Creating Figure 2, Figure S2ab 
###############

# to load the data set descrbing the median values of each biomarker ###
data = read.table("datas/alldata.txt", header = T)
med = data.frame(num = 1:59, med = sapply(1:59, function(x) median(data$obs[data$ytype == x])))

# to save text in Fig2.txt file 
sink("../results/Fig2.txt")

cat("Current scenario", "\n")
give_flowchart(0.25, 0.50, 100) # current scenario

#close the external connection
sink()

sink("../results/FigS2a.txt")

cat("Stringent scenario", "\n")
give_flowchart(0.15, 0.30, 60) # stringent scenario

#close the external connection
sink()

sink("../results/FigS2b.txt")

cat("Flexible scenario", "\n")
give_flowchart(0.35, 0.70, 140) # flexible scenario

#close the external connection
sink()