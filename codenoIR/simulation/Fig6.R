##### To reproduce the upset plot in Figure 6 (first scenario of correlation)

###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Executes the script that simulates datasets under the first scenario of correlation (writes in subfolder "data")
# Warning: some differences in the 13rd digit after the coma point have been reported for some simulated values 
# depending on the operating system and may affect the results 
# Comment the following line if you want to rely on pre-computed intermediate results. Beware that these results are not 
# included by default in this version of the code. They can be downloaded at https://iame.catibiomed.fr/index.php/s/uwuLTLw6ZboYBah
source("simulation.R")

# Executes the script that performs the backward selection process on all simulated datasets under the first scenario of correlation 
# (uses data stored in folder "data" and writes in subfolder "intermediate_results/backward")
# Warning: very time consuming (several days)... to be run on a computing cluster
# Comment the following line if you want to rely on pre-computed intermediate results. Beware that these results are not 
# included by default in this version of the code. They can be downloaded at https://iame.catibiomed.fr/index.php/s/uwuLTLw6ZboYBah
source("backward.R") 

###############
# 2. Creating Figure 6 
###############

# Variable names (for the plot)
vars = c("alpha17_pop", "alpha16_pop", "alpha15_pop", "alpha14_pop", "alpha13_pop", "alpha12_pop", "alpha11_pop")
bmnames = paste0("bm", 1:7)

#creating the data frame with the results of each simulation
mat = as.data.frame(t(sapply(1:p, getdata)))
colnames(mat) = rev(bmnames)

# Creating and saving Figure 6
upset(mat, colnames(mat), sort_sets = F, theme = upset_default_themes(text = element_text(size = 10)), 
	set_sizes = ( upset_set_size() + geom_text(aes(label = ..count..), hjust = 1.2, stat = 'count', size = 3)
		+ expand_limits(y = 110) + ylab('Set size')), 
	base_annotations = list( 'Intersection size' = (intersection_size() + ylab('Biomarker combination size')))
    )

ggsave("Fig6.png", path = "../results", device = "png", width = 15, height = 10, units = "cm", dpi = 1000)

