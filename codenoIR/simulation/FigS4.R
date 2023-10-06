#### Reproduce the upset plot in Figure S4 (strong scenario of correlations)

###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters2.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Warning: the intermediate results are not provided in this version of the code
# Please, refer to the README file for further details

# Executes the script that simulates datasets under the strong scenario of correlation (writes in subfolder "data2") 
# Warning: some differences in the 13rd digit after the coma point have been reported for some simulated values 
# depending on the operated system and may affect the results 
source("simulation2.R")

# Executes the script that performs the backward selection process on all simulated datasets under the strong scenario of correlation 
# (uses data stored in folder "data2" and writes in subfolder "intermediate_results/backward2")
# Warning: very time consuming (several days)... to be run on a computing cluster
source("backward2.R")

###############
# 2. Creating Figure S4
###############

# Variable names (for the plot)
vars = c("alpha17_pop", "alpha16_pop", "alpha15_pop", "alpha14_pop", "alpha13_pop", "alpha12_pop", "alpha11_pop")
bmnames = paste0("bm", 1:7)

# creating the data frame with the results over each simulation
mat = as.data.frame(t(sapply(1:p, getdata_scen2)))
colnames(mat) = rev(bmnames)

# Creating and saving Figure S4
upset(mat, colnames(mat), sort_sets = F, theme = upset_default_themes(text = element_text(size = 10)), set_sizes = (
  upset_set_size()
  + geom_text(aes(label = ..count..), hjust = 1.2, stat = 'count', size = 3)+ expand_limits(y = 110) + ylab('Set size')), 
  base_annotations = list(
    'Intersection size' = (
      intersection_size()
      + ylab('Biomarker combination size')
    )))

ggsave("FigS4.png", path = "../results", device = "png", width = 15, height = 10, units = "cm", dpi = 1000)
