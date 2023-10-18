#### File to summarize biomarker selected for multivariate analysis in Table S4
#### Do not reproduce excatly Table S4 because real data is not shared. 

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
}

###############
# 2. Creating Table S4 
###############

# to load the data set descrbing the median values of each biomarker 
data = read.table("datas/alldata.txt", header = T)

cat = read.table("datas/category.txt", header = T)
load("intermediate_results/uni_fits/selec_lin.RData")
load("intermediate_results/uni_fits/selec_nlin.RData")


sapply(1:8, function (i) assign(paste0("cat", i), data.frame(biomarqueur = NA, N = NA, n = NA, alpah1k = NA, RSEforalpha1k = NA, log10pvalue = NA, longitudinal_submodel = NA), envir = .GlobalEnv))  

for (i in selec_lin){
  catg = cat$category[cat$num == i] 
  trans = get(paste0("cat", catg))
  if (i %in% selec_lin){
    t = read.table(paste0("intermediate_results/uni_fits/lin/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    trans$longitudinal_submodel = "linear"
  }
  trans$biomarqueur = i
  trans$N = length(unique(data$id[data$ytype == i]))
  trans$n = nrow(data[data$ytype == i, ])/length(unique(data$id[data$ytype == i]))
  trans$alpah1k = t$value[t$parameter == "alpha1_pop"]
  trans$RSEforalpha1k = t$rse_sa[t$parameter == "alpha1_pop"]
  stat = t$value[t$parameter == "alpha1_pop"]/t$se_sa[t$parameter == "alpha1_pop"]
  trans$log10pvalue = -log10(2*pnorm(abs(stat), 0, 1, lower.tail = F))
  assign(paste0("cat", catg), rbind(get(paste0("cat", catg)), trans))
}

cat1 = cat1[order(cat1$log10pvalue, decreasing = T, na.last = F), ]
cat2 = cat2[order(cat2$log10pvalue, decreasing = T, na.last = F), ]
cat3 = cat3[order(cat3$log10pvalue, decreasing = T, na.last = F), ]
cat4 = cat4[order(cat4$log10pvalue, decreasing = T, na.last = F), ]
cat5 = cat5[order(cat5$log10pvalue, decreasing = T, na.last = F), ]
cat6 = cat6[order(cat6$log10pvalue, decreasing = T, na.last = F), ]
cat7 = cat7[order(cat7$log10pvalue, decreasing = T, na.last = F), ]
cat8 = cat8[order(cat8$log10pvalue, decreasing = T, na.last = F), ]

cat1[1, ] = c("Category 1", "", "", "", "", "", "")
cat2[1, ] = c("Category 2", "", "", "", "", "", "")
cat3[1, ] = c("Category 3", "", "", "", "", "", "")
cat4[1, ] = c("Category 4", "", "", "", "", "", "")
cat5[1, ] = c("Category 5", "", "", "", "", "", "")
cat6[1, ] = c("Category 6", "", "", "", "", "", "")
cat7[1, ] = c("Category 7", "", "", "", "", "", "")
cat8[1, ] = c("Category 8", "", "", "", "", "", "")



tabS4 = rbind(cat1, cat2, cat3, cat4, cat5, cat6, cat7, cat8)

save(tabS4, file = "../results/TabS4.RData")
