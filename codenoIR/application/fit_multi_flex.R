#### File used to perform the backward on flexible scenario 

###############
# 1. Loading and formating data to keep the three biomarkers involved 
###############

# first step: group the biomarkers following the initial classification 
# and keep only the most associated in each category

data = read.table("datas/alldata.txt", header = T)
cat = read.table("modfiles/category.txt", header = T) # data set giving biomarker classification
#creating an empty data set for each biomarker category (8 in total)
sapply(1:8, function (i) assign(paste0("cat", i), data.frame(biomarqueur = NA, N = NA, n = NA, alpah1k = NA, RSEforalpha1k = NA, log10pvalue = NA, longitudinal_submodel = NA), envir = .GlobalEnv))  

load("intermediate_results/uni_fits/selec_lin_flex.RData")
load("intermediate_results/uni_fits/selec_nlin_flex.RData")

# writing in each category data frame 
for (i in c(selec_lin, selec_nlin)){
  catg = cat$category[cat$num == i] 
  trans = get(paste0("cat", catg))[1, ]
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

# order category data sets by decreasing p-value wald test
cat1 = cat1[order(cat1$log10pvalue, decreasing = T, na.last = F), ]
cat2 = cat2[order(cat2$log10pvalue, decreasing = T, na.last = F), ]
cat3 = cat3[order(cat3$log10pvalue, decreasing = T, na.last = F), ]
cat4 = cat4[order(cat4$log10pvalue, decreasing = T, na.last = F), ]
cat5 = cat5[order(cat5$log10pvalue, decreasing = T, na.last = F), ]
cat6 = cat6[order(cat6$log10pvalue, decreasing = T, na.last = F), ]
cat7 = cat7[order(cat7$log10pvalue, decreasing = T, na.last = F), ]
cat8 = cat8[order(cat8$log10pvalue, decreasing = T, na.last = F), ]

# writing the category name  
cat1[1, ] = c("Category 1", "", "", "", "", "", "")
cat2[1, ] = c("Category 2", "", "", "", "", "", "")
cat3[1, ] = c("Category 3", "", "", "", "", "", "")
cat4[1, ] = c("Category 4", "", "", "", "", "", "")
cat5[1, ] = c("Category 5", "", "", "", "", "", "")
cat6[1, ] = c("Category 6", "", "", "", "", "", "")
cat7[1, ] = c("Category 7", "", "", "", "", "", "")
cat8[1, ] = c("Category 8", "", "", "", "", "", "")

# selecting only the most associated biomarker in each category
sel = as.numeric(na.omit(sapply(1:8, function(i) get(paste0("cat", i))[2, 1])))

# formatting data 
data3biom = data[data$ytype %in% c(sel, 60, 61), ]
data3biom$ytype[data3biom$ytype == 7] = 1
data3biom$ytype[data3biom$ytype == 9] = 2
data3biom$ytype[data3biom$ytype == 24] = 3
data3biom$ytype[data3biom$ytype == 30] = 4
data3biom$ytype[data3biom$ytype == 38] = 5
data3biom$ytype[data3biom$ytype == 46] = 6
data3biom$ytype[data3biom$ytype == 60] = 7
data3biom$ytype[data3biom$ytype == 61] = 8
write.table(data3biom, "intermediate_results/multi_fits/flexible/back6/data.txt", row.names = F)


###############
# 2. Preparing the back6 directory and fit full model
###############

file.copy("modfiles/mod_123456.txt", "intermediate_results/multi_fits/flexible/back6/mod_123456.txt")
file.copy("modfiles/mod123456.mlxtran", "intermediate_results/multi_fits/flexible/back6/mod.mlxtran")

loadProject("intermediate_results/multi_fits/flexible/back6/mod.mlxtran")
runPopulationParameterEstimation()
runStandardErrorEstimation()

###############
# 2. Prepare all possible models in each directory (e.g. /back3 will contain all possible models with 3 biomarkers
###############

sapply(1:5, create_combtxt_flex)

###############
# 3. Backward selection procedure
###############


# Step with 5 biomarkers and next

# create data, model txt and mlxtran files given the results of the backward with 6 biomarkers 
iter_back_flex(5)
# fit the model with 5 biomarkers
fit_all_models_flex(5)

# create data, model txt and mlxtran files given the results of the backward with 5 biomarkers
iter_back_flex(4)
# fit the model with 4 biomarkers
fit_all_models_flex(4)

# create data, model txt and mlxtran files given the results of the backward with 4 biomarkers
iter_back_flex(3)
# fit the model with 3 biomarkers
fit_all_models_flex(3)

# create data, model txt and mlxtran files given the results of the backward with 3 biomarkers
iter_back_flex(2)
# fit the model with 2 biomarkers
fit_all_models_flex(2)

iter_back_flex(1)