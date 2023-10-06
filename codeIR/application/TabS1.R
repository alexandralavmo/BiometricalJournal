# File to summarized the 59 biomarkers in Table S1

#############
# 0. to load required packages and functions used in the analysis 
source("packages.R")
source("functions.R")
#############

#############
# 1. Loading the data sets (biomarkers observations and biomarker categories)
#############
data = read.table("datas/alldata.txt", header = T)
cat = read.table("datas/category.txt", header = T)

#############
# 2. Creating subtables for each category
#############
sapply(1:8, function (i) assign(paste0("cat", i), data.frame(biomarqueur = NA, unit = NA, N = NA, n = NA), envir = .GlobalEnv))  

for (i in 1:59){
  catg = cat$category[cat$num == i] 
  trans = get(paste0("cat", catg))[1, ]
  trans$biomarqueur = i
  trans$N = length(unique(data$id[data$ytype == i]))
  trans$n = nrow(data[data$ytype == i, ])/length(unique(data$id[data$ytype == i]))
  assign(paste0("cat", catg), rbind(get(paste0("cat", catg)), trans))
}


cat1[1, ] = c("Category 1", "", "", "")
cat2[1, ] = c("Category 2", "", "", "")
cat3[1, ] = c("Category 3", "", "", "")
cat4[1, ] = c("Category 4", "", "", "")
cat5[1, ] = c("Category 5", "", "", "")
cat6[1, ] = c("Category 6", "", "", "")
cat7[1, ] = c("Category 7", "", "", "")
cat8[1, ] = c("Category 8", "", "", "")


#############
# 2. Merging the subtables and saving results
#############
tabS1 = rbind(cat1, cat2, cat3, cat4, cat5, cat6, cat7, cat8)
save(tabS1, file = "../results/TabS1.RData")
