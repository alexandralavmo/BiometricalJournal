#### File used to fit the univariate linear joint model for all biomarkers 

###############
# 0. Loading data file 
###############
data = read.table("datas/alldata.txt", header = T)
# list of markers having proportional error
list_prop = c(1, 3, 5, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 22, 23, 25, 27, 28, 29, 31, 32, 33, 34, 37, 38, 39, 40, 41, 42, 45, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59)
# list of markers having constant error
list_cst = c(2, 4, 6, 7, 14, 15, 21, 24, 26, 30, 35, 43, 44, 46, 52)
# list of markers having combined error
list_cmb = c(36)

###############
# 1. Creating the 59 joint model files with linear modeling
###############
sapply(1:59, create_model_lin)

###############
# 2. Estimating the 59 joint models with linear modeling 
###############
sapply(1:59, estimate_lin) 
