#### File used to fit the univariate nonlinear joint model for biomaker not selected with a linear model

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

# median of biomarker observations
med = data.frame(num = 1:59, med = sapply(1:59, function(x) median(data$obs[data$ytype == x])))

### most stringent thresholds to fit nonlinear biomarkers
thrsld_lin = 0.15

# get the vector of biomarker selected with a linear joint model
res = sapply(seq(1, 59), init_lin, thrsld_lin = thrsld_lin)
# nonlinear biomarkers saved in a vector 
nlin_mod = which(res == F)

# creating nonlinear joint model files for the nonlinear biomarkers 
sapply(nlin_mod, create_model_nlin)

# estimating the nonlinear joint models 
sapply(nlin_mod, estimate_nlin)



