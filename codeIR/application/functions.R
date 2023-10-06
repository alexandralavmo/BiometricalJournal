##########
#1. Univariate process - linear modeling
##########

# Function to create all linear models and data files to estimate a joint model on Monolix software involving a biomarker i
# saving data, mlxtran and model txt files 
create_model_lin = function(i){
  # First filter data file in order to keep only three biomarkers
  datab = data[data$ytype %in% c(i, 60, 61), ]
  d = datab
  d$ytype[d$ytype == i] = 1  
  write.table(d, file = paste0("intermediate_results/uni_fits/datas/data", i, ".txt"), row.names = F)
  
  # Then, modify the txt model file for linear modeling
  med = median(datab$obs[datab$ytype == i])
  obj = readLines("modfiles/lin_uni.txt")
  obj[10] = paste0("haz1 = h1 * exp( alpha1*(m-", med, ")+hcov1)")
  obj[11] = paste0("haz2 = h2 * exp( alpha2*(m-", med, ")+hcov2)")
  writeLines(obj, paste0("intermediate_results/uni_fits/lin/mod_", i, ".txt"))
  
  # For nonlinear modeling
  #obj = readLines("modfiles/nonlin_uni.txt")
  #obj[10] = paste0("haz1 = h1 * exp( alpha1*(m-", med, ")+hcov1)")
  #obj[11] = paste0("haz2 = h2 * exp( alpha2*(m-", med, ")+hcov2)")
  #writeLines(obj, paste0("intermediate_results/uni_fits/nonlin/mod_", i, ".txt"))
  
  # Finally, modify the mlxtran parameter file for linear modeling
  # 3 possibilities : constant, proportional or combined error
  # Constant error
  if (i %in% list_cst){
  obj = readLines("modfiles/lin_cst.mlxtran")
  obj[4] = paste0("file = '../datas/data", i, ".txt'")
  obj[36] = paste0("file = 'mod_", i, ".txt'")
  obj[71] = paste0("exportpath = 'mod", i, "'")
  write.table(obj, file = paste0("intermediate_results/uni_fits/lin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
  else if (i%in% list_prop){
    obj = readLines("modfiles/lin_prop.mlxtran")
    obj[4] = paste0("file = '../datas/data", i, ".txt'")
    obj[36] = paste0("file = 'mod_", i, ".txt'")
    obj[71] = paste0("exportpath = 'mod", i, "'")
    write.table(obj, file = paste0("intermediate_results/uni_fits/lin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
  else {
    obj = readLines("modfiles/lin_comb.mlxtran")
    obj[4] = paste0("file = '../datas/data", i, ".txt'")
    obj[36] = paste0("file = 'mod_", i, ".txt'")
    obj[72] = paste0("exportpath = 'mod", i, "'")
    write.table(obj, file = paste0("intermediate_results/uni_fits/lin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
}

# function to estimate a linear joint model involving biomarker i 
# model fit saved in a specific folder
estimate_lin = function(i){
  loadProject(paste0("intermediate_results/uni_fits/lin/mod", i, ".mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}

# function to indicate if a biomarker i is selected for the linear modeling 
# that is to say, if biomarker i has its measurement error < thrsld_lin 
init_lin = function(i, thrsld_lin){
  t = read.table(paste0("intermediate_results/uni_fits/lin/mod", i, "/populationParameters.txt"), header = T, sep = ',')
  names_par = t$parameter
  res_err = ifelse("e_a" %in% names_par & "e_b" %in% names_par, t$value[t$parameter == "e_b"]+t$value[t$parameter == "e_a"]/med$med[med$num == i], ifelse("e_a" %in% names_par, t$value[t$parameter == "e_a"]/med$med[med$num == i], t$value[t$parameter == "e_b"])) 
  abs(res_err) <= thrsld_lin
}

##########
#2. Univariate process - nonlinear modeling
##########

# Function to create all nonlinear models and data files to estimate a joint model on Monolix software involving a biomarker i
# saving data, mlxtran and model txt files
create_model_nlin = function(i){
  # First filter data file in order to keep only the biomarker i
  datab = data[data$ytype %in% c(i, 60, 61), ]
  
  # Then, modify the txt model file for linear modeling
  med = median(datab$obs[datab$ytype == i])
  obj = readLines("modfiles/nonlin_uni.txt")
  obj[10] = paste0("haz1 = h1 * exp( alpha1*(m-", med, ")+hcov1)")
  obj[11] = paste0("haz2 = h2 * exp( alpha2*(m-", med, ")+hcov2)")
  if (med<0){
    obj[10] = paste0("haz1 = h1 * exp( alpha1*(m+", abs(med), ")+hcov1)")
    obj[11] = paste0("haz2 = h2 * exp( alpha2*(m+", abs(med), ")+hcov2)")
  }
  writeLines(obj, paste0("intermediate_results/uni_fits/nonlin/mod_", i, ".txt"))
  
  # Finally, modify the mlxtran parameter file for linear modeling
  # 3 possibilities : constant, proportional or combined error
  # Constant error
  if (i %in% list_cst){
    obj = readLines("modfiles/nonlin_cst.mlxtran")
    obj[4] = paste0("file = '../datas/data", i, ".txt'")
    obj[37] = paste0("file = 'mod_", i, ".txt'")
    obj[76] = paste0("exportpath = 'mod", i, "'")
    write.table(obj, file = paste0("intermediate_results/uni_fits/nonlin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
  else if (i%in% list_prop){
    obj = readLines("modfiles/nonlin_prop.mlxtran")
    obj[4] = paste0("file = '../datas/data", i, ".txt'")
    obj[37] = paste0("file = 'mod_", i, ".txt'")
    obj[76] = paste0("exportpath = 'mod", i, "'")
    write.table(obj, file = paste0("intermediate_results/uni_fits/nonlin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
  else {
    obj = readLines("modfiles/nonlin_comb.mlxtran")
    obj[4] = paste0("file = '../datas/data", i, ".txt'")
    obj[37] = paste0("file = 'mod_", i, ".txt'")
    obj[77] = paste0("exportpath = 'mod", i, "'")
    write.table(obj, file = paste0("intermediate_results/uni_fits/nonlin/mod", i, ".mlxtran"), quote = F, row.names = F, col.names = F)
  }
}

# function to estimate a nonlinear joint model involving biomarker i 
# model fit saved in a specific folder
estimate_nlin = function(i){
  loadProject(paste0("intermediate_results/uni_fits/nonlin/mod", i, ".mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}


##########
#3. Results of univariate process - Flowcharts
##########

# Function that gives the flowchart for univariate analysis depending on: 
# - the threshold of the residual error for linear biomarkers (thrsld_lin)
# - the threshold of the residual error for nonlinear biomarkers (thrsld_nlin)
# - the threshold of the RSE of parameters (thrsld_rse)
# returns a .txt file 

give_flowchart = function(thrsld_lin, thrsld_nlin, thrsld_rse){
  init_lin = function(i, thrsld_lin){
    t = read.table(paste0("intermediate_results/uni_fits/lin/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    names_par = t$parameter
    res_err = ifelse("e_a" %in% names_par & "e_b" %in% names_par, t$value[t$parameter == "e_b"]+t$value[t$parameter == "e_a"]/med$med[med$num == i], ifelse("e_a" %in% names_par, t$value[t$parameter == "e_a"]/med$med[med$num == i], t$value[t$parameter == "e_b"])) 
    abs(res_err) <= thrsld_lin
  }
  
  res = sapply(1:59, init_lin, thrsld_lin = thrsld_lin)
  lin_mod = which(res == T)
  nlin_mod = which(res == F)
  
  ### first step of the flowchart (selection of linear biomarkers or tests for nonlinear modeling)
  cat("Number of biomarkers initially modeled with a linear model: ", length(lin_mod), "\n")
  cat("Number of biomarkers tested for nonlinear modeling: ", length(nlin_mod), "\n\n")
  
  ### for linear biomarkers ###
  
  excl_longi = c()
  excl_surv = c()
  excl_link = c()
  selec_lin = c()
  for (i in lin_mod){
    t = read.table(paste0("intermediate_results/uni_fits/lin/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    names_longi = c("b0_pop", "b1_pop", "omega_b0", "omega_b1", "e_a", "e_b")
    names_survie = c("h1_pop", "h2_pop", "alpha1_pop", "alpha2_pop", "beta_hcov1_score4C", "beta_hcov2_score4C")
    longi_rse = t$rse_sa[t$parameter %in% names_longi]
    surv_rse = t$rse_sa[t$parameter %in% names_survie]
    wald_link = t$rse_sa[t$parameter == "alpha1_pop"]
    if (length(which(longi_rse>thrsld_rse))>0) excl_longi = c(excl_longi, i)
    else if (length(which(surv_rse>thrsld_rse))>0) excl_surv = c(excl_surv, i)
    else if (length(which(wald_link>50))>0) excl_link = c(excl_link, i)
    else selec_lin = c(selec_lin, i)
  }
  
  cat("Number of linear biomarkers excluded because at least one longitudinal RSE is over ", thrsld_rse, "%: ", length(excl_longi), "\n")
  cat("Number of linear biomarkers excluded because at least one survival RSE is over ", thrsld_rse, "%: ", length(excl_surv), "\n")
  cat("Number of linear biomarkers excluded because the Wald test on the link coefficient is not significant: ", length(excl_link), "\n")
  cat("Number of biomarkers selected with a linear model: ", length(selec_lin), "\n\n")
  
  
  ### for nonlinear biomarkers ###
  
  nl_excl_err = c()
  nl_excl_longi = c()
  nl_excl_surv = c()
  nl_excl_link = c()
  selec_nlin = c()
  for (i in nlin_mod){
    t = read.table(paste0("intermediate_results/uni_fits/nonlin/mod", i, "/populationParameters.txt"), header = T, sep = ',')
    if (ncol(t)<4) t$se_sa = NA; t$rse_sa = NA
    t$rse_sa[is.na(t$rse_sa) == T & is.na(t$value) == F] = Inf
    names_par = t$parameter
    names_longi = c("b0_pop", "b1_pop", "b2_pop", "a_pop", "tlag_pop", "omega_b0", "omega_b1", "omega_b2", "omega_a", "omega_tlag", "e_a", "e_b")
    names_survie = c("h1_pop", "h2_pop", "alpha1_pop", "alpha2_pop", "beta_hcov1_score4C", "beta_hcov2_score4C")
    longi_rse = t$rse_sa[t$parameter %in% names_longi]
    surv_rse = t$rse_sa[t$parameter %in% names_survie]
    wald_link = t$rse_sa[t$parameter == "alpha1_pop"]
    res_err = ifelse("e_a" %in% names_par & "e_b" %in% names_par, t$value[t$parameter == "e_b"]+t$value[t$parameter == "e_a"]/med$med[med$num == i], ifelse("e_a" %in% names_par, t$value[t$parameter == "e_a"]/med$med[med$num == i], t$value[t$parameter == "e_b"])) 
    if (abs(res_err)>= thrsld_nlin) nl_excl_err = c(nl_excl_err, i)
    else if (length(which(longi_rse>thrsld_rse))>0) nl_excl_longi = c(nl_excl_longi, i)
    else if (length(which(surv_rse>thrsld_rse))>0) nl_excl_surv = c(nl_excl_surv, i)
    else if (length(which(wald_link>50))>0) nl_excl_link = c(nl_excl_link, i)
    else selec_nlin = c(selec_nlin, i)
  }
  
  cat("Number of nonlinear biomarkers excluded because the threshold for nonlinear residual error is over ", thrsld_nlin, ": ", length(nl_excl_err), "\n")
  cat("Number of nonlinear biomarkers excluded because at least one longitudinal RSE is over ", thrsld_rse, "%: ", length(nl_excl_longi), "\n")
  cat("Number of nonlinear biomarkers excluded because at least one survival RSE is over ", thrsld_rse, "%: ", length(nl_excl_surv), "\n")
  cat("Number of nonlinear biomarkers excluded because the Wald test on the link coefficient is not significant: ", length(nl_excl_link), "\n")
  cat("Number of biomarkers selected with a nonlinear model: ", length(selec_nlin), "\n\n\n")
  
  # save the results of selection
  if (thrsld_lin == 0.25) save(selec_lin, file = "intermediate_results/uni_fits/selec_lin.RData") ; save(selec_nlin, file = "intermediate_results/uni_fits/selec_nlin.RData")
  if (thrsld_lin == 0.15) save(selec_lin, file = "intermediate_results/uni_fits/selec_lin_stringent.RData") ; save(selec_nlin, file = "intermediate_results/uni_fits/selec_nlin_stringent.RData")
  if (thrsld_lin == 0.35) save(selec_lin, file = "intermediate_results/uni_fits/selec_lin_flex.RData") ; save(selec_nlin, file = "intermediate_results/uni_fits/selec_nlin_flex.RData")  
}

##########
#4. Multivariate processes (current, stringent and flexible scenario)
##########

##### For current scenario #####
# function to create model text files with a specified list of biomarkers passed as argument
# create_modtxt(i) creates a txt model including biomarkers i (i can be a vector)

create_modtxt_curr = function(i){  
  num = setdiff(c("1", "2", "3"), i)
  ref = readLines("modfiles/mod_123.txt")
  warn = 0
  l4 = unlist(strsplit(ref[4], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(num)){
    indb = c(grep(paste0("b0", num[j]), l4), grep(paste0("b1", num[j]), l4), grep(paste0("b2", num[j]), l4), grep(paste0("a", num[j], "_pop"), l4), grep(paste0("alpha1", num[j]), l4), grep(paste0("alpha2", num[j]), l4), grep(paste0("a", num[j], ", "), l4))
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l4)-1) warn = 1
  l4 = l4[-c(ind)]
  
  if (warn == 1) {der = paste0(str_sub(rev(l4)[2], 1, 7)) ; l4[c(length(l4)-1)] = der}
  ref[4] = paste(l4, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 8:10){
    for (j in num){
      l = unlist(strsplit(ref[ligne], split = " "))
      if (l[1] == paste0("m", j)) ind_ligne = c(ind_ligne, ligne)  
    }
  }
  
  ind = c()
  indb = c()
  l12 = unlist(strsplit(ref[12], split = " "))
  for (j in 1:length(num)){
    indb = grep(paste0("alpha1", num[j]), l12)
    ind = c(ind, indb)
  }
  l12 = l12[-c(ind)]
  ref[12] = paste(l12, collapse = " ")
  
  l13 = unlist(strsplit(ref[13], split = " "))
  l13 = l13[-c(ind)]
  ref[13] = paste(l13, collapse = " ")
  
  ind = c()
  indb = c()
  l21 = unlist(strsplit(ref[21], split = " "))
  for (j in 1:length(num)){
    indb = grep(num[j], l21)
    ind = c(ind, indb)
  }
  l21 = l21[-c(ind)]
  ref[21] = paste(l21, collapse = " ")
  
  ref = ref[-c(ind_ligne)]
  return(ref)
}


# function create to all possible model text files with k biomarkers
# and write them into the appropriate directory
create_combtxt_curr = function(k){
  create_combination = function(col) {
    txt = create_modtxt_curr(as.character(col))
    path = paste0("intermediate_results/multi_fits/current/back", k, "/mod_", paste(col, collapse = ''), ".txt")
    writeLines(txt, path)
  }
  apply(combn(1:3, k), 2, create_combination)
}


# function to create mlxtran files given a list of biomarkers 
# create_mlx_curr(i) creates a mlxtran file removing biomarkers i (i can be a vector)
create_mlx_curr = function(i){
  mod_ind = 1:3
  mod_ind = paste0(mod_ind[-c(as.integer(i))], collapse = "")
  
  ref3 = readLines("modfiles/mod123.mlxtran")
  ref3[44] = paste0("file = 'mod_", mod_ind, ".txt'")# Set correct model file
  
  l11 = unlist(strsplit(ref3[11], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    ind = grep(i[j], l11)
    indb = c(indb, ind)
  }
  ind = c(indb, which(l11 == "continuous,")[1:length(i)])
  l11 = l11[-c(ind)]
  ref3[11] = paste(l11, collapse = " ") 
  
  l21 = unlist(strsplit(ref3[21], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = c(grep(paste0("b0", i[j]), l21), grep(paste0("b1", i[j]), l21), grep(paste0("b2", i[j]), l21), grep(paste0("a", i[j], "_pop"), l21), grep(paste0("omega_a", i[j]), l21), grep(paste0("alpha1", i[j]), l21), grep(paste0("alpha2", i[j]), l21))
    ind = c(ind, indb)
  }
  l21 = l21[-c(ind)]
  ref3[21] = paste(l21, collapse = " ")
  
  warn = 0
  l42 = unlist(strsplit(ref3[42], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l42)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l42)) warn = 1
  l42 = l42[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l42)[1], 1, 2), "}") ; l42[length(l42)] = der}
  ref3[42] = paste(l42, collapse = " ") 
  
  l52 = unlist(strsplit(ref3[52], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l52)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l52)) warn = 1
  l52 = l52[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l52)[1], 1, 3), "}") ; l52[length(l52)] = der}
  ref3[52] = paste(l52, collapse = " ") 
  
  l53 = unlist(strsplit(ref3[53], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l53)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l53)) warn = 1
  l53 = l53[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l53)[1], 1, 3), "}") ; l53[length(l53)] = der}
  ref3[53] = paste(l53, collapse = " ") 
  
  ind_ligne = c()
  for (ligne in 1:length(ref3)){
    l = unlist(strsplit(ref3[ligne], split = " "))
    for (j in 1:length(i)){
      if (length(grep(paste0("b0", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b1", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b2", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("alpha1", i[j], "_pop"), l)) != 0 | length(grep(paste0("alpha2", i[j], "_pop"), l)) != 0 | length(grep(paste0("b_", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("a_", i[j]), l, fixed = T)) != 0 | length(grep(paste0("omega_a", i[j]), l, fixed = T)) != 0 | length(grep(paste0("a", i[j], "_pop"), l, fixed = T)) != 0) ind_ligne = c(ind_ligne, ligne)
    }
  }
  ref3 = ref3[-c(ind_ligne)]
  return(ref3)
}

# function to create mlxtran files and data files for iteration k of the backward 
# load the results found at the previous iteration and create the new model and data files
iter_back_curr = function(k){
    tab = read.table(paste0("intermediate_results/multi_fits/current/back", k+1, "/mod/populationParameters.txt"), header = T, sep = ',')
    dat = read.table(paste0("intermediate_results/multi_fits/current/back", k+1, "/data.txt"), header = T, sep = ' ')
    stab = tab[c(which(tab$parameter %in% c("alpha11_pop", "alpha12_pop", "alpha13_pop"))), ]
    max = stab$parameter[stab$rse_sa == max(stab$rse_sa[is.nan(stab$rse_sa) == F])&is.nan(stab$rse_sa) == F]
    prev = substr(stab$parameter, 7, 7)
    num = substr(max, 7, 7)
    cur = prev[-which(prev == num)]
    numb = setdiff(c("1", "2", "3"), cur)
    if (tab$rse_sa[tab$parameter == max]>50){
      mlx = create_mlx_curr(numb)
      write.table(mlx, file = paste0("intermediate_results/multi_fits/current/back", k, "/mod.mlxtran"), quote = F, row.names = F, col.names = F)
      dat = dat[-c(which(dat$ytype == as.integer(num))), ]
      write.table(dat, file = paste0("intermediate_results/multi_fits/current/back", k, "/data.txt"), row.names = F)
    }
    else{
      indnan = which(is.nan(stab$rse_sa) == T)
      if(length(indnan)>0){
        min = stab$parameter[stab$value == min(stab$value[indnan])]
        prev = substr(stab$parameter, 7, 7)
        num = substr(min, 7, 7)
        cur = prev[-which(prev == num)]
        numb = setdiff(c("1", "2", "3"), cur)
        mlx = create_mlx_curr(numb)
        write.table(mlx, file = paste0("intermediate_results/multi_fits/current/back", k, "/mod.mlxtran"), quote = F, row.names = F, col.names = F)
        dat = dat[-c(which(dat$ytype == as.integer(num))), ]
        write.table(dat, file = paste0("intermediate_results/multi_fits/current/back", k, "/data.txt"), row.names = F)
      }
    }
}

# function to fit the model involving k biomarkers  
fit_all_models_curr = function(k){
  loadProject(paste0("intermediate_results/multi_fits/current/back", k, "/mod.mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}

##### For flexible scenario #####

# function to create model text files with a specified list of biomarkers passed as argument
# create_modtxt(i) creates a txt model including biomarkers i (i can be a vector)
create_modtxt_flex = function(i){
  
  num = setdiff(c("1", "2", "3", "4", "5", "6"), i)
  ref = readLines("modfiles/mod_123456.txt")
  warn = 0
  l4 = unlist(strsplit(ref[4], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(num)){
    indb = c(grep(paste0("b0", num[j]), l4), grep(paste0("b1", num[j]), l4), grep(paste0("b2", num[j]), l4), grep(paste0("a", num[j], "_pop"), l4), grep(paste0("alpha1", num[j]), l4), grep(paste0("alpha2", num[j]), l4), grep(paste0("a", num[j], ", "), l4))
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l4)-1) warn = 1
  l4 = l4[-c(ind)]
  
  if (warn == 1) {der = paste0(str_sub(rev(l4)[2], 1, 7)) ; l4[c(length(l4)-1)] = der}
  ref[4] = paste(l4, collapse = " ")
  
  ind_ligne = c()
  for (ligne in 8:13){
    for (j in num){
      l = unlist(strsplit(ref[ligne], split = " "))
      if (l[1] == paste0("m", j)) ind_ligne = c(ind_ligne, ligne)  
    }
  }
  
  ind = c()
  indb = c()
  l15 = unlist(strsplit(ref[15], split = " "))
  for (j in 1:length(num)){
    indb = grep(paste0("alpha1", num[j]), l15)
    ind = c(ind, indb)
  }
  l15 = l15[-c(ind)]
  ref[15] = paste(l15, collapse = " ")
  
  l16 = unlist(strsplit(ref[16], split = " "))
  l16 = l16[-c(ind)]
  ref[16] = paste(l16, collapse = " ")
  
  ind = c()
  indb = c()
  l24 = unlist(strsplit(ref[24], split = " "))
  for (j in 1:length(num)){
    indb = grep(num[j], l24)
    ind = c(ind, indb)
  }
  l24 = l24[-c(ind)]
  ref[24] = paste(l24, collapse = " ")
  
  ref = ref[-c(ind_ligne)]
  return(ref)
}


# function to create all possible model text files with k biomarkers
# and write them into the appropriate directory
create_combtxt_flex = function(k){
  create_combination = function(col) {
    txt = create_modtxt_flex(as.character(col))
    path = paste0("intermediate_results/multi_fits/flexible/back", k, "/mod_", paste(col, collapse = ''), ".txt")
    writeLines(txt, path)
  }
  apply(combn(1:6, k), 2, create_combination)
}

# function to create mlxtran files given a list of biomarkers 
# create_mlx_flex(i) creates a mlxtran files removing biomarkers i 
create_mlx_flex = function(i){
  mod_ind = 1:6
  mod_ind = paste0(mod_ind[-c(as.integer(i))], collapse = "")
  
  ref6 = readLines("modfiles/mod123456.mlxtran")
  ref6[56] = paste0("file = 'mod_", mod_ind, ".txt'")# Set correct model file
  
  l11 = unlist(strsplit(ref6[11], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    ind = grep(i[j], l11)
    indb = c(indb, ind)
  }
  ind = c(indb, which(l11 == "continuous,")[1:length(i)])
  l11 = l11[-c(ind)]
  ref6[11] = paste(l11, collapse = " ") 
  
  l21 = unlist(strsplit(ref6[21], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = c(grep(paste0("b0", i[j]), l21), grep(paste0("b1", i[j]), l21), grep(paste0("b2", i[j]), l21), grep(paste0("a", i[j], "_pop"), l21), grep(paste0("omega_a", i[j]), l21), grep(paste0("alpha1", i[j]), l21), grep(paste0("alpha2", i[j]), l21))
    ind = c(ind, indb)
  }
  l21 = l21[-c(ind)]
  ref6[21] = paste(l21, collapse = " ")
  
  warn = 0
  l54 = unlist(strsplit(ref6[54], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l54)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l54)) warn = 1
  l54 = l54[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l54)[1], 1, 2), "}") ; l54[length(l54)] = der}
  ref6[54] = paste(l54, collapse = " ") 
  
  l67 = unlist(strsplit(ref6[67], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l67)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l67)) warn = 1
  l67 = l67[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l67)[1], 1, 3), "}") ; l67[length(l67)] = der}
  ref6[67] = paste(l67, collapse = " ") 
  
  l68 = unlist(strsplit(ref6[68], split = " "))
  ind = c()
  indb = c()
  for (j in 1:length(i)){
    indb = grep(i[j], l68)
    ind = c(ind, indb)
  }
  if (rev(ind)[1] == length(l68)) warn = 1
  l68 = l68[-c(ind)]
  if (warn == 1) {der = paste0(str_sub(rev(l68)[1], 1, 3), "}") ; l68[length(l68)] = der}
  ref6[68] = paste(l68, collapse = " ") 
  
  ind_ligne = c()
  for (ligne in 1:length(ref6)){
    l = unlist(strsplit(ref6[ligne], split = " "))
    for (j in 1:length(i)){
      if (length(grep(paste0("b0", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b1", i[j]), l, fixed = T)) != 0 | length(grep(paste0("b2", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("alpha1", i[j], "_pop"), l)) != 0 | length(grep(paste0("alpha2", i[j], "_pop"), l)) != 0 | length(grep(paste0("b_", i[j]), l, fixed = T)) != 0 | 
          length(grep(paste0("a_", i[j]), l, fixed = T)) != 0 | length(grep(paste0("omega_a", i[j]), l, fixed = T)) != 0 | length(grep(paste0("a", i[j], "_pop"), l, fixed = T)) != 0) ind_ligne = c(ind_ligne, ligne)
    }
  }
  ref6 = ref6[-c(ind_ligne)]
  return(ref6)
}

# function to create mlxtran files and data files for iteration k of the backward 
# load the results found at the previous iteration and create the new model and data files
iter_back_flex = function(k){
    tab = read.table(paste0("intermediate_results/multi_fits/flexible/back", k+1, "/mod/populationParameters.txt"), header = T, sep = ',')
    dat = read.table(paste0("intermediate_results/multi_fits/flexible/back", k+1, "/data.txt"), header = T, sep = ' ')
    stab = tab[c(which(tab$parameter %in% c("alpha11_pop", "alpha12_pop", "alpha13_pop", "alpha14_pop", "alpha15_pop", "alpha16_pop"))), ]
    max = stab$parameter[stab$rse_sa == max(stab$rse_sa[is.nan(stab$rse_sa) == F])&is.nan(stab$rse_sa) == F]
    prev = substr(stab$parameter, 7, 7)
    num = substr(max, 7, 7)
    cur = prev[-which(prev == num)]
    numb = setdiff(c("1", "2", "3", "4", "5", "6"), cur)
    if (tab$rse_sa[tab$parameter == max]>50){
      mlx = create_mlx_flex(numb)
      write.table(mlx, file = paste0("intermediate_results/multi_fits/flexible/back", k, "/mod.mlxtran"), quote = F, row.names = F, col.names = F)
      dat = dat[-c(which(dat$ytype == as.integer(num))), ]
      write.table(dat, file = paste0("intermediate_results/multi_fits/flexible/back", k, "/data.txt"), row.names = F)
    }
    else{
      indnan = which(is.nan(stab$rse_sa) == T)
      if(length(indnan)>0){
        min = stab$parameter[stab$value == min(stab$value[indnan])]
        prev = substr(stab$parameter, 7, 7)
        num = substr(min, 7, 7)
        cur = prev[-which(prev == num)]
        numb = setdiff(c("1", "2", "3", "4", "5", "6"), cur)
        mlx = create_mlx_flex(numb)
        write.table(mlx, file = paste0("intermediate_results/multi_fits/flexible/back", k, "/mod.mlxtran"), quote = F, row.names = F, col.names = F)
        dat = dat[-c(which(dat$ytype == as.integer(num))), ]
        write.table(dat, file = paste0("intermediate_results/multi_fits/flexible/back", k, "/data.txt"), row.names = F)
      }
    }
}

# function to fit the model involving k biomarkers
fit_all_models_flex = function(k){
  loadProject(paste0("intermediate_results/multi_fits/flexible/back", k, "/mod.mlxtran"))
  runPopulationParameterEstimation()
  runStandardErrorEstimation()
}

##########
#5. Results of multivariate process under current scenario
##########

## a. Longitudinal evolution of selected biomarkers

# function that takes a biomarker number bm and returns a ggplot object 
# the ggplot object represents the longitudinal evolution of the biomarker for all patients from time 0 to time 30
plot_mark = function(bm){
  title = paste0("Biomarker ", bm)
  data_mark = data[data$ytype == bm, ]
  data_mark$status = 0
  for (i in surv$id){
    data_mark$status[data_mark$id == i] = surv$obs[surv$id == i]
  }
  data_mark$status = as.factor(data_mark$status)
  levels(data_mark$status) = c("censored", "dead", "discharged")

   gp = ggplot(data_mark, aes(x = time, y = obs, group = id, color = status))+geom_line(lwd = 0.3)+
    xlab("Time (days)")+ylab("")+scale_colour_manual(values = c("grey", "red", "black"))+
    theme_classic()+ggtitle(title)+ylim(min(data_mark$obs), max(data_mark$obs)) + facet_grid(.~status)+
    theme(plot.title = element_text(hjust = 0.5, size = 9), element_line(size = 0.5), axis.text = element_text(size = 9), 
          axis.title = element_text(size = 9), legend.position = "None")
  
  assign(paste0("gp", bm), gp, envir = parent.frame())
}

## b. Goodness of fit tools

# function giving median (with 95% percentiles) empirical individual weighted residuals for biomarker j and time i 
emp_res = function(i, j){
  data = get(paste0("iwres_bm", j))
  emp1 = quantile(data$iwres[data$time<(i+1) & data$time>= i], probs = 0.025)
  emp2 = quantile(data$iwres[data$time<(i+1) & data$time>= i], probs = 0.975)
  emp3 = quantile(data$iwres[data$time<(i+1) & data$time>= i], probs = 0.5)
  c(emp1, emp2, emp3)
}

# function creating a ggplot object representing the parameter "par" distribution 
dist = function(par){
  mn = min(parsim[, par])
  mx = max(parsim[, par])
  bin = (mx-mn)/30
  gp = ggplot(data = parsim, aes(parsim[, par]))+
    geom_histogram(aes(y = ..density..), binwidth = bin, fill = "azure", colour = "lightblue", lwd = 0.6)+
    theme_classic()+
    stat_function(fun = dnorm, colour = "black", args = list(mean = pop$value[pop$parameter == paste0(par, "_pop")], sd = pop$value[pop$parameter == paste0("omega_", par)]))+
    xlab(par)+
    theme(element_line(size = 0.4), axis.text = element_text(size = 5), 
          axis.title = element_text(size = 5))
  assign(paste0("gp_", par), gp, envir = parent.frame())
}

## c. Dynamic predictions 

# function to generate landmark dataset for landmark time lmk under current scenario
gen_landmark_curr = function(lmk) {
	dat = read.table(paste0("intermediate_results/multi_fits/current/back2/data.txt"), header = T)
	
	# Only keep patients still at risk at t = lmk
	ind = sapply(unique(dat$id), function(i) setdiff(dat$time[dat$id == i & dat$ytype %in% 4:5], c(0, 30))<lmk)
	ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
	dat = dat[dat$id %in% ind, ]
	
	# Remove observations after landmark (bayesian prediction at t = lmk)
	ind = (dat$ytype %in% 2:3) & (dat$time>lmk)
	dat = dat[!ind, ]
	
	# Censoring all patients at t = lmk
	dat$obs[dat$ytype %in% 4:5] = 0
	dat$time[(dat$ytype %in% 4:5) & (dat$time>0)] = lmk

	# Export dataset into the appropriate directory
	write.table(dat, file = paste0("intermediate_results/evalAUC/lmk", lmk, "/data_lmk.txt"), row.names = F)
}

# function to generate landmark dataset for landmark time lmk under flexible scenario
gen_landmark_flex = function(lmk) { 
  dat = read.table(paste0("intermediate_results/multi_fits/flexible/back2/data.txt"), header = T)
 
  # Only keep patients still at risk at t = lmk
  ind = sapply(unique(dat$id), function(i) setdiff(dat$time[dat$id == i & dat$ytype %in% 7:8], c(0, 30))<lmk)
  ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
  dat = dat[dat$id %in% ind, ]
  
  # Remove observations after landmark (bayesian prediction at t = lmk)
  ind = (dat$ytype %in% c(3, 6)) & (dat$time>lmk)
  dat = dat[!ind, ]
  
  # Censoring all patients at t = lmk
  dat$obs[dat$ytype %in% 7:8] = 0
  dat$time[(dat$ytype %in% 7:8) & (dat$time>0)] = lmk
  
  # Export dataset into the appropriate directory
  write.table(dat, file = "intermediate_results/compAUC/data_flexible.txt", row.names = F)
}

# function to generate landmark dataset for landmark time lmk under stringent scenario
gen_landmark_strin = function(lmk) {
  dat = read.table(paste0("intermediate_results/multi_fits/stringent/data7.txt"), header = T)
  
  # Only keep patients still at risk at t = lmk
  ind = sapply(unique(dat$id), function(i) setdiff(dat$time[dat$id == i & dat$ytype %in% 60:61], c(0, 30))<lmk)
  ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
  dat = dat[dat$id %in% ind, ]
  
  # Remove observations after landmark (bayesian prediction at t = lmk)
  ind = (dat$ytype %in% c(1)) & (dat$time>lmk)
  dat = dat[!ind, ]
  
  # Censoring all patients at t = lmk
  dat$obs[dat$ytype %in% 60:61] = 0
  dat$time[(dat$ytype %in% 60:61) & (dat$time>0)] = lmk
  
  # Export dataset into the appropriate directory
  write.table(dat, file = "intermediate_results/compAUC/data_stringent.txt", row.names = F)
}

# function to generate survival data in R format for patients that are in the landmark dataset for time = lmk
gen_survdata_curr = function(lmk) {
	dat = read.table(paste0("datas/sim_survival.txt"), header = T)
	
	# Only keep patients still at risk at t = lmk
	ind = sapply(unique(dat$id), function(i) setdiff(dat$time[dat$id == i], c(0, 30))<lmk)
	ind = which(sapply(ind, function(x) length(x) == 0 || !x))	
	dat = dat[dat$id %in% ind, ]
	
	# Export dataset into the appropriate directory
	write.table(dat, file = paste0("intermediate_results/evalAUC/lmk", lmk, "/data_surv.txt"), row.names = F)
}

# generates survival data in R format that are in the landmark dataset for time = lmk
gen_survdata_comp = function(lmk) {
  dat = read.table(paste0("datas/sim_survival.txt"), header = T)
  
  # Only keep patients still at risk at t = lmk
  ind = sapply(unique(dat$id), function(i) setdiff(dat$time[dat$id == i], c(0, 30))<lmk)
  ind = which(sapply(ind, function(x) length(x) == 0 || !x))	# Patient id keep
  dat = dat[dat$id %in% ind, ]
  
  # Export dataset into the appropriate directory
  write.table(dat, file = paste0("intermediate_results/compAUC/data_surv.txt"), row.names = F)
}


# function to simulate individual parameters from a landmark time lmk for current scenario
gen_indv_params_curr = function(lmk) {

	# First duplicate project files and directories
		
	path2 = "intermediate_results/multi_fits/current/back2"
	file.copy(paste0(path2, "/mod_23.txt"), paste0("intermediate_results/evalAUC/lmk", lmk, "/mod_23.txt"))
	file.copy(paste0(path2, "/mod.mlxtran"), paste0("intermediate_results/evalAUC/lmk", lmk, "/mod.mlxtran"))
	file.copy(paste0(path2, "/data.txt"), paste0("intermediate_results/evalAUC/lmk", lmk, "/data.txt"))
	file.copy(paste0(path2, "/mod"), paste0("intermediate_results/evalAUC/lmk", lmk), recursive = T, overwrite = T)

	# load project and fix all population parameters
	loadProject(paste0("intermediate_results/evalAUC/lmk", lmk, "/mod.mlxtran"))
	setInitialEstimatesToLastEstimates()
  	popparams = getPopulationParameterInformation()
  	popparams$method <- "FIXED"
  	setPopulationParameterInformation(popparams)

	# Replace data file
	BaseData = getData()
	setData(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_lmk.txt"), BaseData$headerTypes, BaseData$observationTypes)

	# Simulation under the conditionnal distribution
	runPopulationParameterEstimation() # does nothing since all parameters are fixed (but still mandatory for further steps)
	setConditionalDistributionSamplingSettings(nbsimulatedparameters = 100)
	runConditionalDistributionSampling()
	
	# Export dataset into the appropriate directory
	write.table(getSimulatedIndividualParameters(), file = paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), row.names = F, sep = "\t")
}

# function to simulate individual parameters from a landmark time lmk for the flexible scenario
gen_indv_params_flexible = function(lmk) {
  
  # First duplicate project files and directories
  
  file.copy("intermediate_results/multi_fits/flexible/back2/mod_36.txt", "intermediate_results/compAUC/mod_36.txt", overwrite = T)
  file.copy("intermediate_results/multi_fits/flexible/back2/mod.mlxtran", "intermediate_results/compAUC/mod.mlxtran", overwrite = T)
  file.copy("intermediate_results/multi_fits/flexible/back2/data.txt", "intermediate_results/compAUC/data.txt", overwrite = T)
  file.copy("intermediate_results/multi_fits/flexible/back2/mod", "intermediate_results/compAUC", recursive = T, overwrite = T)
  
  # load project and fix all population parameters
  loadProject(paste0("intermediate_results/compAUC/mod.mlxtran"))
  setInitialEstimatesToLastEstimates()
  popparams = getPopulationParameterInformation()
  popparams$method <- "FIXED"
  setPopulationParameterInformation(popparams)
  
  # Replace data file
  BaseData = getData()
  setData(paste0("intermediate_results/compAUC/data_flexible.txt"), BaseData$headerTypes, BaseData$observationTypes)
  
  # Simulation under the conditionnal distribution
  runPopulationParameterEstimation() # does nothing since all parameters are fixed (but still mandatory for further steps)
  setConditionalDistributionSamplingSettings(nbsimulatedparameters = 10)
  runConditionalDistributionSampling()
  
  # Export dataset into the appropriate directory
  write.table(getSimulatedIndividualParameters(), file = "intermediate_results/compAUC/data_params_flex.txt", row.names = F, sep = "\t")
}

# function to simulate individual parameters from a landmark time lmk for the stringent scenario
gen_indv_params_stringent = function(lmk) {
  
  # First duplicate project files and directories
  
  file.copy("intermediate_results/multi_fits/stringent/mod_7.txt", "intermediate_results/compAUC/mod_7.txt", overwrite = T)
  file.copy("intermediate_results/multi_fits/stringent/mod7.mlxtran", "intermediate_results/compAUC/mod7.mlxtran", overwrite = T)
  file.copy("intermediate_results/multi_fits/stringent/data7.txt", "intermediate_results/compAUC/data7.txt", overwrite = T)
  file.copy("intermediate_results/multi_fits/stringent/mod7", "intermediate_results/compAUC", recursive = T, overwrite = T)
  
  # load project and fix all population parameters
  loadProject(paste0("intermediate_results/compAUC/mod7.mlxtran"))
  setInitialEstimatesToLastEstimates()
  popparams = getPopulationParameterInformation()
  popparams$method <- "FIXED"
  setPopulationParameterInformation(popparams)
  
  # Replace data file
  BaseData = getData()
  setData(paste0("intermediate_results/compAUC/data_stringent.txt"), BaseData$headerTypes, BaseData$observationTypes)
  
  # Simulation under the conditionnal distribution
  runPopulationParameterEstimation() # does nothing since all parameters are fixed (but still mandatory for further steps)
  setConditionalDistributionSamplingSettings(nbsimulatedparameters = 10)
  runConditionalDistributionSampling()
  
  # Export dataset into the appropriate directory
  write.table(getSimulatedIndividualParameters(), file = "intermediate_results/compAUC/data_params_strin.txt", row.names = F, sep = "\t")
}

## d. ROC AUC for all scenario

# function for the computation of time-dependent AUC based on a landmark time lmk, a vector of horizon times hrz 
# for current joint model
getAUC = function(lmk, hrz){
	# Reading data files (survival table and simulated individual parameters)
	surv_table = read.table(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_surv.txt"), h = T)
	sim = read.csv(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), sep = "\t")

	calc_pred = function(x) {
		# Deriving individual biomarker trajectories and hazard functions for both risks 
		bm2 = function(tt) x["b02"]+x["b12"]*tt
		bm3 = function(tt) x["b03"]+x["b13"]*tt
		
		h1 = function(tt) {
			a = x["h1"]*exp(x["alpha12"]*(bm2(tt)-3.6236834))
			a = a*exp(x["alpha13"]*(bm3(tt)-6.2826893))
			a
		}
		h2 = function(tt) {
			a = x["h2"]*exp(x["alpha22"]*(bm2(tt)-3.6236834))
			a = a*exp(x["alpha23"]*(bm3(tt)-6.2826893))
			a
		}
		
		# Cumulative incidence functions
		F1 = function(tt) 1-exp(-(integrate(h1, 0, tt)$value))
		F2 = function(tt) 1-exp(-(integrate(h2, 0, tt)$value))
		F1 = Vectorize(F1)
		F2 = Vectorize(F2)
		
		# Compute survival probabilities for the desired horizon times
		if(is.nan(h1(1000)) | (h1(1000) == Inf) | (h1(1000)>10e+5)){
      		  if(is.nan(h2(1000)) | (h2(1000) == Inf )){
        	  1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      		  }
      		else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    		} else if(is.nan(h2(1000)) | (h2(1000) == Inf) | (h2(1000)>10e+5)){
      			1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    		} else{
      			1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    		}
  	}

  # Making individual predictions for each patient and each required time
  preds = data.frame(id = sim$id, rep = sim$rep, pred = apply(sim, 1, calc_pred))
  preds = aggregate(.~id, preds, mean)
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  # Returning ROC AUC for the required time horizons
  timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)$AUC_2[2] 
}

# function for the computation of time-dependent AUC based on a landmark time lmk, a vector of horizon times hrz 
# for baseline model
getAUC_b = function(lmk, hrz){
  # Reading data files (survival table and simulated individual parameters)
  surv_table = surv_table[surv_table$time>lmk, ]
  pop = read.table("intermediate_results/bsl_fit/mod_bsl/populationParameters.txt", header = T, sep = ',')
  
  beta1 = pop$value[pop$parameter == 'beta_hcov1_score4C']
  beta2 = pop$value[pop$parameter == 'beta_hcov2_score4C']
  h01 = pop$value[pop$parameter == "h1_pop"]
  h02 = pop$value[pop$parameter == "h2_pop"]
  
  calc_pred = function(x) {
    
    h1 = function(t) h01* exp(x["score4C"]*beta1)
    h2 = function(t) h02* exp(x["score4C"]*beta2)
    F1 = function(t) 1-exp(-(h1(t)*t))
    F2 = function(t) 1-exp(-(h2(t)*t))
    
    # Compute survival probabilities for the desired horizon times
    if(is.nan(h1(1000)) | (h1(1000) == Inf)){
      if(is.nan(h2(1000)) | (h2(1000) == Inf)){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      }
      else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    } else if(is.nan(h2(1000)) | (h2(1000) == Inf)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  
  # Making individual predictions for each patient and each required time
  
  preds = data.frame(id = surv_table$id, pred = apply(surv_table, 1, calc_pred))
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)$AUC_2[2]
}

# functions to computate time-dependent AUC for landmark = lmk and horizon = hrz
# returns timeROC object 
# flexible scenario

getAUC_flex = function(lmk, hrz){
  # Reading data files (survival table and simulated individual parameters)
  surv_table = read.table(paste0("intermediate_results/compAUC/data_surv.txt"), h = T)
  sim = read.table("intermediate_results/compAUC/data_params_flex.txt", header = T)
  
  calc_pred = function(x) {
    # Deriving individual biomarker trajectories and hazard functions for both risks 
    bm3 = function(tt) x["b03"]+x["b13"]*tt
    bm6 = function(tt) x["b06"]+x["b16"]*tt
    
    h1 = function(tt) {
      a = x["h1"]*exp( x["alpha13"]*(bm3(tt)-3.6236834) )
      a = a*exp(x["alpha16"]*(bm6(tt)-95.60432)+x["hcov1"])
      a
    }
    h2 = function(tt) {
      a = x["h2"]*exp( x["alpha26"]*(bm3(tt)-3.6236834) )
      a = a*exp(x["alpha26"]*(bm6(tt)-95.60432)+x["hcov2"])
      a
    }
    
    # Cumulative incidence functions
    F1 = function(tt) 1-exp(-(integrate(h1, 0, tt)$value))
    F2 = function(tt) 1-exp(-(integrate(h2, 0, tt)$value))
    F1 = Vectorize(F1)
    F2 = Vectorize(F2)
    
    # Compute survival probabilities for the desired horizon times
    if(is.nan(h1(1000)) | (h1(1000) == Inf) | (h1(1000)>10e+5)){
      if(is.nan(h2(1000)) | (h2(1000) == Inf )){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      }
      else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    } else if(is.nan(h2(1000)) | (h2(1000) == Inf) | (h2(1000)>10e+5)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  # Making individual predictions for each patient and each required time
  preds = data.frame(id = sim$id, rep = sim$rep, pred = apply(sim, 1, calc_pred))
  preds = aggregate(.~id, preds, mean)
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  # Returning ROC AUC for the required time horizons
  roc_flex = timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)
  roc_flex
}

# functions to computate time-dependent AUC for landmark = lmk and horizon = hrz
# returns timeROC object 
# stringent scenario
getAUC_strin = function(lmk, hrz){
  # Reading data files (survival table and simulated individual parameters)
  surv_table = read.table(paste0("intermediate_results/compAUC/data_surv.txt"), h = T)
  sim = read.csv(paste0("intermediate_results/compAUC/data_params_strin.txt"), sep = "\t")
  
  calc_pred = function(x) {
    # Deriving individual biomarker trajectories and hazard functions for both risks 
    bm1 = function(tt) x["b0"]+x["b1"]*tt
    
    h1 = function(tt) {
      a = x["h1"]*exp(x["alpha1"]*(bm1(tt)-95.98391)+x["hcov1"])
      a
    }
    h2 = function(tt) {
      a = x["h2"]*exp(x["alpha2"]*(bm1(tt)-95.98391)+x["hcov2"])
      a
    }
    
    # Cumulative incidence functions
    F1 = function(tt) 1-exp(-(integrate(h1, 0, tt)$value))
    F2 = function(tt) 1-exp(-(integrate(h2, 0, tt)$value))
    F1 = Vectorize(F1)
    F2 = Vectorize(F2)
    
    # Compute survival probabilities for the desired horizon times
    if(is.nan(h1(1000)) | (h1(1000) == Inf) | (h1(1000)>10e+5)){
      if(is.nan(h2(1000)) | (h2(1000) == Inf )){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      }
      else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    } else if(is.nan(h2(1000)) | (h2(1000) == Inf) | (h2(1000)>10e+5)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  # Making individual predictions for each patient and each required time
  preds = data.frame(id = sim$id, rep = sim$rep, pred = apply(sim, 1, calc_pred))
  preds = aggregate(.~id, preds, mean)
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  # Returning ROC AUC for the required time horizons
  roc_strin = timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)
  roc_strin
}

# functions to computate time-dependent AUC for landmark = lmk and horizon = hrz
# returns timeROC object 
# current scenario
getAUC_cur = function(lmk, hrz){
  # Reading data files (survival table and simulated individual parameters)
  surv_table = read.table(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_surv.txt"), h = T)
  sim = read.csv(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), sep = "\t")
  
  calc_pred = function(x) {
    # Deriving individual biomarker trajectories and hazard functions for both risks 
    bm2 = function(tt) x["b02"]+x["b12"]*tt
    bm3 = function(tt) x["b03"]+x["b13"]*tt
    
    h1 = function(tt) {
      a = x["h1"]*exp(x["alpha12"]*(bm2(tt)-3.6236834))
      a = a*exp(x["alpha13"]*(bm3(tt)-6.2826893))
      a
    }
    h2 = function(tt) {
      a = x["h2"]*exp(x["alpha22"]*(bm2(tt)-3.6236834))
      a = a*exp(x["alpha23"]*(bm3(tt)-6.2826893))
      a
    }
    
    # Cumulative incidence functions
    F1 = function(tt) 1-exp(-(integrate(h1, 0, tt)$value))
    F2 = function(tt) 1-exp(-(integrate(h2, 0, tt)$value))
    F1 = Vectorize(F1)
    F2 = Vectorize(F2)
    
    # Compute survival probabilities for the desired horizon times
    if(is.nan(h1(1000)) | (h1(1000) == Inf) | (h1(1000)>10e+5)){
      if(is.nan(h2(1000)) | (h2(1000) == Inf )){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      }
      else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    } else if(is.nan(h2(1000)) | (h2(1000) == Inf) | (h2(1000)>10e+5)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  # Making individual predictions for each patient and each required time
  preds = data.frame(id = sim$id, rep = sim$rep, pred = apply(sim, 1, calc_pred))
  preds = aggregate(.~id, preds, mean)
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  # Returning ROC AUC for the required time horizons
  roc_cur = timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)
  roc_cur
}

# functions to computate time-dependent AUC for landmark = lmk and horizon = hrz
# returns timeROC object 
# baseline model 
getAUC_bsl = function(lmk, hrz){
  # Reading data files (survival table and simulated individual parameters)
  surv_table = surv_table[surv_table$time>lmk, ]
  pop = read.table("intermediate_results/bsl_fit/mod_bsl/populationParameters.txt", header = T, sep = ',')
  
  beta1 = pop$value[pop$parameter == 'beta_hcov1_score4C']
  beta2 = pop$value[pop$parameter == 'beta_hcov2_score4C']
  h01 = pop$value[pop$parameter == "h1_pop"]
  h02 = pop$value[pop$parameter == "h2_pop"]
  
  calc_pred = function(x) {
    
    h1 = function(t) h01* exp(x["score4C"]*beta1)
    h2 = function(t) h02* exp(x["score4C"]*beta2)
    F1 = function(t) 1-exp(-(h1(t)*t))
    F2 = function(t) 1-exp(-(h2(t)*t))
    
    # Compute survival probabilities for the desired horizon times
    if(is.nan(h1(1000)) | (h1(1000) == Inf)){
      if(is.nan(h2(1000)) | (h2(1000) == Inf)){
        1-(1+1-F1(hrz)-F2(lmk))/(1+1-F1(lmk)-F2(lmk))
      }
      else {1-(1+F2(1000)-F1(hrz)-F2(lmk))/(1+F2(1000)-F1(lmk)-F2(lmk)) }
    } else if(is.nan(h2(1000)) | (h2(1000) == Inf)){
      1-(F1(1000)+1-F1(hrz)-F2(lmk))/(F1(1000)+1-F1(lmk)-F2(lmk))  
    } else{
      1-(F1(1000)+F2(1000)-F1(hrz)-F2(lmk))/(F1(1000)+F2(1000)-F1(lmk)-F2(lmk)) 
    }
  }
  
  
  # Making individual predictions for each patient and each required time
  
  preds = data.frame(id = surv_table$id, pred = apply(surv_table, 1, calc_pred))
  
  # Function actually computing AUC from risk individual prediction and the corresponding time horizon
  roc_bsl = timeROC(T = surv_table$time, delta = surv_table$obs, marker = preds$pred, cause = 1, times = hrz, iid = T)
  roc_bsl
}


# function to create a ggplot object representing the dynamic predictions with 95% prediction interval referring to "biomarker 24" (say that it represents neutrophil values)
# for patient id: for time 0 to time lmk, the observed "biomarker 24" values are used, from lmk to time 30 predictions values are used
# returns a ggplot object 

landmark_neutro = function(lmk, id){
  
  t<-0:30
  
  # Load all the simulated individual parameters
  simul = read.table(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), header = T)
  dataobs = dataobs_all[dataobs_all$id == id & dataobs_all$ytype == 24, ]
  # create a indicator variable (1 if time > landmark time and 0 otherwise)
  dataobs$ind = ifelse(dataobs$time<lmk, 0, 1)
  dataobs$ind = factor(dataobs$ind)
  
  # keep individual parameters referring to patient id
  sim<-simul[simul$id == id, ]
  # data frame neutro contains predictions of neutrophil values from time 0 to time 30
  # using the simulated individual parameters
  neutro = data.frame(0:30)
  neutro = as.data.frame(t(neutro))
  neutro = neutro[-1, ]
  for (l in 1:nrow(sim)){
    for (j in 1:31){
      neutro[l, j] = sim$b02[l]+sim$b12[l]*t[j]
    }
  }
  # to build the 95% prediction interval 
  inf_n = apply(neutro, 2, quantile, probs = 0.025)
  sup_n = apply(neutro, 2, quantile, probs = 0.975)
  med_n = apply(neutro, 2, median)
  
  #create the landmark data frame with the neutrophil observations from time 0 to landmmark time, and with neutrophil predictions from landmark time +1 to D30 
  dataland = data.frame(time = 0:30, neut = med_n, low_n = inf_n, upp_n = sup_n, ind = as.factor(c(rep(0, lmk), rep(1, (30-lmk+1)))))
  
  gp<-ggplot(data = dataland, aes(x = time, y = neut))+
    geom_ribbon(aes(x = time, ymin = low_n, ymax = upp_n, group = ind, fill = ind), show.legend = F)+
    scale_fill_manual(values = c("white", "grey"))+
    scale_x_continuous(name = "Time (days)", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
    geom_line(lwd = 0.35, col = '#99FF66')+
    scale_y_continuous(name = bquote("Neutrophils (x"*10^9*"/L)"))+theme_classic()+
    coord_cartesian(ylim = c(0, 20))+
    geom_vline(xintercept = lmk, color = "red", lwd = 0.3)+
    theme(legend.position = "bottom", element_line(size = 0.5), axis.text = element_text(size = 5), axis.title = element_text(size = 5), title = element_text(size = 6))+
    geom_point(data = dataobs, aes(x = time, y = obs, group = ind, col = ind), lwd = 0.5, pch = 16, show.legend = F)+
    scale_color_manual(values = c("black", "#0099FF"))+
    ggtitle("")
  gp
}

# function to create a ggplot object representing the dynamic predictions with 95% prediction interval referring to "biomarker 43" (say that it represents crp values)
# for patient id: for time 0 to time lmk, the observed "biomarker 43" values are used, from lmk to time 30 predictions values are used
# returns a ggplot object 

landmark_crp<-function(lmk, id){
  t = 0:30
  
  ## Load data sets 
  simul = read.table(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), header = T)
  dataobs = dataobs_all[dataobs_all$id == id & dataobs_all$ytype == 43, ]
  dataobs$ind = ifelse(dataobs$time<lmk, 0, 1)
  dataobs$ind = factor(dataobs$ind)
  
  sim<-simul[simul$id == id, ]
  crp = data.frame(0:30)
  crp = as.data.frame(t(crp))
  crp = crp[-1, ]
  for (l in 1:nrow(sim)){
    for (j in 1:31){
      crp[l, j] = sim$b03[l]+sim$b13[l]*t[j]
    }
  }
  inf_c = apply(crp, 2, quantile, probs = 0.025)
  sup_c = apply(crp, 2, quantile, probs = 0.975)
  med_c = apply(crp, 2, median)
  
  
  dataland = data.frame(time = 0:30, crp = med_c, low_c = inf_c, upp_c = sup_c, ind = as.factor(c(rep(0, lmk), rep(1, (30-lmk+1)))))
  
  gp<-ggplot(data = dataland, aes(x = time, y = crp))+
    geom_ribbon(aes(x = time, ymin = low_c, ymax = upp_c, group = ind, fill = ind), show.legend = F)+
    scale_fill_manual(values = c("white", "grey"))+
    scale_x_continuous(name = "Time (days)", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
    geom_line(lwd = 0.35, col = '#9933FF')+
    scale_y_continuous(name = "CRP (log(mg/L))")+theme_classic()+
    coord_cartesian(ylim = c(0, 20))+
    geom_vline(xintercept = lmk, color = "red", lwd = 0.3)+
    theme(legend.position = "bottom", element_line(size = 0.5), axis.text = element_text(size = 5), axis.title = element_text(size = 5), title = element_text(size = 6))+
    geom_point(data = dataobs, aes(x = time, y = obs, group = ind, col = ind), lwd = 0.5, pch = 16, show.legend = F)+
    scale_color_manual(values = c("black", "#0099FF"))+
    ggtitle("")
  gp
}

# function to create a ggplot object representing the dynamic predictions with 95% prediction interval of the survival probability associated with individual id
# returns a ggplot object 

landmark_survie<-function(lmk, id){
  
  sim = read.table(paste0("intermediate_results/evalAUC/lmk", lmk, "/data_params.txt"), header = T)
  sim = sim[sim$id == id, ]
  dataobs = dataobs_all[dataobs_all$id == id, ]
  dataobs$ind = ifelse(dataobs$time<lmk, 0, 1)
  dataobs$ind = factor(dataobs$ind)
  
  t = 0:30
  
  # create data frame surv for the predictions of cumulative incidence function relative to event 1
  # create data frame surv2 for the predictions of cumulative incidence function relative to event 2
  surv<-data.frame(0:30)
  surv<-as.data.frame(t(surv))
  surv = surv2 = surv[-1, ]
  
  # get parameters without variability
  h01 = sim$h1[1]
  h02 = sim$h2[1]
  alpha12 = sim$alpha12[1]
  alpha13 = sim$alpha13[1]
  alpha22 = sim$alpha22[1]
  alpha23 = sim$alpha23[1]
  
  # to perform survival predictions with simulated individual parameters during a replicate l 
  for (l in 1:nrow(sim)){
    # get parameters simulated during a replicate l
    b02 = sim$b02[l]
    b12 = sim$b12[l]
    b03 = sim$b03[l]
    b13 = sim$b13[l]
    h_cov1 = sim$hcov1[l]
    h_cov2 = sim$hcov2[l]
     
    for (j in 1:31){
      neutro = function(t) b02+b12*t
      crp = function(t) b03+b13*t
      h1 = function(t) h01*exp(alpha12*(neutro(t)-3.6236834)+alpha13*(crp(t)-6.2826893)+h_cov1)
      h2 = function(t) h02*exp(alpha22*(neutro(t)-3.6236834)+alpha22*(crp(t)-6.2826893)+h_cov2)
      F1 = function(t) 1-exp(-(integrate(h1, 0, t)$value))
      F2 = function(t) 1-exp(-(integrate(h2, 0, t)$value))
      # save the cumulative incidence function of event 1 and 2 at time j
      # treat the case of infinite values 
      if (h1(30) == Inf){
        surv[l, j] = 1
      }
      else {surv[l, j] = F1(t[j])}
      surv2[l, j] = F2(t[j])
    }
    # compute conditional probability (see equation 8 in the manuscript)
    for (w in (lmk+2):31){
      if (h1(1000) == Inf){
        surv[l, w] = (1+F2(1000)-surv[l, w]-surv2[l, (lmk+1)])/(1+F2(1000)-surv[l, (lmk+1)]-surv2[l, (lmk+1)])
      }
      else if (h2(1000) == Inf){
        surv[l, w] = (F1(1000)+1-surv[l, w]-surv2[l, (lmk+1)])/(F1(1000)+1-surv[l, (lmk+1)]-surv2[l, (lmk+1)])
      }
      else{
        surv[l, w] = (F1(1000)+F2(1000)-surv[l, w]-surv2[l, (lmk+1)])/(F1(1000)+F2(1000)-surv[l, (lmk+1)]-surv2[l, (lmk+1)])
      }
    }
  }
  
  # get the median of predictions, with 95% prediction interval
  inf = apply(surv, 2, quantile, probs = 0.025)
  sup = apply(surv, 2, quantile, probs = 0.975)
  med = apply(surv, 2, median)
  
  med = c(1, med[(lmk+2):31])
  inf = c(1, inf[(lmk+2):31])
  sup = c(1, sup[(lmk+2):31])
  
  # formatting data to plot the results
  datasurv = data.frame(time = lmk:30, pi = med, low = inf, upp = sup)
  gp<-ggplot(data = datasurv, aes(x = time, y = pi))+
    geom_ribbon(aes(x = time, ymin = low, ymax = upp, fill = 'grey'), show.legend = F)+
    scale_fill_manual(values = c("grey"))+
    scale_x_continuous(name = "Time (days)", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
    geom_line(lwd = 0.5, col = 'black')+
    scale_y_continuous(name = "Survival probability", limits = c(0, 1))+theme_classic()+
    geom_vline(xintercept = lmk, color = "red", lwd = 0.3)+
    geom_segment(aes(x = 0, y = 1, xend = lmk, yend = 1), lwd = 0.5)+
    theme(legend.position = "bottom", element_line(size = 0.5), axis.text = element_text(size = 5), axis.title = element_text(size = 5))+
    {if (data.surv$obs[data.surv$id == id] == 1){
      geom_vline(xintercept = data.surv$time[data.surv$id == id], color = 'black', lwd = 0.3)}}+
    {if (data.surv$obs[data.surv$id == id] == 2){
      geom_vline(xintercept = data.surv$time[data.surv$id == id], linetype = 'dashed', lwd = 0.3)}}
  gp
}







