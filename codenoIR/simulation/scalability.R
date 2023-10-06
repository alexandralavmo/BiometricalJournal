# File used to create and run the model to test the scalability of the algorithm depending on the number of patients and the number of biomarkers

###############
# 1. Creating model files
###############
# creating files to test scalability according to the number of patients 
sapply(n_scalable, create_model_scaN)

# creating files to test scalability according to the number of biomarkers
sapply(k_scalable, create_model_scaK)

###############
# 2. Fitting each model with Monolix
###############

sapply(n_scalable, cputime_n) 
sapply(k_scalable, cputime_k) 



