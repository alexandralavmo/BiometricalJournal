
###############
# 1. Creating model files
###############

sapply(1:p, create_model_true)

###############
# 2. Fitting each model with Monolix
###############

sapply(1:p, estimate_true)

