# Load the packages needed to perform the analysis 

library(ggplot2)
library(cowplot)
library(survival)
library(cmprsk)
library(timeROC)
library(R6)
library(RJSONIO)
library(MlxConnectors)
# !!! Change directory to specify correct Monolix localization !!! #
initializeMlxConnectors(software = "monolix", mlxDirectory = "C:/ProgramData/Lixoft/MonolixSuite2018R2")
