#### To derive scalability plots assessing CPU time needed to fit several joint models 

###############
# 0. Loading required functions and parameters used in the simulation
###############

source("parameters.R")
source("functions.R")

###############
# 1. Runing functions performing the analysis
###############

# Executes the script that simulates datasets under the first scenario of correlation (writes in subfolder "data") 
# Warning: some differences in the 13rd digit after the coma point have been reported for some simulated values 
# depending on the operated system and may affect the results 
# Comment the following line if you want to rely on pre-computed intermediate results. Beware that these results are not 
# included by default in this version of the code. They can be downloaded at https://iame.catibiomed.fr/index.php/s/uwuLTLw6ZboYBah
source("simulation.R")

# Executes the script that runs models to test scalability of the approach (uses data stored in subfolder "data_scalability" and 
# writes in subfolder "intermediate_results/fit_scalability") 
# Warning: very time consuming (several days)... to be run on a computing cluster
# Comment the following line if you want to rely on pre-computed intermediate results. Beware that these results are not 
# included by default in this version of the code. They can be downloaded at https://iame.catibiomed.fr/index.php/s/uwuLTLw6ZboYBah
source("scalability.R") 

###############
# 2. Importation of results
###############

# 1.a Scalability according to the number of patients N 

times = sapply(n_scalable, get_cpu)

# divide the times by 3600 to obtain results in hour scale
cpu_par = times[1, ]/3600
cpu_fim = times[2, ]/3600
cpu_time = colSums(times)/3600


# 1.b Scalability according to the number of biomarkers K  

# define 3 vectors for the cpu time spent to estimate parameters, Fisher Information Matrix and both for each model involving 1 to 7 biomarkers
cpu_par_k = c()
cpu_fim_k = c()
cpu_time_k = c()

# incremente the three vectors 
{res = readLines("intermediate_results/fit_scalability/mod1/summary.txt")
l1 = unlist(strsplit(res[31], split = " "))
l2 = unlist(strsplit(res[58], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod2/summary.txt")
l1 = unlist(strsplit(res[39], split = " "))
l2 = unlist(strsplit(res[73], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod3/summary.txt")
l1 = unlist(strsplit(res[47], split = " "))
l2 = unlist(strsplit(res[88], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod4/summary.txt")
l1 = unlist(strsplit(res[55], split = " "))
l2 = unlist(strsplit(res[103], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod5/summary.txt")
l1 = unlist(strsplit(res[63], split = " "))
l2 = unlist(strsplit(res[117], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod6/summary.txt")
l1 = unlist(strsplit(res[71], split = " "))
l2 = unlist(strsplit(res[132], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))

res = readLines("intermediate_results/fit_scalability/mod7/summary.txt")
l1 = unlist(strsplit(res[79], split = " "))
l2 = unlist(strsplit(res[147], split = " "))

cpu_par_k = c(cpu_par_k, as.numeric(rev(l1)[1]))
cpu_fim_k = c(cpu_fim_k, as.numeric(rev(l2)[1]))
cpu_time_k = c(cpu_time_k, as.numeric(rev(l1)[1])+as.numeric(rev(l2)[1]))
}

# divide the times by 3600 to obtain results in hour scale
cpu_par_k = cpu_par_k/3600
cpu_fim_k = cpu_fim_k/3600
cpu_time_k = cpu_time_k/3600

###############
# 3. Results (plots)
############### 

# 2.a Figure S5a: Scalability according to the number of patients N 

data_cpu = data.frame(n = rep(n_scalable, 3), cpu = c(cpu_par, cpu_fim, cpu_time), 
                      grp = c(rep("Parameters", length(n_scalable)), rep("Standard errors", length(n_scalable)), rep("Total", length(n_scalable))))
data_cpu$grp = factor(data_cpu$grp, labels = c("Parameters", "Standard errors", "Total"))

gp = ggplot(data = NULL, aes(x = data_cpu$n, y = data_cpu$cpu, col = data_cpu$grp))+geom_point()+geom_line(linetype = 'dashed')+
  scale_color_manual(values = c("blue3", "cornflowerblue", "lightskyblue1"))+
  theme_classic()+scale_y_continuous(name = "CPU time (hours)")+
  scale_x_continuous(name = "Number of patients", breaks = c(100, 200, 300, 400, 800, 1600))+
  theme(title = element_text(size = 7), legend.title = element_blank(), legend.position = c(0.2, 0.8))
gp
ggsave("FigS5a.png", device = "png", path = "../results", width = 15, height = 8, units = "cm", dpi = 1000)


# 2.b Figure S5b: Scalability according to the number of biomarkers K 

data_cpu_k = data.frame(n = rep(k_scalable, 3), cpu = c(cpu_par_k, cpu_fim_k, cpu_time_k), 
                        grp = c(rep("Parameters", length(k_scalable)), rep("Standard errors", length(k_scalable)), rep("Total", length(k_scalable))))
data_cpu_k$grp = factor(data_cpu_k$grp, labels = c("Parameters", "Standard errors", "Total"))

gpk = ggplot(data = NULL, aes(x = data_cpu_k$n, y = data_cpu_k$cpu, col = data_cpu_k$grp))+geom_point()+geom_line(linetype = 'dashed')+
  scale_color_manual(values = c("blue3", "cornflowerblue", "lightskyblue1"))+
  theme_classic()+scale_y_continuous(name = "CPU time (hours)")+
  scale_x_continuous(name = "Number of biomarkers", breaks = c(1, 2, 3, 4, 5, 6, 7))+
  theme(title = element_text(size = 7), legend.title = element_blank(), legend.position = c(0.2, 0.8))
gpk

ggsave("FigS5b.png", device = "png", path = "../results", width = 15, height = 8, units = "cm", dpi = 1000)

