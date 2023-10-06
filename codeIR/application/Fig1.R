#### This file is used to derive cumulative incidence functions and associated risk table ploted on Figure 1
#### Not on real data (not shared) ####

###############
# 0. Loading required packages and functions used in the analysis
###############

source("packages.R")
source("functions.R")

###############
# 1. Creating Figure 1
###############

# to load survival data set 
CRdata = read.table("datas/sim_survival.txt", header = T)

# get cumulative incidence functions for each event
CIF = cuminc(CRdata$time, CRdata$obs, 1, cencode = 0)

# get the observed time and associated cumulative incidences 
t_death = CIF$`1 1`$time
c_death = CIF$`1 1`$est
t_disch = CIF$`1 2`$time
c_disch = CIF$`1 2`$est

# plot cumulative incidence functions 
ggp = ggplot(data = NULL, aes(x = t_disch, y = c_disch))+geom_line(aes(x = t_death, y = c_death), col = 'red', lwd = 0.7)+
  geom_line(lwd = 0.7, col = 'blue')+ylab("Cumulative Incidence")+theme_classic()+
  ylim(0, 1)+theme(plot.title = element_text(hjust = 0.5), element_line(size = 1), axis.text = element_text(size = 9), panel.grid.major.y = element_line())+
  annotate(geom = "text", x = 25, y = 0.8, label = "CIF of discharge", col = 'blue', size = 3)+
  annotate(geom = "text", x = 25, y = 0.2, label = "CIF of death", col = 'red', size = 3)+
  scale_x_continuous(name = "Time (days)", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))
ggp

# function giving the number of at-risk, died, discharged and censored patients 
at_risk = function(t){
  atrisk = length(CRdata$id[CRdata$time>t])
  ev1 = length(CRdata$id[CRdata$time<= t & CRdata$obs == 1])
  ev2 = length(CRdata$id[CRdata$time<= t & CRdata$obs == 2])
  cens = length(CRdata$id[CRdata$time<= t & CRdata$obs == 0])
  vec = c(atrisk, ev1, ev2, cens)
  names(vec) = c("atrisk", "ev1", "ev2", "cens")
  vec
}

# build the table 
gp_table = ggplot(data = NULL)+theme_classic()+
  scale_x_continuous(name = "", limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_y_discrete(limits = c("Censored", "Discharged", "Died", "At-risk"), name = "") +
  annotate(geom = "text", x = 0, y = "At-risk", label = as.character(at_risk(0)[1]))+
  annotate(geom = "text", x = 5, y = "At-risk", label = as.character(at_risk(5)[1]))+annotate(geom = "text", x = 10, y = "At-risk", label = as.character(at_risk(10)[1]))+
  annotate(geom = "text", x = 15, y = "At-risk", label = as.character(at_risk(15)[1]))+annotate(geom = "text", x = 20, y = "At-risk", label = as.character(at_risk(20)[1]))+
  annotate(geom = "text", x = 25, y = "At-risk", label = as.character(at_risk(25)[1]))+annotate(geom = "text", x = 30, y = "At-risk", label = as.character(at_risk(30)[1]))+
  
  annotate(geom = "text", x = 0, y = "Died", label = as.character(at_risk(0)[2]))+
  annotate(geom = "text", x = 5, y = "Died", label = as.character(at_risk(5)[2]))+annotate(geom = "text", x = 10, y = "Died", label = as.character(at_risk(10)[2]))+
  annotate(geom = "text", x = 15, y = "Died", label = as.character(at_risk(15)[2]))+annotate(geom = "text", x = 20, y = "Died", label = as.character(at_risk(20)[2]))+
  annotate(geom = "text", x = 25, y = "Died", label = as.character(at_risk(25)[2]))+annotate(geom = "text", x = 30, y = "Died", label = as.character(at_risk(30)[2]))+
  
  annotate(geom = "text", x = 0, y = "Discharged", label = as.character(at_risk(0)[3]))+
  annotate(geom = "text", x = 5, y = "Discharged", label = as.character(at_risk(5)[3]))+annotate(geom = "text", x = 10, y = "Discharged", label = as.character(at_risk(10)[3]))+
  annotate(geom = "text", x = 15, y = "Discharged", label = as.character(at_risk(15)[3]))+annotate(geom = "text", x = 20, y = "Discharged", label = as.character(at_risk(20)[3]))+
  annotate(geom = "text", x = 25, y = "Discharged", label = as.character(at_risk(25)[3]))+annotate(geom = "text", x = 30, y = "Discharged", label = as.character(at_risk(30)[3]))+
  
  annotate(geom = "text", x = 0, y = "Censored", label = as.character(at_risk(0)[4]))+
  annotate(geom = "text", x = 5, y = "Censored", label = as.character(at_risk(5)[4]))+annotate(geom = "text", x = 10, y = "Censored", label = as.character(at_risk(10)[4]))+
  annotate(geom = "text", x = 15, y = "Censored", label = as.character(at_risk(15)[4]))+annotate(geom = "text", x = 20, y = "Censored", label = as.character(at_risk(20)[4]))+
  annotate(geom = "text", x = 25, y = "Censored", label = as.character(at_risk(25)[4]))+annotate(geom = "text", x = 30, y = "Censored", label = as.character(at_risk(30)[4]))


plot_grid(ggp, gp_table, nrow = 2, align = "v", rel_heights = c(1, 0.5))


ggsave("Fig1.png", path = "../results", device = "png", width = 14.5, height = 11, units = "cm", dpi = 800)
