#### File used to derive longitudinal evolutions of biomarkers 
#### Do not exactly reproduce Figure 3 as the true data set is not shared 

data = read.table("datas/alldata.txt", header = T)
surv = read.table("datas/sim_survival.txt", header = T)

# get ggplot objects for the evolution of final biomarkers 24 and 43 
plot_mark(24)
plot_mark(43)

# arranging plot 
plot_grid(gp24, gp43, ncol = 2, byrow = T)

# save results
ggsave("Fig3.png", path = "../results", device = "png", width = 20, height = 12, units = "cm", dpi = 800)