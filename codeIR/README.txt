We first provide the Version numbers for the applied R, RStudio and required R-packages:

R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE = French_France.1252  LC_CTYPE = French_France.1252    LC_MONETARY = French_France.1252 LC_NUMERIC = C                  
[5] LC_TIME = French_France.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ComplexUpset_1.3.3   mvtnorm_1.1-3        MlxConnectors_2018.2 RJSONIO_1.3-1.8      R6_2.5.1             cowplot_1.1.1       
 [7] ggplot2_3.3.5        timeROC_0.4          cmprsk_2.2-11        survival_3.3-1      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10         compiler_4.0.3      pillar_1.7.0        pec_2022.03.06      iterators_1.0.14    tools_4.0.3        
 [7] digest_0.6.29       lifecycle_1.0.1     tibble_3.1.6        gtable_0.3.0        lattice_0.20-45     pkgconfig_2.0.3    
[13] rlang_1.0.2         Matrix_1.4-1        foreach_1.5.2       cli_3.3.0           patchwork_1.1.2     parallel_4.0.3     
[19] prodlim_2019.11.13  withr_2.5.0         dplyr_1.0.8         generics_0.1.2      vctrs_0.4.1         globals_0.14.0     
[25] tidyselect_1.1.2    grid_4.0.3          glue_1.6.2          listenv_0.8.0       timereg_2.0.2       future.apply_1.9.0 
[31] fansi_1.0.3         parallelly_1.31.1   lava_1.6.10         purrr_0.3.4         magrittr_2.0.3      scales_1.2.0       
[37] codetools_0.2-18    ellipsis_0.3.2      splines_4.0.3       future_1.25.0       colorspace_2.0-3    numDeriv_2016.8-1.1
[43] utf8_1.2.2          munsell_0.5.0       crayon_1.5.1


------------------------------------------------------------

!!! IMPORTANT NOTE !!!
The code is reproductible using R and Monolix software, version 2018R2 (free for academic users, but you need a valid licence) 
1. Monolix software in the correct version can be download at the following adress: https://lixoft.com/download/win64-monolix-suite-2018r2/
2. Ask for a licence by completing the form at:  https://lixoft.com/downloads/  
3. Connect R with the Monolix API using: 
   > install.packages("C:/ProgramData/Lixoft/MonolixSuite2018R2/mlxConnectors/R/MlxConnectors.tar.gz")
   with the correct path referring to the Monolix location 
4. Change the path in files packages.R (in application folder) and parameters.R (in simulation folder)
   initializeMlxConnectors(software = "monolix", mlxDirectory = "C:/ProgramData/Lixoft/MonolixSuite2018R2")
   give the correct Mlx directory
!!!                !!!

------------------------------------------------------------

The code folder is divived in 3 subfolders: 
- application: contains (1) code needed to reproduce figures and table of the application study, 
			(2) datasets corresponding to jittered versions of the true patient datas, 
			(3) model and mlxtran files used in Monolix, and 
			(4) pre-computed intermediate results
			Therefore, tables and figures might slightly differ from the manuscript because the true data cannot be shared.
			Pre-computed intermediate results are derived from the jittered version of the true patient data.
- simulation: contains  (1) code needed to reproduce figures and tables of the simulation study, 
			(2) simulated datasets that can be used to guarantee the reproductibilty of the results,
			(3) model and mlxtran files used in Monolix, and
			(4) pre-computed intermediate results
- results: contains all the tables and figures that can be reproduced 
	   Figures are saved in png format and Tables as RData objects 
	   
In the following, we detail all the files in "application" and "simulation" folders

########################
# Application folder 
########################

subfolder "datas": contains all datasets relative to the application (similar characteristics than true datasets)
subfolder "modfiles": contains the reference model text files and mlxtran files (for Monolix software)
subfolder "intermediate_results": contains the intermediate results to reproduce figures and tables quickly 

To reproduce a single figure or table, just run the R script whose name corresponds to the figure or table name 
(e.g. To reproduce Figure 1, you must run Fig1.R; to reproduce Table S1, you must run TabS1.R)
By default, these scripts use pre-computed intermediate results. If you do not want to use these intermediate result and compute them
from scratch, you must set the variable interresults = F in the R script files (or use the other version of the code available at 
https://iame.catibiomed.fr/index.php/s/kRoENdS5gG1COuh, in which intermediate results are not provided).

The location where intermediate results are stored is explicitly indicated in each file.
Figure and tables are saved in folder "results".

We do not provide the code for the following tables:  
Table 1 since it only contains descriptive statistics
Table 5 since it only details simulations parameters 
Table S2 since it only details the 4C score components 
Table S3a since it only details thresholds of the sensibility analysis

########################
# Simulation folder 
########################

!! IMPORTANT NOTE !!
Data simulation was performed under Linux. We observed slight differences under Windows.
In particular, some longitudinal observations differs (13rd digit after the decimal point)
Although this phenomenon concerns very few observations, we recommand to use the datasets provided in subfolders "data", "data2" and "data_scalability" 
when using Windows operating system to obtain same results. 
Even if using Linux, we recommand to check if data files generated are exactly the same as the ones provided (for instance, using WinMerge)
!!                !! 

subfolder "modfiles": contains the mlxtran and txt model files needed for Monolix software
subfolder "intermediate_results": contains the intermediate results to reproduce figures and tables quickly as the code may be very time-expensive
subfolder "data": contains the simulated data sets under the first scenario of correlation (using Linux operating system)
subfolder "data2": contains the simulated data sets under the strong scenario of correlation (using Linux operating system)
subfolder "data_scalability": contains the simulated data sets to test scalability of the algorithm (using Linux operating system)

To reproduce a single figure or table, just run the R script whose name corresponds to the figure or table name 
(e.g. To reproduce Figure 6, you must run Fig6.R; to reproduce Table 6, you must run Tab6.R)
By default, these scripts use pre-computed intermediate results. If you do not want to use these intermediate result and compute them
from scratch, you must set the variable interresults = F in the R script files (or use the other version of the code available at 
https://iame.catibiomed.fr/index.php/s/kRoENdS5gG1COuh, in which intermediate results are not provided).

The location where intermediate results are stored is explicitly indicated in each file.
Figure and tables are saved in folder "results".
