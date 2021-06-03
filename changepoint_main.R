pacman::p_load(MASS,
               plyr,
               magrittr,
               knitr,
               tidyverse,
               simstudy, 
               multcomp, 
               lme4, 
               effects,
               pals,
               ggpubr, 
               broom)

#### Set Path and Read in functions ####
# Please adjust to own paths
path_in <- ""
source(str_c(path_in, "functions.R"))

#### Intercepts ####
# starting mean value for each time point for 11 different scenarios
intercepts01 <- rep(25, 12)
intercepts02 <- seq(16,39, 2)
intercepts03 <- c(rep(25, 6), rep(35, 6))
intercepts04 <- c(rep(25, 4), rep(35, 4), rep(45, 4))
intercepts05 <- c(rep(15, 3), rep(25, 3), rep(35, 3), rep(45, 3))
intercepts06 <-  seq(39, 16, -2)
intercepts07 <- c(rep(45, 6), rep(35, 6))
intercepts08 <- c(rep(45, 4), rep(35, 4), rep(25, 4))
intercepts09 <- c(rep(60, 3), rep(50, 3), rep(40, 3), rep(30, 3))
intercepts10 <- c(rep(10, 4), rep(0, 4) , rep(10, 4))
intercepts11 <- c(rep(25, 5), rep(15, 2) , rep(25, 5))

intercepts <- list(intercepts01, intercepts02, intercepts03, intercepts04, 
                    intercepts05, intercepts06, intercepts07, intercepts08, 
                    intercepts09, intercepts10, intercepts11)

# Make Directories for Plots for each Scenario
# Please adjust to own paths
dirs_plot <- 
  str_c(path_in, "plots/course", 
        sprintf("%02d",1:length(intercepts)), "/")
tmp <- llply(dirs_plot, mk_dir)

#### Analysis ####
# Changepoint detection for all selected scenarios
analysis <- llply(1:length(intercepts), 
                  function(i){
                    perform_analysis(intercepts[i][[1]], 
                                     filePath = dirs_plot[i],
                                     fileName = str_c(format(Sys.time(), "%Y%m%d"),
                                                      "TimeSeries"))})
