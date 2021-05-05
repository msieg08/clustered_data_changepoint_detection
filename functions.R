##' Small wrapper for dir.create
##'
##' Small wrapper for dir.create
##' @title Small wrapper for dir.create 
##' @param dir_path 
##' @return NULL
##' @author Jochen Kruppa
##' @export
mk_dir <- function(dir_path){
  if(!file.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

##' generate expression data per timepoint
##'
##' @title generate expression data per timepoint
##' @param i  time point i
##' @param intercept intercept/expected mean value at time point i
##' @param ngen  number of mothers at timepoint i  
##' @return tibble with expression data for all simulated pups at timepoint i 
##' @export

# Simulation of expression Data
generate_timepoint <- function(i, intercept = 50, ngen = 3){
  
  gen.mother <- defData(varname = "mother", dist = "normal", formula = 0, 
                        variance = 5, id = "idMother")
  gen.mother <- defData(gen.mother, varname = "nPups", dist = "noZeroPoisson",
                        formula = 10)
  
  
  dtMother <- genData(ngen, gen.mother)
  
  
  
  dtPups <- genCluster(dtMother, cLevelVar = "idMother", numIndsVar = "nPups",
                       level1ID = "idPups")
  
  gen.pup <- defDataAdd(varname = "gender", dist = "binary", 
                        formula = 0.5)
  
  gen.pup <- defDataAdd(gen.pup, varname = "expression", dist = "normal", 
                        formula = str_c(intercept, " + mother"), 
                        variance = 2)
  
  
  dtMice <- addColumns(gen.pup, dtPups) %>% as_tibble() 
  
  # set zeros to have  no variance
  if(intercept == 0){
    dtMice$expression <- 0
  }
  
  #clean data set
  dtMice <- dtMice %>% 
    select(expression, idMother, gender, group_id = idPups, efMother = mother) %>%
    mutate(idMother = (idMother + ngen*i - ngen), timepoint = i)
  
  return(dtMice)
}


##' generate expression data per timepoint
##'
##' @title generate expression data per timepoint
##' @param intercepts expected mean value for each time point
##' @param seedSet  seed to set to retrieve reproducible results
##' @param filePath path to outputs 
##' @param fileName common filename for outputs 
##' @return list of simulated data, glm results and glht results
##' @export
perform_analysis <- function(intercepts = rep(20, 12), 
                             seedSet = 1308, 
                             filePath = "" ,
                             fileName = "ts_course", ...){
  
  outputAnalysis <- list()
  
  set.seed(seedSet)
  lsMice <- llply(1:length(intercepts), function(i)
    generate_timepoint(i, intercepts[i], ngen = 3))
  
  dtMice <- bind_rows(lsMice) %>% 
    mutate(timepoint = as.factor(timepoint), idMother = as.factor(idMother))
  
  outputAnalysis$dataMice <- dtMice
  
  dtExp <- dtMice %>% group_by(timepoint) %>%
    dplyr::select(timepoint, expression, group_id) %>% 
    spread(key = timepoint, value = expression)%>% 
    dplyr::select(-group_id) %>% 
    t()
  
  outputAnalysis$dataExp <- dtExp
  
  p01 <- 
    ggplot(dtMice) + 
    geom_jitter(mapping = aes(x = timepoint, y = expression, colour = idMother), 
                width = 0.2, size=2.5) +
    scale_colour_manual(values = unname(polychrome(36))) + 
    labs(y = "Gene Expression Activity", x = "Time", title = "Short Time Series") +
    theme_bw() + theme(legend.position = "none", text = element_text(size=20)) 
  
  ggsave(str_c(filePath, fileName, "_Base.png"), p01, width = 10, height = 7)
  
  
  ## linear mixed-effects model with mean parametrization
  lmer_fit <- lmer(expression ~ 0 + timepoint + (1 | idMother), dtMice)
  
  outputAnalysis$linMixMod <- lmer_fit
  
  # multiple testing correction for change point 
  changepoint_contrast <- contrMat(n = as.numeric(table(dtMice$timepoint)),
                                   type = "Changepoint")
  
  glht_fitCP <- glht(lmer_fit, linfct = changepoint_contrast)
  
  outputAnalysis$glhtChange <- glht_fitCP
  
  cp <- confint(glht_fitCP) %>% 
    tidy %>% 
    ggplot(aes(x=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=estimate)) +
    geom_hline(yintercept=0, linetype="11", colour="grey60") +
    geom_segment(aes(xend=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=conf.low, yend=conf.high), size=0.4, 
                 arrow=arrow(ends="both", length=unit(0.09, "inches"), angle=70)) + 
    geom_point() +    
    geom_hline(yintercept = c(-10,10), color = "blue", linetype = "dotted") +
    coord_flip() +
    theme_bw() +
    labs(y = "Estimate", x = "", title = "95% family-wise confidence level\n(with contrast Changepoint)")
  ggsave(str_c(filePath, fileName, "_CP.png"), cp)
  
  # multiple testing correction for Sequen (neigbhor comparison (-1,1))
  sequen_contrast <- contrMat(n = as.numeric(table(dtMice$timepoint)),
                              type = "Sequen")
  
  glht_fitS <- glht(lmer_fit, linfct = sequen_contrast)
  
  outputAnalysis$glhtSequen <- glht_fitS
  
  sequen <- confint(glht_fitS) %>% 
    tidy %>% 
    ggplot(aes(x=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=estimate)) +
    geom_hline(yintercept=0, linetype="11", colour="grey60") +
    geom_segment(aes(xend=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=conf.low, yend=conf.high), size=0.4, 
                 arrow=arrow(ends="both", length=unit(0.09, "inches"), angle=70)) + 
    geom_point() +
    coord_flip() +
    geom_hline(yintercept = c(-10,10), color = "blue", linetype = "dotted") +
    theme_bw() +
    labs(y = "Estimate", x = "", title = "95% family-wise confidence level\n(with contrast Sequen)")
  ggsave(str_c(filePath, fileName, "_Sequen.png"), sequen)
  
  # multiple testing correction for McDermott
  mcdermott_contrast <- contrMat(n = as.numeric(table(dtMice$timepoint)), 
                                 type = "McDermott")
  
  glht_fitMD <- glht(lmer_fit, linfct = mcdermott_contrast) 
  outputAnalysis$glhtMcDermott <- glht_fitMD
  
  mcd <- confint(glht_fitMD) %>% 
    tidy %>% 
    ggplot(aes(x=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=estimate)) +
    geom_hline(yintercept=0, linetype="11", colour="grey60") +
    geom_segment(aes(xend=reorder(fct_inorder(contrast), desc(fct_inorder(contrast))), y=conf.low, yend=conf.high), size=0.4, 
                 arrow=arrow(ends="both", length=unit(0.09, "inches"), angle=70)) + 
    geom_point() +
    coord_flip() +
    geom_hline(yintercept = c(-10,10), color = "blue", linetype = "dotted") +
    theme_bw() +
    labs(y = "Estimate", x = "", title = "95% family-wise confidence level\n(with contrast McDermott)")
  ggsave( str_c(filePath, fileName, "_McDermott.png"), mcd)
  
  # combine all result plots into one
  model_pl <- ggarrange(p01, cp, sequen, mcd, 
                        ncol=2, nrow=2,   legend="none", 
                        labels = c("a)", "b)", "c)", "d)"))
  
  ggsave(str_c(filePath, fileName, "_complete.png"),
         model_pl, width = 10, height = 10)
  
  return(outputAnalysis)
}