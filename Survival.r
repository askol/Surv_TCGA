library(reshape2)
library(dplyr)

source("/gpfs/data/stranger-lab/askol/TCGA2/Survival/Code/Survival_funcs.r")

DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/Expression/"
OutDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Survival/Results/"

## RUN SURVIVAL USING RESIDUALS AFTER ACCOUNTING FOR SVS AND BIRTH YEAR ##
## RUN SURVIVAL USING LCPM PRIOR TO FITTING SVS AND BIRHTYEAR ##
## SURVIVAL PERFORMED ON LCPM SEPERATELY IN MALES AND FEMALES AND THEN COMBINED ##

run.surv <- function(project){
    
    ## GET EXPRESSION AND SURVIVAL DATA ##
    data <- get.data(project, DataDir)
    
    zRmBat <- data[[1]]
    cpm.invnorm <- data[[2]]
    lcpm <- data[[3]]
    gene.info <- data[[4]]
    
    ## PERFORM SURVIVAL ANALYSIS USING EACH GENE ONE AT A TIME ##
    ## PERFORM ONCE USING SVS AND AGAIN USING NONE
    
    ## zRMBat (Batch effects (birthyear, svs) removed and then rank normal transformed
    zRMBat.out <- calc.surv(zRmBat, project, OutDir, prefix="zRMBat", gene.info,
                            use.covs=FALSE)
    
    ## cpm.invnorm (rank inverse normal transformed with and without svs and birthyear ##
    cpm.invnorm.simp <- calc.surv(cpm.invnorm, project, OutDir,
                                  prefix="lcpm.invnorm.simp", gene.info,
                                  use.covs=FALSE)
    cpm.invnorm.covs <- calc.surv(cpm.invnorm, project, OutDir,
                                  prefix = "lcpm.invnorm.covs", gene.info,
                                  use.covs=TRUE)
    
    ## lcpm with and without svs and birthyear ##
    lcpm.simp <- calc.surv(lcpm, project, OutDir, prefix="lcpm.simp", gene.info,
                           use.covs=FALSE)
    lcpm.covs <- calc.surv(lcpm, project, OutDir,
                           prefix = "lcpm.covs", gene.info, use.covs=TRUE)
    
}


args <- commandArgs(TRUE)
project = args[1]
run.surv(project)
print(date())





    


    
