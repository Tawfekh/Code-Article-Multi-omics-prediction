# Abdou Rahmane Wade
# abdou.wade@inrae.fr
# INRAE, UMR BioForA, UMR 0588, F-45075 Orleans, France
# 2021

# Phenotype prediction pipeline using Glmnet
# For the article :
  # eQTLs are key players in the integration of 
  # genomic and transcriptomic data for phenotype prediction

#! R version 4.1.0

# Took ~ 2 days with 50 cores of Intel® Xeon(R) Gold 6150 CPU @ 2.70GHz

# Packages
library(anyLib)
pkg <- c("caret","doSNOW","glmnet")
suppressMessages(anyLib(pkg))

# source the prediction functions
source("Phenotype_Prediction_Functions.R")


# Load Data  --------------------------------------------------------------
  # load poplar Data : PoplarData
  # The object PoplarData is a list containing : 
    # Genomic : Genotyping matrix (241 x 428836) 
    # Trancriptomic : Expression level matrix (241 x 34229) 
    # Phenotype dataframe of the 21 traits (241 x 21) 


# Nested Cross Validation Fold Sampling -----------------------------------
PoplarData$CV.fold <- apply(PoplarData$Phenotype,2, FUN = Ncv_Fold_Sampling)

Ridge.mod <- list()



# Regression Genomic only -------------------------------------------------
Ridge.mod$Genomic <- sapply(colnames(PoplarData$Phenotype), 
                          function(trait){
                            Ncv_RidgeReg(Y=PoplarData$Phenotype[,trait],
                                         First.Omic=PoplarData$Genomic, 
                                         Second.Omic=NULL,
                                         CV.fold=PoplarData$CV.fold[[trait]],
                                         cores = 50)
                            },simplify = FALSE)


# Regression Transcriptomic only ------------------------------------------
Ridge.mod$Expr <- sapply(colnames(PoplarData$Phenotype), 
                          function(trait){
                            Ncv_RidgeReg(Y=PoplarData$Phenotype[,trait],
                                         First.Omic=PoplarData$Transcriptomic, 
                                         Second.Omic=NULL,
                                         CV.fold=PoplarData$CV.fold[[trait]],
                                         cores = 50)
                            },simplify = FALSE)



# Regression Concatenate  Genomic and Transcriptomic ----------------------
Ridge.mod$Comb <- sapply(colnames(PoplarData$Phenotype), 
                          function(trait){
                            Ncv_RidgeReg(Y=PoplarData$Phenotype[,trait],
                                         First.Omic=PoplarData$Genomic, 
                                         Second.Omic=PoplarData$Transcriptomic,
                                         CV.fold=PoplarData$CV.fold[[trait]],
                                         cores = 50)
                            },simplify = FALSE)

# Save Output -------------------------------------------------------------
  # save(Ridge.mod, file = "Pred_MultiOmic_Output.RData")

