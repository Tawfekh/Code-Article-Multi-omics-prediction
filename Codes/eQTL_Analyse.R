# Abdou Rahmane Wade
# abdou.wade@inrae.fr
# INRAE, UMR BioForA, UMR 0588, F-45075 Orleans, France
# 2021

# eQTL Analysis using the Multi-Loci Mixed-Model (MLMM)
# For the article :
  # eQTLs are key players in the integration of 
  # genomic and transcriptomic data for phenotype prediction

#! R version 4.1.0

# With 70 cores of Intel® Xeon(R) Gold 6150 CPU @ 2.70GHz

# Packages
library(anyLib)
pkg <- c("data.table","emma","mlmm","parallel","purrr","pbmcapply","plyr","AGHmatrix")
anyLib(pkg)

# Load Data  --------------------------------------------------------------
 # load poplar Data : PoplarData
 # The object PoplarData is a list containing : 
    # Genomic : Genotyping matrix (241 x 428836) 
    # Trancriptomic : Expression level matrix (241 x 34229)


# Compute genomic relationship matrix (241 x 241) -------------------------
GRM <- Gmatrix(PoplarData$Genomic, method="VanRaden")


# Function Multicore mlmm ---------------------------------------------------
#' @description Function in oder to compute the eQTL analysis
#' @param gene_id Gene name : character
#' @param Geno Genotyping matrix
#' @param GRM Genomic relationship matrix
#' @param nbchunks an integer defining the number of chunks of X to run the analysis, 
#' allows to decrease the memory usage ==> minimum=2, 
#' increase it if you do not have enough memory
#' @param maxsteps maximum number of steps desired in the forward approach
#' @return a list of Significant SNPs for step Bonferonni or Step0 ie No DL
My.mlmm <- function(gene_id, Geno = Geno, GRM = GRM, nbchunks = 6, maxsteps = 10){
  Y <- PoplarData$Transcriptomic[,gene_id]
  mygwas <- mlmm(Y = Y, X = Geno, K = GRM, nbchunks = nbchunks, maxsteps = maxsteps)
  ##output
  ###signif snps @ first step
  if(nrow(subset(mygwas$pval_step[[1]]$out, -log10(pval) >= mygwas$bonf_thresh)) == 0) {
    signifSNPs.step0 <- NULL
  } else {
    signifSNPs.step0 <- subset(mygwas$pval_step[[1]]$out, -log10(pval) >= mygwas$bonf_thresh)
  }
  ###EBIC
  if (is.na(mygwas$opt_extBIC$cof[1])) {
    outEBIC <- NULL
  } else {
    outEBIC <- data.frame("Potri" = gene_id, subset(mygwas$opt_extBIC$out, SNP %in% mygwas$opt_extBIC$cof))
  }
  rownames(outEBIC) <- NULL
  ###mBonf
  if (is.na(mygwas$opt_mbonf$cof[1])) {
    outMBonf <- NULL
  } else {
    outMBonf <- data.frame("Potri" = gene_id, subset(mygwas$opt_mbonf$out, SNP %in% mygwas$opt_mbonf$cof))
  }
  rownames(outMBonf) <- NULL
  ###output
  output <- list("steptable" = mygwas$step_table,
                 "EBIC" = outEBIC,
                 "mBonf" = outMBonf,
                 "signifSNPs.step0" = signifSNPs.step0)
  return(output)
}
My.mlmm_ok <- quietly(possibly(My.mlmm, otherwise = "Error")) # help to identify potential errors



# Run Multicore mlmm ------------------------------------------------------

Gene.list <- colnames(PoplarData$Transcriptomic)

MyeQTLsAnalysis.out <- pbmclapply(Gene.list,
                                  FUN=My.mlmm_ok,
                                  Geno=PoplarData$Genomic, GRM=GRM, 
                                  nbchunks=6, 
                                  maxsteps=10,
                                  mc.cores=70,
                                  mc.style = "txt")

names(MyeQTLsAnalysis.out) <- Gene.list


# Verification that there are no execution errors for certain genes
Error.run <- sapply(MyeQTLsAnalysis.out, function(i) length(i) < 4) # none
# Save eQTL Analysis Output
save(MyeQTLsAnalysis.out,Error.run, file = "MLMM_eQTL.Rdata")


# Extract Step0 result ----------------------------------------------------
eQTLsSignStep0 <- ldply(lapply(MyeQTLsAnalysis.out[!Error.run],
                               function(x){
                                 x$result$signifSNPs.step0
                               }), data.frame)

colnames(eQTLsSignStep0)[1]  <- "Potri"
eQTLsSignStep0$Potri <- substr(eQTLsSignStep0$Potri, 0, 16)



# Extract Step_op result --------------------------------------------------
eQTLsSignStepEnd <- ldply(lapply(MyeQTLsAnalysis.out[!Error.run],
                                 function(x){
                                   x$result$mBonf
                                 }), data.frame)

eQTLsSignStepEnd  <- eQTLsSignStepEnd[-1,]
eQTLsSignStepEnd$Potri <- substr(eQTLsSignStepEnd$Potri, 0, 16)

save(eQTLsSignStep0,eQTLsSignStepEnd, file = "MatriceDesResultas_eQTL.Rdata")