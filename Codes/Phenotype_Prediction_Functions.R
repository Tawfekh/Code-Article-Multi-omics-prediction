# Abdou Rahmane Wade
# abdou.wade@inrae.fr
# INRAE, UMR BioForA, UMR 0588, F-45075 Orleans, France
# 2021

# Phenotype prediction functions  
# For the article :
  # eQTLs are key players in the integration of 
  # genomic and transcriptomic data for phenotype prediction



# Fold Sampling for Nested Cross Validation -------------------------------
#' @description Fold Sampling for Nested Cross Validation
#' @param Y phenotype vector
#' @param outer.fold number of outer fold of the Nested Cross Validation
#' @param inner.fold number of inner fold of the Nested Cross Validation
#' @param iter number of iteration of the Nested Cross Validation
#' @return a list of inner sample for each outer fold of each iteration
Ncv_Fold_Sampling <- function(Y, outer.fold = 5,inner.fold=10, iter = 50){
  sel <- !is.na(Y)
  Y <- Y[sel]
  
  fold <- list()
  for(i in 1:iter)
  {
    
    outer.foldid <- sample(1:outer.fold,size=length(Y),replace=TRUE, 
                           prob = rep(1/outer.fold,outer.fold))
    in.id <- list()
    
    for(outer.t in 1:outer.fold)
    {
      
      outer.train <- which(outer.foldid!=outer.t)
      Y.outer.train <- Y[outer.train]
      in.id[[outer.t]]=sample(1:inner.fold,size=length(Y.outer.train),replace=TRUE, 
                              prob = rep(1/inner.fold,inner.fold))
      
    }
    
    fold[[i]] <- list("out.id"= outer.foldid,
                      "in.id" = in.id)
    
  }
  
  return(fold)
}


# Nested Cross Validation Ridge Regression --------------------------------
#' @description Ridge Regression using the R package glmnet with standard nested cross validation framework
#' @param Y phenotype vector
#' @param First.Omic First omic matrix
#' @param Second.Omic Second omic matrix to concatenate with the First omic matrix
#' @param cores cores number for parallel computing
#' @param alpha glmnet parameter alpha=0 for rigde regression
#' @param CV.fold Fold Sampling for Nested Cross Validation
#' @return list of : variable effects, mu, ridge reg lambda, 
#' @return model accuracies, variation du to genetic and variation du to errors
Ncv_RidgeReg <- function(Y, First.Omic, Second.Omic=NULL, cores = 1, alpha = 0, CV.fold){
  # Package Requirements
  require(doSNOW)
  require(glmnet)
  require(caret)

  #checks
  X <- First.Omic
  Z <- Second.Omic
  stopifnot(length(Y) == nrow(X))
  if(!is.null(Z)){
    stopifnot(length(Y) == nrow(Z))
  }
  sel <- !is.na(Y)
  X <- X[sel,]
  Y <- Y[sel]
  if(!is.null(Z)){
    Z <- Z[sel,]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("glmnet","caret"), .options.snow=opts) %dopar% {
    
    # Output object
    beta <- list(); lambda.best <- mu <- R2 <- RMSE <- Rho <- VarA <- VarE <- numeric()
    
    # Outer-fold
    outer.foldid=CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for(outer.t in 1:outer.fold){
      
      outer.test <- which(outer.foldid==outer.t)
      X.outer.test <- X[outer.test,]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid!=outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train,]
      if(!is.null(Z)){
        Z.outer.test <- Z[outer.test,]
        Z.outer.train <- Z[outer.train,]
      }
      
      # Inner-fold
      inner.foldid=CV.fold[[i]]$in.id[[outer.t]]
      
      # Fit inner model using inner-fold CV: model
      if(is.null(Z)){
        CV.inner.model <- cv.glmnet(x=X.outer.train,y=Y.outer.train,foldid=inner.foldid,alpha=alpha,
                                    keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        VarA[outer.t] <- var(X.outer.test%*%beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train)-VarA[outer.t]
        
      }else{
        
        # Error Msg
        if(is.null(Z)){
          stop("Need Second omic data")
        }
        
        # Fit inner model using inner-fold CV: model
        CV.inner.model <- cv.glmnet(x=cbind(X.outer.train,Z.outer.train),
                                    y=Y.outer.train,foldid=inner.foldid,alpha=alpha,
                                    keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        VarA[outer.t] <- var(cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train)-VarA[outer.t]
        
      }
    }
    
    # Parallel Comp Output
    return(list( 
      "beta"= beta[[which.max(R2)]],
      "lambda" = lambda.best[which.max(R2)],
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)],
      "VarA.list" = VarA[which.max(R2)],
      "VarE.list" = VarE[which.max(R2)],
      "mu" = mu[which.max(R2)]
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list("beta" = do.call("rbind", Y_pred[which(1:length(Y_pred) %% 8 == 1)]), 
              "lambda" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 2)]), 
              "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]),
              "R2" = list("R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)])),
                          "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]))
              ),
              "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]),
              "RMSE" = list("RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)])),
                            "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]))
              ),
              "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]),
              "Rho" = list("Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)])),
                           "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]))
              ),
              "VarA.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 6)]),
              "VarE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 7)]),
              "mu" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 0)]),
              "aplha"=alpha))
}

