#' @title Initial Beta
#' @description Calculate an initial beta value that will be used as the parameter of the optimization function for beta hat, using a least squares equation.
#' @param X An initial \code{matrix} to be used in the function.
#' @param y A \code{vector} of y values to be used in the function.
#' @param B A \code{numeric} (integer) used to denote the number of simulations.
#' @return A \code{matrix} containing the initial beta value(s).
initial_beta<- function(X,y){
  beta = solve(t(X)%*%X)%*%t(X)%*%y
  return(beta)
}

#' @title loss_func
#' @description First, solve for p as a step in finding the beta hat optimization. Second, calculate all of the beta hat function before optimization.
#' @param y A \code{vector} of y values to be used in the function.
#' @param X An initial \code{matrix} to be used in the function.
#' @param beta The \code{matrix} of initial beta values produced in the Initial Beta function.
#' @return A \code{numeric} of a value to be optimized.
loss_func <-function(y,X,beta){

  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
  }

  loss_func = sum(-y*log(p)-(1-y)*log(1-p))
  return(loss_func)
}

#' @title Beta Hat
#' @description Calculate the optimized beta hat value.
#' @param beta The \code{matrix} of initial beta values
#' @param X An initial \code{matrix} to be used in the function.
#' @param y A \code{vector} of y values to be used in the function.
#' @return A \code{matrix} containing the resulting beta hat value(s).
beta_hat <- function(beta, X, y){
  optim(par=beta,fn=loss_func,X=X,y=y)$par
}

#' @title Bootstrap
#' @description Use a bootstrap approach to optimize the logistic regression and begin to find confidence intervals.
#' @param X The initial \code{matrix} to be used in the function.
#' @param y A \code{vector} of y values.
#' @param alpha A \code{numeric} inputted by the user to find desired quantiles.
#' @param B A \code{numeric} describing the desired amount of bootstrap iterations.
#' @return A \code{matrix} containing the calculated quantiles.
bootstrap <- function(X,y,alpha,B){


  beta_mat = matrix(NA, nrow = B, ncol = dim(X)[2])

  for(b in 1:B){

    boot_data = data[sample(1:nrow(data), nrow(data), replace = TRUE), ]

    beta_mat[b, ] = optim(par=beta,fn=loss_func,X=boot_data[,c(2:ncol(boot_data))],y=boot_data[,1])$par

  }
  alpha=0.05
  quantiles <- matrix(NA,nrow=ncol(beta_mat),2)
  for (i in 1:ncol(beta_mat)){
    quantiles[i,]=quantile(beta_mat[,i], c(alpha/2, 1 - alpha/2))
  }
  rownames(quantiles) <- c("B0","B1","B2","B3","B4","B5")
  colnames(quantiles) <- c(alpha/2,1-alpha/2)
  return(quantiles)
}

#' @title Metrics
#' @description Calculate a number of statistical metrics from a confusion matrix comparing predicted and actual data.
#' @param beta The \code{matrix} of initial beta values produced in the Initial Beta function.
#' @param X A \code{matrix} from which data will be pulled to fit the model.
#' @param beta_hat A \code{matrix} of the optimized beta values.
#' @param cutoff A \code{numeric} specified by the user for a cutoff value.
#' @return A \code {list} of the calculated metrics.
metrics <- function(X=X,beta=beta,beta_hat=beta_hat,cutoff=0.5){
  p <- numeric()
  p_predicted <- numeric()

  for (i in 1:10){
    p_predicted[i] <- 1/(1+exp(-t(X[i,])%*%beta_hat))
    if(p_predicted[i]<=cutoff)p_predicted[i]=0
    else p_predicted[i]=1
  }

  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
    if(p[i]<=cutoff)p[i]=0
    else p[i]=1
  }
  TP=0
  TN=0
  FP=0
  FN=0
  confusionmatrix<-rbind(p,p_predicted)
  for (i in 1:ncol(confusionmatrix)){
    if (confusionmatrix["p",i]==1 && confusionmatrix["p_predicted",i]==1) TP=TP+1
    if(confusionmatrix["p",i]==1 && confusionmatrix["p_predicted",i]==0) FN=FN+1
    if(confusionmatrix["p",i]==0 && confusionmatrix["p_predicted",i]==1) FP=FP+1
    if(confusionmatrix["p",i]==0 && confusionmatrix["p_predicted",i]==0) TN=TN+1
  }
  Prevalence <- (FN+TP)/((FN+TP)+(FP+TN))
  Accuracy <- (TP+TN)/(TP+TN+FP+FN)
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  False_Discovery_Rate <- FP/(FP+TP)
  TPR <- Sensitivity
  TNR <- Specificity
  FPR <- 1-TNR
  FNR <- 1-TPR
  LRp <- TPR/FPR
  LRn <- FNR/TNR
  Diagnostic_Odds_Ratio <- LRp/LRn
  list("Prevalence"=Prevalence,"Accuracy"=Accuracy,"Sensitivity"=Sensitivity,
       "Specificity"=Specificity,"False_Discovery_Rate"=False_Discovery_Rate,
       "Diagnostic_Odds_Ratio"=Diagnostic_Odds_Ratio)
}
cutoffs <- seq(0.1,0.9,0.1)
metricsmat <- matrix(NA,nrow = 6,ncol = length(cutoffs))
for (i in 1:length(cutoffs)){
  metricsmat[,i]=metrics(X=X,beta=beta,beta_hat=beta_hat,cutoff = 0.5)
}

#' @title LR Plot
#' @description Plot the logistic regression of the data.
#' @param y The \code{vector} of y values to use for the plot.
#' @param X A \code{matrix} from which data will be pulled to plot.
#' @return A plot depicting a fitted logistic curve to the actual values.
lrplot<- function(X,y){
  plot(y=y,x=X[,3],xlim=c(-3,3))
}

#' @title Get predicted p
#' @description A function that produces the predicted p values for the optimization.
#' @param beta_hat A \code{matrix} of the optimized beta values.
#' @return A \code{vector} of predicted p values.
get_p_predicted <- function(beta_hat){
  p_predicted <- numeric()
  for (i in 1:10){
    p_predicted[i] <- 1/(1+exp(-t(X[i,])%*%beta_hat))
    if(p_predicted[i]<=0.5)p_predicted[i]=0
    else p_predicted[i]=1
  }
}

#' @title Get actual p
#' @description This package produces the actual p values of the data.
#' @param beta The \code{matrix} of the actual, initial beta values.
#' @return A \code{vector} of the actual p values.
#' @author Rukesh Gusain, Michael Zirpoli, Erica Maul
#' @export
get_p_actual <- function(beta){
  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
    if(p[i]<=0.5)p[i]=0
    else p[i]=1
  }
}