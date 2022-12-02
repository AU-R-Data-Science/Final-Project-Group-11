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

#' @title Loss Function
#' @description First, solve for p as a step in finding the beta hat optimization. Second, calculate all of the beta hat function before optimization.
#' @param y A \code{vector} of y values to be used in the function.
#' @param X An initial \code{matrix} to be used in the function.
#' @param beta The \code{matrix} of initial beta values produced in the Initial Beta function.
#' @return A \code{numeric} of a value to be optimized.
loss_func <-function(y,X,beta=initial_beta(X,y)){
  p <- numeric(0)
  for (i in 1:length(y)){
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
################################################################################
beta_hat<-function(){
  beta_hat=optim(par=initial_beta(X,y),fn=loss_func,X=X,y=y)$par
  return(beta_hat)
}

################################################################################

data<-function(y,X){
  data = data.frame(y, X)
  return(data)
}

#' @title Bootstrap
#' @description Use a bootstrap approach to optimize the logistic regression and begin to find confidence intervals.
#' @param X The initial \code{matrix} to be used in the function.
#' @param y A \code{vector} of y values.
#' @param alpha A \code{numeric} inputted by the user to find desired quantiles.
#' @param B A \code{numeric} describing the desired amount of bootstrap iterations.
#' @return A \code{matrix} containing the calculated quantiles.
bootstrap <- function(X,y,alpha,B=20){
  data = data.frame("y" = y, X)
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
  rownames(quantiles) <- c(seq(0,ncol(X)-1,by = 1))
  colnames(quantiles) <- c(alpha/2,1-alpha/2)
  return(quantiles)
}


#' @title Logistic Curve
#' @description Plot the logistic regression of the data.
#' @param data The \code{matrix} of the data to use for the plot.
#' @param predictor A \code{numeric} to assign the predictor variable
#' @param beta_hat A \code{matrix} of the optimized beta values.
#' @return A plot depicting a fitted logistic curve to the actual values.
logistic_curve = function(data, predictor, beta_hat) {

  #find out which beta to use
  index = which(colnames(data) == predictor)

  #define new data frame that contains predictor variable
  newdata = data.frame(X=seq(min(data[,index]), max(data[,index]),len=500))

  #predict new values using beta_hat in p function
  new_p = numeric(0)
  for (i in 1:500){
    new_p[i] = 1/(1+exp(-newdata[i,]*beta_hat[index-1]))
  }

  #plot logistic regression curve
  plot(y ~ data[,index], col="steelblue", xlab=predictor)
  lines(new_p ~ newdata$X, lwd=2)

}

metrics <- function(X=X,beta=initial_beta(X,y),beta_hat=beta_hat(),cutoff=0.5){
  p <- numeric()
  p_predicted <- numeric()

  for (i in 1:length(y)){
    p_predicted[i] <- 1/(1+exp(-t(X[i,])%*%beta_hat))
    if(p_predicted[i]<=cutoff)p_predicted[i]=0
    else p_predicted[i]=1
  }

  for (i in 1:length(y)){
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
  c("Prevalence"=Prevalence,"Accuracy"=Accuracy,"Sensitivity"=Sensitivity,
    "Specificity"=Specificity,"False_Discovery_Rate"=False_Discovery_Rate,
    "Diagnostic_Odds_Ratio"=Diagnostic_Odds_Ratio)
}

################################################################################

metricsplot <- function(metric){
  cutoffs <- seq(0.1,0.9,0.1)
  metricsmat <- as.data.frame(matrix(NA,nrow = 6,ncol = length(cutoffs)))

  for (i in 1:length(cutoffs)){
    metricsmat[,i]<-metrics(X=X,beta=initial_beta(X,y),beta_hat=beta_hat(),cutoff = cutoffs[i])
  }
  rownames(metricsmat) <- names(metrics(X=X,beta=initial_beta(X,y),beta_hat=beta_hat(),cutoff = 0.5))
  colnames(metricsmat) <- cutoffs
  plot(y=metricsmat[metric,],x=cutoffs,xlab="Cut-offs",ylab=metric)
}





#' @author Rukesh Gusain, Michael Zirpoli, Erica Maul
#' @export
