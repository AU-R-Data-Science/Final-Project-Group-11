
X <- cbind(rep(1,10),matrix((rnorm(n=20)), 10, 2))
y <- c(0,0,1,0,0,1,1,1,0,1)

initial_beta<- function(X,y){
  beta = solve(t(X)%*%X)%*%t(X)%*%y
  return(beta)
}
beta <- initial_beta(X,y)

p<-numeric(0)
for (i in 1:10){
  p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
}

loss_func <-function(y,X,beta){
  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
  }
  loss_func = sum(-y*log(p)-(1-y)*log(1-p))
  return(loss_func)
}

beta_hat=optim(par=beta,fn=loss_func,X=X,y=y)$par

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
bootstrap(X=X,y=y,alpha = 0.05)

glm(y~X[,-1],family = binomial)$coefficients
beta_hat

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
  c("Prevalence"=Prevalence,"Accuracy"=Accuracy,"Sensitivity"=Sensitivity,
       "Specificity"=Specificity,"False_Discovery_Rate"=False_Discovery_Rate,
       "Diagnostic_Odds_Ratio"=Diagnostic_Odds_Ratio)
}
cutoffs <- seq(0.1,0.9,0.1)
metricsmat <- as.data.frame(matrix(NA,nrow = 6,ncol = length(cutoffs)))

for (i in 1:length(cutoffs)){
  metricsmat[,i]<-metrics(X=X,beta=beta,beta_hat=beta_hat,cutoff = cutoffs[i])
}
rownames(metricsmat) <- names(metrics(X=X,beta=beta,beta_hat=beta_hat,cutoff = 0.5))
colnames(metricsmat) <- cutoffs


