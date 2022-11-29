
## Create initial beta through the least squares formula
## Function to find initial beta:
initial_beta<- function(X,y){
  beta = solve(t(X)%*%X)%*%t(X)%*%y
  return(beta)
}

## Assign a matrix to "X"
X <- matrix(floor(rnorm(n=90)),10,9)

## Assign vectors to x and y
X <- floor(rnorm(n=10))
y <- c(0,0,1,0,0,1,1,1,0,1)

## Create object "beta" using the initial beta function
beta <- initial_beta(X,y)


## create a function to find bootstrap confidence intervals
bootstrap <- function(alpha,B=20,beta){

  n <- length(y)

  boot_mean <- rep(NA, B)

  for (i in 1:B){

    y_star <- y[sample(1:n, replace = TRUE)]

    boot_mean[i] <- mean(y_star)
  }

  quantile(boot_mean, c(alpha/2, 1 - alpha/2))

}

## Test the bootstrap function using an alpha of 0.05
bootstrap(alpha = 0.05,beta=beta)

beta_hat2 <- 0
p<-numeric(0)


# Solve for p[i]
for (i in 1:10){
  p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
}


# Solve for beta hat
beta_hat <-function(y,X,beta){
  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
  }
  beta_hat = sum(-y*log(p)-(1-y)*log(1-p))
  return(beta_hat)
}

beta_hat=optim(par=beta,fn=beta_hat,X=X,y=y)
#output gives 9 values instead of 1?

# ?
x=t(X)%*%beta

## Working on plotting the linear regression
plot(y~c(1:10))
fit <- glm(y~c(1:10))

## Creating the confusion matrix
y_predicted1 <- beta%*%X
y_predicted2<-numeric()
for (i in 1:length(y_predicted)){
  if (y_predicted[i]<=0.5) y_predicted2[i]<-0
  else y_predicted2[i]<-1
}
confusionmatrix <- rbind(y,y_predicted2)
rownames(confusionmatrix) <- c("actual","predicted")
