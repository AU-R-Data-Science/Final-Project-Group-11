

initial_beta<- function(X,y){
  beta = solve(t(X)%*%X)%*%t(X)%*%y
  return(beta)
}


X <- matrix(floor(rnorm(n=90)),10,9)

X <- floor(rnorm(n=10))
y <- c(0,0,1,0,0,1,1,1,0,1)
beta <- initial_beta(X,y)

bootstrap <- function(alpha,B=20,beta){

  n <- length(y)

  boot_mean <- rep(NA, B)

  for (i in 1:B){

    y_star <- y[sample(1:n, replace = TRUE)]

    boot_mean[i] <- mean(y_star)
  }

  quantile(boot_mean, c(alpha/2, 1 - alpha/2))

}
bootstrap(alpha = 0.05,beta=beta)

beta_hat2 <- 0
p<-numeric(0)


for (i in 1:10){
  p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
}

beta_hat <-function(y,X,beta){
  for (i in 1:10){
    p[i] <- 1/(1+exp(-t(X[i,])%*%beta))
  }
  beta_hat = sum(-y*log(p)-(1-y)*log(1-p))
  return(beta_hat)
}

beta_hat=optim(par=beta,fn=beta_hat,X=X,y=y)

betahat<- beta_hat$value

x=t(X)%*%beta

plot(y~c(1:10))
fit <- glm(y~c(1:10))
