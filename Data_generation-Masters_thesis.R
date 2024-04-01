set.seed(210399)
#Linear model
sim_func1 <- function(n, a=NULL) {
  X1 <- runif(n, -2, 2)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.2)
  #Compute the true value of the intervention specific mean
  if (length(a)>0) {
    A <-a
  } else {
    A <- rbinom(n, 1, prob=plogis(-0.25 + 0.8*X1 + 0.25*X3)) }
  Y <- rbinom(n, 1, prob=plogis(-0.9 + 1.9*X1 + 0.6*X2 + 0.5*A))
  if (length(a)>0) {
    return(mean(Y))
  } else {
    return(data.table(X1=X1,X2=X2,X3=X3,A=A,Y=Y))
  } 
}

#Non-linear model with jumps
sim_func2 <- function(n, a=NULL){
  X1 <- runif(n,0,4)
  X2 <- rnorm(n,0,2)
  X3 <- rnorm(1, 0.5)
  #Compute the true value of the intervention specific mean
  if (length(a)>0) {
    A <- a
  } else {
    A <- rbinom(n, 1, prob=plogis(-as.numeric(X1<= 1/3) + as.numeric(X2 > 4/3 & X2 <= 1/2) - 
                                    min(X3,10)*as.numeric(X1 > 1/2 & X1 <= 3)))}
  Y <- rbinom(n, 1, prob=abs(sin(as.numeric(X2<= 4/3) + 
                                   X2*as.numeric(X1 > 1/2 & X1 <= 3) + (2*A-1)*exp(X1))))
  if (length(a)>0) {
    return(mean(Y))
  } else {
    return(data.table(X1=X1,X2=X2,X3=X3,A=A,Y=Y))
  }
}