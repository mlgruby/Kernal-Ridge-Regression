rm(list = ls())
setwd("~/Desktop/ML 3")

Advertising <- read.csv("~/Desktop/ML 3/Advertising.csv")
Advertising$X <- NULL

set.seed(123)
train <- sample(1:nrow(Advertising), 2*nrow(Advertising)/3)

yTrain <- Advertising$Sales[train]
xTrain <- Advertising[train,-4]
xTest <- Advertising[-train,-4]
yTest <- Advertising$Sales[-train]

#######################################General Functions################################################
rep.row <- function(x,n){
  # returns a matrix in which vector x is repeated n times row wise
  matrix(rep(x,each=n),nrow=n)
}

polykernel <- function(x, xPrime, d, normalize = FALSE){
  # returns polynomial kernel
  if(normalize){
    return (((1 + as.numeric(x) %*% as.numeric(t(xPrime)))^d + 1)/
              (((1 + as.numeric(x) %*% as.numeric(t(x)))^d + 1) * ((1 + as.numeric(xPrime) %*% as.numeric(t(xPrime)))^d + 1))) # intercept addend is done by adding 1 at the end of calculation
  } else {
    return ((1 + as.numeric(x) %*% as.numeric(t(xPrime)))^d + 1) # intercept addend is done by adding 1 at the end of calculation
  }
}

radialkernel <- function(x, xPrime, gamma, normalize = FALSE){
  # returns radial kernel
  if(normalize){
    return((exp(-gamma*as.numeric((sum((x-xPrime)^2)))) + 1)/2) # intercept addend is done by adding 1 at the end of calculation
  } else {
    return(exp(-gamma*as.numeric((sum((x-xPrime)^2)))) + 1) # intercept addend is done by adding 1 at the end of calculation
  }
}


KRR <- function(train, test, y, lambda = 0.1, kernel = NULL, normalize = FALSE, d = 2, gamma = 0.001) {
  ntest <- nrow(test)
  yPred <- 0
  nrw <- nrow(train)
  train <- as.matrix(train)
  KXX <- matrix(data = 0, nrow = nrw, ncol = nrw)
  KX <- matrix(data = 0, nrow = nrw, ncol = 1)
  if (!is.null(kernel) && kernel == "polynomial"){
    for(i in 1:ntest){
      testPoint <- as.numeric(test[i,])
      # train <- as.matrix(train)
      for(j in 1:nrow(train)){
        for(k in 1:nrow(train)){
          KXX[j,k] <- polykernel(train[j,], train[k,], d, normalize)
        }
        KX[j,1] <- polykernel(train[j,], testPoint, d, normalize)
      }
      yPred[i] <- t(y) %*% solve((lambda * diag(nrw) + (KXX))) %*% (KX)
    }
  } else if(!is.null(kernel) && kernel == "radial"){
    for(i in 1:ntest){
      testPoint <- as.numeric(test[i,])
      # train <- as.matrix(train)
      for(j in 1:nrow(train)){
        for(k in 1:nrow(train)){
          KXX[j,k] <- radialkernel(train[j,], train[k,], gamma, normalize)
        }
        KX[j,1] <- radialkernel(train[j,], testPoint, gamma, normalize)
      }
      yPred[i] <- t(y) %*% solve((lambda * diag(nrw) + (KXX+1))) %*% (KX+1)
    }
  } else if (is.null(kernel)){
    print("Ridge Regression without any kernel")
    for(i in 1:ntest){
      testPoint <- as.numeric(test[i,])
      Xx <- as.numeric(apply(data.frame(train * data.frame(rep.row(testPoint, nrow(train)))), 1, sum))
      yPred[i] <- t(y) %*% solve((lambda * diag(nrw) + train %*% t(train))) %*% Xx
    }
  }
  return(yPred)
}

mse <- function(y, yHat){
  return((sum((y-yHat)^2))/length(y))
}

#########################################All functions ends here###############################################


# For gamma best value
gammaSeq <- seq(0.00001, 0.001, 0.00002)
mseList <- 0
for(i in 1:length(gammaSeq)){
  yPred <- KRR(train = xTrain, test = xTest, lambda = 0.1, y = yTrain, kernel = "radial", d = 2, gamma = gammaSeq[i])
  mseList[i] <- mse(yTest, yPred)
  print(i)
}
plot(gammaSeq, mseList, type = "l", xlab = "Gamma", ylab = "MSE", main = "MSE vs Gamma for Radial kernel with lamba = 0.1")

bestGamma <- gammaSeq[which.min(mseList)] # Best Gamma


# For lambda best value
lambdaSeq <- c(0.00000000001, 0.0000000001, 0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1)
lambdaMseList <- 0
for(i in 1:length(lambdaSeq)){
  yPred <- KRR(train = xTrain, test = xTest, lambda = lambdaSeq[i], y = yTrain, kernel = "radial", d = 2, gamma = bestGamma)
  lambdaMseList[i] <- mse(yTest, yPred)
  print(i)
}
plot(lambdaSeq, lambdaMseList, type = "l", xlab = "lambda", ylab = "MSE", main = "MSE vs Lambda for radial kernel with Gamma = 5e-05")

bestLambda <- lambdaSeq[which.min(lambdaMseList)] # Best Lambda

# With normalize Kernel 
# Radial Kenel 
yPred <- KRR(train = xTrain, test = xTest, lambda = 0.1, y = yTrain, kernel = "radial",normalize = TRUE, d = 2, gamma = 0.005)
mse(yTest, yPred)

yPred <- KRR(train = xTrain, test = xTest, lambda = 0.1, y = yTrain, kernel = "polynomial",normalize = TRUE, d = 5, gamma = 0.005)
mse(yTest, yPred)

