#simulando

#X - bidimensional
#Y - tridimensional

simulate_data <- function(N = 50, Q, R, K,Sigma, Lambda, x0){

I <- diag(c(1,1))
I1 <- diag(c(1,0))
I2 <- diag(c(0,1))
Y <- matrix(0, nrow = 3, ncol = N)
origX <- matrix(data = 0, nrow = 2, ncol = N)
Um <- matrix(c(1,1),ncol = 1)
X <- matrix(data = 0, nrow = 2, ncol = N)


X[, 1] <- origX[, 1] <- x0

for( i in 2:N )
{
  
  origX[, i] <- (I + R) %*% origX[, i - 1] - R%*%K%*%( I1%*%origX[, i - 1]%*%t(origX[, i - 1])%*%I1 +
                                                         I2%*%origX[, i - 1]%*%t(origX[, i - 1])%*%I2)%*%Um
  
}


for(i in 1:N )
{
  
  X[, i] <- mvrnorm(n = 1, origX[, i], Lambda)
  
}

for( i in 1:N )
{
  
  Y[, i] <- mvrnorm(n = 1,  Q%*%kronecker(X[, i], I)%*%Um, Sigma)
  
}

return(data.frame(Y = t(Y),X = t(X)))

}


