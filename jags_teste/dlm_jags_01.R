rm(list = ls())
require(rjags)
library(coda)


set.seed(432104)


N <-100

sigma <- 0.5
mu_0 <- 2; r <- 0.4;  tau <- sqrt(sigma^2*1)


x <- rnorm(1 ,mu_0, tau)
y <- rnorm(1, x, sigma)

for(i in 2:N){
  x[i] <- rnorm(1, r*x[i-1], tau)
  y[i] <- rnorm(1, x[i], sigma)
}


plot(y)



model1.string <-"
  model {

  r ~ dlnorm(-1.141,16.50165)
  sigma ~ dgamma(0.001,0.001)
  tau <- pow(sigma, -2)
  
  mu[1] <- 2
  x[1] ~ dnorm(mu[1],tau)

  for(i in 2:N){
  mu[i] <- r*x[i-1]
  x[i] ~ dnorm(mu[i],tau)
  }

    for (i in 1:N){
    y[i] ~ dnorm(x[i], tau)
    }


}
"
model1.spec<-textConnection(model1.string)


jags <- jags.model(model1.spec,
                   data = list('y' = y,
                               'N' = N),
                   n.chains=4,
                   n.adapt=100)



params <- c("sigma","r")

samps <- coda.samples( jags, params, n.iter=10000 )


summary(samps)

