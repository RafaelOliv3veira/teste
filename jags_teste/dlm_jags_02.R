rm(list = ls())
require(rjags)
library(coda)


set.seed(432104)


N <-100

sigma <- 16
kappa <- 0.02
mu_0 <- 20; r <- 0.4; K <- 40; tau <- sqrt(sigma^2*kappa)


x <- rnorm(1 ,mu_0, tau)
y <- rnorm(1, 9*x, sigma)

for(i in 2:N){
  x[i] <- rnorm(1, x[i-1] +  r*x[i-1]*(1-x[i-1]/40) , tau)
  y[i] <- rnorm(1, 9*x[i], sigma)
}


plot(y)



model1.string <-"
model {

r ~ dlnorm(-1.141,16.50165)
K ~ dlnorm(3.5,16.50165)  
sigma ~ dgamma(0.001,0.001)
tau <- pow(sigma, -2)
tau1 <- tau*pow(0.02, -2)
Kw <- 1/K

mu[1] <- K/2
x[1] ~ dnorm(mu[1],tau1)

for(i in 2:N){
mu[i] <- x[i-1] +  r*x[i-1]*(1-Kw*x[i-1])
x[i] ~ dnorm(mu[i],tau1)
}

for (i in 1:N){
ymu[i] <- 9*x[i]
y[i] ~ dnorm(ymu[i], tau)
}


}
"
model1.spec<-textConnection(model1.string)


jags <- jags.model(model1.spec,
                   data = list('y' = y,
                               'N' = N),
                   n.chains=4,
                   n.adapt=100)



params <- c("sigma","r","K")

samps <- coda.samples( jags, params, n.iter=10000 )


summary(samps)
