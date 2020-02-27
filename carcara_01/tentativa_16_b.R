require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(LaplacesDemon)
require(MASS)
require(ggplot2)
require(gridExtra)
require(xtable)

rm(list = ls())

#Estimando os estados, a matriz diagonal R e as matrizes Sigma, Lambda e Q


Tabela <-matrix(0,ncol = 16, nrow = 3)
N <-50
Q <- matrix(c(9,0,0,0,0,5,0,0,0,0,15,0),ncol =4,byrow = T)
R <- diag(c(0.4,0.6))
K <- diag(c(0.025,0.016))
Sigma <- rinvwishart(3, diag(c(40^2,10^2,50^2)))
Lambda <- rinvwishart(2, diag(c(0.1^2,0.1^2)))
x0 <- c(20,30)

source("/home/posmae/leafaros/Allcodigos/teste/carcara_01/simulate_data_function.R")




Dados <- simulate_data(N , Q, R, K,Sigma, Lambda, x0)

Y <- as.matrix(Dados[,1:3])
# X <- Dados[,4:5]

mu <- matrix(c(0, 0), nrow = 2, ncol = 1)
Um <- matrix(c(1,1),ncol = 1)
I <- diag(c(1,1))
I1 <- diag(c(1,0))
I2 <- diag(c(0,1))


initf2 <- function(chain_id = 1) {
  list(q1 =Y[1,1]/rlnorm(1,3.5,0.0606) , q2 =Y[1,2]/rlnorm(1,3.5,0.0606) , q3=Y[1,3]/rlnorm(1,3.5,0.0606) ,
       r1 = rlnorm(1,-1.141,0.0606), r2 = rlnorm(1,-1.141,0.0606), K1 = rlnorm(1,3.5,0.0606), K2 = rlnorm(1,3.5,0.0606), x_tilde = array(rnorm(N),dim = c(N,2)), alpha = chain_id)
} 


n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))


fit <- stan(file = "/home/posmae/leafaros/Allcodigos/teste/carcara_01/modelo.stan", 
            data=list( N=N, y=Y, x0= t(x0),mu=t(mu),I = I, I1 = I1, I2 = I2,Um=t(Um), Lambda = Lambda), 
            iter=10000, chains=n_chains, init=init_ll,control = list(adapt_delta = 0.8,max_treedepth = 20))



Tabela <- summary(fit, pars=c("q1","q2","q3","r1","r2","K1","K2","Sigma"), probs = c(0.025,0.5, 0.975))$summary[,1]




