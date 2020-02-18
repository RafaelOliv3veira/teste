require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(ggplot2)
require(gridExtra)
require(xtable)

rm(list = ls())

N <-50

q <- 10 ; sigma <- 16
mu_0 <- 20; r <- 0.3; K <- 30; tau <- sqrt(16^2*0.012)


x <- rnorm(1 ,mu_0, tau)
y <- rnorm(1, q*x, sigma)

for(i in 2:N){
  x[i] <- rnorm(1, x[i-1] +  r*x[i-1]*(1-x[i-1]/K), tau)
  y[i] <- rnorm(1, q*x[i], sigma)
}



plot(y)


modelo0 <- "
data {
int <lower=0> N;
real y[N];
}

parameters {
vector[N] u_z;
real <lower=0> sigma2;
real <lower=0> q;
real <lower=0.1> K;
real <lower=0> r;
}

transformed parameters {
vector[N] u; //Level
vector[N] u_err = sqrt(sigma2)*0.012* u_z;

u[1] = u_err[1] + 20;
for (t in 2:N) {
u[t] = u[t-1] + r*u[t-1]*(1-u[t-1]/K) + u_err[t] ;
}
}

model {
u_z ~ normal(0,1);
r ~ lognormal(-1.141,0.0606);
K ~ lognormal(3.5,0.0606);
sigma2 ~ inv_gamma(0.001,0.001);
q ~ inv_gamma(0.001,0.001);

for (i in 1:N)
y[i] ~ normal( q*u[i], sqrt(sigma2));

}
"



data <- list(N = N, y = y)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo0, data = data, iter = 40000, thin=10,
               warmup = 10000, chains =3, control = list(adapt_delta = 0.8,max_treedepth = 20))

mu_tau_summary <- summary(ajuste, pars=c( "q","r","K","sigma2"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary


traceplot(ajuste,  pars=c( "K", "r","q","sigma2"))

pairs(ajuste,   pars=c( "K", "r", "q","sigma2"))

stan_ac(ajuste,   pars=c( "K", "r","q","sigma2"))

