require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(ggplot2)
require(gridExtra)
require(xtable)

rm(list = ls())


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


modelo0 <- "
data {
int <lower=0> N;
real y[N];
}

parameters {
vector[N] u_z;
real <lower=0> sigma2;
real <lower=0> r;
}

transformed parameters {
vector[N] u; //Level
vector[N] u_err = sqrt(sigma2)*u_z;

u[1] = u_err[1] + 2;
for (t in 2:N) {
u[t] =  r*u[t-1] + u_err[t] ;
}
}

model {
u_z ~ normal(0,1);
r ~ lognormal(-1.141,0.0606);
sigma2 ~ inv_gamma(0.001,0.001);

target += normal_lpdf(y | u,sqrt(sigma2));

}
"



data <- list(N = N, y = y)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo0, data = data, iter = 40000, thin=10,
               warmup = 10000, chains =3, control = list(adapt_delta = 0.8,max_treedepth = 20))

mu_tau_summary <- summary(ajuste, pars=c( "r", "sigma2"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary


traceplot(ajuste,  pars=c( "r","sigma2"))

pairs(ajuste,   pars=c( "r","sigma2"))

stan_ac(ajuste,   pars=c( "r","sigma2"))
