require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(ggplot2)
require(gridExtra)
require(xtable)


t <- 25; mu_0 <- 1
x <- rnorm(1 ,mu_0, 1)
y <- rnorm(1, x, 1)
a <- 0.5; b <- 1; c <- 1
for(i in 2:t){
  x[i] <- rnorm(1, x[i-1] * a + b, 1)
  y[i] <- rnorm(1, x[i] * c, 1)
}
     

plot(y)


modelo0 <- "
data {
int <lower=0> N;
real y[N];
}

parameters {
vector[N] u_z;
real <lower=0> s_obs;
real <lower=0, upper=1> a;
real b;
}

transformed parameters {
real u[N]; //Level
vector[N] u_err = s_obs * u_z;

u[1] = u_err[1];
for (t in 2:N) {
u[t] = a * u[t-1] + b + u_err[t] ;
}
}

model {
u_z ~ normal(0,1);
a ~ normal(0, 1);
b ~ normal(0, 2.5);
s_obs ~ gamma(2, 1);
y ~ normal(u, s_obs);
}
"



data <- list(N = 25, y = y)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo0, data = data, iter = 10000, thin=10,
               warmup = 5000, chains = 1 )

mu_tau_summary <- summary(ajuste, pars=c( "a","b","s_obs"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary
