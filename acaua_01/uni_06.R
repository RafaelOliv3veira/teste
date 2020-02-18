require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(ggplot2)
require(gridExtra)
require(xtable)

# Simulando

N <- 100
xt <- c()
yt <- c()

q <- 9
sigma <- 12
kappa <- 0.05 

r <- 0.4
K <- 40

xt[1] <- 20

for(t in 2:(N+1)){
  xt[t] <- xt[t-1] + r*xt[t-1]*(1-xt[t-1]/K) + rnorm(1,0,sqrt(kappa*sigma^2))
}

yt  <- q*xt[1:(N)] + rnorm(1,0,sigma)


modelo <-"
data{
int<lower=1> N; // definir N
real<lower=0> y[N]; // definir o vetor y
}
parameters{
//xt
vector[N] u_z;
real<lower=0> K; 
real<lower=0> r;
//yt
real<lower=0> q;
real<lower=0> sigma2;

}
transformed parameters {
real<lower=0> xt[N];
vector[N] xt_err = sqrt(0.05*sigma2)*u_z;

xt[1] = K + xt_err[1];
for (t in 2:N) {
xt[t] = xt[t-1] + r*xt[t-1]*(1-xt[t-1]/K) + xt_err[t] ;
}

}
model{ 

// Prioris
u_z ~ normal(0,1);
K ~ lognormal(3.5,0.1606);
r ~ lognormal(-1.141,0.1606);
sigma2 ~ inv_gamma(0.001,0.001); 
q ~ inv_gamma(0.001,0.001);

y ~ normal( xt, sqrt(sigma2));

}

"


data <- list(N = 100, y = yt )


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo, data = data, iter =4000, thin = 3,
               warmup = 2000, chains = 1)

# ajuste <- stan(model_code = modelo, data = data, iter = 40000, thin = 3,
#                warmup = 30000, chains = 3, control = list(adapt_delta = 0.8,max_treedepth = 20))

print(ajuste, digits = 3)
mu_tau_summary <- summary(ajuste, pars=c( "q","sigma2","K", "r"), probs = c(0.025,0.5, 0.975))$summary

xtable(print(mu_tau_summary, digits = 3), digits = 3)


traceplot(ajuste,  pars=c( "K", "r","q","sigma2"))

pairs(ajuste,   pars=c( "K", "r", "q","sigma2"))

stan_ac(ajuste,   pars=c( "K", "r","q","sigma2"))