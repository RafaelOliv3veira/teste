require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(ggplot2)
require(gridExtra)
require(xtable)

# Simulando

N <- 150
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


plot(yt)

modelo <-"
data{
int<lower=1> N; // definir N
real<lower=0> y[N]; // definir o vetor y
}
parameters{
//xt
real<lower=0> K; 
real<lower=0> r;
real<lower=0> xt[N];
//yt
real<lower=0> q;
real<lower=0> sigma2;

}
transformed parameters {
real<lower=0> tau;

tau = sqrt(0.05*sigma2);
}
model{ 

// Prioris
K ~ lognormal(3.5,0.1606);
r ~ lognormal(-1.141,0.1606);
sigma2 ~ inv_gamma(0.001,0.001); 
q ~ inv_gamma(0.001,0.001);



// Espaco estado
xt[1] ~ normal(K, tau); 


for (i in 2:N)
xt[i] ~ normal( xt[i-1] + r*xt[i-1]*(1-xt[i-1]/K), tau) ;      


for (i in 1:N)
y[i] ~ normal( q*xt[i], sqrt(sigma2));

}

"


data <- list(N = N, y = yt )


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo, data = data, iter = 40000, thin=10,
               warmup = 10000, chains =3, control = list(adapt_delta = 0.8,max_treedepth = 20))

# ajuste <- stan(model_code = modelo, data = data, iter = 100000, thin = 3,
#                warmup = 30000, chains = 3,control = list(adapt_delta = 0.99,max_treedepth = 15))

# ajuste <- stan(model_code = modelo, data = data, iter = 40000, thin = 3,
#                warmup = 30000, chains = 3, control = list(adapt_delta = 0.8,max_treedepth = 20))

# print(ajuste, digits = 3)
mu_tau_summary <- summary(ajuste, pars=c( "q","sigma2","K", "r"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary


xtable(print(mu_tau_summary, digits = 3), digits = 3)


traceplot(ajuste,  pars=c( "K", "r","q","sigma2"))

pairs(ajuste,   pars=c( "K", "r", "q","sigma2"))

stan_ac(ajuste,   pars=c( "K", "r","q","sigma2"))