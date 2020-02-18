require(gridExtra)
require(xtable)

# Simulando

N <- 100
xt <- c()
yt <- c()

q <- 9
sigma <- 2
kappa <- 0.25
tau <- sqrt(0.25*4)

r <- 0.4


xt[1] <- 20

for(t in 2:(N+1)){
  xt[t] <- r*xt[t-1]  + rnorm(1,0,tau)
}

yt  <- q*xt[1:(N)] 


plot(yt)

modelo <-"
data{
int<lower=1> N; // definir N
vector[N] y; // definir o vetor y
}
parameters{
//xt
real<lower=0> r;
vector[N+1] xt;
//yt
real<lower=0> q;
real<lower=0> sigma;

}
transformed parameters {
real<lower=0> taux;
real<lower=0> theta;

theta = sqrt(sigma);
taux = sqrt(sigma)/2;
}
model{ 

r ~ lognormal(-1.141,0.0606);
sigma ~ inv_gamma(0.001,0.001); 
q ~ inv_gamma(0.001,0.001);



// Espaco estado
xt[1] ~ normal(10, theta); 


for (i in 2:N)
xt[i] ~ normal(r*xt[i-1], taux) ;      


for (i in 1:N)
y[i] ~ normal( q*xt[i], sqrt(sigma));

}

"


data <- list(N = 100, y = yt )


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ajuste <- stan(model_code = modelo, data = data, iter = 100000, thin = 10,
               warmup = 70000, chains = 3, control = list(adapt_delta = 0.8,max_treedepth = 20))

# ajuste <- stan(model_code = modelo, data = data, iter = 40000, thin = 3,
#                warmup = 30000, chains = 3, control = list(adapt_delta = 0.8,max_treedepth = 20))

print(ajuste, digits = 3)

mu_tau_summary <- summary(ajuste, pars=c( "q","sigma", "r"), probs = c(0.025,0.5, 0.975))$summary

xtable(print(mu_tau_summary, digits = 3), digits = 3)


traceplot(ajuste,  pars=c( "K", "r","q","sigma"))

pairs(ajuste,   pars=c( "K", "r", "q","sigma"))

stan_ac(ajuste,   pars=c( "K", "r","q","sigma"))