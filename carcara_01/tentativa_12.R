require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(LaplacesDemon)
require(MASS)
require(ggplot2)
require(gridExtra)
require(xtable)


#Estimando os estados, a matriz diagonal R e as matrizes Sigma, Lambda e Q


modelo <-"
data {
int N;
vector<lower=0>[3] y[N];
vector[2] x0[1];
vector[2] Um[1];
vector[2] mu[1];
matrix[2,2] I;
matrix[2,2] I1;
matrix[2,2] I2;
cov_matrix[2] Lambda;
}
parameters {
vector[2] x_tilde[N];
real<lower=0> r1;
real<lower=0> r2;
cov_matrix[3] Sigma;
real<lower=0> K1;
real<lower=0> K2;
real<lower=0> q1;
real<lower=0> q2;
real<lower=0> q3;

}
transformed parameters {
matrix[3,3] Q;
matrix[2,2] R;
matrix[2,2] K;
vector<lower=0>[2] x[N];


Q[1,1] = q1;
Q[1,2] = 0;
Q[1,3] = 0;
Q[2,1] = 0;
Q[2,2] = q2;
Q[2,3] = 0;
Q[3,1] = 0;
Q[3,2] = 0;
Q[3,3] = q3;

R[1,1] = r1;
R[1,2] = 0;
R[2,1] = 0;
R[2,2] = r2;

K[1,1] = 1/K1;
K[1,2] = 0;
K[2,1] = 0;
K[2,2] = 1/K2;


x[1] = x0[1] + x_tilde[1];


for( t in 2:N )
{
  
  x[t] =  (I + R)*x[t-1] - R*K*( I1*x[t-1]*x[t-1]'*I1 + I2*x[t-1]*x[t-1]'*I2)*Um[1] + x_tilde[t];
  
}

}
model {
vector[3] X[N];

r1 ~ lognormal(-1.141,0.0606);
r2 ~ lognormal(-1.141,0.0606);
q1 ~ inv_gamma(0.001,0.001); 
q2 ~ inv_gamma(0.001,0.001);
q3 ~ inv_gamma(0.001,0.001);
K1 ~ lognormal(3.5,0.1606);
K2 ~ lognormal(3.5,0.1606);
Sigma ~ inv_wishart(2.1, diag_matrix(rep_vector(1,3)));

for( t in 1:N )
{
  
  x_tilde[t] ~ multi_normal(mu, Lambda);
   X[t] = [x[t,1],x[t,1],x[t,2]]';
 
}

for( t in 1:N )
{
  
  y[t] ~ multi_normal(Q*X[t], Sigma);
  
}

}
"



N <-50
I <- diag(c(1,1))
I1 <- diag(c(1,0))
I2 <- diag(c(0,1))
Q <- matrix(c(9,0,0,0,0,5,0,0,0,0,15,0),ncol =4,byrow = T)
R <- diag(c(0.4,0.6))
K <- diag(c(0.025,0.016))
Sigma <- rinvwishart(3, diag(c(40^2,10^2,50^2)))
Lambda <- rinvwishart(2, diag(c(0.1^2,0.1^2)))


origX <- matrix(data = 0, nrow = 2, ncol = N)
Um <- matrix(c(1,1),ncol = 1)
X <- matrix(data = 0, nrow = 2, ncol = N)

x0 <- matrix(c(20, 30), nrow = 2, ncol = 1)


X[, 1] <- x0;
origX[, 1] <- x0;

for( i in 2:N )
{
  
  origX[, i] <- (I + R) %*% origX[, i - 1] - R%*%K%*%( I1%*%origX[, i - 1]%*%t(origX[, i - 1])%*%I1 +
                                                         I2%*%origX[, i - 1]%*%t(origX[, i - 1])%*%I2)%*%Um
  
}


for(i in 1:N )
{
  
  X[, i] <- mvrnorm(n = 1, origX[, i], Lambda)
  
}

Y <- matrix(0, nrow = 3, ncol = N)


for( i in 1:N )
{
  
  Y[, i] <- mvrnorm(n = 1,  Q%*%kronecker(X[, i], I)%*%Um, Sigma)
  
}

plot(Y[1,])
plot(Y[2,])
plot(Y[3,])


mu <- matrix(c(0, 0), nrow = 2, ncol = 1)


initf2 <- function(chain_id = 1) {
  list(q1 =Y[1,1]/rlnorm(1,3.5,0.0606) , q2 =Y[2,1]/rlnorm(1,3.5,0.0606) , q3=Y[3,1]/rlnorm(1,3.5,0.0606) ,
       r1 = rlnorm(1,-1.141,0.0606), r2 = rlnorm(1,-1.141,0.0606), K1 = rlnorm(1,3.5,0.0606), K2 = rlnorm(1,3.5,0.0606), x_tilde = array(rnorm(N),dim = c(N,2)), alpha = chain_id)
} 


n_chains <- 1
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))

fit <- stan(model_code = modelo, 
            data=list( N=N, y=t(Y), x0= t(x0),mu=t(mu), I = I, I1 = I1, I2 = I2,Um=t(Um), Lambda = Lambda), 
            iter=10000, chains=n_chains, init=init_ll,control = list(adapt_delta = 0.8,max_treedepth = 20))


fit 


mu_tau_summary <- summary(fit, pars=c("q1","q2","q3","r1","r2","K1","K2","Sigma"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary
