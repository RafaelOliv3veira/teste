require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(MASS)
require(ggplot2)
require(gridExtra)
require(xtable)


#Estimando os estados e a matriz diagonal R as vezes fica bom e outras nao

modelo <-"
data {
int N;
vector[2] y[N];
vector[2] x0[1];
cov_matrix[2] Sigma;
cov_matrix[2] Lambda;
vector[2] Um[1];
matrix[2,2] Q;
matrix[2,2] K;
matrix[2,2] I;
matrix[2,2] I1;
matrix[2,2] I2;
}
transformed data {
matrix[2, 2] Lambda_sqrt;
Lambda_sqrt = cholesky_decompose(Lambda);
}
parameters {
vector[2] x_tilde[N];
real<lower=0> r1;
real<lower=0> r2;

}
transformed parameters {
matrix[2,2] R;
vector[2] x[N];


R[1,1] = r1;
R[1,2] = 0;
R[2,1] = 0;
R[2,2] = r2;



x[1] = x0[1] + Lambda_sqrt*x_tilde[1];

for( t in 2:N )
{
  
  x[t] =  (I + R)*x[t-1] - R*K*( I1*x[t-1]*x[t-1]'*I1 + I2*x[t-1]*x[t-1]'*I2)*Um[1] + Lambda_sqrt*x_tilde[t];
  
}

}
model {
r1 ~ lognormal(-1.141,0.1606);
r2 ~ lognormal(-1.141,0.1606);



for( t in 1:N )
{
  
  x_tilde[t] ~ std_normal();
  
}

for( t in 1:N )
{
  
  y[t] ~ multi_normal(Q*x[t], Sigma);
  
}

}
"



N <- 50;
I <- diag(c(1,1))
I1 <- diag(c(1,0))
I2 <- diag(c(0,1))
Q <- diag(c(9,9))
R <- diag(c(0.4,0.4))
K <- diag(c(0.025,0.025))
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

Lambda <- matrix(c(3, 0, 0,3), nrow = 2, ncol = 2, byrow = TRUE)

for(i in 1:N )
{
  
  X[, i] <- mvrnorm(n = 1, origX[, i], Lambda)
  
}

Y <- matrix(0, nrow = 2, ncol = N)
Sigma <- matrix(c(30^2, 0, 0, 30^2), nrow = 2, ncol = 2)

for( i in 1:N )
{
  
  Y[, i] <- mvrnorm(n = 1,  Q%*%X[, i], Sigma)
  
}

plot(Y[1,])

fit <- stan(model_code = modelo, 
            data=list( N=N, y=t(Y), Sigma=Sigma,Lambda=Lambda, x0= t(x0), Q = Q, K = K, I = I, I1 = I1, I2 = I2,Um=t(Um)), 
            iter=10000, chains=4, control = list(adapt_delta = 0.8));


fit 


mu_tau_summary <- summary(fit, pars=c( "r1","r2"), probs = c(0.025,0.5, 0.975))$summary
mu_tau_summary
