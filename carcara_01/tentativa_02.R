require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(MASS)
require(ggplot2)
require(gridExtra)
require(xtable)


modelo <-"
data {
int N;
vector[2] y[N];
vector[2] x0[1];
cov_matrix[2] Sigma;
cov_matrix[2] Q;
matrix[2,2] R;
}
transformed data {
matrix[2, 2] V_sqrt;
matrix[2, 2] Q_sqrt;
V_sqrt = cholesky_decompose(Sigma);
Q_sqrt = cholesky_decompose(Q);
}
parameters {
vector[2] x_tilde[N];
}
transformed parameters {

vector[2] x[N];


x[1] = x0[1] + Q_sqrt*x_tilde[1];

for( t in 2:N )
{
  
  x[t] = R*x[t-1] + Q_sqrt*x_tilde[t];
  
}

}
model {

for( t in 1:N )
{
  
  x_tilde[t] ~ std_normal();
  
}

for( t in 1:N )
{
  
  y[t] ~ multi_normal(x[t], Sigma);
  
}

}
"



N <- 50;
R <- matrix(c(0.6, 0, 0, 0.3), nrow = 2, ncol = 2, byrow = TRUE);
origX <- matrix(data = 0, nrow = 2, ncol = N);
X <- matrix(data = 0, nrow = 2, ncol = N);

x0 <- matrix(c(20, 10), nrow = 2, ncol = 1);

X[, 1] <- x0;
origX[, 1] <- x0;

for( i in 2:N )
{
  
  origX[, i] <- R %*% origX[, i - 1];
  
}

Q <- matrix(c(3, 0, 0,3), nrow = 2, ncol = 2, byrow = TRUE)

for(i in 1:N )
{
  
  X[, i] <- mvrnorm(n = 1, origX[, i], Q)
  
}

Y <- matrix(0, nrow = 2, ncol = N)
Sigma <- matrix(c(30^2, 0, 0, 30^2), nrow = 2, ncol = 2)

for( i in 1:N )
{
  
  Y[, i] <- mvrnorm(n = 1,  X[, i], Sigma)
  
}

fit <- stan(model_code = modelo, 
            data=list( N=N, y=t(Y), Sigma=Sigma, R=R, Q=Q, x0= t(x0)), 
            iter=10000, chains=4, control = list(adapt_delta = 0.8));


fit 
