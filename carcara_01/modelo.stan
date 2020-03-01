//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
functions {
 vector kronecker_prod(vector A, matrix B, vector D) {
  matrix[2* rows(B),  cols(B)] C;
  vector[4] U;
  int m;
  int n;
  int p;
  int q;
  m = 2;
  n = 1;
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start =  1;
      col_end = q;
      C[row_start:row_end, col_start:col_end] = A[i] * B;
    
  }
  
  U = C*D;
  return U;
}
}

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
matrix[3,4] Q;
matrix[2,2] R;
matrix[2,2] K;
vector<lower=0>[2] x[N];


Q[1,1] = q1;
Q[1,2] = 0;
Q[1,3] = 0;
Q[1,4] = 0;
Q[2,1] = 0;
Q[2,2] = q2;
Q[2,3] = 0;
Q[2,4] = 0;
Q[3,1] = 0;
Q[3,2] = 0;
Q[3,3] = q3;
Q[3,4] = 0;

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

r1 ~ lognormal(-1.141,0.0606);
r2 ~ lognormal(-1.141,0.0606);
q1 ~ inv_gamma(0.001,0.001); 
q2 ~ inv_gamma(0.001,0.001);
q3 ~ inv_gamma(0.001,0.001);
K1 ~ lognormal(3.5,0.1606);
K2 ~ lognormal(3.5,0.1606);
Sigma ~ inv_wishart(4, diag_matrix(rep_vector(1,3)));

for( t in 1:N )
{
  
  x_tilde[t] ~ multi_normal(mu, Lambda);
  
 
}

for( t in 1:N )
{
  
  y[t] ~ multi_normal(Q*kronecker_prod(x[t],I,Um[1]), Sigma);
  
}

}

