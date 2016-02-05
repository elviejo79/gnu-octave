## covariance matrix

v = [1 .5 0; .5 1 0; 0 0 1];

## Cholesky decomposition of
## positive semi-definite covariance
## obtained by multipling the matrix by its transpose
## and must have the Cholesky decomposition applied.
c = chol(v*v');

##time points
n = 1000;
d = randn(n,size(v,1))*c;

## in the code above the matrix de contains a 10,000 X 3 matrix of random numbers
## such that the first two are correlated around 0.8
## while the third  is uncorrelated to the first two.

##note that
cov(d)
## is very similar to
v*v'
