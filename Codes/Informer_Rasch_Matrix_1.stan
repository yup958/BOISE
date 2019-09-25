//
// Posterior sample for informer set with only 1 informer feature. 
// Input some block xb of original x0, and informer data xA,
// Output posterior mean

data {
  int<lower=1> J; //Number of Targets in the Block
  int<lower=1> K; //Number of Compounds
  int<lower=1> A;//Informer Index
  int<lower=0,upper=1> x[J,K];//Original Block Data
  int<lower=0, upper=1> y;//Informer Data
}

// The parameters accepted by the model. Our model
// accepts two parameters 'alpha' (number) and 'bta' (vector).
parameters {
  real alpha;
  vector[K] bta;
}

model {
  alpha ~ std_normal();
  bta ~ std_normal();
  for (j in 1:J)
    x[j] ~ bernoulli_logit(alpha + bta);
  y ~ bernoulli_logit(alpha + bta[A]);
}
