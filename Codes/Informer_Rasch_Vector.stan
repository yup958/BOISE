//
// Posterior sample for informer set. Input some block xb of original x0, and informer data xA,
// Output posterior mean

data {
  int<lower=1> K; //Number of Compounds
  int<lower=1> L; //Number of Informer Compounds (Should greater than 1)
  int<lower=1> A[L];//Informer Index
  int<lower=0,upper=1> x[K];//Original Block Data
  int<lower=0, upper=1> y[L];//Informer Data
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
  x ~ bernoulli_logit(alpha + bta);
  for (l in 1:L)
    y[l] ~bernoulli_logit(alpha + bta[A[l]]);
}
