//(Logistic) Rasch Model with #Rows = 1
data {
  int<lower=1> K; //Number of Compounds
  int<lower=0,upper=1> x[K];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'alpha' and 'bta'.
parameters {
  real alpha;
  vector[K] bta;
}

model {
  alpha ~ std_normal();
  bta ~ std_normal();
    x ~ bernoulli_logit(alpha + bta);
}
