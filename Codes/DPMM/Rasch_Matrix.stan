//(Logistic) Rasch Model with #Rows > 1
data {
  int<lower=1> J; //Number of Targets
  int<lower=1> K; //Number of Compounds
  int<lower=0,upper=1> x[J,K];
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
}
