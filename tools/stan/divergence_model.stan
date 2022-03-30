data {
	int<lower=0> N;
	real k[N];
	real l[N];
}
parameters {
	real<lower=0,upper=0.2> theta;
	real<lower=0,upper=0.2> gamma;
}
model {
target += uniform_lpdf(theta | 0, 0.2);	
target += uniform_lpdf(gamma | 0, 0.2);	
for (n in 1:N) {
	target += k[n]*log(l[n]*theta);
	target += gamma/theta;
	target += -(k[n]+1)*log(l[n]*theta+1);
	target += log(gamma_q(k[n]+1,l[n]*gamma+(gamma/theta)));
}
}
generated quantities {
  real<lower=0> ratio;
  ratio = gamma/theta;
}
