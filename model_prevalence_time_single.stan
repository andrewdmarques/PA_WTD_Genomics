data {
  int<lower=1> nTime;
  int counts[nTime];
  int nCounts[nTime];
}

parameters{
  real meanWeek1;
  vector[nTime-1] changes;
  real<lower=0> changeSigma;
}
transformed parameters{
  vector[nTime] means;
  means[1]=meanWeek1;
  for(ii in 2:nTime) means[ii]=means[ii-1]+changes[ii-1]*changeSigma;
}


model {
  meanWeek1 ~ normal(0,10);
  changes[1] ~ normal(0,1);
  for(ii in 2:(nTime-1)) changes[ii]~normal(changes[ii-1],1);
  counts ~ binomial_logit(nCounts,means);
  changeSigma~gamma(1,2);
}