data {
  int<lower=1> nTime;
  int counts[nTime];
  int group2[nTime];
  int group3[nTime];
  int nCounts[nTime];
  int nGroup2[nTime];
  int nGroup3[nTime];
}

parameters{
  real meanWeek1;
  vector[nTime-1] changes;
  real<lower=0> changeSigma;
  real group2Change;
  real group3Change;
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
  group2 ~ binomial_logit(nGroup2,means+group2Change);
  group3 ~ binomial_logit(nGroup3,means+group3Change);
  changeSigma~gamma(1,2);
  group2Change~double_exponential(0,1);
  group3Change~double_exponential(0,1);
}
