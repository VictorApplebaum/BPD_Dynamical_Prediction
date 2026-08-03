
data {
  int N_individuals;//number of individuals
  int N_d_predictors;//number of predictors
  int max_time;//furthest time point to be predicted
  
  array[N_individuals] int init_times;//time at which baby was born (i.e. gestational age at birth)
  array[N_individuals,max_time,N_d_predictors] int d_predictors;//array of predictors
}

parameters {
  array[N_d_predictors] real<lower=-0.5,upper=0.5> age_decay;//x in paper
  array[N_d_predictors] real<lower=-0.5,upper=0.5> base_decay;//y in paper
  array[N_d_predictors] real<lower=0,upper=5> p;//z in paper
}

model {
  array[N_d_predictors,max_time,N_individuals] real risk;//value representing the probability of a predictor being positive
  for (dpred in 1:N_d_predictors){
    for (time in 1:max_time){
      for (individual in 1:N_individuals){
        //risk as defined in paper
        risk[dpred,time,individual] = tanh(exp(log(p[dpred]) + (time - init_times[individual]) * age_decay[dpred] + time * base_decay[dpred]));
        //fit with a bernoulli distribution
        if (d_predictors[individual,time,dpred] > -1){//-1 signifies missing data, so skip if this value is present
          d_predictors[individual,time,dpred] ~ bernoulli(risk[dpred,time,individual]);
        }
      }
    }
  }
}
