
data {
  int N_individuals;
  int N_d_predictors;
  int max_time;
  
  array[N_individuals] int init_times;
  array[N_individuals,max_time,N_d_predictors] int d_predictors;
}

parameters {
  array[N_d_predictors] real<lower=-0.5,upper=0.5> age_decay;
  array[N_d_predictors] real<lower=-0.5,upper=0.5> base_decay;
  array[N_d_predictors] real<lower=0,upper=5> p;
}

model {
  array[N_d_predictors,max_time,N_individuals] real risk;
  for (dpred in 1:N_d_predictors){
    for (time in 1:max_time){
      for (individual in 1:N_individuals){
        risk[dpred,time,individual] = tanh(exp(log(p[dpred]) + (time - init_times[individual]) * age_decay[dpred] + time * base_decay[dpred]));
        if (d_predictors[individual,time,dpred] > -1){
          d_predictors[individual,time,dpred] ~ bernoulli(risk[dpred,time,individual]);
        }
      }
    }
  }
}
