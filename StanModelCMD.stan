data{
  int N_individuals;
  int N_s_predictors;
  int N_d_predictors;
  int max_time;
    
  array[N_individuals] int init_times;
  array[N_individuals,N_s_predictors] real s_predictors;
  array[N_individuals,max_time,N_d_predictors] int d_predictors;
  
  array[N_individuals,max_time] int clin_outcomes;
}
parameters {
      //normalizer p: contains vector: [prob it eventually deteriorates, prob it remains forever, prob it improves]
  array[4,5] real<lower=0> p;
  array[4,5] real base_par_1;//first parameter in weibull distr
  array[4,5] real base_par_2;//second parameter in weibull distr


  array[N_s_predictors,4,5] real static_coefs;
  array[4,5] real<upper=0> static_decay;

  array[N_d_predictors,4,5] real dynam_coefs;
  array[N_d_predictors,4,5] real<upper=0> dynam_decay;

  array[4,5] real age_coefs;
  array[4,5] real<upper=0> age_decay;

  //array[3,4,5] real prev_state_coefs;


}
model{
   //////////////////////////////////////////////////////First we find parameters for dynamical state predictions
  array[5,max_time,5] real base_propensity;//base deterioration propensity
  //real base_propensity_r[5,max_time];//base remain propensity
  
  array[5,5] real static_multiplier;
  array[5,5] real dynamic_multiplier;
  array[5,5] real age_multiplier;
  array[5,5] real prev_state_multiplier;

  
  //we then combine with the split probabilities to find base transition probs
  array[max_time,5,5] real transition;

  array[max_time,5] real state_prob;
  array[max_time,5] real state_prob_guess;

  //calculate base improvement and deterioration probabilities
  for (time in 1:max_time){
    for (clin_outcome in 1:4){
      for (clin_outcome2 in 1:5){
        base_propensity[clin_outcome,time,clin_outcome2] = p[clin_outcome,clin_outcome2] * exp(-(time/base_par_1[clin_outcome,clin_outcome2])^base_par_2[clin_outcome,clin_outcome2]);
      }
    }
  }
  

  

  

  //run chain
  for (individual in 1:N_individuals){
    
    //initialize chain
    for (clin_outcome2 in 1:5){
      if (clin_outcomes[individual,init_times[individual]] == clin_outcome2){
        state_prob[init_times[individual],clin_outcome2] = 1;
        state_prob_guess[init_times[individual],clin_outcome2] = 1;
      } else {
        state_prob[init_times[individual],clin_outcome2] = 0;
        state_prob_guess[init_times[individual],clin_outcome2] = 0;
      }
      
      for (time in (init_times[individual]+1):max_time){
        state_prob[time,clin_outcome2] = 0;
        state_prob_guess[time,clin_outcome2] = 0;
      }
    }
    
    for (time in (init_times[individual]):max_time){
      

      for (clin_outcome in 1:4){
        for (clin_outcome2 in 1:5){
          
          static_multiplier[clin_outcome,clin_outcome2] = exp(exp((time-init_times[individual])* static_decay[clin_outcome,clin_outcome2]) * dot_product(to_vector(s_predictors[individual,]),to_vector(static_coefs[,clin_outcome,clin_outcome2]) ));// .* to_vector(exp((time-init_times[individual])*to_vector(static_decay[,clin_outcome,clin_outcome2])))));
          dynamic_multiplier[clin_outcome,clin_outcome2] = 1;
          //prev_state_multiplier[clin_outcome,clin_outcome2] = 1;
          for (d_pred in 1:N_d_predictors){
            if (d_predictors[individual,time,d_pred] != -1){
              dynamic_multiplier[clin_outcome,clin_outcome2] = dynamic_multiplier[clin_outcome,clin_outcome2] * exp(dynam_coefs[d_pred,clin_outcome,clin_outcome2]*exp(d_predictors[individual,time,d_pred]*dynam_decay[d_pred,clin_outcome,clin_outcome2]));
            }
          }
          age_multiplier[clin_outcome,clin_outcome2] = exp(age_coefs[clin_outcome,clin_outcome2] * exp((time-init_times[individual])*age_decay[clin_outcome,clin_outcome2]));
        }
      }

   
      //calculate transitions
      for (clin_outcome in 1:1){
        transition[time,clin_outcome,2] = tanh(base_propensity[clin_outcome,time,2] * static_multiplier[clin_outcome,2] * dynamic_multiplier[clin_outcome,2] * age_multiplier[clin_outcome,2]);
        transition[time,clin_outcome,3] = (1 - transition[time,clin_outcome,2]) * tanh(base_propensity[clin_outcome,time,3] * static_multiplier[clin_outcome,3] * dynamic_multiplier[clin_outcome,3] * age_multiplier[clin_outcome,3]);
        transition[time,clin_outcome,4] = (1 - transition[time,clin_outcome,2] - transition[time,clin_outcome,3]) * tanh(base_propensity[clin_outcome,time,4] * static_multiplier[clin_outcome,4] * dynamic_multiplier[clin_outcome,4] * age_multiplier[clin_outcome,4]);
        transition[time,clin_outcome,5] = (1 - transition[time,clin_outcome,2] - transition[time,clin_outcome,3] - transition[time,clin_outcome,4]) * tanh(base_propensity[clin_outcome,time,5] * static_multiplier[clin_outcome,5] * dynamic_multiplier[clin_outcome,5] * age_multiplier[clin_outcome,5]);
        transition[time,clin_outcome,1] = 1 - transition[time,clin_outcome,2]- transition[time,clin_outcome,3]- transition[time,clin_outcome,4]- transition[time,clin_outcome,5];
      }
      
      for (clin_outcome in 2:2){
        transition[time,clin_outcome,1] = tanh(base_propensity[clin_outcome,time,1] * static_multiplier[clin_outcome,1] * dynamic_multiplier[clin_outcome,1] * age_multiplier[clin_outcome,1]);
        transition[time,clin_outcome,3] = (1 - transition[time,clin_outcome,1]) * tanh(base_propensity[clin_outcome,time,3] * static_multiplier[clin_outcome,3] * dynamic_multiplier[clin_outcome,3] * age_multiplier[clin_outcome,3]);
        transition[time,clin_outcome,4] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,3]) * tanh(base_propensity[clin_outcome,time,4] * static_multiplier[clin_outcome,4] * dynamic_multiplier[clin_outcome,4] * age_multiplier[clin_outcome,4]);
        transition[time,clin_outcome,5] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,3] - transition[time,clin_outcome,4]) * tanh(base_propensity[clin_outcome,time,5] * static_multiplier[clin_outcome,5] * dynamic_multiplier[clin_outcome,5] * age_multiplier[clin_outcome,5]);
        transition[time,clin_outcome,2] = 1 - transition[time,clin_outcome,1] - sum(transition[time,clin_outcome,3:5]);
      }
      
      for (clin_outcome in 3:3){
        transition[time,clin_outcome,1] = tanh(base_propensity[clin_outcome,time,1] * static_multiplier[clin_outcome,1] * dynamic_multiplier[clin_outcome,1] * age_multiplier[clin_outcome,1]);
        transition[time,clin_outcome,2] = (1 - transition[time,clin_outcome,1]) * tanh(base_propensity[clin_outcome,time,2] * static_multiplier[clin_outcome,2] * dynamic_multiplier[clin_outcome,2] * age_multiplier[clin_outcome,2]);
        transition[time,clin_outcome,4] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2]) * tanh(base_propensity[clin_outcome,time,4] * static_multiplier[clin_outcome,4] * dynamic_multiplier[clin_outcome,4] * age_multiplier[clin_outcome,4]);
        transition[time,clin_outcome,5] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2] - transition[time,clin_outcome,4]) * tanh(base_propensity[clin_outcome,time,5] * static_multiplier[clin_outcome,5] * dynamic_multiplier[clin_outcome,5] * age_multiplier[clin_outcome,5]);
        transition[time,clin_outcome,3] = 1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2]- transition[time,clin_outcome,4]- transition[time,clin_outcome,5];
      }
      
      for (clin_outcome in 4:4){
        transition[time,clin_outcome,1] = tanh(base_propensity[clin_outcome,time,1] * static_multiplier[clin_outcome,1] * dynamic_multiplier[clin_outcome,1] * age_multiplier[clin_outcome,1]);
        transition[time,clin_outcome,2] = (1 - transition[time,clin_outcome,1]) * tanh(base_propensity[clin_outcome,time,2] * static_multiplier[clin_outcome,2] * dynamic_multiplier[clin_outcome,2] * age_multiplier[clin_outcome,2]);
        transition[time,clin_outcome,3] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2]) * tanh(base_propensity[clin_outcome,time,3] * static_multiplier[clin_outcome,3] * dynamic_multiplier[clin_outcome,3] * age_multiplier[clin_outcome,3]);
        transition[time,clin_outcome,5] = (1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2] - transition[time,clin_outcome,3]) * tanh(base_propensity[clin_outcome,time,5] * static_multiplier[clin_outcome,5] * dynamic_multiplier[clin_outcome,5] * age_multiplier[clin_outcome,5]);
        transition[time,clin_outcome,4] = 1 - transition[time,clin_outcome,1] - transition[time,clin_outcome,2] - transition[time,clin_outcome,3] - transition[time,clin_outcome,5];
      }
      
      transition[time,5,1] = 0;
      transition[time,5,2] = 0;
      transition[time,5,3] = 0;
      transition[time,5,4] = 0;
      transition[time,5,5] = 1;
      
      if (time > init_times[individual]){
        for (clin_outcome in 1:5){
          for (clin_outcome2 in 1:5){
            state_prob[time,clin_outcome2] = state_prob[time,clin_outcome2] + state_prob[time-1,clin_outcome] * transition[time-1,clin_outcome,clin_outcome2];
          }
        }
      }
      
      
      
      //Now if we are at a known point - record guess and correct

        if (clin_outcomes[individual,time] != -1){
          for (clin_outcome in 1:5){
            state_prob_guess[time,clin_outcome] = max([1e-10,state_prob[time,clin_outcome]]);
            //There can be some small computational errors leading to some probabilities being below zero, so we take max
          }
            //Also add to likelihood
            if (time > init_times[individual]){
              clin_outcomes[individual,time] ~ categorical(to_vector(state_prob_guess[time,])/sum(state_prob_guess[time,]));
            }

            
            //Also condition result if we are on a training individual, or if we are below threshold for test individuals
            for (clin_outcome in 1:5){
              if (clin_outcomes[individual,time] == clin_outcome){
                state_prob[time,clin_outcome] = 1;
              } else {
                state_prob[time,clin_outcome] = 0;
              }
            }
        }

      }
  }

}




