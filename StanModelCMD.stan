data{
  int N_individuals;//number of individuals
  int N_s_predictors;//number of static predictors
  int N_d_predictors;//number of dynamic predictors
  int max_time;//largest time to predict forwards
    
  array[N_individuals] int init_times;//gestational ages at birth
  array[N_individuals,N_s_predictors] real s_predictors;//static predictor information
  array[N_individuals,max_time,N_d_predictors] int d_predictors;//dynamical predictor information
  
  array[N_individuals,max_time] int clin_outcomes;//clinical outcomes
}
parameters {
  array[4,5] real<lower=0> p;//a in paper
  array[4,5] real base_par_1;//b in paper
  array[4,5] real base_par_2;//theta in paper


  array[N_s_predictors,4,5] real static_coefs;//lambda in paper
  array[4,5] real<upper=0> static_decay;//psi in paper

  array[N_d_predictors,4,5] real dynam_coefs;//gamma in paper
  array[N_d_predictors,4,5] real<upper=0> dynam_decay;//omega in paper

  array[4,5] real age_coefs;//c in paper
  array[4,5] real<upper=0> age_decay;//d in paper

  //array[3,4,5] real prev_state_coefs;


}
model{
   //////////////////////////////////////////////////////First we find parameters for dynamical state predictions
  array[5,max_time,5] real base_propensity;//base deterioration propensity, alpha in paper

  array[5,5] real static_multiplier;//phi in paper
  array[5,5] real dynamic_multiplier;//kappa in paper
  array[5,5] real age_multiplier;//phi in paper

  
  array[max_time,5,5] real transition;//transition probabilities, P in paper

  array[max_time,5] real state_prob;//probability of belonging to a state, corrected if it is known (D in paper)
  array[max_time,5] real state_prob_guess;//probability of belonging to a state, uncorrected if it is known

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
    for (clin_outcome2 in 1:5){//for each possible clinical outcome, check if the individual was that on the first day
      //then fill in with a 1 if it was
      //fill with a zero otherwise
      if (clin_outcomes[individual,init_times[individual]] == clin_outcome2){
        state_prob[init_times[individual],clin_outcome2] = 1;
        state_prob_guess[init_times[individual],clin_outcome2] = 1;
      } else {
        state_prob[init_times[individual],clin_outcome2] = 0;
        state_prob_guess[init_times[individual],clin_outcome2] = 0;
      }
      //initialise all values after the first day with 0
      for (time in (init_times[individual]+1):max_time){
        state_prob[time,clin_outcome2] = 0;
        state_prob_guess[time,clin_outcome2] = 0;
      }
    }
    
    //now for every future time point
    for (time in (init_times[individual]):max_time){
      
      //calculate the transition probabilities from each state clin_outcome to clin_outcome2
      //apart from state 5, which is death and can't be transitioned out of (this will be done later)
      for (clin_outcome in 1:4){
        for (clin_outcome2 in 1:5){
          
          //calculate multipliers
          static_multiplier[clin_outcome,clin_outcome2] = exp(exp((time-init_times[individual])* static_decay[clin_outcome,clin_outcome2]) * dot_product(to_vector(s_predictors[individual,]),to_vector(static_coefs[,clin_outcome,clin_outcome2]) ));
          dynamic_multiplier[clin_outcome,clin_outcome2] = 1;
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
      
      //also input probabilities of transitioning from death
      transition[time,5,1] = 0;
      transition[time,5,2] = 0;
      transition[time,5,3] = 0;
      transition[time,5,4] = 0;
      transition[time,5,5] = 1;
      
      
      //now calculate next time point forward using transition probabilities
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
            state_prob_guess[time,clin_outcome] = max([1e-10,state_prob[time,clin_outcome]]);//There can be some small computational errors, so we fix this here
            
          }
            //Add to likelihood
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




