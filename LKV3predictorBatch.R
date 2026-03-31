

#for (SeenTimePoints in c(28,1,14,7,21,3)){
fitter <- function(SeenTimePoints){
  #Libraries
  library(geometry)
  library (multiROC)
  library (gtools)
  library(dummies)
  library(ggplot2)
  library(pracma)
  
  
  
  print('Loading Data')
  #select data
  setwd("~/281224")#CHANGE TO RELEVANT DIRECTORY
  #load("~/281224/numbers.RData")
  bootstrap_nums <- 1:200
  #load("Code/DatasetsFull141024.RData")
  #datafull <- datafull[datafull$randomcohort == 'testing',]
  load('Code/DatasetTestFull300125.RData')
  #datafull <- datafull[1:200,]
  data = datafull
  rm(datafull)
  #shunt death time along 1 as we start counting days from 1 rather than zero
  data$deathageday <- data$deathageday + 1
  #remove babies who died before SeenTimePoints
  data <- data[(data$deathageday > SeenTimePoints) | is.na(data$deathageday),]
  max_time <- length(161:252)
  N_individuals <- dim(data)[1]
  pr_bootstrap <- list(None = array(NA,c(N_individuals,max_time)),
                        Oxygen = array(NA,c(N_individuals,max_time)),
                        NIV = array(NA,c(N_individuals,max_time)),
                        Inv = array(NA,c(N_individuals,max_time)),
                        Death = array(NA,c(N_individuals,max_time)))
  
  for (bootstrap in 1:length(bootstrap_nums)){
    set.seed(bootstrap)
    print(bootstrap)
    SaveResults <- TRUE
    
    
    #load("Code/DatasetsFull141024.RData")#CHANGE TO RELEVANT FILE
    #datafull <- datafull[datafull$randomcohort == 'testing',]
    load('Code/DatasetTestFull300125.RData')
    #datafull <- datafull[1:200,]
    data = datafull
    rm(datafull)
    
    
  
    
    
    
    #shunt death time along 1 as we start counting days from 1 rather than zero
    data$deathageday <- data$deathageday + 1
    #remove babies who died before SeenTimePoints
    data <- data[(data$deathageday > SeenTimePoints) | is.na(data$deathageday),]
    
    N_individuals <- dim(data)[1]
    
    #################Extract Static Copredictors######################
    
    #convert gestational age (gest) from weeks to days
    data$gest <- round(7*data$gest)
    
    days_monitored <- 100
    
    data$gender[data$gender == 'male'] <- 1
    data$gender[data$gender == 'female'] <- 0
    data$gender <- as.numeric(data$gender)
    
    data$LevelNNU[data$LevelNNU == 'yes'] <- 1
    data$LevelNNU[data$LevelNNU == 'no'] <- 0
    data$LevelNNU <- as.numeric(data$LevelNNU)
    
    data$surfactant[data$surfactant == 'yes'] <- 1
    data$surfactant[data$surfactant == 'no'] <- 0
    data$surfactant <- as.numeric(data$surfactant)
    
    data$Multip[data$Multip == 'yes'] <- 1
    data$Multip[data$Multip == 'no'] <- 0
    data$Multip <- as.numeric(data$Multip)
    
    data$ANS[data$ANS == 'yes'] <- 1
    data$ANS[data$ANS == 'no'] <- 0
    data$ANS <- as.numeric(data$ANS)
    
    data$PROMchorio[data$PROMchorio == 'yes'] <- 1
    data$PROMchorio[data$PROMchorio == 'no'] <- 0
    data$PROMchorio <- as.numeric(data$PROMchorio)
    
    data$mathighbp[data$mathighbp == 'yes'] <- 1
    data$mathighbp[data$mathighbp == 'no'] <- 0
    data$mathighbp <- as.numeric(data$mathighbp)
    
    data$matbleed[data$matbleed == 'yes'] <- 1
    data$matbleed[data$matbleed == 'no'] <- 0
    data$matbleed <- as.numeric(data$matbleed)
    
    data$matsmoking[data$matsmoking == 'yes'] <- 1
    data$matsmoking[data$matsmoking == 'no'] <- 0
    data$matsmoking <- as.numeric(data$matsmoking)
    
    data$IMDdecile <- (data$IMDdecile-mean(data$IMDdecile,na.rm=TRUE))/sd(data$IMDdecile,na.rm=TRUE)
    
    data$deliverymode[data$deliverymode == 'C-section'] <- 1
    data$deliverymode[data$deliverymode == 'Vaginal'] <- 0
    data$deliverymode <- as.numeric(data$deliverymode)
    
    data$tempbelow[data$tempbelow == 'yes'] <- 1
    data$tempbelow[data$tempbelow == 'no'] <- 0
    data$tempbelow <- as.numeric(data$tempbelow)
    
    
    static_predictors_of_interest <- c('zpreterm',
                                       'gender',
                                       'LevelNNU',
                                       'surfactant',
                                       'Multip',
                                       'ANS',
                                       'PROMchorio',
                                       'mathighbp',
                                       'matbleed',
                                       'matsmoking',
                                       'IMDdecile',
                                       'deliverymode',
                                       'tempbelow')
    N_s_predictors <- length(static_predictors_of_interest)
    static_predictors <- array(NA,c(N_individuals,N_s_predictors))
    load('ave_vals.RData')
    for (s_predictor in 1:N_s_predictors){
      static_predictors[,s_predictor] <- data[[static_predictors_of_interest[s_predictor]]]
      for (individual in 1:N_individuals){#for time being, set missing vals to average of their row
        if(is.na(static_predictors[individual,s_predictor])){
          static_predictors[individual,s_predictor] <- ave_vals[s_predictor]
        }
      }
    }
    
    
    ######################Extract Dynamic Copredictors#######################

    load(paste0('Fits/ClinFit',bootstrap,'.RData'))

    ivhuse <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("IVHuseD",age_day)]) | data[individual,paste0("IVHuseD",age_day)] != 'NA' | age_day > SeenTimePoints){#if value is NA, we must fill in
          if ('no' %in% data[individual,paste0("IVHuseD",age_day:100)]){#if we know that there is a future point with no, then it must be no now
            data[individual,paste0('IVHuseD',age_day)] <- 'no'
          } else {#if we don't know that, then we guess
            data[individual,paste0('IVHuseD',age_day)] <- sample(c('no','yes'),
                                                                 1,
                                                                 prob=c(1-tanh(exp(log(exclinopt$p[1]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[1] + (gest_day - 160) * exclinopt$base_decay[1])),
                                                                        tanh(exp(log(exclinopt$p[1]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[1] + (gest_day - 160) * exclinopt$base_decay[1]))))
          }
        }
        if ('yes' %in% data[individual,paste0("IVHuseD",1:age_day)]){counter <- counter + 1}
        ivhuse[individual,gest_day] <- counter
      }
    }
    
    necsurg <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("necsurgD",age_day)]) | data[individual,paste0("necsurgD",age_day)] != 'NA' | age_day > SeenTimePoints){#if value is NA, we must fill in
          if ('no' %in% data[individual,paste0("necsurgD",age_day:100)]){#if we know that there is a future point with no, then it must be no now
            data[individual,paste0('necsurgD',age_day)] <- 'no'
          } else {#if we don't know that, then we guess
            data[individual,paste0('necsurgD',age_day)] <- sample(c('no','yes'),
                                                                  1,
                                                                  prob=c(1-tanh(exp(log(exclinopt$p[2]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[2] + (gest_day - 160) * exclinopt$base_decay[2])),
                                                                         tanh(exp(log(exclinopt$p[2]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[2] + (gest_day - 160) * exclinopt$base_decay[2]))))
          }
        }
        if ('yes' %in% data[individual,paste0("necsurgD",1:age_day)]){counter <- counter + 1}
        necsurg[individual,gest_day] <- counter
      }
    }
    
    
    
    
    pda <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      pda_used <- FALSE
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("PDAD",age_day)]) | age_day > SeenTimePoints){#if value is NA, we must fill in
          if ('no' %in% data[individual,paste0("PDAD",age_day:100)]){#if we know that there is a future point with no, then it must be no now
            data[individual,paste0('PDAD',age_day)] <- 'no'
          } else {#if we don't know that, then we guess
            data[individual,paste0('PDAD',age_day)] <- sample(c('no','yes'),
                                                              1,
                                                              prob=c(1-tanh(exp(log(exclinopt$p[3]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[3] + (gest_day - 160) * exclinopt$base_decay[3])),
                                                                     tanh(exp(log(exclinopt$p[3]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[3] + (gest_day - 160) * exclinopt$base_decay[3]))))
          }
        }
        if ('yes' %in% data[individual,paste0("PDAD",1:age_day)]){counter <- counter + 1}
        pda[individual,gest_day] <- counter
      }
    }
    
    inotrope <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("inotropeD",age_day)]) | age_day > SeenTimePoints){
          data[individual,paste0('inotropeD',age_day)] <- sample(c('no','yes'),
                                                                 1,
                                                                 prob=c(1-tanh(exp(log(exclinopt$p[4]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[4] + (gest_day - 160) * exclinopt$base_decay[4])),
                                                                        tanh(exp(log(exclinopt$p[4]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[4] + (gest_day - 160) * exclinopt$base_decay[4]))))
        }
        if ('yes' %in% data[individual,paste0("inotropeD",1:age_day)]){counter <- counter + 1}
        if (!is.na(data[individual,paste0("inotropeD",age_day)]) & data[individual,paste0("inotropeD",age_day)] == 'yes'){counter <- 0}
        inotrope[individual,gest_day] <- counter
      }
    }
    
    sepsis <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("sepsisD",age_day)]) | age_day > SeenTimePoints){
          data[individual,paste0('sepsisD',age_day)] <- sample(c('no','yes'),
                                                               1,
                                                               prob=c(1-tanh(exp(log(exclinopt$p[5]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[5] + (gest_day - 160) * exclinopt$base_decay[5])),
                                                                      tanh(exp(log(exclinopt$p[5]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[5] + (gest_day - 160) * exclinopt$base_decay[5]))))
        }
        if ('yes' %in% data[individual,paste0("sepsisD",1:age_day)]){counter <- counter + 1}
        if (!is.na(data[individual,paste0("sepsisD",age_day)]) & data[individual,paste0("sepsisD",age_day)] == 'yes'){counter <- 0}
        sepsis[individual,gest_day] <- counter
      }
    }
    
    #for now, dont use feed
    #feed <- array(-1,c(N_individuals,252))
    #for (individual in 1:N_individuals){
    #  counter <- -1
    #  feed_used <- FALSE
    #  for (gest_day in data$gest[individual]:252){
    #    age_day <- gest_day - data$gest[individual] + 1
    #    if ('feed' %in% data[individual,paste0("feedD",1:age_day)]){counter <- counter + 1}
    #    if (!is.na(data[individual,paste0("feedD",age_day)]) & data[individual,paste0("feedD",age_day)] == 'yes'){counter <- 0}
    #    feed[individual,gest_day] <- counter
    #  }
    #}
    
    pnuse <- array(-1,c(N_individuals,252))
    for (individual in 1:N_individuals){
      counter <- -1
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!is.na(data[individual,paste0("PNuseD",age_day)]) | age_day > SeenTimePoints){
          data[individual,paste0('PNuseD',age_day)] <- sample(c('no','yes'),
                                                              1,
                                                              prob=c(1-tanh(exp(log(exclinopt$p[6]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[6] + (gest_day - 160) * exclinopt$base_decay[6])),
                                                                     tanh(exp(log(exclinopt$p[6]) + (gest_day - data$gest[individual]) * exclinopt$age_decay[6] + (gest_day - 160) * exclinopt$base_decay[6]))))
        }
        if ('yes' %in% data[individual,paste0("PNuseD",1:age_day)]){counter <- counter + 1}
        if (!is.na(data[individual,paste0("PNuseD",age_day)]) & data[individual,paste0("PNuseD",age_day)] == 'yes'){counter <- 0}
        pnuse[individual,gest_day] <- counter
      }
    }
    
    
    d_predictors_of_interest <- c('ivhuse',
                                  'necsurg',
                                  'pda',
                                  'inotrope',
                                  'sepsis',
                                  #'feed',
                                  'pnuse')
    N_d_predictors <- length(d_predictors_of_interest)
    d_predictors <- array(NA,c(N_individuals,252,N_d_predictors))
    d_predictors[,,1] <- ivhuse
    d_predictors[,,2] <- necsurg
    d_predictors[,,3] <- pda
    d_predictors[,,4] <- inotrope
    d_predictors[,,5] <- sepsis
    #d_predictors[,,6] <- feed
    d_predictors[,,6] <- pnuse
    
    
    
    
    ###########################Extract Clinical Outcomes##############################
    #This block is the same when we do model training
    #We remove values past SeenTimePoints later, when creating clin_outcomes_new
    clin_outcomes <- array(-1,c(N_individuals,252,5))
    for (individual in 1:N_individuals){
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        if (!(is.na(data[individual,paste0("CLDD",age_day)])) & data[individual,paste0("CLDD",age_day)] == 'none'){
          clin_outcomes[individual,gest_day,] <- c(1,0,0,0,0)
        } else if (!(is.na(data[individual,paste0("CLDD",age_day)])) &data[individual,paste0("CLDD",age_day)] == 'oxygen'){
          clin_outcomes[individual,gest_day,] <- c(0,1,0,0,0)
        } else if (!(is.na(data[individual,paste0("CLDD",age_day)])) &data[individual,paste0("CLDD",age_day)] == 'NIV'){
          clin_outcomes[individual,gest_day,] <- c(0,0,1,0,0)
        } else if (!(is.na(data[individual,paste0("CLDD",age_day)])) &data[individual,paste0("CLDD",age_day)] == 'vent'){
          clin_outcomes[individual,gest_day,] <- c(0,0,0,1,0)
        } else if ((!is.na(data$deathageday[individual])) & data$deathageday[individual] < age_day){
          clin_outcomes[individual,gest_day,] <- c(0,0,0,0,1)
        }
      }
      #Also use final point, as we don't always have that explicitly:
      final_age_day <- 252 - data$gest[individual] + 1
      if (data$BPDgrade[individual]=='None'){
        clin_outcomes[individual,252,] <- c(1,0,0,0,0)
      } else if (data$BPDgrade[individual] == 'Oxygen'){
        clin_outcomes[individual,252,] <- c(0,1,0,0,0)
      } else if (data$BPDgrade[individual] == 'NIV'){
        clin_outcomes[individual,252,] <- c(0,0,1,0,0)
      } else if (data$BPDgrade[individual] == 'Inv'){
        clin_outcomes[individual,252,] <- c(0,0,0,1,0)
      } else if (data$BPDgrade[individual] == 'NA'){
        clin_outcomes[individual,252,] <- c(0,0,0,0,1)
      }
      #also, if we reach NAs, and last value was none, and final value is none, then baby has been discharged
      for (gest_day in data$gest[individual]:252){
        age_day <- gest_day - data$gest[individual] + 1
        found = FALSE
        found_value = 'init'
        timestepper = age_day
        if ('none' %in% data[individual,paste0("CLDD",1:age_day)] & !(TRUE %in% (c('oxygen','vent','NIV') %in% data[individual,paste0("CLDD",age_day:100)])) & (data$BPDgrade[individual] == 'None')){
          while (found == FALSE & timestepper>1){#look backwards in steps until we find last known value
            timestepper = timestepper-1
            #record the last known value before this time
            if (!is.na(data[individual,paste0("CLDD",timestepper)]) & !(data[individual,paste0("CLDD",timestepper)] == 'NA' )){
              found_value = data[individual,paste0("CLDD",timestepper)]
              found = TRUE
            }
          }
        }
        #if the last known value was none, then assume discharge
        if (data$BPDgrade[individual] == 'None' & found_value == 'none' & found==TRUE){
          clin_outcomes[individual,gest_day,] <- c(1,0,0,0,0)
        }
      }
    }
    
    #Truncate dynamic variables
    d_predictors <- d_predictors[,161:252,]
    clin_outcomes <- clin_outcomes[,161:252,]
    max_time <- length(161:252)
    init_times <- data$gest-160
    
    
    #convert clin_outcomes to contain number of outcome, rather than an array
    clin_outcomes_new <- array(-1,c(N_individuals,max_time))
    for (individual in 1:N_individuals){
      for (time in 1:(init_times[individual] + SeenTimePoints - 1)){#Only use time up to certain value for evaluation
        for (outcome in 1:5){
          if (clin_outcomes[individual,time,outcome] == 1){
            clin_outcomes_new[individual,time] <- outcome
          }
        }
      }
    }
    
    
    ####################Fitting Model###############################################
    ModelFit <- function(data,pars){
      for (parameter in names(pars)){
        do.call("<-",list(parameter, pars[[parameter]]))
      }
      for (parameter in names(data)){
        do.call("<-",list(parameter, data[[parameter]]))
      }
      
      base_propensity <- array(NA,c(5,max_time,5))
      static_multiplier <- array(NA,c(5,5))
      dynamic_multiplier <- array(NA,c(5,5))
      age_multiplier <- array(NA,c(5,5))
      prev_state_multiplier <- array(NA,c(5,5))
      transition <- array(NA,c(max_time,5,5))
      
      state_prob <- array(NA,c(max_time,5))
      state_prob_guess <- array(NA,c(max_time,5))
      state_prob_guess_all <- array(NA,c(N_individuals,max_time,5))
      
      #calculate base improvement and deterioration probabilities
      for (time in 1:max_time){
        for (clin_outcome in 1:4){
          for (clin_outcome2 in 1:5){
            base_propensity[clin_outcome,time,clin_outcome2] = p[clin_outcome,clin_outcome2] * exp(-(time/base_par_1[clin_outcome,clin_outcome2])^base_par_2[clin_outcome,clin_outcome2])
          }
        }
      }
      
      
      
      
      #run chain
      for (individual in 1:N_individuals){
        
        #initialize chain
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
              static_multiplier[clin_outcome,clin_outcome2] = exp(exp((time-init_times[individual]) * static_decay[clin_outcome,clin_outcome2]) * dot(as.vector(s_predictors[individual,]),as.vector(static_coefs[,clin_outcome,clin_outcome2]) ));#.* as.vector(exp(time*as.vector(static_decay[,clin_outcome,clin_outcome2])))));
  dynamic_multiplier[clin_outcome,clin_outcome2] = 1;
  for (d_pred in 1:N_d_predictors){
    if (d_predictors[individual,time,d_pred] != -1){
      dynamic_multiplier[clin_outcome,clin_outcome2] = dynamic_multiplier[clin_outcome,clin_outcome2] * exp(dynam_coefs[d_pred,clin_outcome,clin_outcome2]*exp(d_predictors[individual,time,d_pred]*dynam_decay[d_pred,clin_outcome,clin_outcome2]));
    }
  }

  age_multiplier[clin_outcome,clin_outcome2] = exp(age_coefs[clin_outcome,clin_outcome2] * exp((time-init_times[individual])*age_decay[clin_outcome,clin_outcome2]));
            }
          }
          
          
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
          
        }
        
        for (time in (init_times[individual]+1):max_time){
          
          for (clin_outcome in 1:5){
            for (clin_outcome2 in 1:5){
              state_prob[time,clin_outcome2] = state_prob[time,clin_outcome2] + state_prob[time-1,clin_outcome] * transition[time-1,clin_outcome,clin_outcome2];
            }
          }
          
          #Now if we are at a known point - record guess and correct
          for (clin_outcome in 1:5){
            state_prob_guess[time,clin_outcome] = state_prob[time,clin_outcome];
            #Also condition result if we are on a training individual, or if we are below threshold for test individuals
            
            if (clin_outcomes[individual,time] == clin_outcome){
              state_prob[time,clin_outcome] = 1;
            } else if (clin_outcomes[individual,time] != -1) {
              state_prob[time,clin_outcome] = 0;
            }
            
          }
          
        }
        #Record all guesses for model evaluation
        state_prob_guess_all[individual,,] <- state_prob_guess
      }
      return(state_prob_guess_all)
    }
    
    df <- list(
      N_individuals = N_individuals,
      N_s_predictors = N_s_predictors,
      N_d_predictors = N_d_predictors,
      init_times = init_times,
      s_predictors = static_predictors,
      d_predictors = d_predictors,
      clin_outcomes = clin_outcomes_new,
      max_time = max_time
    )
    
    if (SaveResults == TRUE){save(df,file=paste0('dfDay',SeenTimePoints,'.RData'))}
    
  

    
  
      
      load(paste0('Fits/OptFitBootstrap',bootstrap_nums[bootstrap],'.RData'))
      preds <- ModelFit(df,exopt)
      pr_bootstrap$None[,] <- preds[,,1]
      pr_bootstrap$Oxygen[,] <- preds[,,2]
      pr_bootstrap$NIV[,] <- preds[,,3]
      pr_bootstrap$Inv[,] <- preds[,,4]
      pr_bootstrap$Death[,] <- preds[,,5]
      
      if (SaveResults == TRUE){save(pr_bootstrap,file=paste0('Bootstraps/pr_bootstrapSeed',bootstrap_nums[bootstrap],'Day',SeenTimePoints,'.RData'))}
      }
  



}

library(parallel)
cl <- makeCluster(getOption("cl.cores", 4))
parLapply(cl=cl,c(28,1,21,7),fitter)
cl <- makeCluster(getOption("cl.cores", 2))
parLapply(cl=cl,c(3,14),fitter)


# Required libraries
library(dplyr)
setwd("~/281224")#CHANGE TO RELEVANT DIRECTORY

# Define X values to process
X_values <- c(1, 3, 7, 14, 21, 28)
Y_values <- 1:200  # Adjust to actual range
folders <- c("Bootstraps", "D:/Bootstraps")  # Search both locations

# Function to find the first existing file in the given folders
find_file <- function(filename) {
  for (folder in folders) {
    full_path <- file.path(folder, filename)
    if (file.exists(full_path)) {
      return(full_path)
    }
  }
  return(NULL)  # Return NULL if file is not found
}

# Function to compute quantiles over an array efficiently
compute_quantiles <- function(data, probs = c(0.025, 0.5, 0.975)) {
  apply(data, c(1, 2), function(x) quantile(x, probs = probs, na.rm = TRUE))
}

for (X in X_values) {
  print(paste0('X=',X))
  # Load base bootstrap data
  base_filename <- paste0("pr_bootstrapSeed1Day", X, ".RData")
  base_path <- find_file(base_filename)
  
  if (!is.null(base_path)) {
    load(base_path)  # Loads pr_bootstrap
  } else {
    next  # Skip this X if base file is missing
  }
  
  # Get the number of individuals
  num_individuals <- dim(pr_bootstrap$Death)[1]
  
  # Initialize lists for quantiles
  upper <- list('Death' = array(NA, c(92, num_individuals)),
                'Inv' = array(NA, c(92, num_individuals)),
                'NIV' = array(NA, c(92, num_individuals)),
                'None' = array(NA, c(92, num_individuals)),
                'Oxygen' = array(NA, c(92, num_individuals)))
  
  mid <- upper  # Structure is the same
  lower <- upper  # Structure is the same
  
  # Initialize empty list to store all values for each category
  all_vals <- list('Death' = array(NA, c(num_individuals, 92, length(Y_values))),
                   'Inv' = array(NA, c(num_individuals, 92, length(Y_values))),
                   'NIV' = array(NA, c(num_individuals, 92, length(Y_values))),
                   'None' = array(NA, c(num_individuals, 92, length(Y_values))),
                   'Oxygen' = array(NA, c(num_individuals, 92, length(Y_values))))
  
  # Load all data for Y values and store it in all_vals
  for (Y in Y_values) {
    print(paste0('Y=',Y))
    filename <- paste0("pr_bootstrapSeed", Y, "Day", X, ".RData")
    file_path <- find_file(filename)
    
    if (!is.null(file_path)) {
      load(file_path)  # Loads pr_bootstrap
      
      # Store the data into all_vals
      all_vals$Death[,,Y] <- pr_bootstrap$Death
      all_vals$Inv[,,Y] <- pr_bootstrap$Inv
      all_vals$NIV[,,Y] <- pr_bootstrap$NIV
      all_vals$None[,,Y] <- pr_bootstrap$None
      all_vals$Oxygen[,,Y] <- pr_bootstrap$Oxygen
    }
  }
  
  # Compute quantiles for each category across all time points and individuals
  for (category in names(all_vals)) {
    print(paste0('category=',category))
    category_data <- all_vals[[category]]
    
    # Vectorized quantile computation over entire time/individual dimension
    # Compute quantiles across time points for each individual
    quant_vals <- apply(category_data, c(1, 2), function(x) {
      # Use a direct call to quantile, removing unnecessary apply levels
      quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    })
    
    # Extract quantile values (0.025, 0.5, 0.975)
    lower[[category]] <- quant_vals[1, , ]
    mid[[category]] <- quant_vals[2, , ]
    upper[[category]] <- quant_vals[3, , ]
  }
  
  
  
  # Save results
  save(upper, file = paste0('upperExternalDay', X, '.RData'))
  save(mid, file = paste0('midExternalDay', X, '.RData'))
  save(lower, file = paste0('lowerExternalDay', X, '.RData'))
}
