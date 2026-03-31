#import variables from command line
args = commandArgs(trailingOnly=TRUE)
SEED = args[1]
#SEED = 1
print(paste0('Commencing modelling for seed ',SEED))

print('Importing Libraries')
#import libraries
library(cmdstanr)

print('Loading Data')
#select data
setwd("~/281224")#CHANGE TO RELEVANT DIRECTORY
load("Code/DatasetsFull141024.RData")#CHANGE TO RELEVANT FILE
datafull <- datafull[datafull$randomcohort == 'development',]
N_individuals <- 5000 #use only a section of the data
set.seed(SEED)
selected <- sample(1:dim(datafull)[1],size=N_individuals,replace=TRUE)#select indexes of patients for bootstrapping
data = datafull[selected,]


#shunt death time along 1 as we start counting days from 1 rather than zero
data$deathageday <- data$deathageday + 1


#################Extract Static Copredictors######################

#convert gestational age (gest) from weeks to days
data$gest <- round(7*data$gest)

days_monitored <- 100

#convert some static variables to be in usable number format
data$zpreterm <- as.numeric(data$zpreterm)

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
ave_vals <- array(NA,N_s_predictors)
for (s_predictor in 1:N_s_predictors){
  static_predictors[,s_predictor] <- data[[static_predictors_of_interest[s_predictor]]]
  ave_vals[s_predictor] <- mean(static_predictors[,s_predictor],na.rm=TRUE)
  for (individual in 1:N_individuals){#for time being, set missing vals to average of their row
    if(is.na(static_predictors[individual,s_predictor])){
      static_predictors[individual,s_predictor] <- ave_vals[s_predictor]
    }
  }
}
save(ave_vals,file='ave_vals.RData')
  
######################Extract Dynamic Copredictors#######################
load(paste0('Fits/ClinFit',SEED,'.RData'))
ivhuse <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  counter <- -1
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!is.na(data[individual,paste0("IVHuseD",age_day)]) | data[individual,paste0("IVHuseD",age_day)] != 'NA'){#if value is NA, we must fill in
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
    if (!is.na(data[individual,paste0("necsurgD",age_day)]) | data[individual,paste0("necsurgD",age_day)] != 'NA'){#if value is NA, we must fill in
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
    if (!is.na(data[individual,paste0("PDAD",age_day)])){#if value is NA, we must fill in
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
    if (!is.na(data[individual,paste0("inotropeD",age_day)])){
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
    if (!is.na(data[individual,paste0("sepsisD",age_day)])){
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
    if (!is.na(data[individual,paste0("PNuseD",age_day)])){
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
    } else if (((is.na(data[individual,paste0("CLDD",age_day)])) | data[individual,paste0("CLDD",age_day)] == 'NA') & data$BPDdeath[individual] == 'yes' & !(FALSE %in% (data[individual,paste0("CLDD",99:100)] == array('NA',length(age_day:100))))){
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

#convert clin_outcomes to contain number of outcome, rather than an array
clin_outcomes_new <- array(-1,c(N_individuals,max_time))
for (individual in 1:N_individuals){
  for (time in 1:max_time){
    for (outcome in 1:5){
      if (clin_outcomes[individual,time,outcome] == 1){
        clin_outcomes_new[individual,time] <- outcome
      }
    }
  }
}


####################Fitting Model###############################################
df <- list(
  N_individuals = N_individuals,
  N_s_predictors = N_s_predictors,
  N_d_predictors = N_d_predictors,
  init_times = data$gest-160,
  s_predictors = static_predictors,
  d_predictors = d_predictors,
  clin_outcomes = clin_outcomes_new,
  max_time = max_time
)


print('Loading Model')

model <- cmdstan_model('Code/StanModelCMD.stan')




l <- list()
l$p <- array(0.1,c(4,5))
l$base_par_1 <- array(10,c(4,5))
l$base_par_2 <- array(1,c(4,5))
l$static_coefs <- array(0.5,c(df$N_s_predictors,4,5))
#l$static_decay <- array(0,c(df$N_s_predictors,4,5))
l$static_decay <- array(-0.01,c(4,5))
l$dynam_coefs <- array(0.5,c(df$N_d_predictors,4,5))
l$dynam_decay <- array(-0.1,c(df$N_d_predictors,4,5))
l$age_coefs <- array(0.5,c(4,5))
l$age_decay <- array(-0.1,c(4,5))
#l$prev_state_coefs <- array(0.001,c(3,4,5))
#l$prev_state_decay <- array(-0.1,c(3,4,5))
initf1 <- function() {
  return(l)
}

initf2 <- function(){
  l = exopt

  return(l)
}
#load('OptFitBootstrap1.RData')
load(paste0('Fits/OptFitBootstrap',SEED,'.RData'))

print('Fitting')
opt <- model$optimize(data=df,
                      init=initf2,
                      iter = 20000)
#samps <- model$sample(data=df,init=initf1,chains=2)
exopt <- relist(opt$mle(),skeleton=l)
exopt$static_coefs <- array(exopt$static_coefs,dim=dim(l$static_coefs))
#exopt$static_decay <- array(exopt$static_decay,dim=dim(l$static_decay))
exopt$dynam_coefs <- array(exopt$dynam_coefs,dim=dim(l$dynam_coefs))
exopt$dynam_decay <- array(exopt$dynam_decay,dim=dim(l$dynam_decay))
#exopt$prev_state_coefs <- array(exopt$prev_state_coefs,dim=dim(l$prev_state_coefs))
#exopt$prev_state_decay <- array(exopt$prev_state_decay,dim=dim(l$prev_state_decay))
save(exopt,file=paste0('Fits/OptFitBootstrap',SEED,'.RData'))
print(paste0('Modelling completed'))

quit()
