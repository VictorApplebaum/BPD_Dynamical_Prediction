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

data$matsmoking[data$matsmoking == 'yes'] <- 1
data$matsmoking[data$matsmoking == 'no'] <- 0
data$matsmoking <- as.numeric(data$matsmoking)

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
ivhuse <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  first_point <- 0
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("IVHuseD",age_day)])) & data[individual,paste0("IVHuseD",age_day)] == 'yes' & first_point == 0){
      ivhuse[individual,gest_day] <- 1
      first_point <- 1
    } else if (!(is.na(data[individual,paste0("IVHuseD",age_day)])) & data[individual,paste0("IVHuseD",age_day)] == 'no'){
      ivhuse[individual,gest_day] <- 0
    }
  }
}


necsurg <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  first_point <- 0
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("necsurgD",age_day)])) & data[individual,paste0("necsurgD",age_day)] == 'yes' & first_point == 0){
      necsurg[individual,gest_day] <- 1
      first_point <- 1
    } else if (!(is.na(data[individual,paste0("necsurgD",age_day)])) & data[individual,paste0("necsurgD",age_day)] == 'no'){
      necsurg[individual,gest_day] <- 0
    }
  }
}

pda <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  first_point <- 0
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("PDAD",age_day)])) & data[individual,paste0("PDAD",age_day)] == 'yes' & first_point == 0){
      pda[individual,gest_day] <- 1
      first_point <- 1
    } else if (!(is.na(data[individual,paste0("PDAD",age_day)])) & data[individual,paste0("PDAD",age_day)] == 'no'){
      pda[individual,gest_day] <- 0
    }
  }
}


inotrope <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("inotropeD",age_day)])) & data[individual,paste0("inotropeD",age_day)] == 'yes'){
      inotrope[individual,gest_day] <- 1
    } else if (!(is.na(data[individual,paste0("inotropeD",age_day)])) & data[individual,paste0("inotropeD",age_day)] == 'no'){
      inotrope[individual,gest_day] <- 0
    }
  }
}


sepsis <- array(-1,c(N_individuals,252))
for (individual in 1:N_individuals){
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("sepsisD",age_day)])) & data[individual,paste0("sepsisD",age_day)] == 'yes'){
      sepsis[individual,gest_day] <- 1
    } else if (!(is.na(data[individual,paste0("sepsisD",age_day)])) & data[individual,paste0("sepsisD",age_day)] == 'no'){
      sepsis[individual,gest_day] <- 0
    }
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
  for (gest_day in data$gest[individual]:252){
    age_day <- gest_day - data$gest[individual] + 1
    if (!(is.na(data[individual,paste0("PNuseD",age_day)])) & data[individual,paste0("PNuseD",age_day)] == 'yes'){
      pnuse[individual,gest_day] <- 1
    } else if (!(is.na(data[individual,paste0("PNuseD",age_day)])) & data[individual,paste0("PNuseD",age_day)] == 'no'){
      pnuse[individual,gest_day] <- 0
    }
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




#Truncate dynamic variables
d_predictors <- d_predictors[,161:252,]
max_time <- length(161:252)




####################Fitting Model###############################################
df <- list(
  N_individuals = N_individuals,
  #N_s_predictors = N_s_predictors,
  N_d_predictors = N_d_predictors,
  init_times = data$gest-160,
  #s_predictors = static_predictors,
  d_predictors = d_predictors,
  #clin_outcomes = clin_outcomes_new,
  max_time = max_time
)

clin_model <- cmdstan_model('Code/ClinModel.stan')
l <- list()
l$age_decay <- array(0.01,df$N_d_predictors)
l$base_decay <- array(-0.01,df$N_d_predictors)
l$p <- array(0.01,df$N_d_predictors)
initf <- function(){
  return(l)
}

file1 <- paste0("~/231124/Fits/ClinFit", SEED, ".RData")
file2 <- paste0("~/231124/Fits/ClinFit", 1, ".RData")

if (file.exists(file1)) {
  load(file1)
} else {
  load(file2)
}
opt <- clin_model$optimize(data=df,
                           init=list(exclinopt)
                           )
exclinopt <- relist(opt$mle(),skeleton=l)
save(exclinopt,file=paste0('Fits/ClinFit',SEED,'.RData'))

quit()













