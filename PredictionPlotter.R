# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

SeenTimePoints = 1

#load data
setwd("~/281224")

#load predictions
load(paste0('mid','ExternalDay',SeenTimePoints,'.RData'))

#load data
#load("Code/DatasetsFull141024.RData")#CHANGE TO RELEVANT FILE
#datafull <- datafull[datafull$randomcohort == 'testing',]
load('Code/DatasetTestFull300125.RData')
data = datafull
rm(datafull)
#shunt death time along 1 as we start counting days from 1 rather than zero
data$deathageday <- data$deathageday + 1
#remove babies who died before SeenTimePoints
data <- data[(data$deathageday > SeenTimePoints) | is.na(data$deathageday),]

#synthesize data
max_time = 92
pr_m <- list(None = array(NA,c(dim(mid$None)[1],max_time)),
             Oxygen = array(NA,c(dim(mid$None)[1],max_time)),
             NIV = array(NA,c(dim(mid$None)[1],max_time)),
             Inv = array(NA,c(dim(mid$None)[1],max_time)),
             Death = array(NA,c(dim(mid$None)[1],max_time)))
for (individual in 1:dim(pr_m$None)[1]){
  for (t in 1:max_time){
    pr_m$None[individual,t] <- mid$None[individual,t]
    pr_m$Oxygen[individual,t] <- mid$Oxygen[individual,t]
    pr_m$NIV[individual,t] <- mid$NIV[individual,t]
    pr_m$Inv[individual,t] <- mid$Inv[individual,t]
    pr_m$Death[individual,t] <- mid$Death[individual,t]
    
  }
}



individual = 3
states_labelled <- data[individual,paste0('CLDD',1:length((round(7*data$gest[individual])-160):92))]
states_numbered <- array(NA,length((round(7*data$gest[individual])-160):92))
for (t in 1:length(states_numbered)){
  if (!is.na(states_labelled[t])){
    if (states_labelled[t] == 'none'){
      states_numbered[t] <- 1
    } else if (states_labelled[t] == 'oxygen'){
      states_numbered[t] <- 2
    } else if (states_labelled[t] == 'NIV'){
      states_numbered[t] <- 3
    } else if (states_labelled[t] == 'vent'){
      states_numbered[t] <- 4
    } else if (!is.na(data$deathageday[individual]) & data$deathageday[individual] <= t){
      states_numbered[t] <- 5
    }
  }
}
  
df <- data.frame(
  time = round(7*data$gest[individual]):252,
  state = states_numbered,
  p1 = pr_m$Death[individual,(round(7*data$gest[individual])-160):92],
  p2 = pr_m$Inv[individual,(round(7*data$gest[individual])-160):92],
  p3 = pr_m$NIV[individual,(round(7*data$gest[individual])-160):92],
  p4 = pr_m$Oxygen[individual,(round(7*data$gest[individual])-160):92],
  p5 = pr_m$None[individual,(round(7*data$gest[individual])-160):92]
)

# Normalize probabilities to sum to 1 at each time point
df[, paste0("p", 1:5)] <- df[, paste0("p", 1:5)] / rowSums(df[, paste0("p", 1:5)])

# Calculate cumulative probabilities for ribbons
df <- df %>%
  mutate(cum_p1 = p1,
         cum_p2 = cum_p1 + p2,
         cum_p3 = cum_p2 + p3,
         cum_p4 = cum_p3 + p4,
         cum_p5 = cum_p4 + p5) %>%
  select(time, state, p1, cum_p1, cum_p2, cum_p3, cum_p4, cum_p5)

df$time_index <- seq_len(nrow(df))-1

# Create the plot
ggplot(df, aes(x = time_index)) +#set x = time to count in gestational age rather than age since birth
  geom_ribbon(aes(ymin = 0, ymax = cum_p1, fill = "Death"), alpha = 0.6) +
  geom_ribbon(aes(ymin = cum_p1, ymax = cum_p2, fill = "Inv"), alpha = 0.6) +
  geom_ribbon(aes(ymin = cum_p2, ymax = cum_p3, fill = "NIV"), alpha = 0.6) +
  geom_ribbon(aes(ymin = cum_p3, ymax = cum_p4, fill = "Oxygen"), alpha = 0.6) +
  geom_ribbon(aes(ymin = cum_p4, ymax = cum_p5, fill = "None"), alpha = 0.6) +
  geom_rect(aes(xmin = 0, xmax = SeenTimePoints, ymin = 0, ymax = 1), 
            fill = "white", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(
    data = subset(df, !is.na(state)), # Exclude rows where state is NA
    aes(x = time_index, y = 1, fill = factor(state)), 
    shape = 21, size = 3, stroke = 0.7, color = "black", inherit.aes = FALSE
  ) +
#  geom_line(aes(x=rep(SeenTimePoints,length(time)),y=array(c(0,1),length(time))),size=2,col='black')+
  scale_fill_manual(
    values = c(
      "None" = "green", 
      "Oxygen" = "blue", 
      "NIV" = "orange", 
      "Inv" = "grey", 
      "Death" = "red",
      "1" = "green", 
      "2" = "blue", 
      "3" = "orange", 
      "4" = "grey", 
      "5" = "red"
    ),
    name = "States",
    breaks = c("None", "Oxygen", "NIV", "Inv", "Death"), # Legend order
    guide = guide_legend(order = 1) # Legend customization
  ) +
  labs(
    x = "Days since birth",
    y = "Probability of belonging to state"
  ) +
  theme_minimal()+
  theme(legend.background = element_rect(fill = "white", color = "black"),
        legend.position = c(0.8, 0.7),
        #legend.position = 'None',
        text = element_text(size = 20))
ggsave('D1I3NewStartBlanked.png')
