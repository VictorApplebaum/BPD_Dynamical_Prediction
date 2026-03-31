for (SeenTimePoints in c(1,3,7,14,21,28)){
  for (version in c('mid')){
    for (DaysForward in c(7,14,21,28)){#if only last day is wanted, delete this loop and set FINAL_DAY to true
      FINAL_DAY <- FALSE
      setwd("~/281224")
      
      #load(paste0('pr_bootstrapsDay',SeenTimePoints,'.RData'))
      load(paste0(version,'ExternalDay',SeenTimePoints,'.RData'))
      pr_bootstraps <- get(version)
      load("Code/DatasetsFull141024.RData")#CHANGE TO RELEVANT FILE
      datafull <- datafull[datafull$randomcohort == 'testing',]
      #datafull <- datafull[1:5000,]
      data = datafull
      rm(datafull)
      #shunt death time along 1 as we start counting days from 1 rather than zero
      data$deathageday <- data$deathageday + 1
      #remove babies who died before SeenTimePoints
      data <- data[(data$deathageday > SeenTimePoints) | is.na(data$deathageday),]
      #Also remove babies who we have to predict to outside of the time period we're interested in
      if (FINAL_DAY == FALSE){
        new_indexes = SeenTimePoints+DaysForward+round(data$gest*7)-161<93
        data <- data[new_indexes,]
      }
      
      max_time <- length(161:252)
      N_individuals <- dim(data)[1]
      true_results <- c()
      for (individual in 1:N_individuals){
        true_results <- append(true_results,data[individual,paste0('CLDD',SeenTimePoints+DaysForward)])
        if ((!is.na(data$deathageday[individual]) & data$deathageday[individual] <= SeenTimePoints+DaysForward)){
          true_results[individual] <- 'Death'
        }
      }
      
      
      #Evaluating at final time point
      None_vec <- array(NA,dim=c(N_individuals))
      Ox_vec <- array(NA,dim=c(N_individuals))
      NIV_vec <- array(NA,dim=c(N_individuals))
      Inv_vec <- array(NA,dim=c(N_individuals))
      Death_vec <- array(NA,dim=c(N_individuals))
      for (individual in 1:N_individuals){
        if (FINAL_DAY == TRUE){
          None_vec[individual] <- pr_bootstraps$None[individual,92]
          Ox_vec[individual] <- pr_bootstraps$Oxygen[individual,92]
          NIV_vec[individual] <- pr_bootstraps$NIV[individual,92]
          Inv_vec[individual] <- pr_bootstraps$Inv[individual,92]
          Death_vec[individual] <- pr_bootstraps$Death[individual,92]
        } else {
          None_vec[individual] <- pr_bootstraps$None[which(new_indexes)[individual],round(7*data$gest[individual])+SeenTimePoints+DaysForward-161]
          Ox_vec[individual] <- pr_bootstraps$Oxygen[which(new_indexes)[individual],round(7*data$gest[individual])+SeenTimePoints+DaysForward-161]
          NIV_vec[individual] <- pr_bootstraps$NIV[which(new_indexes)[individual],round(7*data$gest[individual])+SeenTimePoints+DaysForward-161]
          Inv_vec[individual] <- pr_bootstraps$Inv[which(new_indexes)[individual],round(7*data$gest[individual])+SeenTimePoints+DaysForward-161]
          Death_vec[individual] <- pr_bootstraps$Death[which(new_indexes)[individual],round(7*data$gest[individual])+SeenTimePoints+DaysForward-161]
        }
      }
      
      
      if (FINAL_DAY == FALSE){
        None_vec <- None_vec[true_results != 'NA' & !is.na(true_results)]
        Ox_vec <- Ox_vec[true_results != 'NA' & !is.na(true_results)]
        NIV_vec <- NIV_vec[true_results != 'NA' & !is.na(true_results)]
        Inv_vec <- Inv_vec[true_results != 'NA' & !is.na(true_results)]
        Death_vec <- Death_vec[true_results != 'NA' & !is.na(true_results)]
        true_results <- true_results[true_results != 'NA' & !is.na(true_results)]
        
        
        #Ox_vec <- Ox_vec[None_vec != 'NA' & !is.na(None_vec)]
        #NIV_vec <- NIV_vec[None_vec != 'NA' & !is.na(None_vec)]
        #Inv_vec <- Inv_vec[None_vec != 'NA' & !is.na(None_vec)]
        #Death_vec <- Death_vec[None_vec != 'NA' & !is.na(None_vec)]
        #true_results <- true_results[None_vec != 'NA' & !is.na(None_vec)]
        #None_vec <- None_vec[None_vec != 'NA' & !is.na(None_vec)]
      }
      
      if (FINAL_DAY == TRUE){
        true_results = data$BPDgrade
        true_results[true_results == 'NA'] <- 'Death'
        #true_results <- true_results[true_results != 'NA' & !is.na(true_results)]
      } else {
        true_results[true_results == 'none'] <- 'None'
        true_results[true_results == 'oxygen'] <- 'Oxygen'
        true_results[true_results == 'vent'] <- 'Inv'
      }
      
      pr <- list(None = None_vec,
                 Oxygen = Ox_vec,
                 NIV = NIV_vec,
                 Inv = Inv_vec,
                 Death = Death_vec
      )
      
      
      library(geometry)
      library (multiROC)
      library (gtools)
      library(dummies)
      library(ggplot2)
      library(pracma)
      library(gbm)
      performance <- function(data, pr){
        results <- matrix(nrow = 1,ncol = 22)
        pr_df <- data.frame(pr)
        # calculate the performance of the model in the data
        #get AUROC
        colnames(pr_df) <- paste(colnames(pr_df), "_pred_MN")
        true_label <- dummies::dummy(true_results, sep = ".")
        true_label <- data.frame(true_label)
        colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
        colnames(true_label) <- paste(colnames(true_label), "_true")
        final_df <- cbind(true_label, pr_df)
        roc_res <- multi_roc(final_df, force_diag=T)
        #get calibration in the large and calibration slope for each of the BPDgrades and Death
        dat_none <- data.frame(e = pr_df$`None _pred_MN`, o = as.factor(true_label$`None _true`))
        dat_none$e[dat_none$e == 0] = 0.0000000001
        dat_none$e[dat_none$e == 1] = 0.9999999999
        dat_none$logite <- logit(dat_none$e)
        mfit_none = glm(formula = o~I(logite), 
                        family = binomial(link = "logit"), dat_none)
        dat_oxygen <- data.frame(e = pr_df$`Oxygen _pred_MN`, o = as.factor(true_label$`Oxygen _true`))
        dat_oxygen$e[dat_oxygen$e == 0] = 0.0000000001
        dat_oxygen$e[dat_oxygen$e == 1] = 0.9999999999
        dat_oxygen$logite <- logit(dat_oxygen$e)
        mfit_oxygen = glm(formula = o~I(logite), 
                          family = binomial(link = "logit"), dat_oxygen)
        dat_niv <- data.frame(e = pr_df$`NIV _pred_MN`, o = as.factor(true_label$`NIV _true`))
        dat_niv$e[dat_niv$e == 0] = 0.0000000001
        dat_niv$e[dat_niv$e == 1] = 0.9999999999
        dat_niv$logite <- logit(dat_niv$e)
        mfit_niv = glm(formula = o~I(logite), 
                       family = binomial(link = "logit"), dat_niv)
        dat_inv <- data.frame(e = pr_df$`Inv _pred_MN`, o = as.factor(true_label$`Inv _true`))
        dat_inv$e[dat_inv$e == 0] = 0.0000000001
        dat_inv$e[dat_inv$e == 1] = 0.9999999999
        dat_inv$logite <- logit(dat_inv$e)
        mfit_inv = glm(formula = o~I(logite), 
                       family = binomial(link = "logit"), dat_inv)
        dat_Death <- data.frame(e = pr_df$`Death _pred_MN`, o = as.factor(true_label$`Death _true`))
        dat_Death$e[dat_Death$e == 0] = 0.0000000001
        dat_Death$e[dat_Death$e == 1] = 0.9999999999
        dat_Death$logite <- logit(dat_Death$e)
        mfit_Death = glm(formula = o~I(logite), 
                         family = binomial(link = "logit"), dat_Death)
        results[1,1] <- roc_res$AUC$MN$macro
        roc_ci_res <- roc_ci(final_df, conf= 0.95, type='basic', R = 20, index = 6)
        results[1,2] <-  roc_ci_res$basic[[4]] * 2.7004/(nrow (final_df)-1)
        results[1,3] <- mfit_none$coefficients[1]
        results[1,4] <- vcov(mfit_none)[1,1]
        results[1,5] <- mfit_oxygen$coefficients[1]
        results[1,6] <- vcov(mfit_oxygen)[1,1]
        results[1,7] <- mfit_niv$coefficients[1]
        results[1,8] <- vcov(mfit_niv)[1,1]
        results[1,9] <- mfit_inv$coefficients[1]
        results[1,10] <- vcov(mfit_inv)[1,1]
        results[1,11] <- mfit_Death$coefficients[1]
        results[1,12] <- vcov(mfit_Death)[1,1]
        results[1,13] <- mfit_none$coefficients[2]
        results[1,14] <- vcov(mfit_none)[2,2]
        results[1,15] <- mfit_oxygen$coefficients[2]
        results[1,16] <- vcov(mfit_oxygen)[2,2]
        results[1,17] <- mfit_niv$coefficients[2]
        results[1,18] <- vcov(mfit_niv)[2,2]
        results[1,19] <- mfit_inv$coefficients[2]
        results[1,20] <- vcov(mfit_inv)[2,2]
        results[1,21] <- mfit_Death$coefficients[2]
        results[1,22] <- vcov(mfit_Death)[2,2]
        results2 <- as.data.frame(results)
        colnames(results2) <- c("AUROC","var_AUROC", 'citl_None','var_citl_None',"citl_Oxygen", "var_citl_Oxygen", "citl_NIV", "var_citl_NIV", "citl_Inv", "var_citl_Inv", "citl_Death", "var_citl_Death",
                                "cslope_None","var_cslope_None","cslope_Oxygen","var_cslope_Oxygen", "cslope_NIV","var_cslope_NIV", "cslope_Inv","var_cslope_Inv", "cslope_Death","var_cslope_Death")
        return (results2)
      }
      Perf_Results <- performance(true_results, pr)
      Perf_Results
      if (FINAL_DAY==TRUE){
        save(Perf_Results,file=paste0('PlotsAndResults/Perf_Results_Seen',SeenTimePoints,'_ToFinal',version,'.RData'))
      } else {
        save(Perf_Results,file=paste0('PlotsAndResults/Perf_Results_Seen',SeenTimePoints,'_To',SeenTimePoints+DaysForward,version,'.RData'))
      }
      
      
      
      
      
      #Cali plots
      #calibrate plot function
      
      transparent_green <- rgb(0, 255, 0, max = 255, alpha = 50)
      transparent_blue <- rgb(0, 0, 255, max = 255, alpha = 50)
      transparent_orange <- rgb(255 , 165, 0, max = 255, alpha = 50)
      transparent_black <- rgb(0, 0, 0, max = 255, alpha = 50)
      transparent_red <- rgb(255, 0, 0, max = 255, alpha = 50)
      
      calplot<- function (data, pr){
        library(gbm)
        calibrate.plot(ifelse(data=="None",1,0), 
                       pr$None, line.par = list(col = "green", lty = 3), shade.col=transparent_green)
        calibrate.plot(ifelse(data=="Oxygen",1,0),
                       pr$Oxygen, line.par = list(col = "blue", lty = 3), shade.col=transparent_blue, replace = FALSE)
        calibrate.plot(ifelse(data=="NIV",1,0),
                       pr$NIV, line.par = list(col = "orange", lty = 3), shade.col=transparent_orange, replace = FALSE)
        calibrate.plot(ifelse(data=="Inv",1,0),
                       pr$Inv, line.par = list(col = "black", lty = 3), shade.col=transparent_black, replace = FALSE)
        calibrate.plot(ifelse(data=="Death",1,0),
                       pr$Death, line.par = list(col = "red", lty = 3), shade.col=transparent_red, replace = FALSE)
        legend("bottomright", legend=c("No BPD", "Oxygen", "NIV", "Inv", "Death"), col=c("green","blue","orange","black", "red"), lwd=2)
        abline(a=0, b= 1, col = "grey")
      }
      print(calplot(true_results,pr))
      if (FINAL_DAY==TRUE){
        pdf(file=paste0('PlotsAndResults/CPlot_Seen',SeenTimePoints,'_ToFinal',version,'.pdf'))
        calplot(true_results,pr)
        dev.off()
        
      } else {
        pdf(file=paste0('PlotsAndResults/CPlot_Seen',SeenTimePoints,'_To',SeenTimePoints+DaysForward,version,'.pdf'))
        calplot(true_results,pr)
        dev.off()
        
      }
      
      
      
      
      
      
      #####ROC stuff
      
      tp=array(NA,c(9999,5))
      fn=array(NA,c(9999,5))
      fp=array(NA,c(9999,5))
      tn=array(NA,c(9999,5))
      tpr=array(NA,c(9999,5))
      fpr=array(NA,c(9999,5))
      
      
      outcome_asarray <- c('None','Oxygen','NIV','Inv','Death')
      for (percentile in 1:9999){
        for (outcome in 1:5){
          tp[percentile,outcome] <- sum(pr[[outcome_asarray[outcome]]] > (percentile/10000) & true_results == outcome_asarray[outcome])
          fn[percentile,outcome] <- sum(pr[[outcome_asarray[outcome]]] < (percentile/10000) & true_results == outcome_asarray[outcome])
          fp[percentile,outcome] <- sum(pr[[outcome_asarray[outcome]]] > (percentile/10000) & true_results != outcome_asarray[outcome])
          tn[percentile,outcome] <- sum(pr[[outcome_asarray[outcome]]] < (percentile/10000) & true_results != outcome_asarray[outcome])
          
          tpr[percentile,outcome] <- tp[percentile,outcome]/(tp[percentile,outcome]+fn[percentile,outcome])
          fpr[percentile,outcome] <- fp[percentile,outcome]/(fp[percentile,outcome]+tn[percentile,outcome])
        }
      }
      
      
      
      
      
      #ROC plot
      Outcomes <- c("None" = "green",
                    "Oxygen" = "blue",
                    "NIV" = "orange",
                    "INV" = "black",
                    "Death" = 'red')
      ggplot()+
        geom_line(aes(x=c(1,fpr[,1],0),y=c(1,tpr[,1],0),color = 'None'), ,size=2)+
        geom_line(aes(x=c(1,fpr[,2],0),y=c(1,tpr[,2],0),color = 'Oxygen'), ,size=2)+
        geom_line(aes(x=c(1,fpr[,3],0),y=c(1,tpr[,3],0),color = 'NIV'), ,size=2)+
        geom_line(aes(x=c(1,fpr[,4],0),y=c(1,tpr[,4],0),color = 'INV'), ,size=2)+
        geom_line(aes(x=c(1,fpr[,5],0),y=c(1,tpr[,5],0),color = 'Death'), ,size=2)+
        labs(x='False Positive Rate', y='True Positive Rate',color='Outcomes')+
        theme(text = element_text(size = 25),legend.position="right")+
        scale_color_manual(breaks = c('None','Oxygen','NIV','INV','Death'),values = Outcomes)
      if (FINAL_DAY==TRUE){
        ggsave(paste0('PlotsAndResults/AUCCurve_Seen',SeenTimePoints,'_ToFinal',version,'.png'), width=9, height=6.7, dpi=300)
      } else {
        ggsave(paste0('PlotsAndResults/AUCCurve_Seen',SeenTimePoints,'_To',SeenTimePoints+DaysForward,version,'.png'), width=9, height=6.7, dpi=300)
      }
      
      #plot(c(1,fpr[,time],0),c(1,tpr[,time],0),type='l',xlim=c(0,1),ylim=c(0,1))
      auroc = array(NA,6)
      auroc[1] = trapz(rev(c(1,fpr[,1],0)),rev(c(1,tpr[,1],0)))#AUROC score
      auroc[2] = trapz(rev(c(1,fpr[,2],0)),rev(c(1,tpr[,2],0)))#AUROC score
      auroc[3] = trapz(rev(c(1,fpr[,3],0)),rev(c(1,tpr[,3],0)))#AUROC score
      auroc[4] = trapz(rev(c(1,fpr[,4],0)),rev(c(1,tpr[,4],0)))#AUROC score
      auroc[5] = trapz(rev(c(1,fpr[,5],0)),rev(c(1,tpr[,5],0)))#AUROC score
      auroc[6] = mean(auroc[1:5])
      auroc
      if (FINAL_DAY==TRUE){
        save(auroc,file=paste0('PlotsAndResults/Auroc_Seen',SeenTimePoints,'_ToFinal',version,'.RData'))
      } else {
        save(auroc,file=paste0('PlotsAndResults/Auroc_Seen',SeenTimePoints,'_To',SeenTimePoints+DaysForward,version,'.RData'))
      }
      
      library(mlr3measures)
      brierscore = mbrier(as.factor(true_results),do.call(cbind, pr))
      if (FINAL_DAY==TRUE){
        save(brierscore,file=paste0('PlotsAndResults/Brier_Seen',SeenTimePoints,'_ToFinal',version,'.RData'))
      } else {
        save(brierscore,file=paste0('PlotsAndResults/Brier_Seen',SeenTimePoints,'_To',SeenTimePoints+DaysForward,version,'.RData'))
      }
      
    }
  }
}



