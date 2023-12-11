library(lavaan)
library(plyr)


times=8
fits_random <- datagen(model = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
  Eta1 ~ 1
  simuvar1|0*t1', rmsea_cutoff = .05,ID=500,pred_cutoff = .7,rel_cutoff=.6,times=10) 

#############################################################################
#######################       Daten erstellen      ##########################
#############################################################################
times = length(fits_random[[1]])
ID = 500
simu=data.frame()
tms=1:times



for(i in tms) {
  
  dt = as.data.frame(fits_random[["latvars"]][[i]])
  colnames(dt) <- "true_scores"
  dt = cbind(dt,fits_random[["scores"]][[i]])
  dt = cbind(dt,fits_random[["data"]][[i]])
  
  dt$id <- as.integer(1:ID) + (i-1)*ID
  
  dt$sample <- rep(i,ID)
  dt[,paste0("num",i)] <- round(runif(ID,min=1,max=50),2)
  dt[,paste0("cat",i)] <- replicate(ID,sample(c(1,3,5),1))
  dt[,paste0("ord",i)] <- replicate(ID,sample(c(4,5),1))
  
  tmss = tms[-i]
  for(j in tmss){ #erstmal alle Werte fuer andere Vars
    dt[,paste0("num",j)] <- round(runif(ID,min=0,max=200),2)
    dt[,paste0("cat",j)] <- replicate(ID,sample(c(1,2,3,4,5),1))
    dt[,paste0("ord",j)] <- replicate(ID,sample(c(1,2,3,4,5),1))
    
    for(k in 1:ID){
      if( (dt[k,paste0("num",j)] <= 50) &  (dt[k,paste0("cat",j)] %in% c(1,3,5)) ) {dt[k,paste0("ord",j)] = sample(c(1,2,3),1)} 
    }
    
  }
  

  simu <- rbind.fill(simu,dt)
}



##### Prepocessing 

for(i in tms){
  simu[,paste0("num",i)] <- as.numeric(simu[,paste0("num",i)] )
  simu[,paste0("cat",i)] <- as.factor(simu[,paste0("cat",i)] )
  simu[,paste0("ord",i)] <- factor(simu[,paste0("ord",i)], ordered=TRUE, levels=c(1,2,3,4,5) )
}
#for(i in 1:5){
#  simu[,paste0("simuvar",i)] <- factor(simu[,paste0("simuvar",i)], ordered=TRUE, levels=c(1,2,3,4,5,6,7) )
#}

colnames(simu)[which(colnames(simu)=="Eta1")] = "pred_scores"
#shuffle rows
simu = simu[sample(1:nrow(simu)), ]




##### Es hat geklappt
save(simu,fits_random,file="/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/simu_data_univ5.RData")











#############################################################################
##############################   Analyze     ################################
#############################################################################



#correlation whole dataset
cor(simu$true_scores,simu$pred_scores,method="spearman")


#quality of predictions for whole dataset

output <- c("simuvar1","simuvar2","simuvar3","simuvar4","simuvar5")
model= 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
                Eta1 ~ 1
                simuvar1|0*t1'


fit_univ <- lavaan::cfa(model = model, data=simu[,c(3:7)], ordered = output, estimator="WLS", do.fit=T, std.lv=F)
#semTools::compRelSEM(fit_univ)
pred <- lavPredict(fit_univ)
cor(simu$true_scores,pred,method="spearman")






##### individually & whole dataset
#whole <- data.frame()
for(i in 1:times){
  #whole <- rbind(whole,get(paste0("goodfits",i)))
  fit_univ <- lavaan::cfa(model = model, data=get(paste0("goodfits",i)), ordered = output, estimator="WLSMV", do.fit=T, std.lv=F)
  pred <- lavPredict(fit_univ)
  print(mean(pred))
  #print(var(pred))
  #print(fitMeasures(fit_univ)["rmsea"])
  #print(summary(fit_univ, fit.measures=T))
}
#fit_univ <- lavaan::cfa(model = model, data=whole, ordered = output, estimator="WLSMV", do.fit=T, std.lv=F)
#fitMeasures(fit_univ)["rmsea"] #fits NOT well --> good

