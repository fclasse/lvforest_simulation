install.packages("doParallel")
install.packages("doRNG")

library(doRNG)
library(lavaan)
library(plyr)

#dategen einlesen
source("/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/lvforest_simulation/lvforest_univ_simu.R")

ncores <- parallelly::availableCores()
cl <- parallel::makeCluster(spec=ncores) 
doParallel::registerDoParallel(cl)
datasets <- foreach::foreach(j=1:100, .packages=c("lavaan"),  .errorhandling="pass") %do% { #.export=c("datagen"),
  
  
  times = 3
  ID = 500
  items = 8
  schwellen = 4
  
  
  ######## samples erzeugen
  fits_random <- datagen( 
    items = items,
    schwellen = schwellen,
    ID=ID,
    times=times) 
  
  
  #######   Datensatz erstellen    
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
    #dt[,paste0("cat",i)] <- replicate(ID,sample(c(1,3,5),1))
    dt[,paste0("ord",i)] <- replicate(ID,sample(c(4,5),1))
    
    tmss = tms[-i]
    for(j in tmss){ #erstmal alle Werte fuer andere Vars
      dt[,paste0("num",j)] <- round(runif(ID,min=0,max=200),2)
      #dt[,paste0("cat",j)] <- replicate(ID,sample(c(1,2,3,4,5),1))
      dt[,paste0("ord",j)] <- replicate(ID,sample(c(1,2,3,4,5),1))
      
      for(k in 1:ID){
        #if( (dt[k,paste0("num",j)] <= 50) &  (dt[k,paste0("cat",j)] %in% c(1,3,5)) ) {dt[k,paste0("ord",j)] = sample(c(1,2,3),1)} 
        if( dt[k,paste0("num",j)] <= 50 ) {dt[k,paste0("ord",j)] = sample(c(1,2,3),1)} 
      }
      
    }
    
    
    simu <- rbind.fill(simu,dt)
  }
  
  
  ##### Prepocessing 
  for(i in tms){
    simu[,paste0("num",i)] <- as.numeric(simu[,paste0("num",i)] )
    #simu[,paste0("cat",i)] <- as.factor(simu[,paste0("cat",i)] )
    simu[,paste0("ord",i)] <- factor(simu[,paste0("ord",i)], ordered=TRUE, levels=c(1,2,3,4,5) )
  }
  
  colnames(simu)[which(colnames(simu)=="Eta1")] = "pred_scores"
  #shuffle rows
  simu = simu[sample(1:nrow(simu)), ]
  
  
  #####  #quality of predictions for whole dataset
  model = paste(sapply(1:items, function(x) paste0("simuvar",x," +")), collapse = "")
  model = paste0("Eta1 =~ ",substr(model,1,nchar(model)-2))
  model = paste0(model,"\n  Eta1 ~ 1\n  simuvar1|0*t1")
  
  
  fit_univ <- lavaan::cfa(model = model, data=simu, ordered = T, estimator="WLS", do.fit=T, std.lv=F)
  simu$univ_scores <- lavPredict(fit_univ)[,1]
  simu <- simu[,c(ncol(simu),1:(ncol(simu)-1))]
  
  #write
  return(  list(simu, fitMeasures(fit_univ)["rmsea"],semTools::compRelSEM(fit_univ))  )
  
}

save(datasets,file="/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/lvforest_simulation/231213_lvfor_datasets.RData")

################################################################################
###############################  Analyse  ######################################
################################################################################

datasets = datasets[sapply(lapply(datasets , "[[", 2),is.numeric)]
converge_freq = length(datasets)/100

rmseas = sapply(datasets , "[[", 2)
rels = sapply(datasets , "[[", 3)
corrs = sapply(datasets, function(y){cor(  y[[1]]$univ_scores, y[[1]]$pred_scores, method = "spearman")  })
plot(corrs,rels)
cor.test(corrs,rels)  #look at rels to see if lvforest is going to be helpful or not!
