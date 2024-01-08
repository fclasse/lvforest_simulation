setwd("/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional")


datagen <- function(times=3,ID=500,schwellen=4,items=7,rmsea_cutoff=.05){ 
  
  model = paste(sapply(1:items, function(x) paste0("simuvar",x," +")), collapse = "")
  model = paste0("Eta1 =~ ",substr(model,1,nchar(model)-2))
  model = paste0(model,"\n  Eta1 ~ 1\n  simuvar1|0*t1")
  
  randparams <- function(N,num.min,rnd=2,num.max=0.99) {
    random.numbers <- round(runif(min=num.min,max=num.max, N),rnd)
    random.signs   <- sample(c(1,-1), N, replace=TRUE)      
    return(random.numbers * random.signs)
  }
  
  schwellen_probs <- function(schwellen){
    if(schwellen == 6) return(c(0.10,0.25,0.40,0.60,0.75,0.90)  )
    if(schwellen == 5) return(c(0.10,0.30,0.50,0.70,0.90) )
    if(schwellen == 4) return(c(0.15,0.40,0.60,0.85)  )
    if(schwellen == 3) return(c(0.25,0.5,0.75)  )
    if(schwellen == 2) return(c(0.30,0.70)  )
    if(schwellen == 1) return(c(0.50))
    else stop("Error: Function can only produce data with 2 to 7 categories.")
  }
  
  library(plyr)
  library(psych)
  library(lavaan)
  library(resample)
  
  for(i in c("saved_pvalues","saved_rmseas","saved_kappas","saved_betas","saved_data","saved_vars","saved_latvars","saved_predcors","saved_scores","saved_means")){assign(i,list())}
  
  #####################################################################################
  ################# Algorithmus
  #####################################################################################
  h=1;g=1
  while(h<(times+1)){ 
    l <- schwellen #Anzahl an Schwellen
    m <- items #Anzahl an Items pro Level-1 latente Variable
    nlatvar = 1
    
    ## Simu 1: Alle Parameter unterschiedlich, keine random data 
    kappa_shift <- c(0,round(rnorm(n = m-1, mean=0, sd = 1),5)) #"Erwartungswert"-Shift pro Item --> greader sd means greater difference between sc_ind und sc_univ!
    perz_kappa <- matrix(NA,ncol = l, nrow = m)#Perzentile, um die kappa-Parameter zu bestimmen
    for(i in 1:(m)){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.09, l) ,5)}
    beta <- c(1,round(runif(min=0.1,max=2, m-1),5)) #Discrimination parameter
    var <- round(runif(min=0.05,max=2, 1),5)
    
    #Variable erstellen
    Psi <- rnorm(ID,0,sqrt(var))  #psi erstmal mittelwert 0
    
    
    
    
    ## Kappa Matrix erstellen
    kappa <- matrix(ncol = m,nrow = l)  #ncol = L-1 Variablen * Items ; nrow = schwellen
    for(i in 1:m){ # items
      quant <- as.vector(quantile(beta[i]*Psi , probs= perz_kappa[i,])) - kappa_shift[i] 
      kappa[,i] <- quant
    }
    kappa <- round(kappa,2)
    
    
    #Jetzt neues kappa und neues Psi erstellen
    mean = kappa[1,1]*-1
    kappa = kappa + matrix(mean,l,m)
    Psi <- Psi + mean
    
    ## Daten generieren
    probitY <- array(NA,c(ID,l,m)) #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    PY <-  array(NA,c(ID,l,m)) 
    Pis <-  array(NA,c(ID,l+1,m))
    Y <- array(NA,c(ID,m));v=1 #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    
    
    for (i in 1:m) { #Items 
      for (k in 1:l) { #Schwellen
        probitY[,k,i] <-  beta[v]*Psi  - kappa[k,v] 
        PY[,k,i] <- VGAM::probitlink(probitY[,k,i],inverse=T) 
        if (k == 1) Pis[,k,i] <- 1 - PY[,k,i] #P(Y=0)
        if (k > 1) Pis[,k,i] <- PY[,k-1,i] - PY[,k,i] 
        if (k == l) Pis[,k+1,i] <- PY[,k,i]
      }
      for(j in 1:ID){Y[j,i] <- sample(1:(l+1), size=1, replace=TRUE, prob=sapply(1:(l+1), function(x) Pis[j,x,i])  )} 
      v <- v+1
    }
    
    
    
    #### Tabelle erzeugen & Deskriptive Statistik
    table <- matrix(ncol = m,nrow = ID)
    v=1
    
    for (i in 1:m){
      table[,v] <- Y[,i]
      v <- v+1
    }
    
    
    
    ############################           Ueberpruefen            ##########################################
    
    output <- sapply(1:(m),function(y){paste0("simuvar",y)})
    table <- as.data.frame(table)
    colnames(table) <- c(output)
    
    
    
    print(paste("Iteration ",g)) #wieviele Iterationen braucht es?
    g=g+1
    
    
    ## Model schätzen
    fit_ord <- tryCatch({lavaan::cfa(model = model, data=table, ordered = output, estimator="WLS", do.fit=T, std.lv=F, control=list(iter.max=1000))},warning = function(w){NA},error = function(e){NA})
    ##Ergebnisse Zusammenfügen!
    if( suppressWarnings(!is.na(fit_ord))  ){
      
      
      if(fitMeasures(fit_ord)["rmsea"] < rmsea_cutoff){ 
        
        #####Predictions
        scores <- lavaan::lavPredict(fit_ord, method="EBM")
        predcors <- cor(Psi,scores,method="spearman")
        #print(predcors)
        #print(paste("Mean Pred:",round(mean(scores),3),"True:",round(mean(Psi),3)))
        #print(paste("Variance Pred:",round(var(scores),3),"True:",round(var(Psi),3)))
        
        
        saved_pvalues[[h]] <- fitMeasures(fit_ord)["pvalue"]; names(saved_pvalues)[[h]] <- paste0("pvalue",h)
        saved_rmseas[[h]] <- fitMeasures(fit_ord)["rmsea"]; names(saved_rmseas)[[h]] <- paste0("rmsea",h)
        saved_vars[[h]] <- var; names(saved_vars)[[h]] <- paste0("var",h)
        saved_means[[h]] <- mean; names(saved_means)[[h]] <- paste0("means",h)
        saved_kappas[[h]] <- kappa; names(saved_kappas)[[h]] <- paste0("kappa",h)
        saved_betas[[h]] <- beta; names(saved_betas)[[h]] <- paste0("beta",h)
        saved_scores[[h]] <-scores; names(saved_scores)[[h]] <- paste0("pred_scores",h)
        saved_latvars[[h]] <- Psi; names(saved_latvars)[[h]] <- paste0("true_var",h)
        saved_predcors[[h]] <- predcors; names(saved_predcors)[[h]] <- paste0("predcor",h)
        saved_data[[h]] <- table; names(saved_data)[[h]] <- paste0("data",h)
        
        
        h=h+1 
        
      } else {next}
    } else {next}
  } 
  
  fits <- list(saved_pvalues,saved_rmseas,saved_vars,saved_means,saved_kappas,saved_betas,saved_scores,saved_latvars,saved_predcors,saved_data) 
  names(fits) <- c("pvalues","rmseas","vars","means","kappas","betas","scores","latvars","predcors","data")
  return(fits)
}



