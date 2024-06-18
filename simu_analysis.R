################################################################################
################################  SEMTree ######################################
################################################################################
library(devtools)
devtools::install_github("brandmaier/semtree",force=TRUE, build_opts = c()) #aktuellstes Package von github runterladen!
install.packages("doParallel")
install.packages("doRNG")
install.packages("caret")

source("/dss/dsshome1/0B/ra35tik2/paper2/lvforest/lvforest_semtree.R")


################################################################################
#################################  Forest ######################################
################################################################################
library(foreach)
results_forest_part <- foreach(j=1, .errorhandling="pass") %do% { #.export=c("datagen"),
  
  simu = datasets[[j]][[1]]
  input_sampled = sample(c("num1","ord1","num2","ord2","num3","ord3"),3)
  
  start_time <- Sys.time()
  univ_forest <- lvforest(
    data = simu,
    split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8',
    pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8
                Eta1 ~ 1
                simuvar1|0*t1',
    input = input_sampled,  
    idvar = "id",
    ntrees= 20,
    bagging=NULL,
    cutoff_rmsea = .05,
    ordered = TRUE,
    minn = 200,
    split=2,   #sqrt(length(input))
    estmethod = "EBM",
    stab_bonferroni=FALSE
  )
  end_time <- Sys.time()
  runtime_forest <- end_time - start_time
  
  #write result in simu
  simu$forest_scores_part = univ_forest[["Scores"]][["Eta1"]][match(simu$id,univ_forest[["Scores"]][["Eta1"]]$id),"lvscore"]
  simu = simu[,c(ncol(simu),1:(ncol(simu)-1))]
  print(paste("Iteration",j,"finished"))
  return(list(simu,runtime_forest))
  
}




################################################################################
############################  Single SEMTree ###################################
################################################################################
library(lavaan) #semtree funktioniert sonst nicht
library(foreach)

results_trees <- foreach(j=1, .errorhandling="pass") %do% {
  
  simu = results_forest[[j]][[1]]
  
  input = c("num1","ord1","num2","ord2","num3","ord3")
  cols = c("simuvar1","simuvar2","simuvar3","simuvar4","simuvar5","simuvar6","simuvar7","simuvar8",input)
  split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8'
  
  semtree_control=semtree::semtree.control(min.bucket = 200, mtry=NA, method="score", bonferroni=TRUE,progress.bar = FALSE)
  fit_num <- lavaan::cfa(model = split_model, data=simu, std.lv=FALSE, meanstructure=TRUE, estimator = "ML", do.fit=FALSE) 
  
  start_time <- Sys.time()
  tree <- semtree::semtree(model=fit_num, data=simu[,cols], predictors=input, control=semtree_control) 
  end_time <- Sys.time()
  runtime_tree <- end_time - start_time
  
  
  
  #getNumNodes, getTermiNodes, getRules, printFact aus lvforest_semtree.R ziehen
  pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
                Eta1 ~ 1
                simuvar1|0*t1'
  fit_ord <- scores <- list()
  scores <- data.frame()
  ni <- sort(as.vector(getTermiNodes(tree))) 
  rls <- unlist(getRules(tree)[ni]) 
  
  for(i in 1:length(ni)){
    data_refit <- subset(simu,eval(parse(text=rls[i])) )        #richtiger Terminal node
    fit_ord[[i]] <- try(lavaan::cfa(model = pred_model, data = data_refit, ordered = T, estimator = "WLS", std.lv = FALSE,control=list(iter.max=1000)), silent=T  )
    names(fit_ord)[[i]] <- paste0("node",ni[i])
    if(!is.character(fit_ord[[i]])){  
      scc <- cbind(data_refit$id,lavaan::lavPredict(fit_ord[[i]], method='EBM'))
      colnames(scc) <- c("id","Eta1")
      scores <- rbind(scores,scc)  
    } else {
      scna = cbind(data_refit$id,NA)
      colnames(scna) <- c("id","Eta1")
      scores = rbind(scores,scna)
    }
  }
  #colnames(scores) <- c("id","Eta1")
  simu$tree_scores <- scores[match(simu$id,scores$id),"Eta1"]
  simu = simu[,c(ncol(simu),1:(ncol(simu)-1))]
  
  return(list(simu,runtime_tree))
  
}


#plot.semtree aus github kopieren...
nodeFunSemtree<-function(x, labs, digits, varlen)
{
  paste("n=",x$frame$n)
}

plot.semtree(tree) #exported as pdf w=6.27,h=9.69



#get results of models in terminal nodes
for(i in 1:16){
  print(   lavaan::fitMeasures(fit_ord[[i]])["rmsea"] )
  print(fit_ord[[i]]@Data@nobs[[1]])
  }








################################################################################
############################  Joint results ####################################
################################################################################



##plot results
library(ggplot2);library(reshape2);library(latex2exp)

#conv plot
conv_forest_part = sapply(datasets, function(y){sum(is.na(y$forest_part_scores))/1500})
conv_forest = sapply(datasets, function(y){sum(is.na(y$forest_scores))/1500})
conv_univ = sapply(datasets, function(y){sum(is.na(y$univ_scores))/1500})
conv_ind = sapply(datasets, function(y){sum(is.na(y$pred_scores))/1500})
conv_trees = sapply(datasets, function(y){sum(is.na(y$tree_scores))/1500})
dat <- cbind(conv_univ,conv_trees,conv_forest_part,conv_forest,conv_ind)
dat <- melt(dat)
colnames(dat) <- c("counter","mode","value")


ggplot(dat,aes(x=mode, y=value, fill=mode)) + geom_boxplot() +theme(legend.position="none",axis.title.x=element_blank()) +
  scale_x_discrete(labels=c("naive model", "single tree", "LV Forest with\n incomplete\n partitioning variables","LV Forest","distinct models")) +
  scale_y_continuous(labels = scales::percent) + ylab("Non-convergence rate")
ggsave(file="res/conv_plot.pdf", width = 190, height = 120, units = "mm")



#acc plot
acc_forest_part = sapply(datasets, function(y){cor(y$forest_part_scores,y$true_scores, method="spearman", use="pairwise.complete.obs")})
acc_forest = sapply(datasets, function(y){cor(y$forest_scores,y$true_scores, method="spearman", use="pairwise.complete.obs")})
acc_univ = sapply(datasets, function(y){cor(y$univ_scores,y$true_scores, method="spearman", use="pairwise.complete.obs")})
acc_ind = sapply(datasets, function(y){cor(y$pred_scores,y$true_scores, method="spearman", use="pairwise.complete.obs")})
acc_trees = sapply(datasets, function(y){cor(y$tree_scores,y$true_scores, method="spearman", use="pairwise.complete.obs")})
dat <- cbind(acc_univ,acc_trees,acc_forest_part,acc_forest,acc_ind)
dat <- melt(dat)
colnames(dat) <- c("counter","mode","value")


ggplot(dat,aes(x=mode, y=value, fill=mode)) + geom_boxplot() +theme(legend.position="none",axis.title.x=element_blank()) +
  scale_x_discrete(labels=c("naive model", "single tree", "LV Forest with\n incomplete\n partitioning variables","LV Forest","distinct models")) +
  ylab("Correlation latent var. scores w/ true scores")
ggsave(file="res/acc_plot.pdf", width = 190, height = 120, units = "mm")

#ehich scale to use?
#scale_x_discrete(labels=c("naive model", "single tree", "LV Forest with\n incomplete\n partitioning variables","LV Forest","distinct models")) 
#scale_x_discrete(labels=c(TeX("$\\hat{\\eta}_{naive}$"), TeX("$\\hat{\\eta}_{SEMTree}$"), TeX("$\\hat{\\eta}_{part.LVForest}$"),TeX("$\\hat{\\eta}_{LVForest}$"),TeX("$\\hat{\\eta}_{ind.models}$"))) +


################################################################################
##########################  Computation time ###################################
################################################################################

comptime_for = sapply(results_forest,"[[",2);print(mean(comptime_for))
comptime_forpar = sapply(results_forest_part,"[[",2);print(mean(comptime_forpar))

comptime_tree = lapply(results_trees,"[[",2) #;print(mean(comptime_tree))
mean(   unlist(comptime_tree[sapply(comptime_tree, class) != "call"])   )



