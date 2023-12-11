################################################################################
################################  SEMTree ######################################
################################################################################
library(devtools)
library(lavaan) #semtree funktioniert sonst nicht
devtools::install_github("brandmaier/semtree",force=TRUE, build_opts = c()) #aktuellstes Package von github runterladen!
install.packages("doParallel")
install.packages("doRNG")
install.packages("caret")

source("/dss/dsshome1/0B/ra35tik2/paper2/lvforest/lvforest_semtree.R")


start_time <- Sys.time()
univ_forest <- lvforest(
  data = simu,
  split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5',
  pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
                Eta1 ~ 1
                simuvar1|0*t1',
  input = c("num1","cat1","ord1","num2","cat2","ord2","num3","cat3","ord3","num4",
            "cat4","ord4","num5","cat5","ord5","num6","cat6","ord6","num7","cat7",  
            "ord7","num8","cat8","ord8","num9","cat9","ord9","num10","cat10","ord10"),  
  idvar = "id",
  ntrees= 10000,
  bagging=NULL,
  cutoff_rmsea = .05,
  ordered = TRUE,
  minn = 200,
  split=2,   #sqrt(length(input))
  estmethod = "EBM",
  stab_bonferroni=FALSE
)
end_time <- Sys.time()
runtime <- end_time - start_time


save(univ_forest,runtime,file="/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/simu_results_univ10.RData")















################################################################################
############################  Single SEMTree ###################################
################################################################################
input = c("num1","cat1","ord1","num2","cat2","ord2","num3","cat3","ord3","num4",
          "cat4","ord4","num5","cat5","ord5")  #,"num6","cat6","ord6","num7","cat7","ord7","num8","cat8","ord8","num9","cat9","ord9","num10","cat10","ord10")
split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5'

semtree_control=semtree::semtree.control(min.bucket = 200, mtry=NA, method="score", bonferroni=TRUE,progress.bar = TRUE)
fit_num <- lavaan::cfa(model = split_model, data=simu, std.lv=FALSE, meanstructure=TRUE, estimator = "ML", do.fit=FALSE) 

start_time <- Sys.time()
tree <- semtree::semtree(model=fit_num, data=simu, predictors=input, control=semtree_control) 
end_time <- Sys.time()
runtime_tree <- end_time - start_time



#getNumNodes, getTermiNodes, getRules, printFact aus lvforest_semtree.R ziehen
ordered = c("simuvar1","simuvar2","simuvar3","simuvar4","simuvar5")
pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
                Eta1 ~ 1
                simuvar1|0*t1'
fit_ord <- scores <- list()
scores <- data.frame()
ni <- sort(as.vector(getTermiNodes(tree))) 
rls <- unlist(getRules(tree)[ni]) 

for(i in 1:length(ni)){
  data_refit <- subset(simu,eval(parse(text=rls[i])) )        #richtiger Terminal node
  fit_ord[[i]] <- lavaan::cfa(model = pred_model, data = data_refit, ordered = ordered, estimator = "WLS", std.lv = FALSE,control=list(iter.max=1000))
  names(fit_ord)[[i]] <- paste0("node",ni[i])
  scores <- rbind(scores,cbind(data_refit$id,lavaan::lavPredict(fit_ord[[i]], method='EBM')))
}
colnames(scores) <- c("id","Eta1")



save(tree,runtime_tree,scores,file="/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/simu_results_univ_tree5.RData")




#plot.semtree aus github kopieren...
nodeFunSemtree<-function(x, labs, digits, varlen)
{
  paste("n=",x$frame$n)
}
plot.semtree(tree)






################################################################################
######################### Accurate Predictions?  ###############################
################################################################################

#whole model
sc_whole <- lavaan::lavPredict(lavaan::cfa(model = pred_model, data = simu, ordered = ordered, estimator = "WLS", std.lv = FALSE,control=list(iter.max=1000)), method='EBM')
sc_whole <- as.data.frame(cbind(simu$id,sc_whole))
colnames(sc_whole) <- c("id","Eta1")

#####Correlation with true scores
tr <-simu[,c("id","true_scores")]

#individual vs. forest vs. tree vs. whole model
sctr <- tr[,2]
scind <- simu$pred_scores
scfor <- univ_forest[["Scores"]][["Eta1"]][match(tr$id,univ_forest[["Scores"]][["Eta1"]]$id),"lvscore"]
sctre <- scores[match(tr$id,scores$id),"Eta1"]
scwho <- sc_whole[match(tr$id,sc_whole$id),"Eta1"]

scc <- as.data.frame(cbind(sctr,scind,scfor,sctre,scwho))

library(cocor)
cocor.result_for <- cocor(~scind + sctr | scfor + sctr,scc)
cocor.result_tre <- cocor(~scind + sctr | sctre + sctr,scc)
cocor.result_who <- cocor(~scind + sctr | scwho + sctr,scc)

cor(sctr,scind,method='pearson', use="pairwise.complete.obs")







###



##Alle mindestens einmal vorgekommen?
sum(is.na(univ_forest[["Scores"]][["Eta1"]]$lvscore))



####### Analyze NAs
naids <- univ_forest[["Scores"]][["Eta1"]][is.na(univ_forest[["Scores"]][["Eta1"]]$lvscore),"id"]
table(simu[simu$id %in% naids,"sample"]) #in erster Linie sample 1 --> 



####### Analyze subgroups
for(i in 1:10){ #nrow(univ_forest[["Bestnodes"]])
  print(  table( subset(simu,eval(parse(text=univ_forest[["Bestnodes"]]$decision_rule[i])) )$sample  )  )
}

