library(ExtDist)
library(rCausalMGM)

setwd('')

select_continuous_var<-function(dataset){
  maxCat = 4
  continuous_features<-c()
  for (i in 1:ncol(dataset)){
    feature<-dataset[[i]]
    if (length(unique(feature))/length(feature) > 0.05 & length(unique(feature)) > maxCat){
      continuous_features<-c(continuous_features,i)
    } 
  }
  return(continuous_features)
}




simulate_data<-function(data_subset,contunious_var){
  data_subset_con<-data_subset[,contunious_var]
  data_subset_dis<-data_subset[,-contunious_var]
  
  loadings<-rnorm(200*10)
  loadings[abs(loadings)<quantile(abs(loadings),probs = 0.95)]=0
  # hist(loadings)
  loadings<-matrix(loadings,nrow=10)
  
  data_subset_con<-as.matrix(data_subset_con)
  simulated_data<-data_subset_con %*% loadings + rnorm(1000*200)
  
  simulation_list<-list("continuous_data" =data_subset_con, 
                        "data_subset_dis" = data_subset_dis,
                        "loading" = loadings,
                        "simulated_data" = simulated_data,
                        "contunious_var" =colnames(data_subset)[contunious_var] )
  return(simulation_list)
}


dataset<-read.table('./save/data/data.1.txt',header = T)
contunious_var<-select_continuous_var(dataset)
set.seed(666)
folds<-caret::createFolds(1:nrow(dataset), k = 10, list = TRUE, returnTrain = F)
data_subset<-dataset[folds[[1]],]
ig <- rCausalMGM::steps(data_subset, maxDiscrete=5, g=0.1,
            lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=2)), verbose=T)

simulation_1<-simulate_data(data_subset,contunious_var)
saveRDS(simulation_1,'simulation_1.rds')


data_subset<-dataset[folds[[3]],]
ig <- rCausalMGM::steps(data_subset, maxDiscrete=5, g=0.1,
                        lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=2)), verbose=T)
simulation_3<-simulate_data(data_subset,contunious_var)
saveRDS(simulation_3,'simulation_3.rds')


data_subset<-dataset[folds[[4]],]
ig <- rCausalMGM::steps(data_subset, maxDiscrete=5, g=0.1,
                        lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=2)), verbose=T)
simulation_4<-simulate_data(data_subset,contunious_var)
saveRDS(simulation_4,'simulation_4.rds')


data_subset<-dataset[folds[[5]],]
ig <- rCausalMGM::steps(data_subset, maxDiscrete=5, g=0.1,
                        lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=2)), verbose=T)
simulation_5<-simulate_data(data_subset,contunious_var)
saveRDS(simulation_5,'simulation_5.rds')


data_subset<-dataset[folds[[7]],]
ig <- rCausalMGM::steps(data_subset, maxDiscrete=5, g=0.1,
                        lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=2)), verbose=T)
simulation_2<-simulate_data(data_subset,contunious_var)
saveRDS(simulation_2,'simulation_2.rds')
