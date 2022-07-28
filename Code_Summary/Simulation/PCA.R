library(stats)
setwd('')
source('R/MCC.R')

filenames<-c('simulation_1','simulation_2','simulation_3','simulation_4','simulation_5')
info<-matrix(nrow = 5,ncol = 2)

for (i in seq(length(filenames))){
  file_name<-filenames[i]
  simulation<-readRDS(paste0(file_name,".rds" )) 
  dir.create(paste0("./",file_name))
  Y<-simulation$simulated_data
  
  Y_pc<-prcomp(t(Y),scale. = T)
  var_explain<-summary(Y_pc)$importance[2,]
  var_ratio<-c()
  for(j in 1:(length(var_explain)-1)){
    var_ratio<-c(var_ratio,var_explain[j]/var_explain[j+1])
  }
  
  nrank_pc<-which.max(var_ratio[1:200])
  
  info[i,1]<-nrank_pc
  Y_pcs<-Y_pc$rotation[,1:nrank_pc]
  
  dir.create(paste0("./",file_name,"/pc"))
  saveRDS(Y_pc,paste0("./",file_name,"/pc/pc.rds"))
  
  LL<-simulation$continuous_data
  info[i,2]<-MCC(Y_pcs,LL,paste0("./",file_name,"/pc"))

}
