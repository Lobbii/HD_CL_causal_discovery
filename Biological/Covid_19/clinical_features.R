setwd('Covid')
library(GEOquery)
gds <- getGEO("GSE157103")
dim(gds[["GSE157103_series_matrix.txt.gz"]]@phenoData)
clinical_data<-gds[["GSE157103_series_matrix.txt.gz"]]@phenoData@data
saveRDS(clinical_data,"clinical_data.rds")

clinical_data<-readRDS('./clinical_data.rds')
table(clinical_data$batc)


grep('characteristics',colnames(clinical_data))
Cli_DF<-clinical_data[,c(37,grep('characteristics',colnames(clinical_data)))]
rownames(Cli_DF)<-Cli_DF[,1]
Cli_DF<-Cli_DF[,-1]
Cli_DF<-Cli_DF[,1:17]



# remove features with missing value 
Cli_DF<-readRDS('Cli_DF.rds')

index<-Cli_DF[,2] !=':'
Cli_DF<-Cli_DF[index,]

index<-Cli_DF[,3] !='unknown'
Cli_DF<-Cli_DF[index,]

index<-Cli_DF[,9] !=' unknown'
Cli_DF<-Cli_DF[index,]

# 51 unkown apacheii
Cli_DF$apacheii<-NULL

# Ferritin 14 unknwon
Cli_DF$ferritin<-NULL

# crp 17 unkown 
Cli_DF$crp<-NULL

#ddimer 23 unkown 
Cli_DF$ddimer<-NULL

# procalcitonin 22 unkown 
Cli_DF$procalcitonin<-NULL

# lactate 39 unkown 
Cli_DF$lactate<-NULL

# fibrinogen 30 unkown 
Cli_DF$fibrinogen<-NULL

# sofa 50 unkown 
Cli_DF$sofa<-NULL


table(Cli_DF$sofa)
colnames(Cli_DF)[17]

Cli_DF[,2]<-as.numeric(Cli_DF[,2])
Cli_DF[,5]<-as.numeric(Cli_DF[,5])
Cli_DF[,7]<-as.numeric(Cli_DF[,7])
Cli_DF[,9]<-as.numeric(Cli_DF[,9])

colnames(Cli_DF)[8]<-'Diabetes_mellitus'
colnames(Cli_DF)[9]<-'hospital_free_days'

saveRDS(Cli_DF,'Cli_DF_v2.rds')




