setwd('/Users/mjia/Box/EBMF/')
# load metabric sample info
Sample_info<-readRDS('./Sample_info.rds')


setwd('/Users/mjia/Desktop/SSNPA/ssNPA_test/')
# Clinical information for the full dataset.
recurrent_info_1<-read.table('./clinical_feature/41586_2019_1007_MOESM7_ESM.txt',skip = 0,sep="\t", quote="", header=T,stringsAsFactors=F)

# match sample ID
sample_index<-match(Sample_info$SAMPLE_ID,recurrent_info_1$METABRIC.ID)
table(Sample_info$SAMPLE_ID==recurrent_info_1$METABRIC.ID[sample_index])
recurrent_info_1<-recurrent_info_1[sample_index,]


recurrent_info_1$Distant_relapse <-'Survivor'
recurrent_info_1$Distant_relapse[recurrent_info_1$DR ==1 & recurrent_info_1$TDR >= 365*5] <- 'Late'
recurrent_info_1$Distant_relapse[recurrent_info_1$DR ==1 & recurrent_info_1$TDR <= 365*2] <- 'Early'
recurrent_info_1$Distant_relapse[recurrent_info_1$DR ==1 & recurrent_info_1$TDR > 365*2 & recurrent_info_1$TDR < 365*5] <-'Mild'


recurrent_info_1$Loco_regional_relapse <-'Survivor'
recurrent_info_1$Loco_regional_relapse[recurrent_info_1$LR ==1 & recurrent_info_1$TLR >= 365*5] <- 'Late'
recurrent_info_1$Loco_regional_relapse[recurrent_info_1$LR ==1 & recurrent_info_1$TLR <= 365*2] <- 'Early'
recurrent_info_1$Loco_regional_relapse[recurrent_info_1$LR ==1 & recurrent_info_1$TLR > 365*2 & recurrent_info_1$TLR < 365*5] <-'Mild'

saveRDS(recurrent_info_1,'recurrent_info.rds')
