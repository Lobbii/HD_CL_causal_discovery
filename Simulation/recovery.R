source('R/read_causality_graph.R')
source('R/adjacency.R')
source('R/arrowhead.R')
source('R/evaluate_correctness.R')

setwd('')
simulation_5<-readRDS('./simulation_5.rds')
net_true<-read_causality_graph('./Graphs/Graph_0.txt')
graph_true<-Split_graph(graph =net_true ,con_var = simulation_5$contunious_var)

Evaluation<-matrix(NA,nrow = 3,ncol = 16)
colnames(Evaluation)<-apply(expand.grid(c('AP',"AR","AHP","AHR"),c("_Overall","_CC","_CD","_DD")), 1, paste, collapse="")
rownames(Evaluation)<-c("Simulated","EBMF","EMBF_BF")
net_estimated<-read_causality_graph('./simulation_5/orig/ori_mgmpcmax.txt')
Evaluation[1,1:4]<-evaluate_correctness(net_estimated,net_true)
graph_esti<-Split_graph(graph =net_estimated ,con_var = simulation_5$contunious_var)
Evaluation[1,5:16]<-eva_adj(graph_esti,graph_true)

net_estimated<-read_causality_graph('./simulation_5/f/EBMF_mgmpcmax.txt')
Evaluation[2,1:4]<-evaluate_correctness(net_estimated,net_true)
graph_esti<-Split_graph(graph =net_estimated ,con_var = simulation_5$contunious_var)
Evaluation[2,5:16]<-eva_adj(graph_esti,graph_true)


net_estimated<-read_causality_graph('./simulation_5/f_b/EBMF_mgmpcmax.txt')
Evaluation[3,1:4]<-evaluate_correctness(net_estimated,net_true)
graph_esti<-Split_graph(graph =net_estimated ,con_var = simulation_5$contunious_var)
Evaluation[3,5:16]<-eva_adj(graph_esti,graph_true)


write.csv(Evaluation,'simulation_5.csv')
