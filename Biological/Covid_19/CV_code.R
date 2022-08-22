library(rCausalMGM)
library(caret)
library(pROC)
library(plyr)
library(glmnet)
library(edgeR)
library(flashr)
library(permute)
library(Matrix)
library(cowplot)
library(dplyr)

set.seed(42)

count_mat<-read.delim('GSE157103_genes.ec.tsv',row.names = 1)

d0 <- DGEList(count_mat)
d0 <- calcNormFactors(d0)
d0

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(count_mat) # Sample names
snames

cultivar <- substr(snames, 1,  1) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar

group <- interaction(cultivar, time)
group

#plotMDS(d, col = as.numeric(factor(cultivar)))

norm_d <- voom(d, plot = T)
norm_d<-t(norm_d$E)

Y<-as.data.frame(norm_d)
Cli_DF<-readRDS('Cli_DF_v2.rds')

index<-rownames(Y) %in% rownames(Cli_DF)
Y<-Y[index,]

randomize <- shuffle(nrow(Cli_DF))
Y <- Y[randomize,]
Cli_DF <- Cli_DF[randomize,]
clinical_features <- colnames(Cli_DF)

flds <- createFolds(Cli_DF$disease_state, k = 5, list = TRUE, returnTrain = TRUE)

## Results for markov blanket fci
results <- Cli_DF[,c('disease_state', 'ICU')]
results$prob_cov <- 0
results$prob_icu <- 0

## Results for all factors
results2 <- Cli_DF[,c('disease_state', 'ICU')]
results2$prob_cov <- matrix(0, nrow=123, ncol=100)
results2$prob_icu <- matrix(0, nrow=123, ncol=100)

## Results for expression data
results3 <- Cli_DF[,c('disease_state', 'ICU')]
results3$prob_cov <- matrix(0, nrow=123, ncol=100)
results3$prob_icu <- matrix(0, nrow=123, ncol=100)

## Results for markov blanket fges
results4 <- Cli_DF[,c('disease_state', 'ICU')]
results4$prob_cov <- 0
results4$prob_icu <- 0

## Results for markov blanket PC
results5 <- Cli_DF[,c('disease_state', 'ICU')]
results5$prob_cov <- 0
results5$prob_icu <- 0

counter <- 1
for (k in flds) {
  train_expr <- as.matrix(Y[k,])
  train_clin <- as.matrix(Cli_DF[k,])
  test_expr <- as.matrix(Y[-k,])
  test_clin <- as.matrix(Cli_DF[-k,])
  
  #f_ <- readRDS(paste0('k',counter,'.f.rds')) #flash(train_expr,Kmax = 1000)
  #saveRDS(f_,paste0('k',counter,'.f.rds'))
  f.b <- readRDS(paste0('k',counter,'.f.b.rds')) #flash(train_expr, f_init = f_, backfit=TRUE, greedy=FALSE,Kmax = 500)
  #saveRDS(f.b,paste0('k',counter,'.f.b.rds'))
  
  L_train <- f.b$ldf$l
  diag <- Diagonal(x = f.b$ldf$d)
  fact <- t(f.b$ldf$f)
  L_test <- as.matrix(test_expr %*% t(fact) %*% solve((fact %*% t(fact))) %*% solve(diag))
  
  colnames(L_train)<-paste0('LF_',1:ncol(L_train))
  colnames(L_test)<-paste0('LF_',1:ncol(L_test))
  #model_data<-cbind(L, Cli_DF[,c('disease_state', 'ICU')], stringsAsFactors=T)
  train_data<-cbind.data.frame(L_train, train_clin)#, stringsAsFactors=FALSE)
  train_data<-as.data.frame(train_data)
  test_data<-cbind.data.frame(L_test, test_clin)#, stringsAsFactors=FALSE)
  test_data<-as.data.frame(test_data)
  
  train_data$ICU <- as.factor(train_data$ICU)
  train_data$disease_state <- as.factor(train_data$disease_state)
  train_data$mechanical_ventilation <- as.factor(train_data$mechanical_ventilation)
  train_data$Sex <- as.factor(train_data$Sex)

  test_data$ICU <- as.factor(test_data$ICU)
  test_data$disease_state <- as.factor(test_data$disease_state)
  test_data$mechanical_ventilation <- as.factor(test_data$mechanical_ventilation)
  test_data$Sex <- as.factor(test_data$Sex)
  
  ##  FCI Markov Blanket Results
  g <- fciMax(train_data, maxDiscrete = 2, alpha = 0.05, fdr = FALSE, verbose = TRUE)
  icu_mb <- g$markov.blankets$ICU[!(g$markov.blankets$ICU %in% clinical_features)]
  cov_mb <- g$markov.blankets$disease_state[!(g$markov.blankets$disease_state %in% clinical_features)]
  
  icu_data <- train_data[,icu_mb, drop=FALSE]
  icu_data$ICU <- as.factor(train_data$ICU)
  cov_data <- train_data[,cov_mb, drop=FALSE]
  cov_data$disease_state <- as.factor(train_data$disease_state)
  
  icu_model <- glm(ICU~., data = icu_data, family = 'binomial')
  cov_model <- glm(disease_state~., data = cov_data, family = 'binomial')
  
  prob_icu <- predict(icu_model,
                      newdata = test_data[,icu_mb, drop=FALSE],
                      type='response')
  prob_cov <- predict(cov_model,
                      newdata = test_data[,cov_mb, drop=FALSE],
                      type='response')
  
  results$prob_cov[-k] <- prob_cov
  results$prob_icu[-k] <- prob_icu
  
  #### FGES markov blanket
  write.table(train_data, paste0('data/K', counter, '_train.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
  data_file <- paste0('data/K', counter, '_train.txt')
  out_file <- paste0('data/K', counter, '_fges.txt')
  cmd <- paste0('java -jar tetrad.jar -alg FGES -penalty 1 -maxCat 2 -d ', data_file, ' -o ', out_file)
  system(cmd)
  g <- loadGraph(out_file)
  icu_mb <- g$markov.blankets$ICU[!(g$markov.blankets$ICU %in% clinical_features)]
  cov_mb <- g$markov.blankets$disease_state[!(g$markov.blankets$disease_state %in% clinical_features)]
  
  icu_data <- train_data[,icu_mb, drop=FALSE]
  icu_data$ICU <- as.factor(train_data$ICU)
  cov_data <- train_data[,cov_mb, drop=FALSE]
  cov_data$disease_state <- as.factor(train_data$disease_state)
  
  icu_model <- glm(ICU~., data = icu_data, family = 'binomial')
  cov_model <- glm(disease_state~., data = cov_data, family = 'binomial')
  
  prob_icu <- predict(icu_model,
                      newdata = test_data[,icu_mb, drop=FALSE],
                      type='response')
  prob_cov <- predict(cov_model,
                      newdata = test_data[,cov_mb, drop=FALSE],
                      type='response')
  
  results4$prob_cov[-k] <- prob_cov
  results4$prob_icu[-k] <- prob_icu
  
  ##### pcmax Markov Blanket
  g <- pcMax(train_data, maxDiscrete = 2, alpha = 0.05, fdr = FALSE, verbose = TRUE)
  icu_mb <- g$markov.blankets$ICU[!(g$markov.blankets$ICU %in% clinical_features)]
  cov_mb <- g$markov.blankets$disease_state[!(g$markov.blankets$disease_state %in% clinical_features)]
  
  icu_data <- train_data[,icu_mb, drop=FALSE]
  icu_data$ICU <- as.factor(train_data$ICU)
  cov_data <- train_data[,cov_mb, drop=FALSE]
  cov_data$disease_state <- as.factor(train_data$disease_state)
  
  icu_model <- glm(ICU~., data = icu_data, family = 'binomial')
  cov_model <- glm(disease_state~., data = cov_data, family = 'binomial')
  
  prob_icu <- predict(icu_model,
                      newdata = test_data[,icu_mb, drop=FALSE],
                      type='response')
  prob_cov <- predict(cov_model,
                      newdata = test_data[,cov_mb, drop=FALSE],
                      type='response')
  
  results5$prob_cov[-k] <- prob_cov
  results5$prob_icu[-k] <- prob_icu
  
  ## All Factors Results
  icu_data <- L_train
  cov_data <- L_train
  
  cov_full_model <- cv.glmnet(y = as.matrix(as.factor(train_data$disease_state)),
                              x = L_train,
                              family = 'binomial',
                              alpha = 0.5)
  
  icu_full_model <- cv.glmnet(y = as.matrix(as.factor(train_data$ICU)),
                              x = L_train,
                              family = 'binomial',
                              alpha = 0.5)

  l.idx <- 1
  for (l in cov_full_model$lambda) {
    cov_model <- glmnet(x = as.matrix(cov_data),
                        y = train_data$disease_state,
                        family = 'binomial',
                        alpha = 0.5,
                        lambda = l)
    prob_cov <- predict(cov_model, 
                        newx = as.matrix(L_test))#,
    #type='response',
    #s = "lambda.1se")
    results2$prob_cov[-k,l.idx] <- prob_cov
    l.idx <- l.idx + 1
  }
  cov.best <- which.max(apply(results2$prob_cov, 2, function(x) auc(Cli_DF$disease_state, x)))
  
  l.idx <- 1
  for (l in icu_full_model$lambda) {
    icu_model <- glmnet(x = as.matrix(icu_data),
                        y = train_data$ICU,
                        family = 'binomial',
                        alpha = 0.5,
                        lambda = l)
    prob_icu <- predict(icu_model, 
                        newx = as.matrix(L_test))#,
    #type='response',
    #s = "lambda.1se")
    results2$prob_icu[-k,l.idx] <- prob_icu
    l.idx <- l.idx + 1
  }
  icu.best <- which.max(apply(results2$prob_icu, 2, function(x) auc(Cli_DF$ICU, x)))
  
  ## All genes results
  icu_data <- train_expr
  cov_data <- train_expr
  
  cov_full_model <- cv.glmnet(y = as.matrix(as.factor(train_data$disease_state)),
                              x = L_train,
                              family = 'binomial',
                              alpha = 0.5)
  
  icu_full_model <- cv.glmnet(y = as.matrix(as.factor(train_data$ICU)),
                              x = L_train,
                              family = 'binomial',
                              alpha = 0.5)
  
  l.idx <- 1
  for (l in cov_full_model$lambda) {
    cov_model <- glmnet(x = as.matrix(cov_data),
                        y = train_data$disease_state,
                        family = 'binomial',
                        alpha = 0.5,
                        lambda = l)
    prob_cov <- predict(cov_model, 
                        newx = as.matrix(test_expr))#,
    #type='response',
    #s = "lambda.1se")
    results3$prob_cov[-k,l.idx] <- prob_cov
    l.idx <- l.idx + 1
  }
  cov.best3 <- which.max(apply(results3$prob_cov, 2, function(x) auc(Cli_DF$disease_state, x)))
  
  l.idx <- 1
  for (l in icu_full_model$lambda) {
    icu_model <- glmnet(x = as.matrix(icu_data),
                        y = train_data$ICU,
                        family = 'binomial',
                        alpha = 0.5,
                        lambda = l)
    prob_icu <- predict(icu_model, 
                        newx = as.matrix(test_expr))#,
    #type='response',
    #s = "lambda.1se")
    results3$prob_icu[-k,l.idx] <- prob_icu
    l.idx <- l.idx + 1
  }
  icu.best3 <- which.max(apply(results3$prob_icu, 2, function(x) auc(Cli_DF$ICU, x)))
  
  counter <- counter + 1
}

save.image(file = 'covid.RData')
#savehistory(file = 'covid.RHistory')

## Plots for markov blanket
g1 <- plot.roc(ICU ~ prob_icu, data = results,
                print.auc=TRUE,
                ci=TRUE,
                thresholds='best',
                print.thres='best',
                main = "ICU Prediction AUC with FCI Markov Blanket")
g2 <- plot.roc(disease_state ~ prob_cov, data = results,
                print.auc=TRUE,
                ci=TRUE,
                thresholds='best',
                print.thres='best',
                main = "Disease State Prediction AUC with FCI Markov Blanket")

## Plots for all factors
g11 <- plot.roc(ICU ~ prob_icu[,icu.best], data = results2,
                 print.auc=TRUE,
                 ci=TRUE,
                 thresholds='best',
                 print.thres='best',
                 main = "ICU Prediction AUC with All Factors")
g22 <- plot.roc(disease_state ~ prob_cov[,cov.best], data = results2,
                 print.auc=TRUE,
                 ci=TRUE,
                 thresholds='best',
                 print.thres='best',
                 main = "Disease State Prediction AUC with All Factors")

## Plots for all genes
g111 <- plot.roc(ICU ~ prob_icu[,icu.best3], data = results3,
               print.auc=TRUE,
               ci=TRUE,
               thresholds='best',
               print.thres='best',
               main = "ICU Prediction AUC with Gene Expression")
g222 <- plot.roc(disease_state ~ prob_cov[,cov.best3], data = results3,
               print.auc=TRUE,
               ci=TRUE,
               thresholds='best',
               print.thres='best',
               main = "Disease State Prediction AUC with Gene Expression")

#### plots for mb fges
g1111 <- plot.roc(ICU ~ prob_icu, data = results4,
               print.auc=TRUE,
               ci=TRUE,
               thresholds='best',
               print.thres='best',
               main = "ICU Prediction AUC with FGES Markov Blanket")
g2222 <- plot.roc(disease_state ~ prob_cov, data = results4,
               print.auc=TRUE,
               ci=TRUE,
               thresholds='best',
               print.thres='best',
               main = "Disease State Prediction AUC with FGES Markov Blanket")

#### plots for mb pcmax
g11111 <- plot.roc(ICU ~ prob_icu, data = results5,
                  print.auc=TRUE,
                  ci=TRUE,
                  thresholds='best',
                  print.thres='best',
                  main = "ICU Prediction AUC with PC-Max Markov Blanket")
g22222 <- plot.roc(disease_state ~ prob_cov, data = results5,
                  print.auc=TRUE,
                  ci=TRUE,
                  thresholds='best',
                  print.thres='best',
                  main = "Disease State Prediction AUC with PC-Max Markov Blanket")

## Combined ROC curves
jpeg('fig_10.jpg', width = 10, height = 5)

l1 <- plot.roc(ICU ~ prob_icu, data = results,
               print.auc=TRUE,
               print.auc.y=0.2,
               print.auc.x=0.2,
               main = "5-fold CV ICU Prediction AUC",
               col = '#1c61b6')
l11 <- plot.roc(ICU ~ prob_icu[,icu.best], data = results2,
                 print.auc=TRUE,
                 print.auc.y=0.15,
                 print.auc.x=0.2,
                 add = TRUE,
                 col = '#008600')
l111 <- plot.roc(ICU ~ prob_icu[,icu.best3], data = results3,
                  print.auc=TRUE,
                  print.auc.y=0.1,
                  print.auc.x=0.2,
                  add = TRUE,
                  col = '#EE82EE')
l1111 <- plot.roc(ICU ~ prob_icu, data = results4,
                   print.auc=TRUE,
                   print.auc.y=0.05,
                   print.auc.x=0.2,
                   add = TRUE,
                   col = '#FF0000')
l11111 <- plot.roc(ICU ~ prob_icu, data = results5,
                    print.auc=TRUE,
                    print.auc.y=0,
                    print.auc.x=0.2,
                    add = TRUE,
                    col = '#00FFFF')
legend('bottom', legend=c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'),
       col=c('#1c61b6', '#008600', '#EE82EE', '#FF0000', '#00FFFF'), lwd=2)

#dev.off()

#pdf('ROC_Cov.pdf')

l2 <- plot.roc(disease_state ~ prob_cov, data = results,
               print.auc=TRUE,
               print.auc.y=0.2,
               print.auc.x=0.2,
               main = "5-fold CV Disease State Prediction AUC",
               col = '#1c61b6')
l22 <- plot.roc(disease_state ~ prob_cov[,cov.best], data = results2,
                 print.auc=TRUE,
                 print.auc.y=0.15,
                 print.auc.x=0.2,
                 add = TRUE,
                 col = '#008600')
l222 <- plot.roc(disease_state ~ prob_cov[,cov.best3], data = results3,
                  print.auc=TRUE,
                  print.auc.y=0.1,
                  print.auc.x=0.2,
                  add = TRUE,
                  col = '#EE82EE')
l2222 <- plot.roc(disease_state ~ prob_cov, data = results4,
                   print.auc=TRUE,
                   print.auc.y=0.05,
                   print.auc.x=0.2,
                   add = TRUE,
                   col = '#FF0000')
l22222 <- plot.roc(disease_state ~ prob_cov, data = results5,
                    print.auc=TRUE,
                    print.auc.y=0,
                    print.auc.x=0.2,
                    add = TRUE,
                    col = '#00FFFF') +
legend('bottom', legend=c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'),
       col=c('#1c61b6', '#008600', '#EE82EE', '#FF0000', '#00FFFF'), lwd=2)

dev.off()
