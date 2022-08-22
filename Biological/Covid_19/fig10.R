icu.list <- list()

icu.list[[1]] <- roc(ICU ~ prob_icu, data = results,
               print.auc=TRUE,
               print.auc.y=0.2,
               print.auc.x=0.2,
               col = '#1c61b6')
icu.list[[2]] <- roc(ICU ~ prob_icu[,icu.best], data = results2,
                print.auc=TRUE,
                print.auc.y=0.15,
                print.auc.x=0.2,
                add = TRUE,
                col = '#008600')
icu.list[[3]] <- roc(ICU ~ prob_icu[,icu.best3], data = results3,
                 print.auc=TRUE,
                 print.auc.y=0.1,
                 print.auc.x=0.2,
                 add = TRUE,
                 col = '#EE82EE')
icu.list[[4]] <- roc(ICU ~ prob_icu, data = results4,
                  print.auc=TRUE,
                  print.auc.y=0.05,
                  print.auc.x=0.2,
                  add = TRUE,
                  col = '#FF0000')
icu.list[[5]] <- roc(ICU ~ prob_icu, data = results5,
                   print.auc=TRUE,
                   print.auc.y=0,
                   print.auc.x=0.2,
                   add = TRUE,
                   col = '#00FFFF')

g1 <- ggroc(icu.list) +
  labs(title = '5-fold CV ICU Prediction AUC') +
  scale_color_hue(name = 'Method',
    labels = c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'))

#+ legend('bottom', legend=c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'),
       #col=c('#1c61b6', '#008600', '#EE82EE', '#FF0000', '#00FFFF'), lwd=2)

cov.list <- list()

cov.list[[1]] <- roc(disease_state ~ prob_cov, data = results,
               print.auc=TRUE,
               print.auc.y=0.2,
               print.auc.x=0.2,
               plot = TRUE,
               main = "5-fold CV Disease State Prediction AUC",
               col = '#1c61b6')
cov.list[[2]] <- roc(disease_state ~ prob_cov[,cov.best], data = results2,
                print.auc=TRUE,
                print.auc.y=0.15,
                print.auc.x=0.2,
                plot = TRUE,
                add = TRUE,
                col = '#008600')
cov.list[[3]] <- roc(disease_state ~ prob_cov[,cov.best3], data = results3,
                 print.auc=TRUE,
                 print.auc.y=0.1,
                 print.auc.x=0.2,
                 plot = TRUE,
                 add = TRUE,
                 col = '#EE82EE')
cov.list[[4]] <- roc(disease_state ~ prob_cov, data = results4,
                  print.auc=TRUE,
                  print.auc.y=0.05,
                  print.auc.x=0.2,
                  plot = TRUE,
                  add = TRUE,
                  col = '#FF0000')
cov.list[[5]] <- roc(disease_state ~ prob_cov, data = results5,
                   print.auc=TRUE,
                   print.auc.y=0,
                   print.auc.x=0.2,
                   plot = TRUE,
                   add = TRUE,
                   col = '#00FFFF')
  #legend('bottom', legend=c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'),
  #       col=c('#1c61b6', '#008600', '#EE82EE', '#FF0000', '#00FFFF'), lwd=2)

g2 <- ggroc(cov.list) + 
  labs(title = '5-fold CV Disease state Prediction AUC') +
  scale_color_hue(name = 'Method',
                  labels = c('FCI-Max Markov Blanket', 'All EBMF Factors', 'Gene Expression', 'FGES Markov Blanket', 'PC-Max Markov Blanket'))

library(tidyverse)
data.auc1 <- icu.list %>%
  map(~tibble(AUC = .x$auc)) %>%
  bind_rows(.id = "name")
data.auc2 <- cov.list %>%
  map(~tibble(AUC = .x$auc)) %>%
  bind_rows(.id = "name")

prow <- plot_grid(g1 + theme(legend.position="none"), g2 + theme(legend.position="none"), labels = 'AUTO')

legend <- get_legend(
  g1 #+ theme(legend.box.margin = margin(0, 0, 0, 50))
)

pall <- plot_grid(prow, legend, rel_widths = c(2, 0.6))

