library(dplyr)
library(survminer)
library(survival)

surv <- read.table('41586_2019_1007_MOESM7_ESM.txt', sep = '\t', header = TRUE)
rownames(surv) <- surv$METABRIC.ID
surv_status <- surv %>% select(DeathBreast, Death, T, TLR, LR, TDR, DR)
surv_status <- as.data.frame(na.omit(surv_status))

surv_status$event <- (surv_status$Death & surv_status$DeathBreast) | surv_status$LR | surv_status$DR
surv_status$time <- apply(surv_status[, c('T', 'TLR', 'TDR')], 1, FUN = min)
surv_status <- surv_status %>% select(time, event)

latents_df <- as.data.frame(latents)
l_mb <- latents_df %>% select(Z1, Z6, Z10, Z59)
colnames(l_mb) <- c('LF1', 'LF6', 'LF10', 'LF59')

surv_status <- surv_status %>% filter(row.names(surv_status) %in% row.names(l_mb))
surv_data <- merge(l_mb, surv_status, by=0, all=TRUE)
rownames(surv_data) <- surv_data$Row.names

res.cut <- surv_cutpoint(surv_data, time = 'time', event = 'event',
                         variables = c('LF1', 'LF6', 'LF10', 'LF59'))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)

##### Factor 1 #####
#plot(res.cut, 'V1', palette = 'npg')
fit <- survfit(Surv(time, event) ~LF1, data = res.cat)
ggsurvplot(fit, data = res.cat, conf.int = TRUE, title = 'LF1', pval = T)#, risk.table = TRUE)

##### Factor 6 #####
#plot(res.cut, 'V6', palette = 'npg')
fit <- survfit(Surv(time, event) ~LF6, data = res.cat)
ggsurvplot(fit, data = res.cat, conf.int = TRUE, title = 'LF6', pval = T)#, risk.table = TRUE)

##### Factor 10 #####
#plot(res.cut, 'V10', palette = 'npg')
fit <- survfit(Surv(time, event) ~LF10, data = res.cat)
ggsurvplot(fit, data = res.cat, conf.int = TRUE, title = 'LF10', pval = T)#, risk.table = TRUE

##### Factor 59 #####
#plot(res.cut, 'V59', palette = 'npg')
fit <- survfit(Surv(time, event) ~LF59, data = res.cat)
ggsurvplot(fit, data = res.cat, conf.int = TRUE, title = 'LF59', pval = T)#, risk.table = TRUE

#####
fit2 <- coxph(Surv(time, event) ~ LF1 + LF6 + LF10 + LF59, data = surv_data)
ggadjustedcurves(fit2, data = surv_data, method = 'average')
ggforest(fit2, fontsize = 1.0, cpositions = c(0.02, 0.075, 0.25))

# Check whether follow cox regression model
cox.zph(fit2)

# Use cutoffs for multiple survival curve plots (based on combination of low/high of factors)
ggsurvplot(survfit(Surv(time, event) ~ 
                     (LF1>-0.0006932175)
                   + (LF6>0.0138811038)
                   + (LF10>-0.0092665728), 
                   surv_data),
           pval = T,
           title = 'Survival based on Factor Stratification Combination',
           legend.labs = c('All Low', 'LF1-, LF6-, LF10+', 'LF1-, LF6+, LF10-', 'LF1-, LF6+, LF10+', 'LF1+, LF6-, LF10-', 'LF1+, LF6-, LF10+', 'LF1+, LF6+, LF10-', 'All High'),
           legend = 'right')
