load("PCa_TCGA.Rdata")
up_gene <- read.delim('Query_Genes.txt',stringsAsFactors = F)
up_gene <- up_gene$Up.[1:9]

library(survival)
library(survminer)
library(ggplot2)
library(gplots)
library(gridExtra)

#分四组
  a <- quantile(merge$NT5E)
  merge$append <- ifelse(merge$NT5E<a[2],'NT5E:lowest',
                         ifelse(merge$NT5E<a[3],'NT5E:lower',
                                ifelse(merge$NT5E<a[4],'NT5E:higher','NT5E:higheset')))
  su <- Surv(merge$OS.time,merge$OS)
  sfit <- survfit(su~append, data=merge)
  sfit
  summary(sfit)
  su2 <- Surv(merge$RFS.time,merge$RFS)
  sfit_2 <- survfit(su2~append, data=merge)
  sfit_2
  summary(sfit_2)
  oa <- ggsurvplot(sfit, conf.int=F, pval=TRUE,data = merge,
             title='Overall Survival',combine = T)
  rf <- ggsurvplot(sfit_2, conf.int=F, pval=TRUE,data = merge,
                   title="Recurrence/Disease Free Survival",combine = T)
  arrange_ggsurvplots(x = list(oa,rf),title = 'gene',nrow = 1,ncol = 2)
  print(paste0('Is ',gene,' in the dataset?  ',a))
  survival_data <- merge[,c('OS.time','OS','RFS.time','RFS',gene)]
  write.table(survival_data,file = paste0(gene,"_Survival_data.txt"),
              quote = F,sep = "\t",row.names = F)
