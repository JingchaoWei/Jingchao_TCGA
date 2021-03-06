# Thu Nov 08 15:13:40 2018 ------------------------------
rm(list = ls())
load('PCa_TCGA.Rdata')
marker <- c('ADRA1A','CXCL12','PPP1R1A')
data <- merge[,c('sampleID',marker)]
attach(data)
high_group <- subset(data,
                     subset = ADRA1A>median(ADRA1A) & CXCL12>median(CXCL12) &
                       PPP1R1A>median(PPP1R1A),select = 1)
high_group$group <- rep('high',139)
low_group <- subset(data,
                    subset = ADRA1A<=median(ADRA1A) & CXCL12<=median(CXCL12) &
                    PPP1R1A<=median(PPP1R1A),select = 1)
low_group$group <- rep('low',136)
group <- rbind(high_group,low_group)
detach(data)
write.table(high_group,'High_group_3_markers.txt',quote = F,sep = '\t',row.names = F,col.names = T)
write.table(low_group,'Low_group_3_markers.txt',quote = F,sep = '\t',row.names = F,col.names = T)
write.table(marker,'3_markers.txt',quote = F,sep = '\t',row.names = F,col.names = F)

#gene expression between high and low groups (t test)
MyTtest <- function(gene){
  data1 <- merge[match(group$sampleID,merge$sampleID),c('sampleID',gene)]
  data1 <- merge(data1,group,by='sampleID',all=F)
  group1 <- subset(data1,subset = data1$group=='high',select = gene)[[gene]]
  group2 <- subset(data1,subset = data1$group=='low',select = gene)[[gene]]
  #var equal test
  var(group1)
  var(group2)
  if (var(group1)==0&var(group2)==0) {
    print(gene)
  } else{
    var.test(group1,group2)
    tmp <- ifelse(var.test(group1,group2)$p.value>=0.05,"TRUE","FALSE")
    t.test(group1,group2,paired = F,var.equal = as.logical(tmp))
    pvalue_ttest <- t.test(group1,group2,paired = F)$p.value
    mean <- aggregate(data1[[gene]],by=list(data1$group),FUN=mean)
    mean$type <- rep('mean',2)
    sd <- aggregate(data1[[gene]],by=list(data1$group),FUN=sd)
    sd$type <- rep('sd',2)
    result <- rbind(mean,sd)
    result$gene <- rep(gene,4)
    result$p.value <- rep(pvalue_ttest,4)
    return(result)
  }
}

gene <- rownames(exp)

#检查基因名是否存在于数据库
colnames <- colnames(merge)
istrue <- vector()
for (i in gene) {
  findit <- isTRUE(i %in% colnames)
  istrue <- c(istrue,findit)
}
unique(istrue)

#利用循环，将每个基因的结果汇总
sum <- data.frame(Group.1=character(),
                      x=numeric(),
                      type=character(),
                      gene=character(),
                      p.value=numeric())
for (j in gene) {
  tmp <- MyTtest(j) 
  sum <- rbind(sum,tmp)
}
sum

sig <- sum[sum$p.value<0.05,]
sig.genes <- unique(sig$gene)
sig.genes


#figure out whether it's higher in Group.1 or Group.2 for each gene in sig.genes
MyComparison <- function(a){
  b <- sig[sig$gene==a,]
  value_groupHigh <- subset(b,subset = b$type=='mean'&b$Group.1=='high',select = 'x')[1,1]
  value_groupLow <- subset(b,subset = b$type=='mean'&b$Group.1=='low',select = 'x')[1,1]
  tmp <- as.numeric(value_groupHigh)-as.numeric(value_groupLow)
  comparison <- ifelse(tmp > 0,'higher in Group_High','higher in Group_Low')
  c <- data.frame(gene=a,append=comparison)
  return(c)
}

sig.genes_new <- data.frame(gene=character(),append=character(),stringsAsFactors = F)
for (a in sig.genes) {
  tmp <- MyComparison(a)
  sig.genes_new <- rbind(sig.genes_new,tmp)
}
sig.genes_new <- sig.genes_new[order(sig.genes_new$append),]
sig.genes_new

save(sig.genes_new,sig.genes,sum,group,file = 'Sig.genes.Rdata')


#other features between high and low groups
load('Sig.genes.Rdata')
load('PCa_TCGA.Rdata')
phe_group <- merge(group,phe,by= 'sampleID',all=F)
group_low <- phe_group[phe_group$group=='low',]
group_high <- phe_group[phe_group$group=='high',]

a <- table(group_high$gleason_score)
b <- table(group_low$gleason_score)
c <- rbind(a,b)
c
par(mfrow=c(2,1))
barplot(a)
barplot(b)
chisq.test(c)#if got any warning, better to use fishertest
fisher.test(c)
p <- fisher.test(c)$p.value
format(p,scientific = F)# use non-scientific format, easier to read


a <- table(group_high$biochemical_recurrence)
b <- table(group_low$biochemical_recurrence)
c <- rbind(a,b)
c
par(mfrow=c(2,1))
barplot(a)
barplot(b)
chisq.test(c)#if got any warning, better to use fishertest
fisher.test(c)
p <- fisher.test(c)$p.value
format(p,scientific = F)# use non-scientific format, easier to read

a <- table(group_high$followup_treatment_success)
b <- table(group_low$followup_treatment_success)
c <- rbind(a,b)
c
par(mfrow=c(2,1))
barplot(a)
barplot(b)
chisq.test(c)#if got any warning, better to use fishertest
fisher.test(c)
p <- fisher.test(c)$p.value
format(p,scientific = F)# use non-scientific format, easier to read






