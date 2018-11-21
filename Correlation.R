# Wed Nov 21 13:41:25 2018 ------------------------------
#ref:https://ytlogos.github.io/2017/06/25/R%E8%AF%AD%E8%A8%80%E5%AD%A6%E4%B9%A0%E7%AC%94%E8%AE%B0%E4%B9%8B%E7%9B%B8%E5%85%B3%E6%80%A7%E7%9F%A9%E9%98%B5%E5%88%86%E6%9E%90%E5%8F%8A%E5%85%B6%E5%8F%AF%E8%A7%86%E5%8C%96/


rm(list=ls())
load("PCa_TCGA.Rdata")

data_demo <- merge[1:100,8888:8896]
colnames(data_demo)

library(Hmisc)
library(tibble)
library(reshape2)

res <- rcorr(as.matrix(data_demo),type = 'pearson')
res#总览
res$r#提取相关系数
res$P#提取p value

#combine r and p into a matrix
res_r <- as.data.frame(res$r)
res_r <- add_column(res_r,gene = rownames(res_r),.before = 1)
res_r <- melt(res_r,id.vars=1,measure.vars=2:10,
              variable.name="gene2",value.name="r")

res_p <- as.data.frame(res$P)
res_p <- add_column(res_p,gene = rownames(res_p),.before = 1)
res_p <- melt(res_p,id.vars=1,measure.vars=2:10,
              variable.name="gene2",value.name="p")

res2 <- cbind(res_r,res_p)
res2 <- res2[,-c(4,5)]
res2 <- res2[res2$gene!=res2$gene2,]
res2$gene2 <- as.character(res2$gene2)

#visualize result
#1.
symnum(res$r)
#2.
library(corrplot)
corrplot(res$r,type='full',order='hclust',tl.col = "black", 
         diag=T,outline = T,bg = 'white',addgrid.col = 'black',
         addCoef.col = 'black',tl.cex = 1,
         p.mat = res$P,sig.level = 0.01,insig = 'blank',
         number.digits = 2)
res$r
#3.可以直接输入矩阵，很方便！
library(PerformanceAnalytics)
chart.Correlation(data_demo,histogram = T,method = 'pearson')
#4.热图
library(pheatmap)
pheatmap(res$r,
         col = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = F,cluster_cols=F,
         display_numbers = T,number_color = 'black',fontsize_number = 10,
         border_color = NA,legend=TRUE,
         show_colnames=T,show_rownames=T,
         fontsize_col = 10, 
         fontsize_row = 10,fontsize = 10,
         main = 'Correlation')


