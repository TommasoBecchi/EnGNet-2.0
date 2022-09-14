library("GEOquery")
info=getGEO(filename="GSE157194_series_matrix.txt")
gene_expression=read.delim(file="Raw_gene_counts_matrix.txt",row.names = 1)
gene_expression<-gene_expression[,1:111]

library(dplyr)
sampleinfo<-pData(info)
sampleinfo<-sampleinfo[1:111,]
sampleinfo2<-select(sampleinfo, source_name_ch1 )
rownames(sampleinfo2)<-colnames(gene_expression)
sampleinfo2$source_name_ch1[which(sampleinfo2$source_name_ch1 == "m0_AL") ] <- "Lesional"
sampleinfo2$source_name_ch1[which(sampleinfo2$source_name_ch1 == "m0_AN") ] <- "Non_Lesional"

BiocManager::install("pheatmap")
library("pheatmap")
corMatrix<-cor(gene_expression, use="c")
pheatmap(corMatrix)
pheatmap(corMatrix,annotation_col=sampleinfo2,show_rownames=FALSE,show_colnames=FALSE,annotation_legend=TRUE)

library("ggplot2")
pca<-prcomp(t(gene_expression))
colnames(sampleinfo2)<-c("group")
cbind(sampleinfo2,pca$x) %>%
ggplot(aes(x = PC1, y=PC2, col=group)) + geom_point()

newc<-factor(sampleinfo2$group)
sampleinfo2$group<-newc
levels(sampleinfo2$group)
colors<-c("red","blue")[sampleinfo2$group]
plotMDS(gene_expression,col=colors,pch=19) 
legend("topleft",fill=c("red","blue"),legend=levels(sampleinfo2$group),cex=0.8)
##############
##############

x<-DGEList(gene_expression)
dim(x)

library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
biomart <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position', 'percentage_gene_gc_content' ,'description'), filters = 'ensembl_gene_id', values =rownames(gene_expression), mart = hsmart)


#REMOVE LOW EXPRESSED DATA
cpm<-cpm(x)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
plot(density(lcpm[,1]), lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
keep.exprs <- filterByExpr(x, group=sampleinfo2$group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")

######
par(mfrow=c(1,2))
lcpm <- cpm(x, log=TRUE)
boxplot(lcpm[,1:50], las=2, main="")
x2 <- calcNormFactors(x, method="TMM")  
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm[,1:50], las=2, main="")
#######
#######
design <- model.matrix(~ 0 + sampleinfo2$group)
colnames(design) <- c("Lesional","Non_Lesional")
cont.matrix <- makeContrasts(LesVsNonLes= Lesional - Non_Lesional,levels=design)
par(mfrow=c(1,2))
v <- voom(x2, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
#####
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
results<-data.frame(dt)
results$id<-rownames(results)
DEGs<-results[results$LesVsNonLes!=0,]

library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
biomart <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position', 'percentage_gene_gc_content' ,'description'), filters = 'ensembl_gene_id', values =rownames(DEGs), mart = hsmart)
######
les.vs.nles <- topTreat(tfit, coef=1, n=Inf)
head(les.vs.nles)
par(mfrow=c(1,1))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-4,13))


les.vs.nles %>% 
  mutate(Significant = adj.P.Val < 0.05, abs(logFC) > 0.5) %>% 
  ggplot(aes(x = logFC, y = -adj.P.Val, col=Significant)) + geom_point()+geom_hline(aes(yintercept=-0.05), colour="#990000", linetype="dashed")+ geom_vline(aes(xintercept=1), colour="#990000", linetype="dashed")+geom_vline(aes(xintercept=-1), colour="#990000", linetype="dashed")

les.vs.nles%>% 
  mutate(Significant = adj.P.Val < 0.05, abs(logFC) > 1) %>% 
  ggplot(aes(x = logFC, y = -adj.P.Val, col=Significant)) + geom_point()+geom_hline(aes(yintercept=-0.05), colour="#990000", linetype="dashed")+ geom_vline(aes(xintercept=1), colour="#990000", linetype="dashed")+geom_vline(aes(xintercept=-1), colour="#990000", linetype="dashed")+ylim(-0.2,0)


###############################
###############################
dim(DEGs)
dim(x$counts)
gene_final<-data.frame(x$counts)
count_AL<-gene_final[rownames(gene_final)%in% DEGs$id,sampleinfo2$group=="Lesional"]
count_NL<-gene_final[rownames(gene_final)%in% DEGs$id,sampleinfo2$group=="Non_Lesional"]
count_AL<-count_AL[,1:54]
####################
####################
cpmAL<-cpm(count_AL)
plot(cpmAL[,1],count_AL[,1],ylim=c(0,50),xlim=c(0,1))
abline(h=10,v=0.4)
thresh<-cpmAL>0.5
head(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >= 54
summary(keep)
countAL_final<-count_AL[keep==TRUE,]
##################
##################
library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
biomart2 <- getBM(attributes=c('ensembl_gene_id','start_position','end_position', 'percentage_gene_gc_content'), filters = 'ensembl_gene_id', values =rownames(count_AL), mart = hsmart)   
biomart2$length<-biomart2$end_position-biomart2$start_position
gene_length_gc<-biomart2[,c(1,4,5)]
dim(gene_length_gc)

firstids<-rownames(gene_final)
both <-firstids %in% gene_length_gc$ensembl_gene_id
summary(both)
gene_final$present<-both
gene_final_AL<-gene_final[gene_final$present==TRUE,sampleinfo2$group=="Lesional"]
gene_final_AL$present<-NULL
dim(gene_final_AL)
gene_final_NL<-gene_final[gene_final$present==TRUE,sampleinfo2$group=="Lesional"]
gene_final_NL$present<-NULL
dim(gene_final_AL)

library("cqn")
cqn.subset<-cqn(count_AL,lengths=gene_length_gc$length,x=gene_length_gc$percentage_gene_gc_content)
normalized_count_AL<-data.frame(cqn.subset$y+cqn.subset$offset)
dim(normalized_count_AL)

par(mfrow=c(1,2))
lcpm2 <- cpm(count_AL, log=TRUE)
boxplot(lcpm2[,20:40], las=2, main="")
boxplot(normalized_count_AL[,20:40])

normalized_count_AL<-normalized_count[,sampleinfo2$group=="Lesional"]
normalized_count_NL<-normalized_count[,sampleinfo2$group=="Non_Lesional"]

library("psych")
correlation_AL<-corr.test(t(normalized_count_AL),method="kendall",adjust="BH",ci=FALSE)
cormat<-correlation_AL$r
head(cormat,n=c(6,6))
pmat<-correlation_AL$p
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
cormat_flat<-flattenCorrMatrix(cormat,pmat)
significant<-cormat_flat$p<0.05
summary(significant)
final_data<-cormat_flat[significant==TRUE,]
write.csv(final_data,file="data_deg.csv")

library("igraph")
BiocManager::install("rgl")
library("rgl")
g<-graph_from_data_frame(final_data)

par(mfrow=c(1,2))
thr<-c()
score<-c()
for (i in seq (0,1,0.01)){
  data<-final_data[abs(final_data$cor)>=i,]
  g<-graph_from_data_frame(data,directed=F)
  dist<-degree.distribution(g)
  power_law<-fit_power_law(dist)
  thr<-c(thr,i)
  score<-c(score,power_law$KS.stat)
}
result=data.frame(threshold=thr,score=score)
best<-min(score)
result[result$score==best,]
plot(result,type="l",main="Lesional")

thr<-c()
score<-c()
for (i in seq (0,1,0.01)){
  data<-final_data_NL[abs(final_data_NL$cor)>=i,]
  g<-graph_from_data_frame(data,directed=F)
  dist<-degree.distribution(g)
  power_law<-fit_power_law(dist)
  thr<-c(thr,i)
  score<-c(score,power_law$KS.stat)
}
result=data.frame(threshold=thr,score=score)
best<-min(score)
result[result$score==best,]
plot(result,type="l",main="Non_Lesional")

score

BiocManager::install("pROC")
