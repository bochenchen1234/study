library(R.utils)
library(Matrix) 
library(Seurat)
library(data.table)
library(stringr)
library(tidyverse)
library(GEOquery)
library(org.Hs.eg.db)
library(data.table)
library(clusterProfiler)
library(affy)
library(hgu133plus2.db)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rio)
library(limma)
library(future)
library(future.apply)
library(mlr3verse)
library(dplyr)
library(tidyr)
library(glmnet)
library(caret)
library(broom)
library(car) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(ggrepel)
library(sva)
library(parallel)
library(monocle)
library(CellChat)
library(oligo)
library(tinyarray)
library(FactoMineR)
library(factoextra)
library(ggridges)
library(forcats) 
library(ImmuCellAI)
library(ImmuCellAImouse)
library(scCustomize)
library(DGEobj.utils)
library(GenomicFeatures)
library(txdbmaker)
library(patchwork)
library(rms)
library(survival)

gse_GSE14323 <- getGEO('GSE14323', destdir = ".", AnnotGPL = FALSE ,getGPL = F)
pdata1 = gse_GSE14323$`GSE14323-GPL571_series_matrix.txt.gz`@phenoData@data
pdata2 = gse_GSE14323$`GSE14323-GPL96_series_matrix.txt.gz`@phenoData@data
pdata <- filter(pdata1, 
                (str_detect(source_name_ch1, 'cirrhosis') | 
                   str_detect(source_name_ch1, 'Normal')) &
                  !str_detect(source_name_ch1, 'HCC'))
expression_data <- exprs(gse_GSE14323[[1]])
expression_data<-as.data.frame(expression_data)
expression_data<-expression_data[,rownames(pdata)]
#save(pdata,file='GSE14323_pdata.csv')
#save(expression_data,file='GSE14323_expression_data.csv')
soft=getGEO(filename ="GSE14323_family.soft/GSE14323_family.soft")
gpl=soft@gpls[["GPL571"]]@dataTable@table
gpl<-gpl[,c(1,12)]
test_function <- apply(gpl,1,function(x){
  paste(x[1],str_split(x[2],'///',simplify=T),sep='...')})
X = tibble(unlist(test_function))
colnames(X) <- "ttt"
gpl<- separate(X,ttt,c("ID","ENTREZID"),sep = "\\...")
gpl=distinct(gpl,ID,.keep_all = T)
expression_data<-merge(gpl,expression_data,by.x=1,by.y=0)
expression_data<-expression_data[,-1]
expression_data<-na.omit(expression_data)
x1=expression_data$ENTREZID
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")
expression_data$SYMBOL<-x2$SYMBOL
expression_data<-na.omit(expression_data)
symbol <- x2[match(expression_data$ENTREZID,x2$ENTREZID),"SYMBOL"]
table(duplicated(symbol))
expression_data<- aggregate(expression_data, by=list(symbol), FUN=mean)
table(duplicated(expression_data$Group.1)) 
expression_data<- column_to_rownames(expression_data,'Group.1')
expression_data<-expression_data[,-c(1,2)]
write.csv(expression_data,file='GSE14323_expression_data.csv')
write.csv(pdata,file='GSE14323_pdata.csv')

gse_GSE164760 <- getGEO('GSE164760', destdir = ".", AnnotGPL = FALSE ,getGPL = F)
pdata = gse_GSE164760$GSE164760_series_matrix.txt.gz@phenoData@data
pdata<-filter(pdata,str_detect(pdata$characteristics_ch1,'Healthy')|str_detect(pdata$characteristics_ch1,'Cirrhotic'))
expression_data <- exprs(gse_GSE164760[[1]])
expression_data<-as.data.frame(expression_data)
expression_data<-expression_data[,rownames(pdata)]
soft=getGEO(filename ="GSE164760_family.soft/GSE164760_family.soft")
gpl=soft@gpls[["GPL13667"]]@dataTable@table
gpl<-gpl[,c('ID','Entrez Gene')]
test_function <- apply(gpl,1,function(x){
  paste(x[1],str_split(x[2],'///',simplify=T),sep='...')})
X = tibble(unlist(test_function))
colnames(X) <- "ttt"
gpl<- separate(X,ttt,c("ID","ENTREZID"),sep = "\\...")
gpl=distinct(gpl,ID,.keep_all = T)
expression_data<-merge(gpl,expression_data,by.x=1,by.y=0)
expression_data<-expression_data[,-1]
expression_data<-na.omit(expression_data)
x1=expression_data$ENTREZID
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")
expression_data$SYMBOL<-x2$SYMBOL
expression_data<-na.omit(expression_data)
symbol <- x2[match(expression_data$ENTREZID,x2$ENTREZID),"SYMBOL"]
table(duplicated(symbol))
expression_data<- aggregate(expression_data, by=list(symbol), FUN=mean)
table(duplicated(expression_data$Group.1)) 
expression_data<- column_to_rownames(expression_data,'Group.1')
expression_data<-expression_data[,-c(1,16)]
#write.csv(expression_data,file='GSE164760_expression_data.csv')
#write.csv(pdata,file='GSE164760_pdata.csv')

gse_GSE49541 <- getGEO('GSE49541', destdir = ".", AnnotGPL = FALSE ,getGPL = F,GSEMatrix = TRUE)
pdata = gse_GSE49541$GSE49541_series_matrix.txt.gz@phenoData@data
expression_data <- exprs(gse_GSE49541[[1]])
expression_data<-as.data.frame(expression_data)
expression_data<-expression_data[,rownames(pdata)]
soft=getGEO(filename ="GSE49541_family.soft.gz")
gpl=soft@gpls[["GPL570"]]@dataTable@table
gpl<-gpl[,c('ID','Gene Symbol')]
test_function <- apply(gpl,1,function(x){
  paste(x[1],str_split(x[2],'///',simplify=T),sep='...')})
X = tibble(unlist(test_function))
colnames(X) <- "ttt"
gpl<- separate(X,ttt,c("ID","Gene Symbol"),sep = "\\...")
gpl=distinct(gpl,ID,.keep_all = T)
expression_data<-merge(gpl,expression_data,by.x=1,by.y=0)
expression_data<-expression_data[,-1]
expression_data<-na.omit(expression_data)
expression_data1 <- aggregate(expression_data, 
                              by = list(`Gene Symbol` = expression_data$`Gene Symbol`), 
                              FUN = mean)
table(duplicated(expression_data1$`Gene Symbol`)) 
expression_data <- expression_data1[-1,] 
rownames(expression_data)<-expression_data$`Gene Symbol`
expression_data$`Gene Symbol`<-NULL
expression_data[,"Gene Symbol"]<-NULL
expression_data<-expression_data[,rownames(pdata)]
pdata<-arrange(pdata,desc(pdata$characteristics_ch1))
expression_data<-expression_data[,rownames(pdata)]
#write.csv(expression_data,file = 'GSE49541_expression_data.csv')
#write.csv(pdata,file = 'GSE49541_pdata.csv')

setwd("肝纤维化/肝纤维化转录组")
GSE49541_expression_data<-read.csv('GSE49541_expression_data.csv')
GSE49541_pdata<-read.csv('GSE49541_pdata.csv')
GSE49541_expression_data<-column_to_rownames(GSE49541_expression_data,'X')
GSE49541_pdata<-column_to_rownames(GSE49541_pdata,'X')
group_list=c(rep('Mild fibrosis',40),rep('Cirrhosis',32))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Mild fibrosis")
boxplot(GSE49541_expression_data,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
legend('topright',
       legend = levels(group_list), 
       fill = unique(as.numeric(group_list)), 
       title = "group",
       cex = 1,
       pt.cex  = 1.2,
       x.intersp  = 0.5,
       y.intersp  = 0.8,
       title.adj  = 0.5,
       text.font=2)
GSE49541_expression_data1=normalizeBetweenArrays(GSE49541_expression_data)
boxplot(GSE49541_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
legend('topright',
       legend = levels(group_list), 
       fill = unique(as.numeric(group_list)), 
       title = "group",
       cex = 1,
       pt.cex  = 1.2,
       x.intersp  = 0.5,
       y.intersp  = 0.8,
       title.adj  = 0.5,
       text.font = 2)
library(limma) 
library(sva)
design=model.matrix(~ group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(GSE49541_expression_data1)
fit=lmFit(GSE49541_expression_data1,design)
fit=eBayes(fit) 
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf)
getwd()
write.csv(allDiff,file = "deglog_GSE49541.csv")
#png(filename = "boxplot.png", width = 9, height = 6, units = "in", res = 800)
#par(mfrow=c(1,2))
#boxplot(GSE49541_expression_data,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
#legend('topright',
#       legend = levels(group_list), 
#       fill = unique(as.numeric(group_list)), 
#       title = "group",
#       cex = 2,
#       pt.cex  = 1.2,
#       x.intersp  = 0.5,
#       y.intersp  = 0.8,
#       title.adj  = 0.5,
#       text.font = 2)
#boxplot(GSE49541_expression_data1,outline=FALSE, notch=F,col=group_list, las=2,xaxt="n")
#legend('topright',
#       legend = levels(group_list), 
#       fill = unique(as.numeric(group_list)), 
#       title = "group",
#       cex = 2,
#       text.font = 2,
#       pt.cex  = 1.2,
#       x.intersp  = 0.5,
#       y.intersp  = 0.8,
#       title.adj  = 0.5)
#dev.off()

setwd("肝纤维化/肝纤维化转录组")
library(WGCNA)
GSE164760_expression_data<-read.csv('GSE164760_expression_data.csv')
GSE164760_pdata<-read.csv('GSE164760_pdata.csv')
GSE164760_expression_data<-column_to_rownames(GSE164760_expression_data,'X')
GSE164760_pdata<-column_to_rownames(GSE164760_pdata,'X')
GSE14323_expression_data<-read.csv('GSE14323_expression_data.csv')
GSE14323_pdata<-read.csv('GSE14323_pdata.csv')
GSE14323_expression_data<-column_to_rownames(GSE14323_expression_data,'X')
GSE14323_pdata<-column_to_rownames(GSE14323_pdata,'X')
new_expression_data<-merge(GSE164760_expression_data,GSE14323_expression_data,by.x=0,by.y=0)
new_expression_data<-column_to_rownames(new_expression_data,'Row.names')
GSE164760_expression_data<-GSE164760_expression_data[rownames(new_expression_data),]
expr<-GSE164760_expression_data
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { expr[which(expr <= 0)] <- NaN
expr <- log2(expr)
print('log2')}
GSE164760_expression_data1<-log2(GSE164760_expression_data+1)
GSE164760_pdata$source<-ifelse(str_detect(GSE164760_pdata$title,'Cirrhotic'),'Cirrhosis','Healthy')
p.n=filter(GSE164760_pdata,source=="Healthy")
p.t=filter(GSE164760_pdata,source=="Cirrhosis")
GSE164760_expression_data1=GSE164760_expression_data1[,c(rownames(p.n),rownames(p.t))]
group_list=c(rep('Healthy',6),rep('Cirrhosis',8))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Healthy")
dev.new()
boxplot(GSE164760_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
GSE164760_expression_data2<-normalizeBetweenArrays(GSE164760_expression_data1)
boxplot(GSE164760_expression_data2,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")

GSE14323_expression_data<-read.csv('GSE14323_expression_data.csv')
GSE14323_pdata<-read.csv('GSE14323_pdata.csv')
GSE14323_expression_data<-column_to_rownames(GSE14323_expression_data,'X')
GSE14323_pdata<-column_to_rownames(GSE14323_pdata,'X')
GSE14323_expression_data<-GSE14323_expression_data[rownames(new_expression_data),]
expr<-GSE14323_expression_data
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { expr[which(expr <= 0)] <- NaN
expr <- log2(expr)
print('log2')}
GSE14323_expression_data1<-GSE14323_expression_data
GSE14323_pdata$source<-ifelse(str_detect(GSE14323_pdata$source_name_ch1,'cirrhosis'),'Cirrhosis','Healthy')
p.n=filter(GSE14323_pdata,source=="Healthy")
p.t=filter(GSE14323_pdata,source=="Cirrhosis")
GSE14323_expression_data1=GSE14323_expression_data1[,c(rownames(p.n),rownames(p.t))]
group_list=c(rep('Healthy',19),rep('Cirrhosis',41))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Healthy")
boxplot(GSE14323_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
GSE14323_expression_data2<-normalizeBetweenArrays(GSE14323_expression_data1)
boxplot(GSE14323_expression_data2,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")

new_expression_data<-merge(GSE164760_expression_data2,GSE14323_expression_data2,by.x=0,by.y=0)
new_expression_data<-column_to_rownames(new_expression_data,'Row.names')
GSE164760_pdata<-data.frame(GEO=rownames(GSE164760_pdata),source=GSE164760_pdata$source,geo_accession='GSE164760')
GSE14323_pdata<-data.frame(GEO=rownames(GSE14323_pdata),source=GSE14323_pdata$source,geo_accession='GSE14323')
new_pdata<-rbind(GSE164760_pdata,GSE14323_pdata)
new_pdata<-arrange(new_pdata,desc(new_pdata$source),new_pdata$GEO)
new_pdata<-column_to_rownames(new_pdata,'GEO')
setwd("肝纤维化/肝纤维化转录组/WGCNA")
#write.csv(new_expression_data,file = 'new_expression_data1.csv')
#write.csv(new_pdata,'new_pdata1.csv')

setwd("肝纤维化/肝纤维化转录组/WGCNA")
new_expression_data<-read.csv('new_expression_data1.csv')
new_pdata<-read.csv('new_pdata1.csv')
new_expression_data<-column_to_rownames(new_expression_data,'X')
new_pdata<-column_to_rownames(new_pdata,'X')
p.n=filter(new_pdata,source=="Healthy")
p.t=filter(new_pdata,source=="Cirrhosis")
new_expression_data=new_expression_data[,c(rownames(p.n),rownames(p.t))]
new_pdata$source<-as.factor(new_pdata$source)
new_pdata$source <- relevel(new_pdata$source, ref="Healthy")
group_list=c(rep('Healthy',25),rep('Cirrhosis',49))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Healthy")
GEO<-as.factor(new_pdata$geo_accession)
library(tinyarray)
library(FactoMineR)
library(factoextra)
draw_pca(exp = new_expression_data, group_list = GEO)
new_expression_data1 <- new_expression_data 
batch <-GEO
design=model.matrix(~group_list)
rownames(design)<-colnames(new_expression_data1)
boxplot(new_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
new_expression_data2 <- ComBat(new_expression_data1,
                               batch=batch,
                               mod=design,
                               par.prior=TRUE, prior.plots=FALSE)
dev.new()
boxplot(new_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
boxplot(new_expression_data2,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
draw_pca(exp = new_expression_data2, group_list = GEO)
#png(filename = "boxplot.png", width = 9, height = 6, units = "in", res = 800)
#par(mfrow=c(1,2))
#
#boxplot(new_expression_data1,outline=FALSE, notch=T,col=group_list, las=2,xaxt="n")
#legend('topright',
#       legend = levels(group_list), 
#       fill = unique(as.numeric(group_list)), 
#       title = "group",
#       cex = 2,
#       pt.cex  = 1.2,
#       x.intersp  = 0.5,
#       y.intersp  = 0.8,
#       title.adj  = 0.5,
#       text.font = 2)
#
#
#boxplot(new_expression_data2,outline=FALSE, notch=F,col=group_list, las=2,xaxt="n")
#
#legend('topright',
#       legend = levels(group_list), 
#       fill = unique(as.numeric(group_list)), 
#       title = "group",
#       cex = 2,
#       text.font = 2,
#       pt.cex  = 1.2,
#       x.intersp  = 0.5,
#       y.intersp  = 0.8,
#       title.adj  = 0.5)
#dev.off()
#png(filename = "PCA.png", width = 8, height = 4, units = "in", res = 800)
#draw_pca(exp = new_expression_data1, group_list = batch)+theme(legend.position  = "none")+
#draw_pca(exp = new_expression_data2, group_list =  batch)&
# theme(legend.text  = element_text(size = 14,face = 'bold'),
#       legend.title = element_text(face='bold',size = 13))
#
#dev.off()

new_expression_data2 = t(new_expression_data2[order(apply(new_expression_data2, 1, mad), 
                                                    decreasing = T)[1:5000],])
new_expression_data2[1:4,1:4]
gene=new_expression_data2
gsg = goodSamplesGenes(gene, verbose = 3)
gsg$allOK
gsg$goodSamples
table(gsg$goodSamples)

powers <- 1:20
sft <- pickSoftThreshold(gene, powerVector = powers, verbose = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     type = 'n', 
     xlab = 'Soft Threshold (power)', 
     ylab = 'Scale Free Topology Model Fit, signed R^2', 
     main = paste('Scale Independence'),
     cex.lab  = 1.2,
     cex.main  = 1.5,
     font.lab  = 2,
     font.main  = 2)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, 
     col = 'red',
     cex = 1.2,
     font = 2)
abline(h = 0.90, col = 'red')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = 'Soft Threshold (power)', 
     ylab = 'Mean Connectivity', 
     type = 'n',
     main = paste('Mean Connectivity'),
     cex.lab  = 1.2,
     cex.main  = 1.5,
     font.lab  = 2,
     font.main  = 2) 
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels = powers, 
     col = 'red',
     cex = 1.2,
     font = 2)
sft$powerEstimate
#png(filename = "power.png", width = 11, height = 6, units = "in", res = 900)
#par(mfrow=c(1,2))
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
#     type = 'n', 
#     xlab = 'Soft Threshold (power)', 
#     ylab = 'Scale Free Topology Model Fit, signed R^2', 
#     main = paste('Scale Independence'),
#     cex.lab  = 1.2,
#     cex.main  = 1.5,
#     font.lab  = 2,
#     font.main  = 2)
#
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels = powers, 
#     col = 'red',
#     cex = 1.2,
#     font = 2)
#
#abline(h = 0.90, col = 'red')
#
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
#     xlab = 'Soft Threshold (power)', 
#     ylab = 'Mean Connectivity', 
#     type = 'n',
#     main = paste('Mean Connectivity'),
#     cex.lab  = 1.2,
#     cex.main  = 1.5,
#     font.lab  = 2,
#     font.main  = 2) 
#
#text(sft$fitIndices[,1], sft$fitIndices[,5], 
#     labels = powers, 
#     col = 'red',
#     cex = 1.2,
#     font = 2)
#dev.off()

power = sft$powerEstimate
power

powers <- sft$powerEstimate
adjacency <- adjacency(gene, power = powers)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]
#write.table(tom_sim, 'TOMsimilarity.txt', sep = '\t', col.names = NA, quote = FALSE)

tom_dis  <- 1 - tom_sim
par(mfrow = c(1, 1))
geneTree <- hclust(as.dist(tom_dis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)

minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(
  geneTree, 
  dynamicColors, 
  groupLabels='Dynamic Tree Cut',
  dendroLabels = FALSE, 
  addGuide = TRUE, 
  hang = 0.03, 
  guideHang = 0.05,
  main = 'Gene dendrogram and module colors',
  cex.colorLabels = 1.5,
  cex.main  = 2,
  cex.lab  = 1.8,
  cex.axis  = 1.5,
  cex.rowText  = 1.7,
  font.main  = 2,
  font.lab  = 2,
  font.axis  = 2
)
#png(filename = "plotDendroAndColors.png", width = 14, height =8, units = "in", res = 900)
#par(font=2,mar=c(5,6,4,1))
#plotDendroAndColors(
#  geneTree, 
#  dynamicColors, 
#  groupLabels='Dynamic Tree Cut',
#  dendroLabels = FALSE, 
#  addGuide = TRUE, 
#  hang = 0.03, 
#  guideHang = 0.05,
#  main = 'Gene dendrogram and module colors',
#  cex.main  = 2,
#  cex.lab  = 1.8,
#  cex.axis  = 1.5,
#  cex.rowText  = 1.7,
#  font.main  = 2,
#  font.lab  = 2,
#  font.axis  = 2,
#  cex.colorLabels = 1
#)
#dev.off()

plot_sim <- -(1-tom_sim)
diag(plot_sim) <- NA
TOMplot(plot_sim, geneTree, dynamicColors,
        main = 'Network heatmap plot, selected genes')

MEList <- moduleEigengenes(gene, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]
#write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)

ME_cor <- cor(MEs)
ME_cor[1:6,1:6]
METree <- hclust(as.dist(1-ME_cor), method = 'average')
dev.new()
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')
abline(h = 0.4, col = 'blue')
abline(h = 0.25, col = 'red')
plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90,cex.adjacency  = 1.5,
                      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))
#png(filename = "plotEigengeneNetworks.png", width = 12, height = 8, units = "in", res = 900)
#par(font.lab  = 2, font.axis  = 2, cex.lab  = 1.2, cex.main  = 1.3,mar=c(5,8,2,1))
#
#plotEigengeneNetworks(
#  MEs, 
#  "", 
#  cex.lab  = 1.2,
#  xLabelsAngle = 90,
#  cex.adjacency  = 1.5,
#  marDendro = c(0, 4, 1, 2), 
#  marHeatmap = c(3, 4, 1, 2),
#  plotAdjacency = TRUE,
#  plotDendrograms = TRUE 
#)
#dev.off()

merge_module <- mergeCloseModules(gene, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge_module$colors
table(mergedColors)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)
#png(filename = "plotDendroAndColors2.png", width = 14, height = 8, units = "in", res = 900)
#par(font=2)
#
#plotDendroAndColors(geneTree, 
#                    cbind(dynamicColors, mergedColors),
#                    c('Dynamic Tree Cut', 'Merged dynamic'),
#                      dendroLabels = FALSE, 
#                      addGuide = TRUE, 
#                      hang = 0.03, 
#                      guideHang = 0.05,
#                      main = 'Gene dendrogram and module colors',
#                      cex.main  = 2,
#                      cex.lab  = 1.8,
#                      cex.axis  = 1.5,
#                      cex.rowText  = 1.7,
#                      font.main  = 2,
#                      font.lab  = 2,
#                      font.axis  = 2,
#                      cex.colorLabels = 1
#                    )
#         
#dev.off()

module <- merge_module$newMEs
new_pdata$source<-as.factor(new_pdata$source)
new_pdata$source <- relevel(new_pdata$source, ref="Healthy")
trait=new_pdata
nGenes = ncol(gene)
nSamples = nrow(gene)
x=as.factor(new_pdata$source) 
design=model.matrix(~0+ new_pdata$source)
colnames(design)=levels(x)
rownames(design)=colnames(new_expression_data1)
trait=design
module <- merge_module$newMEs
moduleTraitCor <- cor(module, trait, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))
#write.table(moduleTraitCor, 'moduleTraitCor.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(moduleTraitPvalue, 'moduleTraitPvalue.txt', sep = '\t', col.names = NA, quote = FALSE)
textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 7, 3, 1))
dev.new()
labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait relationships'),
               xLabels =colnames(design), yLabels = names(module), ySymbols = names(module),
               colorLabels = FALSE, colors = greenWhiteRed(50), cex.text = 1.2, zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = FALSE)
dev.off()
getwd()
#png(filename = "labeledHeatmap.png", width = 8, height = 8, units = "in", res = 800)
#pdf(file = "labeledHeatmap.pdf", width = 8, height = 8)
#par(font=2,mar=c(4,7,3,1))
#labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait relationships'),
#               xLabels =colnames(design), yLabels = names(module), ySymbols = names(module),
#               colorLabels = FALSE, colors = greenWhiteRed(50), cex.text = 1.2, zlim = c(-1,1),
#               textMatrix = textMatrix, setStdMargins = FALSE)
#dev.off()

table(mergedColors)
gene_module <- data.frame(gene_name = colnames(gene), module = mergedColors, stringsAsFactors = FALSE)
head(gene_module)
table(gene_module$module)
gene_module_select_brown <- subset(gene_module, module == c('brown'))$gene_name
#write.csv(gene_module_select_brown,'WGCNA_brown_gene.csv')
gene_module_select_blue <- subset(gene_module, module == c('blue'))$gene_name
#write.csv(gene_module_select_blue,'WGCNA_blue_gene.csv')

library(VennDiagram)
setwd("肝纤维化/肝纤维化转录组")
blue<-read.csv('WGCNA/WGCNA_blue_gene.csv')
RNA_seq_deg<-read.csv('差异分析/deglog_GSE49541.csv')
WGCNA_deg<-blue
RNA_seq_deg<-column_to_rownames(RNA_seq_deg,'X')
WGCNA_deg<-column_to_rownames(WGCNA_deg,'X')
gene_list<-read.table(file = '韦恩图/GeneList (1).txt',header = T,sep = '\t')
RNA_seq_deg1<-filter(RNA_seq_deg,RNA_seq_deg$logFC>=0.25&RNA_seq_deg$adj.P.Val<=0.05)
a1<-rownames(RNA_seq_deg1)
b1<-WGCNA_deg$x
c1<-gene_list$Symbol
library(VennDiagram)
venn_list <- list(a = a1, b = b1, c = c1)  
venn.diagram(  
  x = venn_list,  
  filename = 'venn1.png',
  imagetype = 'png',
  height = 2000,
  width = 2000,
  resolution = 300,
  fill = c("#ffb2b2","#b2e7cb","#b2d4ec"),
  alpha = 0.50,
  category.names = c('RNA_seq', 'WGCNA', 'IRG'),
  cat.col = rep('black', 3),
  col = 'black',
  cex = 1.2,
  fontfamily = 'serif',
  cat.cex = 1.2,
  cat.fontfamily = 'serif',
  cat.dist = c(0.1, 0.1, 0.1),
  margin = 0.1,
  main='UP'
)
brown<-read.csv('WGCNA/WGCNA_brown_gene.csv')
WGCNA_deg<-brown
WGCNA_deg<-column_to_rownames(WGCNA_deg,'X')
RNA_seq_deg2<-filter(RNA_seq_deg,RNA_seq_deg$logFC<=-0.25&RNA_seq_deg$adj.P.Val<=0.05)
a2<-rownames(RNA_seq_deg2)
b2<-WGCNA_deg$x
c2<-gene_list$Symbol
venn_list <- list(a = a2, b = b2, c = c2)  
venn.diagram(  
  x = venn_list,  
  filename = 'venn2.png',
  imagetype = 'png',
  height = 2000,
  width = 2000,
  resolution = 300,
  fill = c("#ffb2b2","#b2e7cb","#b2d4ec"),
  alpha = 0.50,
  category.names = c('RNA_seq', 'WGCNA', 'IRG'),
  cat.col = rep('black', 3),
  col = 'black',
  cex = 1.2,
  fontfamily = 'serif',
  cat.cex = 1.2,
  cat.fontfamily = 'serif',
  cat.dist = c(0.1, 0.1, 0.1),
  margin = 0.1,
  main='DOWN'
)
Core_gene_up<- Reduce(intersect, list(a1, b1, c1))
Core_gene_down<- Reduce(intersect, list(a2, b2, c2))
Core_gene<-c(Core_gene_up,Core_gene_down)
#write.csv(Core_gene_up,file = 'Core_gene_up.csv')
#write.csv(Core_gene_down,file = 'Core_gene_down.csv')
#write.csv(Core_gene,file = 'Core_gene_25.csv')
#write.csv(RNA_seq_deg_core,file ='RNA_seq_deg_core.csv' )

GSE49541_expression_data<-read.csv('GSE49541_expression_data.csv')
GSE49541_expression_data<-column_to_rownames(GSE49541_expression_data,'X')
GSE49541_pdata<-read.csv('GSE49541_pdata.csv')
GSE49541_pdata<-column_to_rownames(GSE49541_pdata,'X')
group_list=factor(c(rep('Mild fibrosis',40),rep('Cirrhosis',32)))
data<-GSE49541_expression_data[Core_gene,]
data1<-normalizeBetweenArrays(data)
data2<-t(data1)
data2<-as.data.frame(data2)
table(rownames(data2)==colnames(GSE49541_expression_data))
table(rownames(data2)==rownames(GSE49541_pdata))
data2$class=group_list
data2$class <- as.factor(ifelse(data2$class == 'Mild fibrosis', 0, 1))
x <- as.matrix(data2[, -which(names(data2) == "class")])
y <- data2$class
cv_control <- trainControl(method = "cv", number = 10)
set.seed(1234)
lasso_model <- train(x, y, method = "glmnet",
                     trControl = cv_control,
                     tuneLength = 10)
best_lambda <- lasso_model$bestTune$lambda
y_numeric <- as.numeric(as.character(y)) - 1
library(glmnet)
final_model <- glmnet(x, y_numeric, alpha = 1, lambda = best_lambda, family = "binomial")
final_model1 <- glmnet(x, y_numeric, alpha = 1,  family = "binomial")
print(coef(final_model))
plot(final_model1, xvar = "lambda", label = TRUE)
abline(v = log(best_lambda), lty = 2, col = "red")

coef.min<-coef(final_model)
coef.min1<-data.frame(coef.min)
active.min <- which(coef.min@i != 0)
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]
lasso_geneids <- lasso_geneids[-1] %>% as.data.frame()
colnames(lasso_geneids) <- 'symbol'
lasso<-data.frame(symbol=lasso_geneids$symbol,s0=coef.min1[lasso_geneids$symbol,])

coefs <- as.matrix(coef(final_model, s = best_lambda))
nonzero_coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]
nonzero_coefs <- nonzero_coefs[-1, , drop = FALSE]
coef_df <- data.frame(
  feature = factor(rownames(nonzero_coefs), levels = rownames(nonzero_coefs)[order(nonzero_coefs)]),
  coefficient = as.numeric(nonzero_coefs)
)
ggplot(coef_df, aes(x = feature, y = coefficient, fill = coefficient > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("tomato", "dodgerblue"), guide = FALSE) +
  coord_flip() +
  labs(
    title = "Lasso Regression Coefficients",
    subtitle = paste("Lambda =", round(best_lambda, 5)),
    x = "Features",
    y = "Coefficient Value"
  ) +
  theme_classic() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle  = element_text(hjust = 0.5, size = 14, face = "italic"),
    axis.text.y  = element_text(face = "bold", size = 17),
    axis.text.x  = element_text(face = "bold", size = 17),
    axis.title.x  = element_text(size = 16, face = "bold"),
    axis.title.y  = element_text(size = 16, face = "bold"),
    panel.background  = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    legend.position  = "none"
  )
#ggsave(filename = 'bar.pdf',width = 12,height = 9,units = 'in')

par(mar = c(5, 5, 4, 2) )
plot(final_model1, 
     xvar = "lambda", 
     main = "",
     xlab = "",
     ylab = "",
     cex.axis  = 1.5,
     font.axis  = 2
)
title(xlab = "Log Lambda",
      ylab = "Coefficients",
      cex.lab  = 1.8,
      font.lab  = 2,
      col.lab  = "black"
)
title(main = "Regularization Path", 
      line = 2.5, 
      cex.main  = 2.0, 
      font.main  = 2)
abline(v = log(best_lambda), 
       lty = 2, 
       col = "red", 
       lwd = 2)
text(x = log(best_lambda), 
     y = max(coef(final_model1)[-1, ]), 
     labels = "Optimal Lambda", 
     pos = 4, 
     col = "red", 
     cex = 1.5, 
     font = 2)
#tiff(width = 12,height = 9,units = 'in',res = 900,filename ='Regularization Path.tiff')
#pdf(file = "Regularization Path.pdf",  width = 12, height = 9)
#plot(final_model1, xvar = "lambda", main = "")
#mtext("Regularization Path", side = 3, line = 2, cex = 1)
#abline(v = log(best_lambda), lty = 2, col = "red")
#text(x = log(best_lambda), y = max(coef(final_model)[-1, ]), 
#     labels = "Optimal Lambda", pos = 4, col = "red")
#dev.off()

coef_path <- as.matrix(final_model1$beta)
coef_path <- as.data.frame(t(coef_path))
coef_path$lambda <- final_model1$lambda
library(tidyr)
coef_path_long <- pivot_longer(
  coef_path,
  cols = -lambda,
  names_to = "variable",
  values_to = "coefficient"
)
library(ggplot2)
ggplot(coef_path_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line() +
  geom_vline(xintercept = log(best_lambda), linetype = "dashed", color = "red") +
  labs(x = "Log Lambda", y = "Coefficients", title = "Lasso Coefficient Path") +
  theme_classic() +
  theme(legend.position = "none")
#ggsave(filename = 'Lasso Coefficient Path.pdf',width = 12,height = 9,units = 'in')

cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
cv_lasso 
par(mar = c(5, 5, 5.5, 2),
    font.axis  = 2,
    cex.axis  = 1.3,
    mgp = c(3, 1, 0)
)
plot(cv_lasso, 
     main = "Cross-Validation Curve",
     xlab = "log(λ)", 
     ylab = "Binomial Deviance",
     col = "blue", 
     cex.lab  = 1.3,
     cex.main  = 1.5,
     font.main  = 2
)
abline(v = log(best_lambda), lty = 2, col = "red", lwd = 2)
text(log(best_lambda), max(cv_lasso$cvm), 
     labels = paste("λ =", round(best_lambda, 5)), 
     pos = 4, col = "red",cex = 2)
lasso<-lasso[lasso$s0>0,]
#write.csv(lasso,file = 'lasso.csv')

library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(GSVA)
library(IOBR)
library(linkET)
setwd("肝纤维化/肝纤维化转录组/免疫浸润")
load('geneset.RData')
setwd("肝纤维化/肝纤维化转录组")
lasso<-read.csv('肝纤维化/肝纤维化转录组/机器学习/X/lasso.csv')
new_expression_data<-read.csv('new_expression_data1.csv')
group_list=c(rep('Healthy',25),rep('Cirrhosis',49))
new_expression_data<-column_to_rownames(new_expression_data,'X')
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Healthy")
GSE49541_expression_data<-read.csv('GSE49541_expression_data.csv')
GSE49541_pdata<-read.csv('GSE49541_pdata.csv')
GSE49541_expression_data<-column_to_rownames(GSE49541_expression_data,'X')
GSE49541_pdata<-column_to_rownames(GSE49541_pdata,'X')
group_list1=c(rep('Mild fibrosis',40),rep('Cirrhosis',32))
group_list1=factor(group_list1)
group_list1 <- relevel(group_list1, ref="Mild fibrosis")

library(IOBR)
tme_deconvolution_methods
set.seed(1234)
im_cibersort <- deconvo_tme(eset = new_expression_data,
                            method = "cibersort",
                            arrays = T,
                            perm = 1000
)
library(tidyr)
#write.csv(im_cibersort,file = 'im_cibersort.csv')
im_cibersort<-read.csv('im_cibersort_new.csv')
im_cibersort<-column_to_rownames(im_cibersort,'X')

im_cibersort$group<-group_list
df_long <- melt(im_cibersort, id.vars  =colnames(im_cibersort)[c(1,24:27)])
ggplot(df_long, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ group, scales = "free_x") +
  labs(fill = "celltype") +
  guides(fill = guide_legend(ncol = 1,
                             label.theme  = element_text(size = 12,face = 'bold'))) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    strip.text  = element_text(size = 13,
                               face = "bold"),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(size = 10)
  )
#ggsave(filename = 'barplot_new_expression.tiff',width = 9,height = 9,dpi = 1600,units = 'in')
#ggsave(filename = 'barplot_new_expression.pdf',width = 9,height = 9)

library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(GSVA)
library(IOBR)
library(linkET)
setwd("肝纤维化/肝纤维化转录组/免疫浸润")
lasso<-read.csv('肝纤维化/肝纤维化转录组/机器学习/X/lasso.csv')
load('geneset.RData')
new_expression_data<-read.csv('new_expression_data_14_16.csv')
new_expression_data<-column_to_rownames(new_expression_data,'X')

lasso<-arrange(lasso,desc(lasso$s0))
geneset_lasso<-list(lasso_symbol=lasso$symbol[1:8])
new_expression_data<- as.matrix(new_expression_data) 
param <- ssgseaParam(exprData = new_expression_data,
                     geneSets = geneset_lasso)
im_ssgsea_lasso  <- gsva(param)
im_ssgsea_lasso<-as.data.frame(im_ssgsea_lasso)
param_gse4954 <- ssgseaParam(exprData = GSE49541_expression_data1,
                             geneSets = geneset_lasso)
im_ssgsea_lasso_gse4954  <- gsva(param_gse4954)
im_ssgsea_lasso_gse4954<-as.data.frame(im_ssgsea_lasso_gse4954)
im_ssgsea_lasso<-t(im_ssgsea_lasso)
im_ssgsea_lasso<-as.data.frame(im_ssgsea_lasso)
im_ssgsea_lasso$group<-group_list
im_ssgsea_lasso$GSE<-'GSE164760&GSE14323'
group_list1<-as.data.frame(group_list1)
im_ssgsea_lasso_gse4954<-t(im_ssgsea_lasso_gse4954)
im_ssgsea_lasso_gse4954<-as.data.frame(im_ssgsea_lasso_gse4954)
im_ssgsea_lasso_gse4954$group<-group_list1$group_list1
im_ssgsea_lasso_gse4954$GSE<-'GSE49541'
im_ssgsea_lasso<-rbind(im_ssgsea_lasso,im_ssgsea_lasso_gse4954)
#write.csv(im_ssgsea_lasso,file = 'im_ssgsea_lasso.csv')
im_ssgsea_lasso<-read.csv('im_ssgsea_lasso.csv')
im_ssgsea_lasso<-column_to_rownames(im_ssgsea_lasso,'X')
p1<-ggplot(im_ssgsea_lasso[im_ssgsea_lasso$GSE==names(table(im_ssgsea_lasso$GSE))[1],], aes(x = GSE, y = lasso_symbol, fill = group)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name = "Group"
  ) +
  labs(y = 'ssGSEA Score') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "left",
    legend.text  = element_text(size = 8, face = 'bold'),
    legend.title  = element_blank(),
    axis.title.y  = element_text(size = 8),
    axis.text.y  = element_text(size = 6),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(size = 6, color = "black"),
  )
p1
p2<-ggplot(im_ssgsea_lasso[im_ssgsea_lasso$GSE==names(table(im_ssgsea_lasso$GSE))[2],], aes(x = GSE, y = lasso_symbol, fill = group)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name = "Group",
    labels = c("Cirrhosis","Mild fibrosis" ) 
  ) +
  labs(y = 'ssGSEA Score') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "left",
    legend.text  = element_text(size = 8, face = 'bold'),
    legend.title  = element_blank(),
    axis.title.y  = element_text(size = 8),
    axis.text.y  = element_text(size = 6),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(size = 6, color = "black"),
  )
p2
p2+p1

new_expression_data_lasso<-new_expression_data[lasso$symbol,]
GSE49541_expression_data1_lasso<-GSE49541_expression_data1[lasso$symbol,]
new_expression_data_lasso<-t(new_expression_data_lasso)
new_expression_data_lasso<-as.data.frame(new_expression_data_lasso)
new_expression_data_lasso$group<-group_list
new_expression_data_lasso$GSE<-'GSE164760&GSE14323'
GSE49541_expression_data1_lasso<-t(GSE49541_expression_data1_lasso)
GSE49541_expression_data1_lasso<-as.data.frame(GSE49541_expression_data1_lasso)
GSE49541_expression_data1_lasso$group<-group_list1$group_list1
GSE49541_expression_data1_lasso$GSE<-'GSE49541'
new_expression_data_lasso<-rbind(new_expression_data_lasso,GSE49541_expression_data1_lasso)
new_expression_data_lasso$id<-rownames(new_expression_data_lasso)
new_expression_data_lasso <- new_expression_data_lasso %>%
  pivot_longer(
    cols = -c(id,group, GSE),
    names_to = "gene",
    values_to = "expression_value"
  )
#write.csv(new_expression_data_lasso,file = 'new_expression_data_lasso.csv')
new_expression_data_lasso<-read.csv('new_expression_data_lasso.csv')
new_expression_data_lasso<-column_to_rownames(new_expression_data_lasso,'X')
new_expression_data_lasso1<-new_expression_data_lasso[new_expression_data_lasso$GSE=='GSE164760&GSE14323',]
new_expression_data_lasso2<-new_expression_data_lasso[new_expression_data_lasso$GSE=='GSE49541',]
ggplot(new_expression_data_lasso1, aes(x = gene, y = expression_value, fill = group)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name =  names(table(new_expression_data_lasso1$GSE))
  ) +
  labs(y = 'Expression Value') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "center",
    legend.text  = element_text(size = 12, face = 'bold'),
    legend.title  = element_text(size = 12, face = 'bold'),
    axis.title.y  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(size = 10, color = "black"),
  )
ggplot(new_expression_data_lasso2, aes(x = gene, y = expression_value, fill = group)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name =  names(table(new_expression_data_lasso2$GSE)),
    labels = c("Cirrhosis","Mild fibrosis" ) 
  ) +
  labs(y = 'Expression Value') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "center",
    legend.text  = element_text(size = 12, face = 'bold'),
    legend.title  = element_text(size = 12, face = 'bold'),
    axis.title.y  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(size = 10, color = "black"),
  )
#ggsave(filename = 'boxplot.tiff',height = 8,width = 10)

new_expression_data<- as.matrix(new_expression_data) 
param <- ssgseaParam(exprData = new_expression_data,
                     geneSets = geneset)
im_ssgsea  <- gsva(param)
im_ssgsea<-as.data.frame(im_ssgsea)
im_ssgsea<-t(im_ssgsea)
im_ssgsea<-cbind(im_ssgsea,im_ssgsea_lasso[im_ssgsea_lasso$GSE=="GSE164760&GSE14323",])
im_ssgsea<-as.data.frame(im_ssgsea)
im_ssgsea$GSE<-NULL
im_ssgsea$id<-rownames(im_ssgsea)
im_ssgsea_long <- im_ssgsea %>%
  pivot_longer(
    cols = -c(id,group,lasso_symbol),
    names_to = "cell",
    values_to = "ssGSEA Score"
  )
#write.csv(im_ssgsea_long,file='im_ssgsea_long.csv')
im_ssgsea_long<-read.csv('im_ssgsea_long.csv')
im_ssgsea_long<-column_to_rownames(im_ssgsea_long,'X')
ggplot(im_ssgsea_long, aes(x = cell, y = ssGSEA.Score, fill = group)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name = "Group"
  ) +
  labs(y = 'ssGSEA Score') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "center",
    legend.text  = element_text(size = 12, face = 'bold'),
    legend.title  = element_blank(),
    axis.title.y  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(
      size = 10, 
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.margin  = margin(10, 10, 10, 10)
  )
ggplot(im_ssgsea_long[im_ssgsea_long$score!='none',], aes(x = cell, y = ssGSEA.Score, fill = score)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.fill  = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = score),
    method = "wilcox.test", 
    label = "p.signif",
    size = 4.5,
    show.legend  = FALSE
  ) +
  scale_fill_manual(
    values = c("#1d4a9b", "#e5171a"),
    name = "Group"
  ) +
  labs(y = 'ssGSEA Score') +
  theme_bw()+
  theme(
    legend.position  = "top",
    legend.direction  = "horizontal",
    legend.justification  = "center",
    legend.text  = element_text(size = 12, face = 'bold'),
    legend.title  = element_blank(),
    axis.title.y  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    axis.title.x  = element_blank(),
    axis.text.x  = element_text(
      size = 10, 
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.margin  = margin(10, 10, 10, 10)
  )

mygene <- 'lasso_symbol'
new_expression_data<- as.matrix(new_expression_data) 
param <- ssgseaParam(exprData = new_expression_data,
                     geneSets = geneset)
im_ssgsea  <- gsva(param)
im_ssgsea<-as.data.frame(im_ssgsea)
nc = t(rbind(im_ssgsea,t(im_ssgsea_lasso[im_ssgsea_lasso$GSE=='GSE164760&GSE14323',])[mygene,]))
m = rcorr(nc)$r[1:nrow(im_ssgsea),(ncol(nc)-length(mygene)+1):ncol(nc)]
p = rcorr(nc)$P[1:nrow(im_ssgsea),(ncol(nc)-length(mygene)+1):ncol(nc)]
head(p)
p<-as.matrix(p)
colnames(p)<-mygene
tmp <- matrix(case_when(as.vector(p) <= 0.01 ~ "**",
                        as.vector(p) <= 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col ='45',
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               cluster_rows=F,
               cluster_cols = F)
p1
#ggsave(filename = 'pheatmap_lasso.tiff',plot = p1,height = 5,width = 10)

setwd("肝纤维化/肝纤维化转录组/单细胞分析")
human_ctrl1<-Read10X('GSM5741404_human_ctrl1/human_ctrl1/')
human_ctrl2<-Read10X('GSM5741405_human_ctrl2/human_ctrl2/')
human_cir1<-Read10X('GSM5741406_human_cir1/human_cir1/')
human_cir2<-Read10X('GSM5741407_human_cir2/human_cir2/')

human_ctrl1_object<- CreateSeuratObject(counts = human_ctrl1,
                                        project = "human_ctrl1", 
                                        min.cells = 0,
                                        min.features = 200
)
human_ctrl1_object[["percent.mt"]] <- PercentageFeatureSet(human_ctrl1_object, 
                                                           pattern = "^MT-")
hist(human_ctrl1_object[["percent.mt"]]$percent.mt)
VlnPlot(human_ctrl1_object, features = c("nFeature_RNA",
                                         "nCount_RNA",
                                         "percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(human_ctrl1_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_ctrl1_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

human_ctrl1_val<- subset(human_ctrl1_object, 
                         subset = nFeature_RNA > 200 &
                           nFeature_RNA < 7000 & 
                           percent.mt < 15)
VlnPlot(human_ctrl1_val, features = c("nFeature_RNA",
                                      "nCount_RNA",
                                      "percent.mt"),
        ncol = 3)
colnames(human_ctrl1_val)
rownames(human_ctrl1_val@meta.data)
table(colnames(human_ctrl1_val)==rownames(human_ctrl1_val@meta.data))
colname<-paste("human_ctrl1_",colnames(human_ctrl1_val),sep="")
colname
human_ctrl1_val<-RenameCells(object = human_ctrl1_val,colname)
plot3 <- FeatureScatter(human_ctrl1_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(human_ctrl1_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3,plot4))

human_ctrl2_object<- CreateSeuratObject(counts = human_ctrl2,
                                        project = "human_ctrl2", 
                                        min.cells = 0,
                                        min.features = 200
)
human_ctrl2_object[["percent.mt"]] <- PercentageFeatureSet(human_ctrl2_object, 
                                                           pattern = "^MT-")
hist(human_ctrl2_object[["percent.mt"]]$percent.mt)
VlnPlot(human_ctrl2_object, features = c("nFeature_RNA",
                                         "nCount_RNA",
                                         "percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(human_ctrl2_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_ctrl2_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

human_ctrl2_val<- subset(human_ctrl2_object, 
                         subset = nFeature_RNA > 200 &
                           nFeature_RNA < 7000 & 
                           percent.mt < 15)
VlnPlot(human_ctrl2_val, features = c("nFeature_RNA",
                                      "nCount_RNA",
                                      "percent.mt"),
        ncol = 3)
colnames(human_ctrl2_val)
rownames(human_ctrl2_val@meta.data)
table(colnames(human_ctrl2_val)==rownames(human_ctrl2_val@meta.data))
colname<-paste("human_ctrl2_",colnames(human_ctrl2_val),sep="")
colname
human_ctrl2_val<-RenameCells(object = human_ctrl2_val,colname)
plot3 <- FeatureScatter(human_ctrl2_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(human_ctrl2_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3,plot4))

human_cir1_object<- CreateSeuratObject(counts = human_cir1,
                                       project = "human_cir1", 
                                       min.cells = 0,
                                       min.features = 200
)
human_cir1_object[["percent.mt"]] <- PercentageFeatureSet(human_cir1_object, 
                                                          pattern = "^MT-")
hist(human_cir1_object[["percent.mt"]]$percent.mt)
VlnPlot(human_cir1_object, features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(human_cir1_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_cir1_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

human_cir1_val<- subset(human_cir1_object, 
                        subset = nFeature_RNA > 200 &
                          nFeature_RNA < 7000 & 
                          percent.mt < 15)
VlnPlot(human_cir1_val, features = c("nFeature_RNA",
                                     "nCount_RNA",
                                     "percent.mt"),
        ncol = 3)
colnames(human_cir1_val)
rownames(human_cir1_val@meta.data)
table(colnames(human_cir1_val)==rownames(human_cir1_val@meta.data))
colname<-paste("human_cir1_",colnames(human_cir1_val),sep="")
colname
human_cir1_val<-RenameCells(object = human_cir1_val,colname)
plot3 <- FeatureScatter(human_cir1_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(human_cir1_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3,plot4))

human_cir2_object<- CreateSeuratObject(counts = human_cir2,
                                       project = "human_cir2", 
                                       min.cells = 0,
                                       min.features = 200
)
human_cir2_object[["percent.mt"]] <- PercentageFeatureSet(human_cir2_object, 
                                                          pattern = "^MT-")
hist(human_cir2_object[["percent.mt"]]$percent.mt)
VlnPlot(human_cir2_object, features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(human_cir2_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_cir2_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

human_cir2_val<- subset(human_cir2_object, 
                        subset = nFeature_RNA > 200 &
                          nFeature_RNA < 7000 & 
                          percent.mt < 15)
VlnPlot(human_cir2_val, features = c("nFeature_RNA",
                                     "nCount_RNA",
                                     "percent.mt"),
        ncol = 3)
colnames(human_cir2_val)
rownames(human_cir2_val@meta.data)
table(colnames(human_cir2_val)==rownames(human_cir2_val@meta.data))
colname<-paste("human_cir2_",colnames(human_cir2_val),sep="")
colname
human_cir2_val<-RenameCells(object = human_cir2_val,colname)
plot3 <- FeatureScatter(human_cir2_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(human_cir2_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3,plot4))

sce_mergeTEN<- merge(human_ctrl1_val,
                     y=c(human_ctrl2_val,
                         human_cir1_val,
                         human_cir2_val),
                     project = "scTEN")
Idents(sce_mergeTEN)<-sce_mergeTEN$orig.ident
sce_mergeTEN <- NormalizeData(sce_mergeTEN,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
sce_mergeTEN <- FindVariableFeatures(sce_mergeTEN,
                                     selection.method = "vst", 
                                     nfeatures = 2000)
sce_mergeTEN <- ScaleData(sce_mergeTEN,vars.to.regress = c('percent.mt'))
sce_mergeTEN <- RunPCA(sce_mergeTEN)
library(harmony)
str(sce_mergeTEN@meta.data)
sce_mergeTEN_harmony <- IntegrateLayers(object = sce_mergeTEN,
                                        method = HarmonyIntegration,
                                        orig.reduction = 'pca',
                                        new.reduction = 'harmony',
                                        group.by.vars="orig.ident",verbose=F)
sce_mergeTEN_harmony <- RunHarmony(object = sce_mergeTEN,
                                   group.by.vars  = "orig.ident", 
                                   reduction.save  = "harmony",
                                   orig.reduction = 'pca'
)
sce_mergeTEN_harmony[['RNA']] <- JoinLayers(sce_mergeTEN_harmony[['RNA']])
data.use <- Stdev(object = sce_mergeTEN_harmony, reduction = 'pca')
a=0
b=0
for (i in data.use) {
  a<-a+i
  b=b+1
  if(a/sum(data.use)>=0.85){
    print(b)
    break
  }
}
ElbowPlot(sce_mergeTEN_harmony,ndims = 50)
sce_mergeTEN_harmony <- FindNeighbors(sce_mergeTEN_harmony,reduction = 'harmony',dims = 1:40)
sce_mergeTEN_harmony <- FindClusters(sce_mergeTEN_harmony,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
sce_mergeTEN_harmony <- RunUMAP(sce_mergeTEN_harmony,dims = 1:40,reduction = 'harmony')
sce_mergeTEN_harmony$source<-ifelse(sce_mergeTEN_harmony@meta.data$orig.ident=='human_cir1'|sce_mergeTEN_harmony@meta.data$orig.ident=='human_cir2','Cirrhosis','Healthy')
library(clustree)
clustree(sce_mergeTEN_harmony) 
DimPlot(sce_mergeTEN_harmony,reduction = 'umap',group.by = 'orig.ident') +
  ggtitle('harmony')

marker_cluster=c(
  'CD3D',
  'FCGR3A','NCAM1','NKG7','NCR1',"KLRD1","GNLY",
  'IGHG1', 'MZB1',
  'CD79A','CD79B',"MS4A1",
  'VWF','PECAM1','PLPP1',
  'EPCAM','KRT19','SOX9',
  "TTR","ALB","APOE",
  'CD14','CD68','CD163',
  'MARCO','TIMD4',
  'FCGR3B', 'CSF3R',
  'CD1C', 'FCER1A','CLEC10A',
  'CLEC9A', 'XCR1', 'CADM1',
  'LILRA4','TCF4','CLEC4C',
  'PDGFRB','ACTA2','COL1A1',
  'KIT','CPA3')
DotPlot(sce_mergeTEN_harmony,features = marker_cluster,group.by = 'RNA_snn_res.0.4')+
  theme_classic() +
  RotatedAxis() +
  theme(
    axis.title.x  = element_text(size = 16, face = "bold"), 
    axis.title.y  = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14,)
  )
sce_mergeTEN_harmony$celltype <- ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4 %in% c(0,1,5,6), "T cell",
                                        ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4%in%c(2,4,13),"NK cell",
                                               ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==14,"Plasma cell",
                                                      ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==8,"B cell",
                                                             ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==9,"Endothelial cell",
                                                                    ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==15,"Epithelial cell",
                                                                           ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==16,"Hepatocyte",
                                                                                  ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==10,"MoMFs",
                                                                                         ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==7,"Kupffer cell",
                                                                                                ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==3,"Neutrophils",
                                                                                                       ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==11,"DC2",
                                                                                                              ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==12,"DC1",
                                                                                                                     ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==18,"pDC",
                                                                                                                            ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==19,"Hepatic stellate cell",
                                                                                                                                   ifelse(sce_mergeTEN_harmony$RNA_snn_res.0.4==17,"Mast cell","Unknown"
                                                                                                                                   )))))))))))))))
sce_mergeTEN_harmony<-subset(sce_mergeTEN_harmony,celltype!='Unknown')
celltype_order <- c(
  "T cell", 
  "NK cell", 
  "Plasma cell", 
  "B cell", 
  "Endothelial cell", 
  "Epithelial cell", 
  "Hepatocyte", 
  "MoMFs", 
  "Kupffer cell", 
  "Neutrophils", 
  "DC2", 
  "DC1", 
  "pDC", 
  "Hepatic stellate cell", 
  "Mast cell"
)
sce_mergeTEN_harmony$celltype <- factor(sce_mergeTEN_harmony$celltype, levels = celltype_order)
table(sce_mergeTEN_harmony$celltype)
DotPlot(sce_mergeTEN_harmony, features = marker_cluster, group.by  = 'celltype') +
  theme_classic() +
  RotatedAxis() +
  theme(
    axis.title.x  = element_blank(), 
    axis.title.y  = element_blank(),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 12)
  )
DotPlot(sce_mergeTEN_harmony,features = marker_cluster,group.by = 'RNA_snn_res.0.4')+
  theme_classic() +
  RotatedAxis() +
  theme(
    axis.title.x  = element_text(size = 16, face = "bold"), 
    axis.title.y  = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14,)
  )
p1<-DimPlot(sce_mergeTEN_harmony,
            reduction = 'umap',
            group.by = 'RNA_snn_res.0.4',
            label = F)+
  labs(title = 'Cluster')+
  theme(
    legend.text  = element_text(size = 12),
    legend.title  = element_text(size = 14, face = "bold")
  )
p2<-DimPlot(sce_mergeTEN_harmony, 
            reduction = 'umap', 
            group.by  = 'celltype', 
            label = F, 
            label.size  = 4
)+
  theme(
    legend.text  = element_text(size = 12),
    legend.title  = element_text(size = 14, face = "bold")
  )
p1|p2

immune_celltypes <- c(
  "T cell", "NK cell", "Plasma cell", "B cell",
  "MoMFs", "Kupffer cell", "Neutrophils",
  "DC2", "DC1", "pDC", "Mast cell"
)
sce_mergeTEN_harmony <- subset(
  x = sce_mergeTEN_harmony,
  subset = celltype %in% immune_celltypes 
)
sce_mergeTEN_harmony$celltype <- factor(sce_mergeTEN_harmony$celltype, levels = immune_celltypes)
table(sce_mergeTEN_harmony$celltype)
p3<-DimPlot(sce_mergeTEN_harmony, 
            reduction = 'umap', 
            group.by  = 'celltype', 
            label = F, 
            label.size  = 4
)+
  theme(
    legend.text  = element_text(size = 12),
    legend.title  = element_text(size = 14, face = "bold")
  )
p1|p2|p3
#ggsave(filename = 'DotPlot.pdf',plot = p,width = 16,height = 9)
#ggsave(filename = 'celltype.pdf',plot = p,width = 10.5,height = 8)
#ggsave(filename = 'Cluster.tiff',plot = p,width = 9,height = 8)

sce_mergeTEN_harmony_T<-subset(sce_mergeTEN_harmony,celltype=='T cell')
sce_mergeTEN_harmony_T <- NormalizeData(sce_mergeTEN_harmony_T,
                                        normalization.method = "LogNormalize",
                                        scale.factor = 10000)
sce_mergeTEN_harmony_T <- FindVariableFeatures(sce_mergeTEN_harmony_T,
                                               selection.method = "vst", 
                                               nfeatures = 2000)
sce_mergeTEN_harmony_T <- ScaleData(sce_mergeTEN_harmony_T,vars.to.regress = c('percent.mt'))
sce_mergeTEN_harmony_T <- RunPCA(sce_mergeTEN_harmony_T)
data.use <- Stdev(object = sce_mergeTEN_harmony_T, reduction = 'pca')
a=0
b=0
for (i in data.use) {
  a<-a+i
  b=b+1
  if(a/sum(data.use)>=0.85){
    print(b)
    break
  }
}
ElbowPlot(sce_mergeTEN_harmony_T,ndims = 50)
sce_mergeTEN_harmony_T <- FindNeighbors(sce_mergeTEN_harmony_T,reduction = 'harmony',dims = 1:40)
library(harmony)
sce_mergeTEN_harmony_T<-RunHarmony(sce_mergeTEN_harmony_T,"orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce_mergeTEN_harmony_T, 'harmony')
dim(harmony_embeddings)
sce_mergeTEN_harmony_T <- sce_mergeTEN_harmony_T %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) 
sce_mergeTEN_harmony_T<-FindClusters(sce_mergeTEN_harmony_T,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
library(clustree)
clustree(sce_mergeTEN_harmony_T)
Idents(sce_mergeTEN_harmony_T)<-sce_mergeTEN_harmony_T$RNA_snn_res.0.7
DimPlot(sce_mergeTEN_harmony_T,reduction = "umap",label = T)
marker<-c('CD4','CD8A')
DotPlot(sce_mergeTEN_harmony_T,features = marker,group.by = 'RNA_snn_res.0.7')
sce_mergeTEN_harmony_T$celltype2<-ifelse(sce_mergeTEN_harmony_T$RNA_snn_res.0.7=='6'
                                         |sce_mergeTEN_harmony_T$RNA_snn_res.0.7=='7',
                                         'CD4+ T cell','CD8+ T cell')
DotPlot(sce_mergeTEN_harmony_T,features = marker,group.by = 'celltype2')
sce_mergeTEN_harmony$celltype2<-sce_mergeTEN_harmony$celltype
sce_mergeTEN_harmony$celltype2<-as.character(sce_mergeTEN_harmony$celltype2)
T_cell_indices <- sce_mergeTEN_harmony$celltype == 'T cell'
sce_mergeTEN_harmony$celltype2[T_cell_indices] <- sce_mergeTEN_harmony_T$celltype2
table(sce_mergeTEN_harmony$celltype2)
table(sce_mergeTEN_harmony_T$celltype2)

setwd("肝纤维化/肝纤维化转录组/单细胞分析")
lasso<-read.csv("肝纤维化/肝纤维化转录组/机器学习/X/lasso.csv")
lasso$symbol
table(sce_mergeTEN_harmony$celltype)
geneset<-list(ssGSEA_lasso=lasso$symbol)
sce_mergeTEN_harmony<-subset(sce_mergeTEN_harmony,
                             celltype!='Endothelial cell'&
                               celltype!='Epithelial cells'&
                               celltype!='Hepatocyte'&
                               celltype!='Hepatic stellate cell')
expr <- as.matrix(sce_mergeTEN_harmony@assays$RNA@layers$data)
rownames(expr)<-rownames(sce_mergeTEN_harmony)
colnames(expr)<-colnames(sce_mergeTEN_harmony)
param <- ssgseaParam(exprData = expr,
                     geneSets = geneset)
geneset_gsva  <- gsva(param)
geneset_gsva<-as.data.frame(geneset_gsva)
geneset_gsva<-t(geneset_gsva)
geneset_gsva<-as.data.frame(geneset_gsva)
names(geneset_gsva)='ssGSEA_lasso'
table(rownames(geneset_gsva)==rownames(sce_mergeTEN_harmony@meta.data))
sce_mergeTEN_harmony$ssGSEA_lasso=geneset_gsva$ssGSEA_lasso
#write.csv(geneset_gsva,file = 'geneset_gsva.csv')

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
sce_mergeTEN_harmony$ssGSEA_lasso=geneset_gsva$ssGSEA_lasso
plot_data<-sce_mergeTEN_harmony@meta.data
plot_data$celltype <- factor(plot_data$celltype, levels = unique(plot_data$celltype))
plot_data$source <- factor(plot_data$source, levels = unique(plot_data$source))
plot_data$celltype2 <- factor(plot_data$celltype2, levels = unique(plot_data$celltype2))

p <- ggplot(plot_data, 
            aes(x = celltype2, 
                y = ssGSEA_lasso, 
                fill = source)) +
  geom_boxplot(
    width = 0.7,                
    alpha = 0.6,                
    position = position_dodge(0.8),  
    color = "black",            
    outlier.shape   = 16,         
    outlier.size   = 1.5,         
    outlier.alpha   = 0.6         
  ) +  
  stat_compare_means(
    aes(group = source), 
    method = "wilcox.test",   
    label = "p.signif"  
  ) +  
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  
  theme_classic() +
  labs( 
    y = "ssGSEA Score", 
    x = "Cell Type"
  ) +
  theme(
    axis.title.x  = element_blank(),
    axis.title.y  = element_text(size = 12),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1,size = 12),
    axis.text.y   = element_text( size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 15, face = "bold"),
    plot.title  = element_blank()
  ) +
  guides(fill = guide_legend(override.aes  = list(size = 3)))
p
ggsave(filename = paste('ssGSEA_lasso_gene_CD8','tiff',sep = '.'),plot = p,width = 12,height = 9,dpi = 1600)
ggsave(filename = paste('ssGSEA_lasso_gene_CD8','pdf',sep = '.'),plot = p,width = 12,height = 9)

sce_mergeTEN_harmony_T_CD8<-subset(sce_mergeTEN_harmony,celltype2=='CD8+ T cell')
sce_mergeTEN_harmony_T_CD8 <- NormalizeData(sce_mergeTEN_harmony_T_CD8,
                                            normalization.method = "LogNormalize",
                                            scale.factor = 10000)
sce_mergeTEN_harmony_T_CD8 <- FindVariableFeatures(sce_mergeTEN_harmony_T_CD8,
                                                   selection.method = "vst", 
                                                   nfeatures = 2000)
sce_mergeTEN_harmony_T_CD8 <- ScaleData(sce_mergeTEN_harmony_T_CD8,vars.to.regress = c('percent.mt'))
sce_mergeTEN_harmony_T_CD8 <- RunPCA(sce_mergeTEN_harmony_T_CD8)
data.use <- Stdev(object = sce_mergeTEN_harmony_T_CD8, reduction = 'pca')
a=0
b=0
for (i in data.use) {
  a<-a+i
  b=b+1
  if(a/sum(data.use)>=0.85){
    print(b)
    break
  }
}
ElbowPlot(sce_mergeTEN_harmony_T_CD8,ndims = 50)
sce_mergeTEN_harmony_T_CD8 <- FindNeighbors(sce_mergeTEN_harmony_T_CD8,reduction = 'harmony',dims = 1:40)
library(harmony)
sce_mergeTEN_harmony_T_CD8<-RunHarmony(sce_mergeTEN_harmony_T_CD8,"orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce_mergeTEN_harmony_T_CD8, 'harmony')
dim(harmony_embeddings)
sce_mergeTEN_harmony_T_CD8 <- sce_mergeTEN_harmony_T_CD8 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) 
sce_mergeTEN_harmony_T_CD8<-FindClusters(sce_mergeTEN_harmony_T_CD8,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
library(clustree)
clustree(sce_mergeTEN_harmony_T_CD8)
Idents(sce_mergeTEN_harmony_T_CD8)<-sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6
DimPlot(sce_mergeTEN_harmony_T_CD8,reduction = "umap",label = T)
Idents(sce_mergeTEN_harmony_T_CD8)<-sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6
markers<-FindAllMarkers(sce_mergeTEN_harmony_T_CD8,
                        min.pct = 0.25,
                        logfc.threshold = 0.25,
                        only.pos=T
)
gene_list <- markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = ",")) %>%
  deframe() %>%
  as.list()
gene_list
marker<-c(
  "SLC4A10","ZBTB16","NCR3",
  "S1PR1",'GPR183',
  "GNLY","GZMB",'FGFBP2',
  "PDCD1",'TNFRSF9',
  'CD68','CD14'
)
sce_mergeTEN_harmony_T_CD8$celltype3<-ifelse(sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6%in%c(1,2,8,10),'MAIT',
                                             ifelse(sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6=='4','CD8+ TEX',
                                                    ifelse(sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6%in%c(0,3,9),'CD8+ TEF',
                                                           ifelse(sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6=='11','MoMFs',
                                                                  ifelse(sce_mergeTEN_harmony_T_CD8$RNA_snn_res.0.6=='6','CD8+ TCM','UNknow')))))
DotPlot(sce_mergeTEN_harmony_T_CD8,features = marker,group.by ='RNA_snn_res.0.6')
DotPlot(sce_mergeTEN_harmony_T_CD8,features = marker,group.by = 'celltype3')
T_cell_indices <- sce_mergeTEN_harmony$celltype2 == 'CD8+ T cell'
sce_mergeTEN_harmony$celltype3<-sce_mergeTEN_harmony$celltype2
sce_mergeTEN_harmony$celltype3[T_cell_indices] <- sce_mergeTEN_harmony_T_CD8$celltype3
table(sce_mergeTEN_harmony$celltype3)
table(sce_mergeTEN_harmony_T_CD8$celltype3)
sce_mergeTEN_harmony_T_CD8<-subset(sce_mergeTEN_harmony_T_CD8,celltype3!='UNknow'&celltype3!='MoMFs')
levels=c('MAIT','CD8+ TCM','CD8+ TEF','CD8+ TEX')
sce_mergeTEN_harmony_T_CD8$celltype3<-factor(sce_mergeTEN_harmony_T_CD8$celltype3,levels = levels)
marker<-c(
  "SLC4A10","ZBTB16","NCR3",
  "S1PR1",'GPR183',
  "GNLY","GZMB",'FGFBP2',
  "PDCD1",'TNFRSF9'
)
DotPlot(sce_mergeTEN_harmony_T_CD8,features = marker,group.by = 'RNA_snn_res.0.6')
DotPlot(sce_mergeTEN_harmony_T_CD8,features = marker,group.by = 'celltype3')+
  theme_classic() +
  RotatedAxis() +
  theme(
    axis.title.x  = element_blank(), 
    axis.title.y  = element_blank(),
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8)
  )
p1<-DimPlot(sce_mergeTEN_harmony_T_CD8,group.by = 'celltype3')+labs(title = 'CD8+ T cell')+theme(legend.position = 'none')
p1
DimPlot(sce_mergeTEN_harmony_T_CD8,group.by = 'RNA_snn_res.0.6')
FeaturePlot(sce_mergeTEN_harmony_T_CD8, 
            features = 'CCL20', 
            split.by  = 'source') +
  theme(legend.position  = "right") 
p2<-FeaturePlot_scCustom(
  seurat_object = sce_mergeTEN_harmony_T_CD8,
  features = "CCL20",
  split.by  = "source",
  colors_use = viridis::viridis(11)
)
p2
p1<-DimPlot(sce_mergeTEN_harmony_T_CD8,group.by = 'celltype3',
            label = F,split.by = 'source')+
  labs(title = 'celltype')+
  theme(
    legend.text  = element_text(size = 12,
                                face = "bold")
  )
p1

sce_mergeTEN_harmony_T_CD8_Cirrhosis<-subset(sce_mergeTEN_harmony_T_CD8,source=='Cirrhosis')
sce_mergeTEN_harmony_T_CD8_Cirrhosis$celltype_lasso_ssGESA_CD8_Cirrhosis<-ifelse(sce_mergeTEN_harmony_T_CD8_Cirrhosis$ssGSEA_lasso>=median(sce_mergeTEN_harmony_T_CD8_Cirrhosis$ssGSEA_lasso),'high score CD8+ T cell','low score CD8+ T cell')
Idents(sce_mergeTEN_harmony_T_CD8_Cirrhosis)<-sce_mergeTEN_harmony_T_CD8_Cirrhosis$celltype_lasso_ssGESA_CD8_Cirrhosis
markers_source<-FindAllMarkers(sce_mergeTEN_harmony_T_CD8_Cirrhosis, 
                               only.pos = T, 
                               min.pct = 0.25
)
high_markers_source<-markers_source[markers_source$p_val_adj<0.05&
                                      markers_source$cluster=='high score CD8+ T cell'&markers_source$avg_log2FC>0.25,]

erich_go = enrichGO(gene =high_markers_source$gene,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    pAdjustMethod = 'BH',
                    ont = "ALL",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 1)
erich_go_data<-as.data.frame(erich_go)
p1_GO<-barplot(erich_go,showCategory =20,font.size=12,label_format=60, split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale="free",space="free_y")
p2_GO<-clusterProfiler::dotplot(erich_go,showCategory =10,font.size=14,label_format=50,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale="free",space="free_y")+
  ggtitle('GO')+
  theme(plot.title = element_text(face = 'bold',size = 15))
#ggsave(filename = 'GO_barplot_high.tiff',plot =p1_GO,width = 9,height = 10)
#ggsave(filename = 'GO_dotplot_high.tiff',plot =p2_GO,width = 11,height = 10.8)
library(org.Hs.eg.db)
library(AnnotationDbi)
k=keys(org.Hs.eg.db,keytype = "ENTREZID")
list=AnnotationDbi::select(org.Hs.eg.db,keys=k,columns = c("SYMBOL","ENTREZID"), keytype="ENTREZID")
ID_list=list[match(high_markers_source$gene,list[,"SYMBOL"]),]
erich_kegg<-enrichKEGG(ID_list$ENTREZID,
                       organism = 'hsa',
                       pAdjustMethod='BH',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 1)
KEGG<-as.data.frame(erich_kegg)
p1_KEGG<-barplot(erich_kegg,showCategory =20,font.size=15,label_format=20)
p2_KEGG<-clusterProfiler::dotplot(erich_kegg,showCategory =20,font.size=13,label_format=40)+
  ggtitle('KEGG')+
  theme(plot.title = element_text(face = 'bold',size = 14))
#ggsave(filename = 'KEGG_barplot_high.tiff',plot =p1_KEGG,width = 10,height = 9)
#ggsave(filename = 'KEGG_dotplot_high.tiff',plot =p2_KEGG,width = 10,height = 9)

Idents(sce_mergeTEN_harmony_T_CD8)<-sce_mergeTEN_harmony_T_CD8$celltype3
mmt<-sce_mergeTEN_harmony_T_CD8
mmt<-FindVariableFeatures(mmt)
matrix<-as(as.matrix(mmt@assays$RNA$counts),'sparseMatrix')
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix)
)
sample_ann <- mmt@meta.data
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,
                           phenoData = pd,
                           featureData =fd,
                           expressionFamily = negbinomial.size(),
                           lowerDetectionLimit=0.1
)
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sc_cds_2),
                                    num_cells_expressed >= 10))
fData(sc_cds_2)[1:5,]
ordering_genes<-VariableFeatures(mmt)
length(ordering_genes)
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~celltype3",cores = 32,
                                      verbose = T)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
plot_ordering_genes(sc_cds2)
plot_pc_variance_explained(sc_cds2) + geom_vline(xintercept = 6) 
sc_cds2<- reduceDimension(sc_cds2, max_components = 2,
                          num_dim=18,
                          reduction_method  = "DDRTree",
                          residualModelFormulaStr = "~orig.ident",
                          verbose=T
)
sc_cds2 <- orderCells(sc_cds2,root_state = 4)
#save(sc_cds2,file='sc_cds2_dim18_state4.RData')
library(monocle)
load('sc_cds2_dim18_state4.RData')
lasso<-read.csv('肝纤维化/肝纤维化转录组/机器学习/X/lasso.csv')
lasso<-arrange(lasso,desc(lasso$s0))
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T)
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T)+
  facet_wrap(~State)
p11<-plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = T)+
  theme(
    legend.title  = element_text(size = 11),
    legend.text  = element_text(size = 9)
  ) 
p22<-plot_cell_trajectory(sc_cds2, color_by = "celltype3",show_branch_points = T)+
  labs(col='celltype')
p33 <- plot_cell_trajectory(sc_cds2, markers = c("CCL20"), use_color_gradient = TRUE, show_branch_points = FALSE) +
  theme(
    legend.title  = element_text(size = 11),
    legend.text  = element_text(size = 9)
  )
(p11+p22+p33)&
  theme(
    axis.title = element_text(size=20),
    axis.text = element_text(size=18))
#ggsave(filename = 'plot_Pseudotime.tiff',plot = p11,width = 6,height = 6,dpi = 900)
#ggsave(filename = 'plot_celltype.tiff',plot = p22,width = 8,height = 8,dpi = 900)
#ggsave(filename = 'plot_CLL20.tiff',plot = p33,width = 6,height = 6,dpi = 900)
#ggsave(filename = 'plot_Pseudotime.pdf',plot = p11,width = 6,height = 6)
#ggsave(filename = 'plot_celltype.pdf',plot = p22,width = 8,height = 8)
#ggsave(filename = 'plot_CLL20.pdf',plot = p33,width = 6,height = 6)

to_be_tested <- row.names(subset(fData(sc_cds2), gene_short_name %in% lasso$symbol[1:8]))
cds_subset <- sc_cds2[to_be_tested,]
p11 <- plot_genes_jitter(
  cds_subset, 
  grouping = "celltype3", 
  color_by = "celltype3"
) +
  xlab(label = "CD8+ T cell")+
  labs(color = "CD8+ T cell") +
  theme(
    legend.position  ='none',
    axis.text.x  = element_blank()
  ) +
  guides(color = guide_legend(override.aes  = list(size = 3)))
p11
p22 <- plot_genes_violin(cds_subset, grouping = "celltype3", color_by = "celltype3")
p33 <- plot_genes_in_pseudotime(
  cds_subset,
  color_by = "celltype3",
  trend_formula = "~ sm.ns(Pseudotime, df=2)",
  relative_expr = T,
  vertical_jitter = 0.1
)+
  labs(color = "CD8+ T cell") +
  theme(
    legend.text  = element_text( size = 12),  
    legend.title  = element_text( size = 14),
    axis.text.x  = element_blank()
  ) +
  guides(color = guide_legend(override.aes  = list(size = 3)))
plotc1 <- p11|p33
plotc1 

library(ggplot2)
library(patchwork)
p11_modified <- p11 + theme(legend.position  = "none")
combined_plot <- p11_modified + p33 + plot_layout(ncol = 2, guides = "collect")
#ggsave(filename = 'p11.tiff',plot = combined_plot,width = 6,height = 6,dpi = 900)
#ggsave(filename = 'p11.pdf',plot = combined_plot,width = 6,height = 6)

library(TCellSI)
library(DGEobj.utils)
sce_mergeTEN_harmony_T_CD8_MAIT<-subset(sce_mergeTEN_harmony,celltype3=='MAIT')
sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis<-subset(sce_mergeTEN_harmony_T_CD8_MAIT,source=='Cirrhosis')
sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis$score <- ifelse(
  GetAssayData(sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis)["CCL20", ] >= 
    median(GetAssayData(sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis)["CCL20", ]),
  'high CCL20 MAIT', 'low CCL20 MAIT'
)
scRNA<-sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis
sample_scRNA <- as.matrix(scRNA@assays$RNA@layers$counts)
rownames(sample_scRNA)<-rownames(scRNA)
colnames(sample_scRNA)<-colnames(scRNA)

pseudobulk <- AggregateExpression(
  sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis, 
  assays = "RNA", 
  slot = "counts",
  group.by  = c("orig.ident",'score'),
  return.seurat  = FALSE,
  aggregation.fun  = "sum"
)
pseudobulk_matrix<-as.matrix(pseudobulk$RNA)
pseudobulk_matrix_TPM<-generate_scRNA_TPM(txdb, pseudobulk_matrix)
pseudobulk_matrix_TPM_log<-log(pseudobulk_matrix_TPM+1)
batch=c('human_cir1','human_cir1','human_cir2','human_cir2')
group_list<-c('high CCL20 MAIT','low CCL20 MAIT','high CCL20 MAIT','low CCL20 MAIT')
design=model.matrix(~group_list)
pseudobulk_matrix_TPM_log1 <- ComBat(pseudobulk_matrix_TPM_log,
                                     batch=batch,
                                     mod=design,
                                     par.prior=TRUE, prior.plots=FALSE)
scRNA_scores <- TCellSI::TCSS_Calculate(pseudobulk_matrix_TPM_log1)
scRNA_scores <- scRNA_scores %>%
  t() %>%
  as.data.frame()
scRNA_scores$score<-group_list 
names(scRNA_scores)[1]<-'ID'
scRNA_scores_long<-scRNA_scores %>%
  pivot_longer(
    cols = -c(ID,score),
    names_to = "type",
    values_to = "value"
  )

generate_scRNA_TPM <- function(txdb, sample_scRNA) {
  library(DGEobj.utils)
  library(GenomicFeatures)
  library(txdbmaker)
  library(org.Hs.eg.db) 
  library(AnnotationDbi)
  exons.list.per.gene  <- exonsBy(txdb, by = "gene")
  exonic.gene.sizes  <- sum(width(GenomicRanges::reduce(exons.list.per.gene))) 
  gfe <- data.frame( 
    gene_id = names(exonic.gene.sizes), 
    length = exonic.gene.sizes  
  )
  gfe$SYMBOL <- mapIds(
    org.Hs.eg.db,   
    keys = sub("\\.[0-9]+$", "", gfe$gene_id),  
    column = "SYMBOL", 
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  gfe_na <- na.omit(gfe) 
  gfe_na$gene_id <- NULL
  match_genes <- intersect(rownames(sample_scRNA), gfe_na$SYMBOL)
  sample_scRNA_match <- sample_scRNA[match_genes, ]
  gfe_na_matched <- gfe_na[gfe_na$SYMBOL %in% match_genes, ]
  gfe_na_matched <- gfe_na_matched %>% arrange(desc(length)) %>% distinct(SYMBOL, .keep_all = TRUE)
  gfe_na_matched <- gfe_na_matched[order(match(gfe_na_matched$SYMBOL, match_genes)), ]
  sample_scRNA_TPM <- convertCounts(
    counts = sample_scRNA_match,
    unit = "TPM",
    geneLength = gfe_na_matched$length,
    log = FALSE,
    normalize = "none"
  )
  return(sample_scRNA_TPM)
}

txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf.gz", format = "gtf")
sample_scRNA_TPM<-generate_scRNA_TPM(txdb, sample_scRNA)
library(DGEobj.utils)
library(GenomicFeatures)
library(txdbmaker)
library(org.Hs.eg.db) 
library(AnnotationDbi)
sample_scRNA_TPM1 <- log(sample_scRNA_TPM + 1)
sample_scRNA[1:5, 1:5]
batch=sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis$orig.ident
group_list<-sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis$score
design=model.matrix(~group_list)
sample_scRNA_TPM2  <- ComBat(sample_scRNA_TPM1 ,
                             batch=batch,
                             mod=design,
                             par.prior=TRUE, prior.plots=FALSE)
scRNA_scores <- TCellSI::TCSS_scRNAseqCalculate(sample_scRNA_TPM1, core= 7)
scRNA_scores <- scRNA_scores %>%
  t() %>%
  as.data.frame()
scRNA_scores$score<-sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis$score
write.csv(scRNA_scores,file = 'scRNA_scores_MAIT_TPM1.csv')

sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis<-subset(sce_mergeTEN_harmony_T_CD8_MAIT,source=='Cirrhosis')
ccl20_expression <- GetAssayData(
  sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis, 
  assay = "RNA", 
  slot = "data")["CCL20", ]
sce_mergeTEN_harmony_T_CD8_MAIT_cirrhosis$score<-ifelse(ccl20_expression >= median(ccl20_expression),
                                                        "high CCL20 MAIT",
                                                        "low CCL20 MAIT"
)
scRNA_scores$score<-scRNA$score
names(scRNA_scores)[1]<-'ID'
sce_mergeTEN_harmony_T_CD8_MAIT<-subset(sce_mergeTEN_harmony,celltype3=='MAIT')
scRNA_scores_long<-scRNA_scores %>%
  pivot_longer(
    cols = -c(ID,score),
    names_to = "type",
    values_to = "value"
  )
scRNA_scores_long$type<-factor(scRNA_scores_long$type,levels=c('Quiescence','Regulating','Proliferation','Helper','Cytotoxicity','Progenitor_exhaustion','Terminal_exhaustion','Senescence'))

scRNA@meta.data$Quiescence <- scRNA_scores$Quiescence
scRNA@meta.data$Regulating <- scRNA_scores$Regulating
scRNA@meta.data$Proliferation <- scRNA_scores$Proliferation
scRNA@meta.data$Helper <- scRNA_scores$Helper
scRNA@meta.data$Cytotoxicity <- scRNA_scores$Cytotoxicity
scRNA@meta.data$Progenitor_exhaustion <- scRNA_scores$Progenitor_exhaustion
scRNA@meta.data$Terminal_exhaustion <- scRNA_scores$Terminal_exhaustion
scRNA@meta.data$Senescence <- scRNA_scores$Senescence
VlnPlot(scRNA, features = "Quiescence")
FeaturePlot(object = scRNA, features = "Quiescence")  
RidgePlot(scRNA, features = "Quiescence", ncol = 1)

library(ggplot2)
library(ggpubr)
library(dplyr)
ggplot(scRNA_scores_long, aes(x = type, y = value, fill = score)) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(0.75),
    outlier.shape  = 21,
    outlier.size  = 2.5,
    outlier.stroke  = 0.4,
    outlier.color  = "black",
    alpha = 0.8
  ) +
  stat_compare_means(
    aes(group = score),
    method = "wilcox.test", 
    label = "p.signif", 
    show.legend  = FALSE,
    size = 4.5,
    label.y = max(scRNA_scores_long$value) * 1.12,
    vjust = 0.5
  ) +  
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  
  theme_classic() +
  labs( 
    y = "Score"
  ) +
  theme(
    axis.title.x  = element_blank(),
    axis.title.y  = element_text(size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1,size = 10),
    axis.text.y   = element_text( size = 10),
    legend.title = element_blank(),
    legend.text  = element_text(size = 10),
    plot.title  = element_blank()
  ) +
  guides(fill = guide_legend(override.aes  = list(size = 3)))
#ggsave(filename = 'low_and_high_MAIT.png',width =12 ,height = 6,dpi = 900)

sce_mergeTEN_harmony_Cirrhotic<-subset(sce_mergeTEN_harmony,
                                       source=='Cirrhosis')
Idents(sce_mergeTEN_harmony_Cirrhotic)<-sce_mergeTEN_harmony_Cirrhotic$celltype3
library(CellChat)
sce_mergeTEN_harmony_Cirrhotic_input <- GetAssayData(sce_mergeTEN_harmony_Cirrhotic, assay = "RNA", slot = "data")
Idents(sce_mergeTEN_harmony_Cirrhotic)<-sce_mergeTEN_harmony_Cirrhotic$celltype2
labels <- factor(sce_mergeTEN_harmony_Cirrhotic$celltype3,levels=levels(Idents(sce_mergeTEN_harmony_Cirrhotic)))
labels
meta <- data.frame(group = labels, row.names = rownames(sce_mergeTEN_harmony_Cirrhotic@meta.data))
Cirrhotic_cellchat <- createCellChat(object = sce_mergeTEN_harmony_Cirrhotic_input, meta = meta, group.by = "group")
Cirrhotic_cellchat@DB <- CellChatDB.human
gc()
Cirrhotic_cellchat <- subsetData(Cirrhotic_cellchat) 
Cirrhotic_cellchat <- identifyOverExpressedGenes(Cirrhotic_cellchat)
Cirrhotic_cellchat <- identifyOverExpressedInteractions(Cirrhotic_cellchat)
Cirrhotic_cellchat <- projectData(Cirrhotic_cellchat, PPI.human)
Cirrhotic_cellchat <- computeCommunProb(Cirrhotic_cellchat, raw.use = TRUE)
Cirrhotic_cellchat <- filterCommunication(Cirrhotic_cellchat, min.cells = 10)
Cirrhotic_cellchat <- computeCommunProbPathway(Cirrhotic_cellchat)
Cirrhotic_cellchat <- netAnalysis_computeCentrality(Cirrhotic_cellchat, slot.name = "netP")
Cirrhotic_cellchat <- aggregateNet(Cirrhotic_cellchat)

groupSize.Cirrhotic <- as.numeric(table(Cirrhotic_cellchat@idents))
groupSize.Cirrhotic
table(Cirrhotic_cellchat@idents)
par(mfrow = c(1,2), xpd=TRUE)
#save(Cirrhotic_cellchat,file='Cirrhotic_cellchat_MAIT.RData')

setwd("肝纤维化/肝纤维化转录组/单细胞分析/细胞通讯IRG/更新/")
load('Cirrhotic_cellchat_MAIT.RData')
cellchat<-Cirrhotic_cellchat
df<-subsetCommunication(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

par(mfrow=c(2,2))
par(oma = c(0, 0, 0, 0), mar = c(0, 0., 0, 0))
netVisual_circle(cellchat@net$count,arrow.size = 0.005, vertex.weight = groupSize, weight.scale = T, label.edge= F,vertex.label.cex = 1,
                 targets.use =c(3) ,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, arrow.size = 0.005,vertex.weight = groupSize, weight.scale = T, label.edge= F,vertex.label.cex = 1,
                 targets.use = c(3),
                 title.name = "Interaction weights/strength")

bubble.sender<-netVisual_bubble(cellchat, 
                                sources.use = c(3),
                                targets.use = c(1:15),
                                remove.isolate = FALSE,title.name = 'Sender',
                                signaling = c('CCL'),
                                angle.x = 45
)+ theme(
  plot.title = element_text(size = 14,hjust = 0.5),
  axis.text.x  = element_text(size = 12),
  axis.text.y  = element_text(size = 12),
  legend.text  = element_text(size = 11),  
  legend.title  = element_text( size = 12)
)
bubble.sender  
bubble.reciver<-netVisual_bubble(cellchat, 
                                 sources.use = c(1:2,4:15), 
                                 targets.use = c(3), 
                                 remove.isolate = FALSE,
                                 title.name = 'Reciver',
                                 signaling = c('CCL'),
                                 angle.x = 45
)+ theme(
  plot.title = element_text(size = 14,hjust = 0.5),
  axis.text.x  = element_text(size = 12),
  axis.text.y  = element_text(size = 12),
  legend.text  = element_text(size = 11),  
  legend.title  = element_text( size = 12)
)
bubble.sender
bubble.reciver
#ggsave(filename = 'sender.png',plot = bubble.sender,height = 4,width = 8)
#ggsave(filename = 'reciver.png',plot = bubble.reciver,height = 4,width = 8)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height=16,font.size=11,font.size.title = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height=16,font.size=11,font.size.title = 20)
ht1 + ht2

netAnalysis_signalingRole_network(cellchat,
                                  signaling = c("CCL"),
                                  width = 8, height = 4, font.size = 8)
#pdf(file = 'signalingRole_network.pdf',height = 3,width = 6)
#dev.off()

netVisual_heatmap(cellchat, 
                  signaling = 'CCL',
                  color.heatmap = 'Reds',
                  title.name = 'CCL Pathway'
)
#pdf(file = 'heatmap.pdf',width = 8,height = 8)
#dev.off()

plotGeneExpression(cellchat, signaling = "CCL"
)
#pdf(file = 'plotGeneExpression.pdf',width = 5,height = 5)
#dev.off()

sce_mergeTEN_harmony_MoMFs<-subset(sce_mergeTEN_harmony,celltype3=='MoMFs')
sce_mergeTEN_harmony_MoMFs <- NormalizeData(sce_mergeTEN_harmony_MoMFs,
                                            normalization.method = "LogNormalize",
                                            scale.factor = 10000)
sce_mergeTEN_harmony_MoMFs <- FindVariableFeatures(sce_mergeTEN_harmony_MoMFs,
                                                   selection.method = "vst", 
                                                   nfeatures = 2000)
sce_mergeTEN_harmony_MoMFs <- ScaleData(sce_mergeTEN_harmony_MoMFs,vars.to.regress = c('percent.mt'))
sce_mergeTEN_harmony_MoMFs <- RunPCA(sce_mergeTEN_harmony_MoMFs)
data.use <- Stdev(object = sce_mergeTEN_harmony_MoMFs, reduction = 'pca')
a=0
b=0
for (i in data.use) {
  a<-a+i
  b=b+1
  if(a/sum(data.use)>=0.85){
    print(b)
    break
  }
}
ElbowPlot(sce_mergeTEN_harmony_MoMFs,ndims = 50)
sce_mergeTEN_harmony_MoMFs <- FindNeighbors(sce_mergeTEN_harmony_MoMFs,reduction = 'harmony',dims = 1:40)
library(harmony)
sce_mergeTEN_harmony_MoMFs<-RunHarmony(sce_mergeTEN_harmony_MoMFs,"orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce_mergeTEN_harmony_MoMFs, 'harmony')
dim(harmony_embeddings)
sce_mergeTEN_harmony_MoMFs <- sce_mergeTEN_harmony_MoMFs %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) 
sce_mergeTEN_harmony_MoMFs<-FindClusters(sce_mergeTEN_harmony_MoMFs,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
library(clustree)
clustree(sce_mergeTEN_harmony_MoMFs)
Idents(sce_mergeTEN_harmony_MoMFs)<-sce_mergeTEN_harmony_MoMFs$RNA_snn_res.1
DimPlot(sce_mergeTEN_harmony_MoMFs,reduction = "umap",label = T)
marker<-c('CD86','TNF',
          'CD163','TGFB1')
DotPlot(sce_mergeTEN_harmony_MoMFs,features = marker,group.by = 'RNA_snn_res.0.3')
sce_mergeTEN_harmony_MoMFs$celltype4<-ifelse(sce_mergeTEN_harmony_MoMFs$RNA_snn_res.0.1=='1',
                                             'M1','M2')
DotPlot(sce_mergeTEN_harmony_MoMFs,features = marker,group.by = 'celltype4')
table(sce_mergeTEN_harmony_MoMFs$celltype4)
sce_mergeTEN_harmony$celltype4<-sce_mergeTEN_harmony$celltype3
T_cell_indices <- sce_mergeTEN_harmony$celltype3 == 'MoMFs'
sce_mergeTEN_harmony$celltype4[T_cell_indices] <- sce_mergeTEN_harmony_MoMFs$celltype4
table(sce_mergeTEN_harmony$celltype4)

sce_mergeTEN_harmony_Cirrhotic<-subset(sce_mergeTEN_harmony,
                                       source=='Cirrhosis')
sce_mergeTEN_harmony_Healthy<-subset(sce_mergeTEN_harmony,
                                     source=='Healthy')
sce_mergeTEN_harmony_Cirrhotic<-subset(sce_mergeTEN_harmony_Cirrhotic,
                                       celltype4=='MAIT'|celltype4=='M1'|celltype4=='M2')
sce_mergeTEN_harmony_Healthy<-subset(sce_mergeTEN_harmony_Healthy,
                                     celltype4=='MAIT'|celltype4=='M1'|celltype4=='M2')
library(CellChat)
sce_mergeTEN_harmony_Cirrhotic_input <- GetAssayData(sce_mergeTEN_harmony_Cirrhotic, assay = "RNA", slot = "data")
Idents(sce_mergeTEN_harmony_Cirrhotic)<-sce_mergeTEN_harmony_Cirrhotic$celltype4
labels <- factor(sce_mergeTEN_harmony_Cirrhotic$celltype4,levels=levels(Idents(sce_mergeTEN_harmony_Cirrhotic)))
labels
meta <- data.frame(group = labels, row.names = rownames(sce_mergeTEN_harmony_Cirrhotic@meta.data))
Cirrhotic_cellchat <- createCellChat(object = sce_mergeTEN_harmony_Cirrhotic_input, meta = meta, group.by = "group")
Cirrhotic_cellchat@DB <- CellChatDB.human
gc()
Cirrhotic_cellchat <- subsetData(Cirrhotic_cellchat) 
Cirrhotic_cellchat <- identifyOverExpressedGenes(Cirrhotic_cellchat)
Cirrhotic_cellchat <- identifyOverExpressedInteractions(Cirrhotic_cellchat)
Cirrhotic_cellchat <- projectData(Cirrhotic_cellchat, PPI.human)
Cirrhotic_cellchat <- computeCommunProb(Cirrhotic_cellchat, raw.use = TRUE)
Cirrhotic_cellchat <- filterCommunication(Cirrhotic_cellchat, min.cells = 10)
Cirrhotic_cellchat <- computeCommunProbPathway(Cirrhotic_cellchat)
Cirrhotic_cellchat <- netAnalysis_computeCentrality(Cirrhotic_cellchat, slot.name = "netP")
Cirrhotic_cellchat <- aggregateNet(Cirrhotic_cellchat)

groupSize.Cirrhotic <- as.numeric(table(Cirrhotic_cellchat@idents))
groupSize.Cirrhotic
table(Cirrhotic_cellchat@idents)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Cirrhotic_cellchat@net$count,sources.use = c(1),
                 arrow.size = 0.01, vertex.weight = groupSize.Cirrhotic, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Cirrhotic_cellchat@net$weight,sources.use = c(1), 
                 arrow.size = 0.01,vertex.weight = groupSize.Cirrhotic, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

vertex.receiver=c(2,3)
pairLR.CCL <- extractEnrichedLR(Cirrhotic_cellchat, signaling = "CCL", geneLR.return  = FALSE)
LR.show  <- pairLR.CCL[4,]
netVisual_individual(Cirrhotic_cellchat, signaling = c('CCL'), 
                     vertex.receiver = c(1),
                     pairLR.use = LR.show,
                     arrow.size=0.1,
                     layout = 'hierarchy')
netVisual_individual(Cirrhotic_cellchat, signaling = c('CCL'),
                     layout = 'circle',
                     pairLR.use = LR.show,
                     arrow.size = 0.01,
                     label.edge = F)
netVisual_individual(Cirrhotic_cellchat, signaling = c('CCL'),
                     layout = 'chord',
                     pairLR.use = LR.show,
                     arrow.size = 0.01,
                     vertex.label.cex = 1,
                     label.edge = TRUE  ,
                     title.cex=2)
bubble.sender<-netVisual_bubble(Cirrhotic_cellchat, 
                                sources.use = c(1),
                                targets.use = c(1:3),
                                remove.isolate = FALSE,title.name = 'Sender'
                                ,signaling = c('CCL'),angle.x = 45
)+ theme(
  plot.title = element_text(size = 16,hjust = 0.5),
  axis.text.x  = element_text(size = 14),
  axis.text.y  = element_text(size = 14),
  legend.text  = element_text( size = 13),  
  legend.title  = element_text( size = 14),
  legend.position  = "right"
)
bubble.sender 
#ggsave('bubble_sender.png',width = 5,height = 4,dpi = 900,plot = bubble.sender)
#ggsave('bubble_sender.pdf',width = 5,height = 4,plot = bubble.sender)
bubble.reciver<-netVisual_bubble(Cirrhotic_cellchat, 
                                 sources.use = c(2:3), 
                                 targets.use = c(1), 
                                 remove.isolate = FALSE,
                                 title.name = 'Reciver'
                                 ,signaling = c('CCL'),angle.x = 45
)+ theme(
  plot.title = element_text(size = 16,hjust = 0.5),
  axis.text.x  = element_text(size = 14),
  axis.text.y  = element_text(size = 14),
  legend.text  = element_text(size = 13),  
  legend.title  = element_text(size = 14)
)
bubble.reciver<-bubble.reciver+ guides(
  colour = guide_colourbar(title = "Commun. Prob.", order = 1),
  size = guide_legend(title = "p-value", order = 2)
)
#ggsave('bubble_reciver.png',width = 5,height = 4,dpi = 900,plot = bubble.reciver)
#ggsave('bubble_reciver.pdf',width = 5,height = 4,plot = bubble.reciver)
netVisual_heatmap(Cirrhotic_cellchat, 
                  signaling = 'CCL',
                  color.heatmap = 'Reds',
                  title.name = 'CCL Pathway')
plotGeneExpression(Cirrhotic_cellchat, signaling = "CCL"
                   ,type = 'violin'
)

sce_mergeTEN_harmony_Healthy
sce_mergeTEN_harmony_Healthy_input <- GetAssayData(sce_mergeTEN_harmony_Healthy, assay = "RNA", slot = "data")
Idents(sce_mergeTEN_harmony_Healthy)<-sce_mergeTEN_harmony_Healthy$celltype4
labels <- factor(sce_mergeTEN_harmony_Healthy$celltype4,levels=levels(Idents(sce_mergeTEN_harmony_Healthy)))
labels
meta <- data.frame(group = labels, row.names = rownames(sce_mergeTEN_harmony_Healthy@meta.data))
Healthy_cellchat <- createCellChat(object = sce_mergeTEN_harmony_Healthy_input, meta = meta, group.by = "group")
Healthy_cellchat@DB <- CellChatDB.human
gc()
Healthy_cellchat <- subsetData(Healthy_cellchat) 
Healthy_cellchat <- identifyOverExpressedGenes(Healthy_cellchat)
Healthy_cellchat <- identifyOverExpressedInteractions(Healthy_cellchat)
Healthy_cellchat <- projectData(Healthy_cellchat, PPI.human)
Healthy_cellchat <- computeCommunProb(Healthy_cellchat, raw.use = TRUE)
Healthy_cellchat <- filterCommunication(Healthy_cellchat, min.cells = 10)
Healthy_cellchat <- computeCommunProbPathway(Healthy_cellchat)
Healthy_cellchat <- netAnalysis_computeCentrality(Healthy_cellchat, slot.name = "netP")
Healthy_cellchat <- aggregateNet(Healthy_cellchat)

groupSize.Healthy <- as.numeric(table(Healthy_cellchat@idents))
groupSize.Healthy
table(Healthy_cellchat@idents)
#save(Cirrhotic_cellchat,file='Cirrhotic_cellchat_MAIT_M1.RData')
#save(Healthy_cellchat,file='Healthy_cellchat_MAIT_M1.RData')

load('Cirrhotic_cellchat_MAIT_M1.RData')
load('Healthy_cellchat_MAIT_M1.RData')
object.list <- list(Healthy = Healthy_cellchat, Cirrhotic = Cirrhotic_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T ,arrow.size = 0.1)
netVisual_diffInteraction(cellchat, weight.scale = T, arrow.size = 0.1,measure = "weight")

gg1 <- netVisual_heatmap(cellchat,signaling = c('CCL'),title.name = 'CCL Pathway')
gg2 <- netVisual_heatmap(cellchat, measure = "weight",signaling = c('CCL'))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,signaling = 'CCL')
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,signaling = 'CCL')
combined_plot <- gg1 + gg2 + 
  plot_layout(guides = "collect")
final_plot <- combined_plot &
  theme(
    axis.title.x   = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 14,face = 'bold'),
    axis.text.y   = element_text(size = 18,face = 'bold'),
    legend.title  = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14,face = 'bold'),
    axis.line  = element_line(linewidth = 1),
    legend.key.size  = unit(1, "lines")
  )
final_plot
#ggsave(filename = 'final_plot.png',width = 8,height = 6,dpi = 900)

library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6,font.size=14,font.size.title = 11.5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6,font.size=14,font.size.title = 11.5)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(16, "cm"),width=unit(16, "cm"))
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu",font.size=14,font.size.title = 11.5)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu",font.size=14,font.size.title = 11.5)
draw(ht1 + ht2+ht3+ht4, ht_gap = unit(0.5, "cm"),height=unit(16, "cm"),width=unit(16, "cm"))
#pdf(file = "outgoing_incoming_patterns.pdf",  width = 20, height = 10)
#draw(ht1 + ht2+ht3+ht4, ht_gap = unit(0.7, "cm"),height=unit(16, "cm"))
#dev.off() 
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(16, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(16, "cm"))
par(mfrow = c(1,1), xpd=TRUE)

vertex.receiver=c(2,3)
pairLR.CCL <- extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return  = FALSE)
LR.show  <- pairLR.CCL[4,]
netVisual_individual(cellchat, signaling = c('CCL'),
                     layout = 'chord',
                     pairLR.use = LR.show,
                     arrow.size = 0.01,
                     vertex.label.cex = 2,
                     vertex.label.font = 2,
                     label.edge = TRUE  ,
                     title.cex=2)
par(mfrow = c(1,1), xpd=TRUE)
#pdf(file = paste('CCL_', names(object.list)[i+1],'_circle','.pdf',sep = ''),width = 6,height = 6)
#png(filename  = paste('CCL_', names(object.list)[i+1],'_chord','.png',sep = ''),width = 2000,height = 2000,res = 300)
#dev.off()
p1 <- netVisual_aggregate(object.list[[i]], signaling = 'CCL',pairLR.use=LR.show,arrow.size = 0.01,vertex.label.cex = 1.2,
                          layout = "chord", signaling.name = paste('CCL', names(object.list)[i]))
p2 <- netVisual_aggregate(object.list[[i+1]], signaling = 'CCL',pairLR.use=LR.show,arrow.size = 0.01,vertex.label.cex = 1.2,
                          layout = "chord", signaling.name = paste('CCL', names(object.list)[i+1]))
p1 <- netVisual_aggregate(object.list[[i]], signaling = 'CCL',vertex.label.cex   = 2,
                          layout = "circle", signaling.name = paste('CCL', names(object.list)[i]))
p2 <- netVisual_aggregate(object.list[[i+1]], signaling = 'CCL',vertex.label.cex   = 2,
                          layout = "circle", signaling.name = paste('CCL', names(object.list)[i+1]))
p1+p2
cellchat@idents

bubble_sender<-netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(1:3), 
                                signaling = c("CCL"), 
                                title.name = 'Sender',
                                comparison = c(1, 2), angle.x = 45)+ theme(
                                  plot.title = element_text(size = 16,hjust = 0.5),
                                  axis.text.x  = element_text(size = 14),
                                  axis.text.y  = element_text(size = 14),
                                  legend.text  = element_text(size = 13),  
                                  legend.title  = element_text(size = 14)
                                )+ guides(
                                  colour = guide_colourbar(title = "Commun. Prob.", order = 1),
                                  size = guide_legend(title = "p-value", order = 2)
                                )
bubble_reciver<-netVisual_bubble(cellchat, sources.use = c(2:3), targets.use = c(1), 
                                 signaling = c("CCL"), 
                                 title.name = 'Reciver',
                                 comparison = c(1, 2), angle.x = 45)+ theme(
                                   plot.title = element_text(size = 16,hjust = 0.5),
                                   axis.text.x  = element_text(size = 14),
                                   axis.text.y  = element_text(size = 14),
                                   legend.text  = element_text(size = 13),  
                                   legend.title  = element_text(size = 14)
                                 )+ guides(
                                   colour = guide_colourbar(title = "Commun. Prob.", order = 1),
                                   size = guide_legend(title = "p-value", order = 2)
                                 )
#ggsave('bubble_sender_comparison.png',width = 5,height = 4,dpi = 900,plot = bubble_sender)
#ggsave('bubble_sender_comparison.pdf',width = 5,height = 4,plot = bubble_sender)
#ggsave('bubble_reciver_comparison.png',width = 5,height = 4,dpi = 900,plot = bubble_reciver)
#ggsave('bubble_reciver_comparison.pdf',width = 5,height = 4,plot = bubble_reciver)

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Healthy", "Cirrhotic"))
#png(filename = 'plotGeneExpression_comparison.png',width = 3000,height = 3000,res=900)
pdf(file = 'plotGeneExpression_comparison.pdf',width = 5,height = 5)
plotGeneExpression(cellchat, signaling = "CCL", split.by = "datasets", colors.ggplot = T)
dev.off()
#save(sce_mergeTEN_harmony,file='sce_mergeTEN_harmony_Immunity.RData')

library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
set.seed(1234)
target_gene="CCL20"
load("sce_mergeTEN_harmony_Immunity.RData")
sce_mergeTEN_harmony <- subset(sce_mergeTEN_harmony, source == "Cirrhosis")
sce_mergeTEN_harmony <- subset(sce_mergeTEN_harmony, celltype3 == "MAIT")
countMat = GetAssayData(sce_mergeTEN_harmony, layer = "counts")
sce_mergeTEN_harmony <- FindVariableFeatures(object=sce_mergeTEN_harmony, selection.method="vst", nfeatures=3000)
hvgs <- VariableFeatures(sce_mergeTEN_harmony)
data=as.data.frame(countMat[unique(c(target_gene,hvgs)),])
rm(sce_mergeTEN_harmony, countMat);gc()
result <- scTenifoldKnk(countMatrix = data, 
                        gKO = target_gene,
                        qc_mtThreshold = 0.1,
                        qc_minLSize = 500,
                        nc_nNet = 10,
                        nc_nCells = 500,
                        nCores = 64 
)