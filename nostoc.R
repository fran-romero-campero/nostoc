library(FactoMineR)
library(factoextra)
library(ballgown)
library(VennDiagram)

experimental.design <- read.csv("experimental_design.csv",as.is=T)
experimental.design

bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)
bg.data
sampleNames(bg.data)

gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)

colnames(gene.expression) <- c("rep1_0","rep1_3","rep1_9","rep1_24",
                               "rep2_0","rep2_3","rep2_9","rep2_24",
                               "rep3_0","rep3_3","rep3_9","rep3_24")


gene.names.correspondence <- read.table(file="gene_names_correspondence.tsv",sep="\t",header=F,as.is=T)
gene.names <- gene.names.correspondence$V2
names(gene.names) <- gene.names.correspondence$V1

rownames(gene.expression)[which(is.na(gene.names[rownames(gene.expression)]))]

df.gene.expression <- data.frame(gene_id=rownames(gene.expression),gene_name=gene.names[rownames(gene.expression)],gene.expression,stringsAsFactors = F)
df.gene.expression$gene_name[which(df.gene.expression$gene_name == "")] <- df.gene.expression$gene_id[which(df.gene.expression$gene_name == "")]

df.gene.expression <- df.gene.expression[!is.na(df.gene.expression$gene_name),]

write.table(x=df.gene.expression,file = "gene_expression_complete.tsv",sep = "\t",quote = F,row.names = F)
df.gene.expression <- read.table(file = "gene_expression_complete.tsv",sep="\t",as.is=T,header=T)


gene.expression <- as.matrix(df.gene.expression[,3:ncol(df.gene.expression)])
is.matrix(gene.expression)
head(gene.expression)
rownames(gene.expression) <- df.gene.expression$gene_name

h0 <- (df.gene.expression$rep1_0 + df.gene.expression$rep2_0 + df.gene.expression$rep3_0)/3
h3 <- (df.gene.expression$rep1_3 + df.gene.expression$rep2_3 + df.gene.expression$rep3_3)/3
h9 <- (df.gene.expression$rep1_9 + df.gene.expression$rep2_9 + df.gene.expression$rep3_9)/3
h24 <- (df.gene.expression$rep1_24 + df.gene.expression$rep2_24 + df.gene.expression$rep3_24)/3

mean.gene.expression <- matrix(c(h0,h3,h9,h24),ncol=4)
colnames(mean.gene.expression) <- c("h0","h3","h9","h24")
rownames(mean.gene.expression) <- df.gene.expression$gene_name

pca.gene.expression <- data.frame(colnames(gene.expression),t(gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"

res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70),main = "")

## Hierarchical clustering
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 3)    

png(filename = "figures/hierarchical_tree.png",width = 800,height = 400)
fviz_dend(res.hcpc,k=4,
          cex = 1,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 400      # Augment the room for labels
)
dev.off()

nrow(gene.expression)
gene.expression <- gene.expression[apply(X = gene.expression,MARGIN = 1,FUN = sum) != 0,]

nrow(gene.expression)

log2.gene.expression <- log2(gene.expression + 1)


replicate.scatter.plot <- function(rep1,rep2,gexpress)
{
  plot(x=gexpress[,rep1],y = gexpress[,rep2],pch=19,cex=0.8,col="darkgrey",xlab=paste(rep1,"h"),ylab=paste(rep2,"h"),cex.lab=1.5)
  lines(x=c(-1,20),y=c(-1,20),lwd=2,col="red")
  text(x = 4,y=14,labels = paste("cor =",paste(round(100*cor(gexpress[,rep1],gexpress[,rep2]),digits = 2),"%")),cex = 1.2)
}

replicate.scatter.plot(rep1 = "rep1_0", rep2 = "rep2_0", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep1_0", rep2 = "rep3_0", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep2_0", rep2 = "rep3_0", gexpress = log2.gene.expression)

replicate.scatter.plot(rep1 = "rep1_3", rep2 = "rep2_3", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep1_3", rep2 = "rep3_3", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep2_3", rep2 = "rep3_3", gexpress = log2.gene.expression)

replicate.scatter.plot(rep1 = "rep1_9", rep2 = "rep2_9", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep1_9", rep2 = "rep3_9", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep2_9", rep2 = "rep3_9", gexpress = log2.gene.expression)

replicate.scatter.plot(rep1 = "rep1_24", rep2 = "rep2_24", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep1_24", rep2 = "rep3_24", gexpress = log2.gene.expression)
replicate.scatter.plot(rep1 = "rep2_24", rep2 = "rep3_24", gexpress = log2.gene.expression)

## Differential expression analysis
library(limma)

## Specification of the experimental design
factor.experimental.design <- c(1,2,3,4,1,2,3,4,1,2,3,4)
limma.experimental.design <- model.matrix(~ -1+factor(factor.experimental.design))
colnames(limma.experimental.design) <- c("h0","h3","h9","h24")

## Linear model fit
linear.fit <- lmFit(log2.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(h3-h0, h9-h0, h24-h0,
                                 h9-h3, h24-h3,
                                 h24-h9,
                                 levels=c("h0","h3","h9","h24"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract results
fc.threshold <- 2
q.val.threshold <- 0.05
final.degs <- c()
for(i in 1:6)
{
  degs.results.i <- topTable(contrast.results, number=nrow(log2.gene.expression),coef=i,sort.by="logFC")
  head(degs.results.i)
  
  fold.change.i <- degs.results.i$logFC
  p.values.i <- degs.results.i$P.Value
  q.values.i <- degs.results.i$adj.P.Val
  genes.ids.i <- rownames(degs.results.i)
  
  activated.genes.i <- genes.ids.i[fold.change.i > log2(fc.threshold) & q.values.i < q.val.threshold]
  repressed.genes.i <- genes.ids.i[fold.change.i < - log2(fc.threshold) & q.values.i < q.val.threshold]
  
  final.degs <- c(final.degs,activated.genes.i,repressed.genes.i)
}

final.degs <- unique(final.degs)
length(final.degs)

degs.expression <- gene.expression[final.degs,]

library(WGCNA)
library(cluster)
allowWGCNAThreads()

gene.correlation <- cor(t(degs.expression))

distance.matrix <- 1 - gene.correlation

hierarchical.clustering <- hclust(as.dist(distance.matrix),method="average")

hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)

pam.2 <- pam(as.dist(distance.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(distance.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(distance.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(distance.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(distance.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(distance.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(distance.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(distance.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(distance.matrix),k=10,diss=TRUE)

sil2 <- silhouette(hclust.2,dist=distance.matrix)
sil3 <- silhouette(hclust.3,dist=distance.matrix)
sil4 <- silhouette(hclust.4,dist=distance.matrix)
sil5 <- silhouette(hclust.5,dist=distance.matrix)
sil6 <- silhouette(hclust.6,dist=distance.matrix)
sil7 <- silhouette(hclust.7,dist=distance.matrix)
sil8 <- silhouette(hclust.8,dist=distance.matrix)
sil9 <- silhouette(hclust.9,dist=distance.matrix)
sil10 <- silhouette(hclust.10,dist=distance.matrix)

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

plot(2:10,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.3,0.8),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:10,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)


plot(2:6,pam.sil.values[1:5],type="o",col="blue",pch=0,ylim=c(0.3,0.8),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:6,hclust.sil.values[1:5],type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)

## Three clustres
cluster1.pam3 <- names(which(pam.3[["clustering"]] == 1))
length(cluster1.pam3)
cluster2.pam3 <- names(which(pam.3[["clustering"]] == 2))
length(cluster2.pam3)
cluster3.pam3 <- names(which(pam.3[["clustering"]] == 3))
length(cluster3.pam3)

expr.cluster1.pam3 <- mean.gene.expression[cluster1.pam3,]
expr.cluster2.pam3 <- mean.gene.expression[cluster2.pam3,]
expr.cluster3.pam3 <- mean.gene.expression[cluster3.pam3,]

mean.profile.cluster1.pam3 <- colMeans(expr.cluster1.pam3)
mean.profile.cluster2.pam3 <- colMeans(expr.cluster2.pam3)
mean.profile.cluster3.pam3 <- colMeans(expr.cluster3.pam3)

plot(x = c(0,3,9,24),mean.profile.cluster1.pam3,type="o",col="blue",xlab="Time (h)",ylab="Expression Level (FPKM)",lwd=3,pch=0,main="Cluster Profiles",ylim=c(0,800))
lines(x = c(0,3,9,24),mean.profile.cluster2.pam3,type="o",col="red",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster3.pam3,type="o",col="darkgreen",lwd=3,pch=1)

write.table(x = cluster1.pam3, file = "cluster_azul_respuesta_rapida.txt",quote = F,row.names = F,col.names = F)
write.table(x = cluster2.pam3, file = "cluster_red_respuesta_progresiva.txt",quote = F,row.names = F,col.names = F)
write.table(x = cluster3.pam3, file = "cluster_verde_represion_rapida.txt",quote = F,row.names = F,col.names = F)


kk.1 <- enrichKEGG(gene         = cluster1.pam3,
                 organism     = 'ana',
                 pvalueCutoff = 0.05)

as.data.frame(kk.1)

kk.2 <- enrichKEGG(gene         = cluster2.pam3,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.2)

kk.3 <- enrichKEGG(gene         = cluster3.pam3,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.3)

## heatmap

head(mean.gene.expression)
dim(mean.gene.expression)
clusters.mean.expression <- mean.gene.expression[c(cluster1.pam3,cluster2.pam3,cluster3.pam3),]

genes <- rownames(clusters.mean.expression)

scaled.gene.expression <- matrix(nrow=length(genes),ncol=97)

for(i in 1:length(genes))
{
  current.gene <- genes[i]
  
  expression.current.gene.i <- clusters.mean.expression[current.gene,]/max(clusters.mean.expression[current.gene,])
  
  scaled.gene.expression[i,] <- approx(x=c(0,3,9,24),expression.current.gene.i,xout=seq(from=0,to=24,by=0.25))[[2]]
}

library(gplots)
colfunc <- colorRampPalette(c("black","blue","yellow"))

png(filename = "figures/heatmap.png",width = 600)
heatmap.2(scaled.gene.expression,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = "")
dev.off()


scaled.gene.expression <- matrix(nrow=length(genes),ncol=4)

for(i in 1:length(genes))
{
  current.gene <- genes[i]
  
  scaled.gene.expression[i,] <- clusters.mean.expression[current.gene,]/max(clusters.mean.expression[current.gene,])
}

library(gplots)
colfunc <- colorRampPalette(c("black","blue","yellow"))

png(filename = "figures/heatmap_sin_interpolacion.png",width = 600)
heatmap.2(scaled.gene.expression,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = "")
dev.off()






##########################
## EXPLORATORY ANALYSIS ##
##########################


## Four clustres
cluster1.pam4 <- names(which(pam.4[["clustering"]] == 1))
length(cluster1.pam4)
cluster2.pam4 <- names(which(pam.4[["clustering"]] == 2))
length(cluster2.pam4)
cluster3.pam4 <- names(which(pam.4[["clustering"]] == 3))
length(cluster3.pam4)
cluster4.pam4 <- names(which(pam.4[["clustering"]] == 4))
length(cluster4.pam4)

expr.cluster1.pam4 <- mean.gene.expression[cluster1.pam4,]
expr.cluster2.pam4 <- mean.gene.expression[cluster2.pam4,]
expr.cluster3.pam4 <- mean.gene.expression[cluster3.pam4,]
expr.cluster4.pam4 <- mean.gene.expression[cluster4.pam4,]

mean.profile.cluster1.pam4 <- colMeans(expr.cluster1.pam4)
mean.profile.cluster2.pam4 <- colMeans(expr.cluster2.pam4)
mean.profile.cluster3.pam4 <- colMeans(expr.cluster3.pam4)
mean.profile.cluster4.pam4 <- colMeans(expr.cluster4.pam4)

plot(x = c(0,3,9,24),mean.profile.cluster1.pam4,type="o",col="blue",xlab="Time (h)",ylab="Expression Level (FPKM)",lwd=3,pch=0,main="Cluster Profiles",ylim=c(0,800))
lines(x = c(0,3,9,24),mean.profile.cluster2.pam4,type="o",col="red",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster3.pam4,type="o",col="darkgreen",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster4.pam4,type="o",col="black",lwd=3,pch=1)


kk.1 <- enrichKEGG(gene         = cluster1.pam4,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.1)

kk.2 <- enrichKEGG(gene         = cluster2.pam4,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.2)

kk.3 <- enrichKEGG(gene         = cluster3.pam4,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.3)

kk.4 <- enrichKEGG(gene         = cluster4.pam4,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.4)


## Five clustres
cluster1.pam5 <- names(which(pam.5[["clustering"]] == 1))
length(cluster1.pam5)
cluster2.pam5 <- names(which(pam.5[["clustering"]] == 2))
length(cluster2.pam5)
cluster3.pam5 <- names(which(pam.5[["clustering"]] == 3))
length(cluster3.pam5)
cluster4.pam5 <- names(which(pam.5[["clustering"]] == 4))
length(cluster4.pam5)
cluster5.pam5 <- names(which(pam.5[["clustering"]] == 5))
length(cluster5.pam5)

expr.cluster1.pam5 <- mean.gene.expression[cluster1.pam5,]
expr.cluster2.pam5 <- mean.gene.expression[cluster2.pam5,]
expr.cluster3.pam5 <- mean.gene.expression[cluster3.pam5,]
expr.cluster4.pam5 <- mean.gene.expression[cluster4.pam5,]
expr.cluster5.pam5 <- mean.gene.expression[cluster5.pam5,]

mean.profile.cluster1.pam5 <- colMeans(expr.cluster1.pam5)
mean.profile.cluster2.pam5 <- colMeans(expr.cluster2.pam5)
mean.profile.cluster3.pam5 <- colMeans(expr.cluster3.pam5)
mean.profile.cluster4.pam5 <- colMeans(expr.cluster4.pam5)
mean.profile.cluster5.pam5 <- colMeans(expr.cluster5.pam5)

plot(x = c(0,3,9,24),mean.profile.cluster1.pam5,type="o",col="blue",xlab="Time (h)",ylab="Expression Level (FPKM)",lwd=3,pch=0,main="Cluster Profiles",ylim=c(0,800))
lines(x = c(0,3,9,24),mean.profile.cluster2.pam5,type="o",col="red",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster3.pam5,type="o",col="darkgreen",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster4.pam5,type="o",col="black",lwd=3,pch=1)
lines(x = c(0,3,9,24),mean.profile.cluster5.pam5,type="o",col="orange",lwd=3,pch=1)


kk.1 <- enrichKEGG(gene         = cluster1.pam5,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.1)

kk.2 <- enrichKEGG(gene         = cluster2.pam5,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.2)

kk.3 <- enrichKEGG(gene         = cluster3.pam5,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.3)

kk.4 <- enrichKEGG(gene         = cluster4.pam5,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.4)

kk.5 <- enrichKEGG(gene         = cluster5.pam5,
                   organism     = 'ana',
                   pvalueCutoff = 0.05)
as.data.frame(kk.5)

