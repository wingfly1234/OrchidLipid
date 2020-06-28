setwd("~/orchid_lipid/")
# DEG analysis use DESeq2
o <- read.csv("gene_count_matrix.csv", row.names ="gene_id" )   #gene_count_matrix.csv got from stringtie analysis
library("DESeq2")#install.packages( "DESeq2" )
countData <- as.matrix(o)
colData <- data.frame("condition" = rep(c("Flower","Leaf"),each=3))
rownames(colData) <- colnames(countData)
all(rownames(colData) %in% colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
#dds$condition
dds$condition <- factor(dds$condition, levels = c("Leaf","Flower"))
dds$condition
dds <- DESeq(dds)
res <- results(dds)
res
hist( res$pvalue, breaks=20, col="grey" )
plotDispEsts( dds , ylim = c(1e-6, 1e1))
summary(res)
#vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds)
sampleDists <- dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:3, 1:3 ]
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- c("Flower rep1","Flower rep2","Flower rep3","Leaf rep1","Leaf rep2","Leaf rep3")
colnames(sampleDistMatrix)<- c("Flower rep1","Flower rep2","Flower rep3","Leaf rep1","Leaf rep2","Leaf rep3")
library( "gplots" )#install.packages( "gplots" )
library( "pheatmap" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap( sampleDistMatrix,
          clustering_distance_rows = sampleDists, 
          clustering_distance_cols = sampleDists,
          col = colours, border_color = NA )
plotPCA(rld)
#plotPCA(vsd)
#library(scales)
#show_col(hue_pal()(2)) 
library(ggplot2)
plotmyPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE,main="") 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = group)) + 
    geom_point(size = 2,shape=1) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed()
}
plotmyPCA(rld)
group <- factor(group, levels = c("Flower","Leaf"))
require("ggrepel")
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = group)) + 
     geom_point(size = 3) + geom_text_repel(label=d$name)+
 xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
   ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+ ylim(-15 ,15)+xlim(-65,65)
plotMA(res,alpha=0.05)
sum(res$padj < 0.05, na.rm=TRUE)
library("ggplot2")
table(res$Significant)
resOrdered <- res[order(res$pvalue),]
re <- as.data.frame(res)
re$Significant  <- "No"
re[is.na(res$padj),]$Significant  <- "Not test"
re[re$padj < 0.05 &re$log2FoldChange >=1 &!is.na(res$padj),]$Significant <-"Up"
re[re$padj < 0.05 &re$log2FoldChange<=-1&!is.na(res$padj),]$Significant  <-"Down"
re$Significant  <- factor(re$Significant,levels = c("Up","Down","No"))
v <- ggplot(re[!is.na(res$padj),],aes(log2FoldChange,-log10(padj)))
v + geom_point(aes(colour=Significant))+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)

write.table(re,"allDEGs.txt",sep="\t",quote=F)

#select <- re$Significant=="Up" | re$Significant=="Down"
#selectgene <- !is.na(rownames(re)[select])
pheatmap(assay(rld)[rownames(assay(rld))
                    %in%
                    rownames(re[re$Significant=="Up" | re$Significant=="Down",])
                    ,],
         border_color = NA,show_rownames = F,scale="row" ,labels_col = c("Flower rep1","Flower rep2","Flower rep3","Leaf rep1","Leaf rep2","Leaf rep3")
)

colnames(assay(rld)) <- c("Flower rep1","Flower rep2","Flower rep3","Leaf rep1","Leaf rep2","Leaf rep3")


#GO and KEGG Enrichment use clusterProfiler

#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("clusterProfiler")
library("clusterProfiler")
peq_go <- read.delim("peq_go_gene04",header=F)
pgo <- build_Anno(peq_go,peq_go_dis)
peq_go_dis <- read.delim("go_dis04",header=T)
up<-as.data.frame(rownames(re[re$Significant=="Up",]))
down<-as.data.frame(rownames(re[re$Significant=="Down",]))
colnames(up)<-c("V1")
colnames(down)<-c("V1")
ego <- enricher(up$V1, pvalueCutoff = 0.05, pAdjustMethod = "BH",minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = peq_go,TERM2NAME = peq_go_dis)
dego <- enricher(down$V1, pvalueCutoff = 0.05, pAdjustMethod = "BH",minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = peq_go,TERM2NAME = peq_go_dis)
head(ego)
barplot(ego, showCategory=20)
barplot(dego, showCategory=20)
emapplot(ego, showCategory = 20)
emapplot(dego,showCategory = 20)
goplot(ego)
dim(ego)
dim(ego[ego$ONTOLOGY=='BP',])
write.table(ego,"leaf_up_regulated.txt",sep="\t",quote=FALSE)
write.table(dego,"leaf_down_regulated.txt",sep="\t",quote=FALSE)


pca2ncbi <- read.delim("pha_plaza2ncbi",head=F)
head(pca2ncbi)
length(up)
up<-as.data.frame(rownames(re[re$Significant=="Up",]))
down<-as.data.frame(rownames(re[re$Significant=="Down",]))
colnames(up)<-c("V1")
colnames(down)<-c("V1")
fl_up_ncbi <- merge(up,pca2ncbi,by="V1")
fl_down_ncbi <- merge(down,pca2ncbi,by="V1")
fl_up_kegg <- bitr_kegg(fl_up_ncbi$V2, fromType='ncbi-proteinid', toType='kegg', organism='peq')
fl_down_kegg <- bitr_kegg(fl_down_ncbi$V2, fromType='ncbi-proteinid', toType='kegg', organism='peq')
fl_up<- merge(fl_up_ncbi,fl_up_kegg,by.x="V2",by.y="ncbi-proteinid")
fl_up<- fl_up[,c(2,1,3)]
colnames(fl_up)[1:2] <- c("ID","ncbi-proteinid")
write.table(fl_up,"flower_up_regulated_ID_ncbi_pid_keggid.txt",sep="\t",quote=FALSE,row.names = FALSE)
fl_down<- merge(fl_down_ncbi,fl_down_kegg,by.x="V2",by.y="ncbi-proteinid")
fl_down<- fl_down[,c(2,1,3)]
colnames(fl_down)[1:2] <- c("ID","ncbi-proteinid")
write.table(fl_down,"flower_down_regulated_ID_ncbi_pid_keggid.txt",sep="\t",quote=FALSE,row.names = FALSE)
kk <- enrichKEGG(gene         = fl_up_kegg$kegg,
                 organism     = 'peq',
                 pvalueCutoff = 0.05 )
head(kk)
browseKEGG(kk, 'peq00062')
browseKEGG(kk, 'peq00945')
browseKEGG(kk,'peq00941')
browseKEGG(kk, 'peq00945')
browseKEGG(kk,'peq00073')
browseKEGG(kk,'peq01040')
browseKEGG(kk,'peq00904')
browseKEGG(kk,'peq04626')
browseKEGG(kk,'peq04075')
browseKEGG(kk,'peq00940')
browseKEGG(kk,'peq00480')
browseKEGG(kk,'peq00360')
browseKEGG(kk,'peq00500')
write.table(kk,"flower_up_regulated_kegg_enrichment.txt",sep="\t",quote=FALSE)
mkk <- enrichMKEGG(gene = fl_up_kegg$kegg,
                   organism = 'peq')
head(mkk)
write.table(mkk,"flower_up_regulated_kegg_module_enrichment.txt",sep="\t",quote=FALSE)
browseKEGG(mkk, 'M00415')

dkk <- enrichKEGG(gene         = fl_down_kegg$kegg,
                 organism     = 'peq',
                 pvalueCutoff = 0.05 )
head(dkk)
write.table(dkk,"flower_down_regulated_kegg_enrichment.txt",sep="\t",quote=FALSE)
browseKEGG(dkk,'peq00195')
browseKEGG(dkk,'peq00710')
browseKEGG(dkk,'peq00196')
browseKEGG(dkk,'peq00860')
browseKEGG(dkk,'peq00970')
browseKEGG(dkk,'peq01200')
browseKEGG(dkk,'peq00630')
browseKEGG(dkk,'peq00220')
browseKEGG(dkk,'peq00010')
 browseKEGG(dkk,'peq01230')
browseKEGG(dkk,'peq00260')
browseKEGG(dkk,'peq00051')
browseKEGG(dkk,'peq00030')
 browseKEGG(dkk,'peq00910')
browseKEGG(dkk,'peq00500')
browseKEGG(dkk,'peq00052')
browsemKEGG(dmkk,'M00611')
dmkk <- enrichMKEGG(gene = fl_down_kegg$kegg,
                   organism = 'peq')
write.table(dmkk,"flower_down_regulated_kegg_module_enrichment.txt",sep="\t",quote=FALSE)

