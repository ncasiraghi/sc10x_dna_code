setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

load('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/intcellids_AllCells_GeneModel_Transcript.RData')

dim(mat)

colnames(mat) <- GeneModel$GeneID

ids <- gsub(basename(cellsList),pattern = '-genemodel.sort.bed',replacement = '')
ids <- gsub(ids,pattern = 'intcellids_',replacement = 'cell_')

rownames(mat) <- ids

# clean dna data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample1.Rdata")

keep <- as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)])
keep <- paste0('cell_',keep)

mat <- mat[keep,]

dim(mat)

sum(is.na(mat))

# NA per cell
hist(rowSums(is.na(mat)),50)
summary(rowSums(is.na(mat)))

# NA per gene
hist(colSums(is.na(mat)))
summary(colSums(is.na(mat)))
        
mat <- mat[,which(colSums(is.na(mat)) == 0)]

dim(mat)

# PCA
pca <- prcomp(x = mat,scale. = FALSE)
pca.var <- (pca$sdev^2)/sum(pca$sdev^2)
loading_scores_pca1 <- sort(abs(pca$rotation[,1]),decreasing = T)
loading_scores_pca2 <- sort(abs(pca$rotation[,2]),decreasing = T)
pca.data <- data.frame(cell_id=rownames(mat), PCA_1=pca$x[,1], PCA_2=pca$x[,2],stringsAsFactors = F)

top_ls_pca1 <- loading_scores_pca1[1:1000]
top_ls_pca2 <- loading_scores_pca2[1:1000]

par(pty="s")
plot(abs(pca$rotation[,1]),abs(pca$rotation[,2]),ylim=c(0,0.05),xlim=c(0,0.03))

top_loading_scores <- intersect(names(top_ls_pca1),names(top_ls_pca2))


# tSNE
library( Rtsne )
set.seed(9)
tsne = Rtsne(mat, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=.5, dims=2)
tsne.data = data.frame(cell_id=rownames(mat), TSNE_1=tsne$Y[,1], TSNE_2=tsne$Y[,2], stringsAsFactors = F)

# plot
#par(pty='s',mfrow=c(1,2))
plot(x = pca.data$PCA_1,pca.data$PCA_2,pch=16,col=grey(.3,.6),main = paste('n. cells: ',nrow(tsne.data)))
plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=16,col=grey(.3,.6),main = paste('n. cells: ',nrow(tsne.data)))

# scDNA clustering 

# hierarchical cluster model
fit_cluster_hierarchical <- hclust(dist(scale(tsne.data[,2:3])))
tsne.data$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3)) 

# assign cell-barcode to cell-id
per_cell_summary_metrics <- '/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv'

barcodes_tab <- read.delim(per_cell_summary_metrics,header = T,as.is = T,stringsAsFactors = F,sep = ',')
barcodes_tab$cell_id <- paste0('cell_',barcodes_tab$cell_id)

scDNA_clustering_info <- merge(x = tsne.data,y = barcodes_tab,by = 'cell_id',all.x = T)

write.table(scDNA_clustering_info,file = 'outs/Sample1_scDNA_clustering_info_Transcripts.tsv',col.names = T,row.names = F,quote = F,sep = '\t')

# get median cn per gene by cluster

clusters <- as.character(unique(tsne.data$cl_hierarchical))

agg_mat <- matrix(nrow = ncol(mat),
                  ncol = length(clusters),
                  data = NA,
                  dimnames = list(colnames(mat), clusters))

dim(agg_mat)

for(i in clusters){
  cluster_barcodes <- tsne.data$cell_id[grep(tsne.data$cl_hierarchical,pattern = i)]
  agg_mat[,i] <- apply(mat[cluster_barcodes,],MARGIN = 2,FUN = median)
}


layout.matrix <- matrix(c(1,1,1,2,
                          1,1,1,2,
                          1,1,1,2), nrow = 3, ncol = 4,byrow = T)

layout(layout.matrix)

colDNA <- c("#3182bd","#2ca25f","#feb24c")
names(colDNA) <- 1:3

plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=21,col='white',cex=3,las=1,cex.lab=1.2,ylab="TSNE_2",xlab="TSNE_1",
     bg=colDNA[as.character(tsne.data$cl_hierarchical)],main = paste('n. cells: ',nrow(tsne.data)))
grid(lty = 1)
legend("bottomleft",fill = colDNA,legend = names(colDNA))

boxplot(agg_mat[top_loading_scores,],outline=F,ylab="Median Gene Copy Number",las=1,xlab="Clusters",cex.lab=1.2,col=colDNA)
# add info on pca loading scores

agg_mat_table <- as.data.frame(agg_mat,stringsAsFactors = F)
agg_mat_table <- cbind(rownames(agg_mat),agg_mat_table)
colnames(agg_mat_table)[1] <- 'gene_id'

pca_loading_scores <- as.data.frame(pca$rotation[,1:2],stringsAsFactors = F)
pca_loading_scores$gene_id <- rownames(pca$rotation)
  
agg_mat_table <- merge(x = agg_mat_table,y = pca_loading_scores,by = 'gene_id',all.x = T)

write.table(agg_mat_table,file = 'outs/Sample1_scDNA_medianCN_Transcript_clusters.tsv',col.names = T,row.names = F,quote = F,sep = '\t')
