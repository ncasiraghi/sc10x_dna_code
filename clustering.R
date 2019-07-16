# sample1
# wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/"
# igvname <- "igvsif.sample1.txt"

# sample2
# wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample2.genemodel/"
# igvname <- "igvsif.sample2.txt"

# sample1 + sample2
wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1sample2.genemodel/"
igvname <- "igvsif.sample1sample2.txt"

setwd(wd)
load("distmat.genebased.RData")
class(distmat)

distmat <- as.matrix(distmat)

# remove noisy cells as labelled from 10x cnv pipeline
metrics_s1 <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv",header = T,sep=",",stringsAsFactors = F)
metrics_s2 <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample2-CNV/outs/per_cell_summary_metrics.csv",header = T,sep=",",stringsAsFactors = F)

keep.cell_id.s1 <- paste0(metrics_s1$cell_id[which(metrics_s1$is_noisy==0)],"-s1")
keep.cell_id.s2 <- paste0(metrics_s2$cell_id[which(metrics_s2$is_noisy==0)],"-s2")

distmat <- distmat[c(keep.cell_id.s1,keep.cell_id.s2),c(keep.cell_id.s1,keep.cell_id.s2)]
# distmat <- distmat[gsub(keep.cell_id.s1,pattern = "-s1",replacement = ""),gsub(keep.cell_id.s1,pattern = "-s1",replacement = "")]
# distmat <- distmat[gsub(keep.cell_id.s2,pattern = "-s2",replacement = ""),gsub(keep.cell_id.s2,pattern = "-s2",replacement = "")]

# remove noisy cells based on copy number segs [ ad-hoc filter]

masterDFgenemodel_s1 <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample1.Rdata"
masterDFgenemodel_s2 <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample2.Rdata"

clean_noisy_cells <- function(masterDFgenemodel, sample_label){
  load(masterDFgenemodel)
  keep.cells.customfilter <- master[which(master$flag10x_is.noisy == 0 ),]
  keep.cells.customfilter <- keep.cells.customfilter[which(keep.cells.customfilter$loss < .5 ),]
  keep.cells.customfilter <- paste0(as.character(keep.cells.customfilter$id),sample_label)
  return(keep.cells.customfilter)
}

keep.cells.customfilter_s1 <- clean_noisy_cells(masterDFgenemodel = masterDFgenemodel_s1,sample_label = "-s1")
keep.cells.customfilter_s2 <- clean_noisy_cells(masterDFgenemodel = masterDFgenemodel_s2,sample_label = "-s2")

keep.clean.cells <- unique(c(keep.cells.customfilter_s1,keep.cells.customfilter_s2))

distmat <- distmat[keep.clean.cells,keep.clean.cells]

distmat <- as.dist(distmat)

# Hierarchical clustering
hc<-hclust(distmat,method="ward.D2")
dend <- as.dendrogram(hc)
plot(hc,hang = -1, cex = 0.06)

ids <- rownames(as.matrix(distmat))

#sif <- data.frame(id=ids,dendrogram.index=match(ids,labels(dend)))
sif <- data.frame(id=ids,dendrogram.index=match(ids,labels(dend)),sample=NA)
sif$sample[grep(sif$id,pattern = "-s1")] <- "primary.nuclei"
sif$sample[grep(sif$id,pattern = "-s2")] <- "relapse.nuclei"

#write.table(x = sif,file = igvname,col.names = T,row.names = F,quote = F,sep = "\t")
write.table(x = sif,file = paste0("no_noisy_cells_",igvname),col.names = T,row.names = F,quote = F,sep = "\t")

# correct seg
seg <- read.delim("sample1sample2.10xcn.seg",as.is = T,stringsAsFactors = F,header = F)
seg <- seg[which(seg[,1] %in% as.character(sif$id)),]
write.table(x = seg,file = "no_noisy_cells_sample1sample2.10xcn.seg",col.names = F,row.names = F,quote = F,sep = "\t")
