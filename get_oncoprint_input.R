load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/intcellids_AllCells_GeneModel.RData")

dim(mat)

drivers <- unique(readLines("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sc10x_dna_code/mb_driver_genes.txt"))

# clean dna data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample1.Rdata")

keep <- paste0("intcellids_",as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)]))

rownames_mat <- gsub(basename(cellsList),pattern = "-genemodel.sort.bed",replacement = "")

mat <- mat[which(rownames_mat %in% keep),]

dim(mat)

# Sample, Gene, Alteration (AMP), Type (CNA).

goi_index <- which(GeneModel$GeneName %in% drivers)

mat <- mat[,goi_index]
goi <- GeneModel[goi_index,]

dim(mat)

oncoprint <- c()

for(i in seq(1,length(keep))){
  message(keep[i])
  this <- data.frame(Sample = keep[i],
                     Gene = goi$GeneName,
                     Alteration = NA,
                     Type = "CNA",
                     cn = mat[i,],
                     stringsAsFactors = F)
  
  this$Alteration[which(this$cn > 2)] <- "GAIN"
  this$Alteration[which(this$cn > 3)] <- "AMP"
  #this$Alteration[which(this$cn < 2)] <- "HETLOSS"
  #this$Alteration[which(this$cn < 1)] <- "HOMDEL"
  this <- this[which(!is.na(this$Alteration)),]
  
  oncoprint <- rbind(oncoprint, this)
}

unaltered_samples <- setdiff(keep, unique(oncoprint$Sample))

file.name <- 'oncoprint_scdna_tumor_GAIN.tsv'
write.table(oncoprint[,1:4],file = file.name,col.names = F,row.names = F,quote = F,sep = '\t')
cat(unaltered_samples,file = file.name,sep = "\n",append = T)
