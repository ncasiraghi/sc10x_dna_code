library( data.table )

#GeneModel <- read.delim("/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed",stringsAsFactors = F,header = T)
GeneModel <- read.delim("/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.bed",stringsAsFactors = F,header = T)
colnames(GeneModel)[1] <- 'chr'

xxx <- which(duplicated(GeneModel$GeneID))
length(xxx)

# which cells to be merged
#folder_cells_int_genemodel <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/tmp.parallel"
folder_cells_int_genemodel <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/tmp.parallel"

cellsList <- list.files(folder_cells_int_genemodel,pattern = "intcellids_",full.names = T)

get_and_filter <- function(bed){
  this <- fread(bed,data.table = FALSE)
  colnames(this)[5] <- "GeneID"
  m <- merge(x = GeneModel,y = this,by = "GeneID",all.x = TRUE)
  rownames(m) <- m$GeneID
  m <- m[GeneModel$GeneID,]
  return(round(m$V9,2))
}

cells <- lapply(cellsList, get_and_filter)

mat <- do.call(rbind,cells)

# save(mat,cellsList,GeneModel,file = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/intcellids_AllCells_GeneModel.RData",compress = TRUE)
save(mat,cellsList,GeneModel,file = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/intcellids_AllCells_GeneModel_Transcript.RData",compress = TRUE)

colnames(mat) <- GeneModel$GeneID

ids <- gsub(basename(cellsList),pattern = '-genemodel.sort.bed',replacement = '')
ids <- gsub(ids,pattern = 'intcellids_',replacement = 'cell_')

rownames(mat) <- ids

write.table(mat,file = '/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/intcellids_AllCells_GeneModel_Transcript.tsv',col.names = T,row.names = T,sep = '\t')

