load('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/intcellids_AllCells_GeneModel.RData')

dim(mat)

colnames(mat) <- GeneModel$GeneID

ids <- gsub(basename(cellsList),pattern = '-genemodel.sort.bed',replacement = '')
ids <- gsub(ids,pattern = 'intcellids_',replacement = 'cell_')

rownames(mat) <- ids

write.table(mat,file = '/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/intcellids_AllCells_GeneModel.tsv',col.names = T,row.names = T,sep = '\t')
