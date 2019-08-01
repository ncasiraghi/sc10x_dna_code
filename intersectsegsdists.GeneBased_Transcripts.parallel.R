library(parallel)
mc.cores = 5

# sample1
wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript"
bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F,col.names = c("chrom","start","end","cell_id","copy_number","event_confidence"))
info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

setwd(wd)

# get single cells only
dd <- merge(x = bed,y = info.percell,by = "cell_id",all.y = T)
dd <- dd[with(dd,order(cell_id,chrom,start,end)),]

# get all combinations
AllCellIDs <- unique(dd$cell_id)
#combList <- t(combn(AllCellIDs,m = 2))

#dir.create(file.path(wd,"tmp.parallel"))

# write bed for each cell
writebed <- function(i,AllCellIDs){
  cell <- AllCellIDs[i]
  cell.dd <- dd[which(dd$cell_id == cell),c("chrom","start","end","copy_number","cell_id")]
  cell.bed <- file.path(wd,"tmp.parallel",paste0("cellid_",cell,".bed"))
  write.table(cell.dd,file = cell.bed,quote = F,col.names = F,row.names = F,sep = "\t")
}

mclapply(seq(1,length(AllCellIDs),1),writebed,AllCellIDs=AllCellIDs,mc.cores = mc.cores)

# Gene Model
GeneModel <- "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.bed"

# keep one-seg per gene
selseg <- function(int.sort.bed){
  sb <- read.delim(int.sort.bed,as.is = T,stringsAsFactors = F,header = F)
  # use ensmbl IDs
  genes.with.multiple.segs <- unique(sb[,5][which(duplicated(sb[,5]))])
  rows.to.remove <- c()
  for(enid in genes.with.multiple.segs){
    k <- which(sb[,5]==enid)
    if(length(unique(sb[k,1]))==1){
      #sb[k[1],9] <- mean(sb[k[1:length(k)],9])
      sb[k[1],9] <- weighted.mean(x = sb[k[1:length(k)],9],w = sb[k[1:length(k)],3]-sb[k[1:length(k)],2])
      sb[k[1],3] <- sb[k[length(k)],3]
      sb[k[1],6] <- paste(sb[k[1:length(k)],6],collapse=",")
      sb[k[1],7] <- paste(sb[k[1:length(k)],7],collapse=",")
      sb[k[1],8] <- paste(sb[k[1:length(k)],8],collapse=",")
      rows.to.remove <- c(rows.to.remove,k[2:length(k)])
    } else {
      rows.to.remove <- c(rows.to.remove,k)
    }
  }
  sb <- sb[-rows.to.remove,]
  return(sb)
}

# compute intersections between sc-bed and GeneModel
intSingleCellGeneModel <- function(i,AllCellIDs,GeneModel){
  cell <- AllCellIDs[i]
  cell.bed <- file.path(wd,"tmp.parallel",paste0("cellid_",cell,".bed"))
  # intersect sc-bed a with GeneModel
  int.bed <- file.path(wd,"tmp.parallel",paste0("intcellids_",cell,"-","genemodel.bed"))
  cmd.int <- paste("intersectBed","-a",GeneModel,"-b",cell.bed,"-wb >",int.bed)
  system(cmd.int)
  int.sort.bed <- gsub(int.bed,pattern = "\\.bed",replacement = ".sort.bed")
  cmd.sort <- paste("sortBed -i",int.bed,">",int.sort.bed)
  system(cmd.sort)
  # keep only one-seg per gene
  gm.bed <- selseg(int.sort.bed=int.sort.bed)
  gm.bed.file <- int.sort.bed
  write.table(gm.bed,file = gm.bed.file,col.names = F,row.names = F,quote = F,sep = "\t")
  system(paste("rm",int.bed))
}

mclapply(seq(1,length(AllCellIDs),1),intSingleCellGeneModel,AllCellIDs=AllCellIDs,GeneModel=GeneModel,mc.cores = mc.cores)
