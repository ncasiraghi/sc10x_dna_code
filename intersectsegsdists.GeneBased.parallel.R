library(parallel)
mc.cores = 50

# sample1
#wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel"
#bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F)
#info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

# sample2
# wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample2.genemodel"
# bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample2-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F)
# colnames(bed) <- c("chrom","start","end","cell_id","copy_number","event_confidence")
# info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample2-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

# sample3
# wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample3.genemodel"
# bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample3-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F)
# colnames(bed) <- c("chrom","start","end","cell_id","copy_number","event_confidence")
# info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample3-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

# sample4
wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample4.genemodel"
bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample4-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F)
colnames(bed) <- c("chrom","start","end","cell_id","copy_number","event_confidence")
info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample4-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

setwd(wd)

# get single cells only
dd <- merge(x = bed,y = info.percell,by = "cell_id",all.y = T)
dd <- dd[with(dd,order(cell_id,chrom,start,end)),]

# get all combinations
AllCellIDs <- unique(dd$cell_id)
combList <- t(combn(AllCellIDs,m = 2))

dir.create(file.path(wd,"tmp.parallel"))

# write bed for each cell
writebed <- function(i,AllCellIDs){
  cell <- AllCellIDs[i]
  cell.dd <- dd[which(dd$cell_id == cell),c("chrom","start","end","copy_number","cell_id")]
  cell.bed <- file.path(wd,"tmp.parallel",paste0("cellid_",cell,".bed"))
  write.table(cell.dd,file = cell.bed,quote = F,col.names = F,row.names = F,sep = "\t")
}

mclapply(seq(1,length(AllCellIDs),1),writebed,AllCellIDs=AllCellIDs,mc.cores = mc.cores)

# Gene Model
GeneModel <- "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed"

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

# compare gene models
compareBedGeneBased <- function(i,combList){
  p <- combList[i,]
  message(paste(p,collapse = "|"))
  header <- c("chr.gene","start.gene","end.gene","Gene.name","Gene.ID","chr.seg","start.seg","end.seg","cn","cellid")
  a.bed.file <- file.path(wd,"tmp.parallel",paste0("intcellids_",p[1],"-","genemodel.sort.bed"))
  b.bed.file <- file.path(wd,"tmp.parallel",paste0("intcellids_",p[2],"-","genemodel.sort.bed"))
  # compare a.gm.bed and b.gm.bed
  a.bed <- read.delim(a.bed.file,as.is = T,header = F,stringsAsFactors = F)
  colnames(a.bed) <- header
  b.bed <- read.delim(b.bed.file,as.is = T,header = F,stringsAsFactors = F)
  colnames(b.bed) <- header
  ib <- merge(x = a.bed,y = b.bed,by = c("Gene.ID","chr.gene","start.gene","end.gene"),all = F,suffixes = c("_a","_b"))
  out <- data.frame(a.cell_id=ib$cellid_a[1],b.cell_id=ib$cellid_b[1],dist=as.numeric(dist(rbind(ib$cn_a,ib$cn_b),method = "euclidean")),stringsAsFactors = F)
  return(out)
}

intbedout <- mclapply(seq(1,nrow(combList),1),compareBedGeneBased,combList=combList,mc.cores = mc.cores)
save(intbedout,AllCellIDs,file = "intbedout.genebased.RData",compress = T)

if(TRUE){

  distmat <- matrix(ncol = length(AllCellIDs),nrow = length(AllCellIDs),data = NA,dimnames = list(AllCellIDs,AllCellIDs))
  for(i in 1:length(intbedout)){
    d <- intbedout[[i]]
    distmat[as.character(d$a.cell_id),as.character(d$b.cell_id)] <- d$dist
    distmat[as.character(d$b.cell_id),as.character(d$a.cell_id)] <- d$dist
  }
  
  distmat <- as.dist(distmat)
  save(distmat,AllCellIDs,file = "distmat.genebased.RData",compress = T)

}
