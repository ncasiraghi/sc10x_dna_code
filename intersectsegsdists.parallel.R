library(parallel)
mc.cores = 20

wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna"
setwd(wd)

bed <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/node_cnv_calls.bed",comment.char = "#",header = F)
colnames(bed) <- c("chrom","start","end","cell_id","copy_number","event_confidence")

# how many cells?
info.percell <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv",header = T,sep = ",",stringsAsFactors = F)

# get single cells only
dd <- merge(x = bed,y = info.percell,by = "cell_id",all.y = T)
dd <- dd[with(dd,order(cell_id,chrom,start,end)),]

# get all combinations
AllCellIDs <- unique(dd$cell_id)
combList <- t(combn(AllCellIDs,m = 2))

# set.seed(88)
# combList <- combList[sample(nrow(combList),size = 5000),]

# write bed for each cell

writebed <- function(i,AllCellIDs){
  cell <- AllCellIDs[i]
  cell.dd <- dd[which(dd$cell_id == cell),c("chrom","start","end","copy_number","cell_id")]
  cell.bed <- file.path(wd,"paralleltest",paste0("cellid_",cell,".bed"))
  write.table(cell.dd,file = cell.bed,quote = F,col.names = F,row.names = F,sep = "\t")
}

mclapply(seq(1,length(AllCellIDs),1),writebed,AllCellIDs=AllCellIDs,mc.cores = mc.cores)

# compute intersections

intbeds <- function(i,combList){
  p <- combList[i,]
  # bed a
  a <- dd[which(dd$cell_id == p[1]),c("chrom","start","end","copy_number","cell_id")]
  a.bed <- file.path(wd,"paralleltest",paste0("cellid_",p[1],".bed"))
  # bed b
  b <- dd[which(dd$cell_id == p[2]),c("chrom","start","end","copy_number","cell_id")]
  b.bed <- file.path(wd,"paralleltest",paste0("cellid_",p[2],".bed"))
  # intersect bed
  int.bed <- file.path(wd,"paralleltest",paste0("intcellids_",p[1],"-",p[2],".bed"))
  cmd.int <- paste("intersectBed","-a",a.bed,"-b",b.bed,"-wb >",int.bed)
  system(cmd.int)
  int.sort.bed <- gsub(int.bed,pattern = "\\.bed",replacement = ".sort.bed")
  cmd.sort <- paste("sortBed -i",int.bed,">",int.sort.bed)
  system(cmd.sort)
  system(paste("rm",int.bed))
  ib <- read.delim(file = int.sort.bed,as.is = T,header = F,stringsAsFactors = F)
  ib <- ib[,c(1:5,9,10)]
  colnames(ib) <- c("chr","start","end","a.cn","a.cell_id","b.cn","b.cell_id")
  out <- data.frame(a.cell_id=ib$a.cell_id[1],b.cell_id=ib$b.cell_id[1],dist=as.numeric(dist(rbind(ib$a.cn,ib$b.cn),method = "euclidean")),stringsAsFactors = F)
  system(paste("rm",int.sort.bed))
  return(out)
}

intbedout <- mclapply(seq(1,nrow(combList),1),intbeds,combList=combList,mc.cores = mc.cores)
save(intbedout,file = "intbedout.RData",compress = T)

# fill distance matrix
distmat <- matrix(ncol = length(AllCellIDs),nrow = length(AllCellIDs),data = NA,dimnames = list(0:(length(AllCellIDs)-1),0:(length(AllCellIDs)-1)))

for(i in 1:length(intbedout)){
  d <- intbedout[[i]]
  distmat[as.character(d$a.cell_id),as.character(d$b.cell_id)] <- d$dist
  distmat[as.character(d$b.cell_id),as.character(d$a.cell_id)] <- d$dist
}

distmat <- as.dist(distmat)
save(distmat,file = "distmat.RData",compress = T)

