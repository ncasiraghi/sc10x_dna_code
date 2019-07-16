# compare gene models
library(parallel)
wd <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna"
mc.cores <- 50

cells.s1 <- list.files(file.path(wd,"sample1.genemodel","tmp.parallel"),pattern = "intcellids_",full.names = T)
df.s1 <- data.frame(bedFile=cells.s1,id=gsub(basename(cells.s1),pattern = "\\.bed",replacement = ".s1"),stringsAsFactors = F)

cells.s2 <- list.files(file.path(wd,"sample2.genemodel","tmp.parallel"),pattern = "intcellids_",full.names = T)
df.s2 <- data.frame(bedFile=cells.s2,id=gsub(basename(cells.s2),pattern = "\\.bed",replacement = ".s2"),stringsAsFactors = F)

df <- rbind(df.s1,df.s2)
rownames(df) <- df$id

AllCellIDs <- unique(df$id)
combList <- t(combn(AllCellIDs,m = 2))

compareBedGeneBased <- function(i,combList){
  p <- combList[i,]
  message(paste(p,collapse = "|"))
  header <- c("chr.gene","start.gene","end.gene","Gene.name","Gene.ID","chr.seg","start.seg","end.seg","cn","cellid")
  a.bed.file <- df[p[1],1]
  b.bed.file <- df[p[2],1]
  # compare a.gm.bed and b.gm.bed
  a.bed <- read.delim(a.bed.file,as.is = T,header = F,stringsAsFactors = F)
  colnames(a.bed) <- header
  b.bed <- read.delim(b.bed.file,as.is = T,header = F,stringsAsFactors = F)
  colnames(b.bed) <- header
  ib <- merge(x = a.bed,y = b.bed,by = c("Gene.ID","chr.gene","start.gene","end.gene"),all = F,suffixes = c("_a","_b"))
  out <- data.frame(a.cell_id=gsub(p[1],pattern = "intcellids_|genemodel.sort.",replacement = ""),
                    b.cell_id=gsub(p[2],pattern = "intcellids_|genemodel.sort.",replacement = ""),
                    dist=as.numeric(dist(rbind(ib$cn_a,ib$cn_b),method = "euclidean")),stringsAsFactors = F)
  return(out)
}

setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1sample2.genemodel")
intbedout <- mclapply(seq(1,nrow(combList),1),compareBedGeneBased,combList=combList,mc.cores = mc.cores)
save(intbedout,df,AllCellIDs,file = "intbedout.genebased.RData",compress = T)

if(FALSE){
  dimnames <- gsub(AllCellIDs,pattern = "intcellids_|genemodel.sort.",replacement = "")
  distmat <- matrix(ncol = length(AllCellIDs),nrow = length(AllCellIDs),data = NA,dimnames = list(dimnames,dimnames))
  
  for(i in 1:length(intbedout)){
    d <- intbedout[[i]]
    distmat[as.character(d$a.cell_id),as.character(d$b.cell_id)] <- d$dist
    distmat[as.character(d$b.cell_id),as.character(d$a.cell_id)] <- d$dist
  }
  
  distmat <- as.dist(distmat)
  save(distmat,AllCellIDs,file = "distmat.genebased.RData",compress = T)
  
}
