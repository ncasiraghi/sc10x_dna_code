SEG_FILE <- '/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/no_noisy_cells_sample1.10xcn.seg'

seg <- read.delim(SEG_FILE,as.is = T,header = F)
colnames(seg) <- c('cell_id','chr','start','end','cn')

seg <- seg[which(seg$chr == 7 ),]

cell_id <- unique(seg$cell_id)

bs <- sort(unique(seg$start))

m <- matrix(data = 0,nrow = length(cell_id),ncol = length(bs))

for(i in seq_len(nrow(seg))){
  idx_row <- which(cell_id==seg$cell_id[i])
  idx_col <- which(bs==seg$start[i])
  m[idx_row,idx_col] <- 1
}

library(gplots)

barplot(colSums(m),las=2,names.arg = NA)

x <- c()
for(i in 2:length(bs)){
  x <- c(x,bs[i]-bs[i-1])
}
barplot(x)
summary(x)

borders <- as.numeric(which(colSums(m)==length(cell_id)))

mp <- m[,(borders[1]+1):(borders[2]-1)]

j <- which.max(colSums(mp,na.rm = TRUE))

mp <- mp[order(mp[,j]),]

