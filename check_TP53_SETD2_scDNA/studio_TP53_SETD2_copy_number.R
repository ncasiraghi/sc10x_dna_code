setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sc10x_dna_code/check_TP53_SETD2_scDNA")

# 3	47057919	47205457	SETD2	ENSG00000181555
# 17	7565097	7590856	TP53	ENSG00000141510

scbed <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/tmp.parallel/",pattern = 'intcellids_',full.names = TRUE)

tab <- c()
for(cell in scbed){
  m <- read.delim(cell,stringsAsFactors = F,as.is = T,header = F)
  m <- m[which(m$V4 %in% c("SETD2","TP53")),]
  
  this <- data.frame(CELL=m$V10[1],
                     cn_TP53=m[which(m$V4 == "TP53"),9],
                     cn_SETD2=m[which(m$V4 == "SETD2"),9])
  
  tab <- rbind(tab,this)
}

# clean scDNA data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample1.Rdata")
keep <- as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)])
tab$NOISY_CELL <- 1
tab$NOISY_CELL[which(tab$CELL %in% keep)] <- 0

# compute freq loss

tab$cn_TP53 <- round(tab$cn_TP53)
tab$cn_SETD2 <- round(tab$cn_SETD2)

sum(table(tab$cn_TP53))
sum(table(tab$cn_SETD2))

# all cells 
out1 <- data.frame(tot.cells=nrow(tab),
                  n.cells_lossTP53=length(which(tab$cn_TP53 < 2)),
                  n.cells_lossSETD2=length(which(tab$cn_SETD2 < 2)),
                  n.cells_lossTP53_and_SETD2=length(which(tab$cn_TP53 < 2 & tab$cn_SETD2 < 2)),
                  fraction.cells_lossTP53=round(length(which(tab$cn_TP53 < 2))/nrow(tab),2),
                  fraction.cells_lossSETD2=round(length(which(tab$cn_TP53 < 2))/nrow(tab),2),
                  fraction.cells_lossTP53_and_SETD2=round(length(which(tab$cn_TP53 < 2 & tab$cn_SETD2 < 2))/nrow(tab),2),
                  stringsAsFactors = F)

# not noisy cells 
clean <- tab[which(tab$NOISY_CELL==0),]
out2 <- data.frame(tot.cells=nrow(clean),
                  n.cells_lossTP53=length(which(clean$cn_TP53 < 2)),
                  n.cells_lossSETD2=length(which(clean$cn_SETD2 < 2)),
                  n.cells_lossTP53_and_SETD2=length(which(clean$cn_TP53 < 2 & clean$cn_SETD2 < 2)),
                  fraction.cells_lossTP53=round(length(which(clean$cn_TP53 < 2))/nrow(clean),2),
                  fraction.cells_lossSETD2=round(length(which(clean$cn_TP53 < 2))/nrow(clean),2),
                  fraction.cells_lossTP53_and_SETD2=round(length(which(clean$cn_TP53 < 2 & clean$cn_SETD2 < 2))/nrow(clean),2),
                  stringsAsFactors = F)


out <- rbind(out1,out2)
write.table(out,file = "TP53_SETD2_lossfreq_scDNA_sample1.tsv",sep = "\t",col.names = T,row.names = F,quote = F)

# summary plot
par(pty="s",mfrow=c(2,2))
x <- barplot(table(tab$cn_TP53),ylab = "n. cells",xlab="copy number status",main="TP53 (all cells)",ylim=c(0,410))
text(x = x,y = table(tab$cn_TP53),labels = table(tab$cn_TP53),pos = 3)

x <- barplot(table(tab$cn_SETD2),ylab = "n. cells",xlab="copy number status",main="SETD2 (all cells)",ylim=c(0,410))
text(x = x,y = table(tab$cn_SETD2),labels = table(tab$cn_SETD2),pos = 3)

x <- barplot(table(clean$cn_TP53),ylab = "n. cells",xlab="copy number status",main="TP53 (not-noisy cells)",ylim=c(0,410))
text(x = x,y = table(clean$cn_TP53),labels = table(clean$cn_TP53),pos = 3)

x <- barplot(table(clean$cn_SETD2),ylab = "n. cells",xlab="copy number status",main="SETD2 (not-noisy cells)",ylim=c(0,410))
text(x = x,y = table(clean$cn_SETD2),labels = table(clean$cn_SETD2),pos = 3)

