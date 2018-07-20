library(segmentSeq)
library(doParallel)
registerDoParallel(cores = 6)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

inDir <- "/projects/ajt200/Heidi_RAC1/"
covDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/"
covDir2 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/"
covDir3 <- "/projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/log2ChIPinput/"
covDir4 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/"
covDir5 <- "/projects/ajt200/BAM_masters/RNAseq/wt/multi_unique/bam_all/coverage/"
plotDir <- "/projects/ajt200/Heidi_RAC1/plots/"
gff.data <- read.csv(file = paste0(inDir, "genes_gff.csv"),header=F)
exon.gff <- gff.data[which(as.character(gff.data[,3]) == "exon"),]
gene.gff <- gff.data[which(as.character(gff.data[,3]) == "gene"),]
tir <- c(11293676, 11293257)
walkera <- c(11292988, 11292959)
walkerb <- c(11292733, 11292707)
arc1glpl <- c(11292472, 11292458)
arc2mhd <- c(11292163, 11292149)
lrr1 <- c(11291617, 11291898)
lrr2 <- c(11291096, 11291506)
lrr3 <- c(11289346, 11290527)
lrr4 <- c(11289244, 11289261)

all.tes <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand.txt", header = T)
chr.tes <- all.tes[which(all.tes[,1] == "Chr1"),]
rac1.tes <- chr.tes[which(chr.tes[,2] >= 11288000 & chr.tes[,3] <= 11300000),]

################################################################
# read crossover maps for RAC1 from Col/Ler, Col/Mh and Col/Wi #
################################################################

col.ler <- read.csv(file = paste0(inDir, "RAC1_Col_Ler.csv"),header = T)
col.mh <- read.csv(file = paste0(inDir, "RAC1_Col_Mh.csv"),header = T)
col.wl <- read.csv(file = paste0(inDir, "RAC1_Col_Wl.csv"), header = T)

########################
# polymorphism density #
########################

coller.coords <- col.ler[,1]
colmh.coords <- col.mh[,1]
colwl.coords <- col.wl[,1]
start <- coller.coords[1]
end <- coller.coords[length(coller.coords)]
wins <- seq(start, end, by=100)
wins <- c(wins, end)
coller.win <- NULL
colmh.win <- NULL
colwl.win <- NULL
for(k in 1:length(wins)-1) {
	print(k)
	coller.win <- c(coller.win, length(coller.coords[which(coller.coords >= wins[k] & coller.coords < wins[k+1])]))
	colmh.win <- c(colmh.win, length(colmh.coords[which(colmh.coords >= wins[k] & colmh.coords < wins[k+1])]))
	colwl.win <- c(colwl.win, length(colwl.coords[which(colwl.coords >= wins[k] & colwl.coords < wins[k+1])]))
}

# Calculate mean coverage values within windows of defined size
inNames <- c("REC8_log2ChIPinput", "WT_log2nucNakedDNA", "WT_H3K4me3_log2ChIPinput", "WT_H3K9me2_log2ChIPinput", "MSH4_log2ChIPinput")
inNames2 <- c("log2SPO11oligoNakedDNA")
inNames3 <- c("WT_log2SPO11ChIP4REC8input")
inNames4 <- c("H2A_log2ChIPinput", "H2AW_log2ChIPinput", "H2AX_log2ChIPinput", "H2AZ_log2ChIPinput")
inNames5 <- c("wt_RNAseq_ATCACG")
libNames <- c("REC8", "MNase", "H3K4me3", "H3K9me2", "MSH4")
libNames2 <- c("SPO11-1-oligos")
libNames3 <- c("SPO11-1 ChIP")
libNames4 <- c("H2A", "H2A.W", "H2A.X", "H2A.Z")
libNames5 <- c("RNA-seq")

windows <- c(100)
winNames <- c("100bp")
for(s in 1:length(windows)) {
  print(s)
  covWins <- windows[s]
  print(winNames[s])
  foreach(k = 1:length(inNames), .combine = 'c') %dopar% {
    print(k)
    dat <- read.table(file = paste0(covDir, inNames[k], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
    winDat <- NULL
    chrDat <- dat[dat$V1 == "Chr1",]
    cov <- dat[,4]
    covCoords <- seq(1, length(cov), by = 1)
    covIRcoords <- IRanges(start = covCoords, width = 1)
    covGRcoords <- GRanges(seqnames = "Chr1", strand = "+", ranges = covIRcoords)
    seqWindows <- seq(1, length(cov), by = covWins)
    seqWindows <- c(seqWindows, length(cov))
    windowsIRanges <- IRanges(start = seqWindows, width = covWins)
    windowsGRanges <- GRanges(seqnames = "Chr1", strand = "+", ranges = windowsIRanges)
    overlaps <- getOverlaps(windowsGRanges, covGRcoords, whichOverlaps = T)
    covWinVals <- sapply(overlaps, function(x) mean(cov[x]))
    winDat <- cbind(seqWindows, covWinVals)
    write.table(winDat, file = paste0(inDir, inNames, "_RAC1_norm_coverage_", winNames[s], ".txt"))
  }
}

#REC8 <- read.table(file = paste0(inDir, "REC8_log2ChIPinput_RAC1_norm_coverage_", winNames[1], ".txt"))
#REC8_RAC1 <- REC8[REC8[,1] >= start & REC8[,1] <= end,]

# Or plot coverage at single base-pair resolution
dat1 <- mclapply(seq_along(inNames), function(x) {
  read.table(file = paste0(covDir, inNames[x], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
}, mc.cores = 5)
dat2 <- lapply(seq_along(inNames2), function(x) {
  read.table(file = paste0(covDir2, inNames2[x], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
})
dat3 <- lapply(seq_along(inNames3), function(x) {
  read.table(file = paste0(covDir3, inNames3[x], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
})
dat4 <- mclapply(seq_along(inNames4), function(x) {
  read.table(file = paste0(covDir4, inNames4[x], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
}, mc.cores = 4)
dat5 <- lapply(seq_along(inNames5), function(x) {
  read.table(file = paste0(covDir5, inNames5[x], "_norm_allchrs_coverage_coord_tab_11300000.bed"))
})
dat <- c(dat1, dat2, dat3, dat4, dat5)
dat_RAC1 <- lapply(seq_along(dat), function(x) {
  dat[[x]][dat[[x]][,2] >= start & dat[[x]][,2] <= end,]
})


###########
# Col/Ler #
###########


plot
plot(as.numeric(as.character(col.ler[,1])), as.numeric(col.ler[,6]), type = "l", ylim=c(-5,150), main="Col vs Ler")
rug(as.numeric(col.ler[,1]))
for(k in 1:length(exon.gff[,1])){
        y.plot <- 150
        segments(exon.gff[k,4],y.plot,exon.gff[k,5],y.plot,col=1)
}
y.plot <- 145
segments(tir[1],y.plot,tir[2],y.plot,col="green")
segments(walkera[1],y.plot,walkera[2],y.plot,col=2)
segments(walkerb[1],y.plot,walkerb[2],y.plot,col=2)
segments(arc1glpl[1],y.plot,arc1glpl[2],y.plot,col=2)
segments(arc2mhd[1],y.plot,arc2mhd[2],y.plot,col=2)
segments(lrr1[1],y.plot,lrr1[2],y.plot,col=4)
segments(lrr2[1],y.plot,lrr2[2],y.plot,col=4)
segments(lrr3[1],y.plot,lrr3[2],y.plot,col=4)
segments(lrr4[1],y.plot,lrr4[2],y.plot,col=4)
for(k in 1:length(gene.gff[,1])){
        abline(v=gene.gff[k,4],col=1,lwd=0.5)
        abline(v=gene.gff[k,5],col=1,lwd=0.5)
}
for(k in 1:length(rac1.tes[,1])){
        y.plot <- 140
        segments(rac1.tes[k,2],y.plot,rac1.tes[k,3],y.plot,col=2,lwd=1)
}
par(new=T)
plot(wins,coller.win,type="l",col=4,lwd=0.6,lty=2,yaxt="n")
axis(side=4)



plotCovRAC1 <- function(dat, ylab) { 
  plot(dat[,3], dat[,4], type = "l", ylim = c(min(dat[,4]), max(dat[,4])), ylab = ylab)
  for(k in 1:length(exon.gff[,1])) {
    y.plot <- max(dat[,4])
    segments(x0 = exon.gff[k, 4], y0 = y.plot, x1 = exon.gff[k, 5], y1 = y.plot, col = 1)
  }
  y.plot <- max(dat[,4])-0.2
  segments(x0 = tir[1], y0 = y.plot, x1 = tir[2], y1 = y.plot, col = "green")
  segments(x0 = walkera[1], y0 = y.plot, x1 = walkera[2], y1 = y.plot, col = 2)
  segments(x0 = walkerb[1], y0 = y.plot, x1 = walkerb[2], y1 = y.plot, col = 2)
  segments(x0 = arc1glpl[1], y0 = y.plot, x1 = arc1glpl[2], y1 = y.plot, col = 2)
  segments(x0 = arc2mhd[1], y0 = y.plot, x1 = arc2mhd[2], y1 = y.plot, col = 2)
  segments(x0 = lrr1[1], y0 = y.plot, x1 = lrr1[2], y1 = y.plot, col = 4)
  segments(x0 = lrr2[1], y0 = y.plot, x1 = lrr2[2], y1 = y.plot, col = 4)
  segments(x0 = lrr3[1], y0 = y.plot, x1 = lrr3[2], y1 = y.plot, col = 4)
  segments(x0 = lrr4[1], y0 = y.plot, x1 = lrr4[2], y1 = y.plot, col = 4)
  for(k in 1:length(gene.gff[,1])) {
    abline(v = gene.gff[k, 4], col = 4, lwd = 1) 
    abline(v = gene.gff[k, 5], col = 4, lwd = 1)
  }
  for(k in 1:length(rac1.tes[,1])) {
    y.plot <- max(dat[,4])-0.4
    segments(x0 = rac1.tes[k, 2], y0 = y.plot, x1 = rac1.tes[k, 3], y1 = y.plot, col = 2, lwd = 1)
  }
  box(lwd = 1.5)
}

pdf(paste0(plotDir, "RAC1_log2-transformed_coverage_profiles_v080917.pdf"), height = 26, width = 7)
par(mfcol=c(12,1))
par(mar=c(2.1,5.1,2.1,5.1))
plotCovRAC1(dat = dat_RAC1[[1]], ylab = libNames[1])
plotCovRAC1(dat = dat_RAC1[[2]], ylab = libNames[2])
plotCovRAC1(dat = dat_RAC1[[6]], ylab = libNames2[1])
plotCovRAC1(dat = dat_RAC1[[7]], ylab = libNames3[1])
plotCovRAC1(dat = dat_RAC1[[3]], ylab = libNames[3])
plotCovRAC1(dat = dat_RAC1[[4]], ylab = libNames[4])
plotCovRAC1(dat = dat_RAC1[[5]], ylab = libNames[5])
plotCovRAC1(dat = dat_RAC1[[8]], ylab = libNames4[1])
plotCovRAC1(dat = dat_RAC1[[9]], ylab = libNames4[2])
plotCovRAC1(dat = dat_RAC1[[10]], ylab = libNames4[3])
plotCovRAC1(dat = dat_RAC1[[11]], ylab = libNames4[4])
plotCovRAC1(dat = dat_RAC1[[12]], ylab = libNames5[1])
dev.off()


