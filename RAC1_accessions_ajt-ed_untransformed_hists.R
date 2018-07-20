library(segmentSeq)
library(ggplot2)
library(doParallel)
registerDoParallel(cores = 6)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

inDir <- "/projects/ajt200/Heidi_RAC1/"
covDir1 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/"
covDir2 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/"
covDir3 <- "/projects/ajt200/BAM_masters/H3K4me3/WT/coverage/"
covDir4 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/"
covDir5 <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/"
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
rac1.tes <- chr.tes[which(chr.tes[,2] >= 11288000 & chr.tes[,3] <= 11310368),]

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
end <- 11310400
#end <- coller.coords[length(coller.coords)]
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
inNames1 <- c("WT_SPO11-oligo_RPI1_nakedDNA_R1")
inNames2 <- c("WT_nuc_nakedDNA")
inNames3 <- c("WT_H3K4me3_ChIP_WT_H3K9me2_input")
inNames4 <- c("WT_H3K9me2_ChIP_input")
inNames5 <- c("REC8_ChIP_input", "MSH4_ChIP_input")
libNamesCombined <- c(inNames1, inNames2, inNames3, inNames4, inNames5)

windows <- c(100)
winNames <- c("100bp")
for(s in 1:length(windows)) {
  print(s)
  covWins <- windows[s]
  print(winNames[s])
  foreach(k = 1:length(inNames1), .combine = 'c') %dopar% {
    print(k)
    dat <- read.table(file = paste0(covDir1, inNames1[k], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
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
    write.table(winDat, file = paste0(inDir, inNames1, "_RAC1_norm_coverage_", winNames[s], ".txt"))
  }
}

#REC8 <- read.table(file = paste0(inDir, "REC8_log2ChIPinput_RAC1_norm_coverage_", winNames[1], ".txt"))
#REC8_RAC1 <- REC8[REC8[,1] >= start & REC8[,1] <= end,]

# Or plot coverage at single base-pair resolution
dat1 <- lapply(seq_along(inNames1), function(x) {
  read.table(file = paste0(covDir1, inNames1[x], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
})
dat2 <- lapply(seq_along(inNames2), function(x) {
  read.table(file = paste0(covDir2, inNames2[x], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
})
dat3 <- lapply(seq_along(inNames3), function(x) {
  read.table(file = paste0(covDir3, inNames3[x], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
})
dat4 <- lapply(seq_along(inNames4), function(x) {
  read.table(file = paste0(covDir4, inNames4[x], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
})
dat5 <- lapply(seq_along(inNames5), function(x) {
  read.table(file = paste0(covDir5, inNames5[x], "_norm_allchrs_coverage_coord_tab_11320000.bed"))
})

dat <- c(dat1, dat2, dat3, dat4, dat5)
dat_RAC1 <- lapply(seq_along(dat), function(x) {
  dat[[x]][dat[[x]][,2] >= start & dat[[x]][,2] <= end,]
})
names(dat_RAC1) <- libNamesCombined

ggHistCovRAC1_ChIP <- function(dat, ...) {

y.seg <- max(dat_RAC1[[1]]$V4)-0.5
df <- data.frame(x0 = tir[1], y0 = y.seg, x1 = tir[2], y1 = y.seg)
plotRAC1 <- ggplot(dat_RAC1[[1]], aes(x = V2, y = V4)) +
            geom_bar(stat = "identity", fill = "red") +
            labs(x = "", y = "SPO11-1") + 
            geom_segment(data = df, aes(x = x0, y = y0, xend = x1, yend = y1), colour = "green", size = 2) +   
            theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
            #geom_smooth(method = lm)
ggsave(plotRAC1, file = paste0(plotDir, "plotRAC1.pdf"))

 
ChIPplotCovRAC1 <- function(dat, ylab) { 
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

pdf(paste0(plotDir, "RAC1_coverage_profiles_untransformed_SPO11-1-oligos_and_MNase_v110917.pdf"), height = 26, width = 7)
par(mfcol=c(12,1))
par(mar=c(2.1,5.1,2.1,5.1))
plotCovRAC1(dat = dat_RAC1[[1]], ylab = names(dat_RAC1[1]))
plotCovRAC1(dat = dat_RAC1[[12]], ylab = names(dat_RAC1[12]))
plotCovRAC1(dat = dat_RAC1[[5]], ylab = names(dat_RAC1[5]))
plotCovRAC1(dat = dat_RAC1[[6]], ylab = names(dat_RAC1[6]))
plotCovRAC1(dat = dat_RAC1[[2]], ylab = names(dat_RAC1[2]))
plotCovRAC1(dat = dat_RAC1[[3]], ylab = names(dat_RAC1[3]))
plotCovRAC1(dat = dat_RAC1[[4]], ylab = names(dat_RAC1[4]))
plotCovRAC1(dat = dat_RAC1[[7]], ylab = names(dat_RAC1[7]))
plotCovRAC1(dat = dat_RAC1[[8]], ylab = names(dat_RAC1[8]))
plotCovRAC1(dat = dat_RAC1[[9]], ylab = names(dat_RAC1[9]))
plotCovRAC1(dat = dat_RAC1[[10]], ylab = names(dat_RAC1[10]))
plotCovRAC1(dat = dat_RAC1[[11]], ylab = names(dat_RAC1[11]))
dev.off()


