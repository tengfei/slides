### R code from vignette source 'main.Rnw'

###################################################
### code chunk number 1: up (eval = FALSE)
###################################################
## library(BiocInstaller)
## useDevel(TRUE)
## biocLite("ggbio")


###################################################
### code chunk number 2: rs (eval = FALSE)
###################################################
## source("~/.Rprofile")


###################################################
### code chunk number 3: chip-sample
###################################################
library(chipseq)
data(cstest)
## estimate fragment length
fraglen <- estimate.mean.fraglen(cstest$ctcf)
fraglen[!is.na(fraglen)]
ctcf.ext <- resize(cstest$ctcf, width = 200)
cov.ctcf <- coverage(ctcf.ext)
gfp.ext <- resize(cstest$gfp, width = 200)
cov.gfp <- coverage(gfp.ext)


###################################################
### code chunk number 4: peak-diff
###################################################
peakCutoff(cov.ctcf, fdr = 0.0001)
peaks.ctcf <- slice(cov.ctcf, lower = 8)
peaks.gfp <- slice(cov.gfp, lower = 8)
peakSummary <- diffPeakSummary(peaks.gfp, peaks.ctcf)
peakSummary <-within(peakSummary,
    {
      diffs <- asinh(sums2) - asinh(sums1)
      resids <- (diffs - median(diffs)) / mad(diffs)
      up <- resids > 2
      down <- resids < -2
      change <- ifelse(up, "up", 
                       ifelse(down, "down", "flat"))
    })
pks <- as(peakSummary, 'GRanges')


###################################################
### code chunk number 5: gr
###################################################
head(pks)


###################################################
### code chunk number 6: chip-seqnames
###################################################
chrs <- unique(as.character(seqnames(pks)))
library(GenomicRanges)
pks <- keepSeqlevels(pks, chrs)
seqlengths(pks) <- seqlengths(cstest)[chrs]
head(pks)


###################################################
### code chunk number 7: overview
###################################################
pks <- pks[order(values(pks)$diffs, 
                 decreasing = TRUE)][1:50]
autoplot(pks, layout = "karyogram", 
         aes(color = diffs, fill = diffs))


###################################################
### code chunk number 8: chip-gr
###################################################
wh <- GRanges("chr11", IRanges(3062594, 3092593))


###################################################
### code chunk number 9: gf-default
###################################################
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
library(ggbio)
p.gene <- autoplot(txdb, which = wh)    
p.gene


###################################################
### code chunk number 10: gf-chevron
###################################################
p.gene.c <- autoplot(txdb, which = wh, 
                     gap.geom = "chevron")        
p.gene.c


###################################################
### code chunk number 11: gf-reduce
###################################################
p.gene.r <- autoplot(txdb, which = wh, 
                     stat = "reduce")        
library(grid)
print(p.gene.r, vp = viewport(height = 0.15, width = 1))


###################################################
### code chunk number 12: gr-instance
###################################################
## shorter name
cstest.s <- stack(cstest)
cstest.s <- resize(cstest.s, width = 200)
chipseq <- subsetByOverlaps(cstest.s, wh)
head(chipseq)


###################################################
### code chunk number 13: default-chipseq
###################################################
autoplot(chipseq)


###################################################
### code chunk number 14: chipseq-facets
###################################################
autoplot(chipseq, facets = strand ~ .)
## equivalent to 
## autoplot(gr, geom = "chevron") +  
##              facet_grid(strand ~ seqnames)


###################################################
### code chunk number 15: prac-facets
###################################################
autoplot(chipseq, facets = sample ~ .)


###################################################
### code chunk number 16: chipseq-cov
###################################################
autoplot(chipseq, stat = "coverage", 
         facets = sample ~ .)


###################################################
### code chunk number 17: chipseq-cov-area
###################################################
autoplot(chipseq, stat = "coverage", 
         geom = "area", facets = sample ~ .)


###################################################
### code chunk number 18: geom-stat-gr
###################################################
ggbio:::.ggbio.geom
ggbio:::.ggbio.stat


###################################################
### code chunk number 19: so-cov
###################################################
p.cov <- autoplot(chipseq, stat = "coverage", 
                  facets = sample ~ ., geom = "area")
p.cov


###################################################
### code chunk number 20: chipseq-arb
###################################################
autoplot(chipseq, color = "blue", fill = "red")


###################################################
### code chunk number 21: chipseq-aes
###################################################
autoplot(chipseq, aes(fill = strand))


###################################################
### code chunk number 22: prac-aes
###################################################
autoplot(chipseq, color = "black", aes(fill = strand), 
         facets = sample ~.)


###################################################
### code chunk number 23: gr-scale
###################################################
p.cov + theme_bw() + 
  geom_hline(yintercept = 8, color = "red")


###################################################
### code chunk number 24: chip-tracks
###################################################
tks <-tracks("coverage" = p.cov, 
             "gene" = p.gene, xlim = wh,
              xlim.change = c(TRUE, TRUE)) 
tks


###################################################
### code chunk number 25: tracks-label
###################################################
tks <- tks + coord_cartesian(xlim  = c(3.085e6, 3.09e6)) + theme_alignment()
tks


###################################################
### code chunk number 26: own-data
###################################################
simul <- seq(from = start(wh), to = end(wh), by = 1)
df <- data.frame(x = simul, y = rnorm(length(simul)))
p <- qplot(data = df, x = x, y = y, geom = "line")
tracks("coverage" = p.cov, "gene" = p.gene, 
       "user" = p, xlim = wh)


###################################################
### code chunk number 27: ideogram
###################################################
library(rtracklayer)
## require network
head(ucscGenomes())
library(biovizBase)
## getIdeogram()
obj <- getIdeogram("hg19")
obj <- getIdeogram("hg19", cytoband = FALSE)


###################################################
### code chunk number 28: seqinfo
###################################################
autoplot(seqinfo(obj)[paste0("chr", c(1:22, 'X'))])


###################################################
### code chunk number 29: ideogram
###################################################
library(biovizBase)
## mm9 <- getIdeogram("mm9")
cyto.def <- getOption("biovizBase")$cytobandColor
cyto.new <- c(cyto.def, c(gpos33 = "grey80", gpos66 = "grey60"))
p.ideo <- plotIdeogram(mm9, "chr10", zoom = c(start(wh),end(wh)))  + 
  scale_fill_manual(values = cyto.new) 
print(p.ideo)


###################################################
### code chunk number 30: tracks-ideo
###################################################
tks <- tracks("ideogram" = p.ideo, "coverage" = p.cov, 
              "gene" = p.gene, 
              xlim = wh, 
              xlim.change = c(FALSE, TRUE, TRUE), 
              heights = c(1, 5, 5))
tks


###################################################
### code chunk number 31: processing
###################################################
crc1 <- system.file("extdata", "crc1-missense.csv", package = "biovizBase")
crc1 <- read.csv(crc1)
library(GenomicRanges)
mut.gr <- with(crc1,GRanges(Chromosome, IRanges(Start_position, End_position),
                            strand = Strand))
values(mut.gr) <- subset(crc1, select = -c(Start_position, End_position, Chromosome))
data("hg19Ideogram", package = "biovizBase")
seqs <- seqlengths(hg19Ideogram)
## subset_chr
chr.sub <- paste("chr", 1:22, sep = "")
## levels tweak
seqlevels(mut.gr) <- c(chr.sub, "chrX")
mut.gr <- keepSeqlevels(mut.gr, chr.sub)
seqs.sub <- seqs[chr.sub]
## remove wrong position
bidx <- end(mut.gr) <= seqs.sub[match(as.character(seqnames(mut.gr)),
              names(seqs.sub))]
mut.gr <- mut.gr[which(bidx)]
## assign_seqlengths
seqlengths(mut.gr) <- seqs.sub
## reanme to shorter names
new.names <- as.character(1:22)
names(new.names) <- paste("chr", new.names, sep = "")
mut.gr.new <- renameSeqlevels(mut.gr, new.names)


###################################################
### code chunk number 32: ideo-track
###################################################
hg19Ideo <- hg19Ideogram
hg19Ideo <- keepSeqlevels(hg19Ideogram, chr.sub)
hg19Ideo <- renameSeqlevels(hg19Ideo, new.names)
p <- ggplot() + layout_circle(hg19Ideo, geom = "ideo", 
          fill = "gray70", radius = 30, trackWidth = 4)
p


###################################################
### code chunk number 33: circle-scale
###################################################
p <- p + layout_circle(hg19Ideo, geom = "scale", 
              size = 2, radius = 35, trackWidth = 2)
p


###################################################
### code chunk number 34: circle-text
###################################################
p <- p + layout_circle(hg19Ideo, geom = "text", 
                       aes(label = seqnames), vjust = 0,
                       radius = 38, trackWidth = 7)
p


###################################################
### code chunk number 35: circle-mut
###################################################
p <- p + layout_circle(mut.gr, geom = "rect", 
                       color = "steelblue",
                       radius = 23 ,trackWidth = 6)
p


###################################################
### code chunk number 36: circle-link-process
###################################################
rearr  <- read.csv(system.file("extdata", 
                               "crc-rearrangment.csv", package = "biovizBase"))
## start position
gr1 <- with(rearr, GRanges(chr1, IRanges(pos1, width = 1)))
## end position
gr2 <- with(rearr, GRanges(chr2, IRanges(pos2, width = 1)))
## add extra column
nms <- colnames(rearr)
.extra.nms <- setdiff(nms, c("chr1", "chr2", "pos1", "pos2"))
values(gr1) <- rearr[,.extra.nms]
## remove out-of-limits data
seqs <- as.character(seqnames(gr1))
.mx <- seqlengths(hg19Ideo)[seqs]
idx1 <- start(gr1) > .mx
seqs <- as.character(seqnames(gr2))
.mx <- seqlengths(hg19Ideo)[seqs]
idx2 <- start(gr2) > .mx
idx <- !idx1 & !idx2
gr1 <- gr1[idx]
seqlengths(gr1) <- seqlengths(hg19Ideo)
gr2 <- gr2[idx]
seqlengths(gr2) <- seqlengths(hg19Ideo)
values(gr1)$to.gr <- gr2
## rename to gr
gr <- gr1
values(gr)$rearrangements <- ifelse(as.character(seqnames(gr))
              == as.character(seqnames((values(gr)$to.gr))),
             "intrachromosomal", "interchromosomal")
gr.crc1 <- gr[values(gr)$individual == "CRC-1"]


###################################################
### code chunk number 37: circle-point
###################################################
p <- p + layout_circle(gr.crc1, geom = "point", 
         aes(y = score, size = tumreads), color = "red",
         radius = 12 ,trackWidth = 10, grid = TRUE) +
  scale_size(range = c(1, 2.5))
p


###################################################
### code chunk number 38: circle-link
###################################################
p <- p + layout_circle(gr.crc1, geom = "link", 
                       linked.to = "to.gr", 
                       aes(color = rearrangements),
                       radius = 10 ,trackWidth = 1)
p


