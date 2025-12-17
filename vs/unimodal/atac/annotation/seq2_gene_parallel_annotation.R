# Run this bit in R/4.3.0
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # dm6 genome annotation
library(GenomicRanges)
library(dplyr)
library(doParallel)
options(warn=-1)

region_file <- "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/Annotations/peaks.txt"
regions <- readLines(region_file)

# --- Step 2: Convert peak strings to GRanges ---
parse_peak <- function(peak) {
  chrom <- sub(":.*", "", peak)
  pos <- sub(".*:", "", peak)
  start <- as.integer(sub("-.*", "", pos))
  end <- as.integer(sub(".*-", "", pos))
  GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
}

gr_list <- lapply(regions, parse_peak)
gr <- do.call(c, gr_list)
seqlevels(gr) <- paste0("chr", seqlevels(gr))


# --- Step 3: Annotate with ChIPseeker using ±1 kb window around TSS --- (modified this to increase upstream to 5kb)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
peakAnno <- annotatePeak(
  peak = gr,
  TxDb = txdb,
  tssRegion = c(-5000, 1000),
  verbose = FALSE
)

cluster <- makeCluster(20)

registerDoParallel(cluster)

x <- foreach(i = 1:length(gr), .combine = c) %dopar% {
  
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # dm6 genome annotation
  library(GenomicRanges)
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  peakAnno <- annotatePeak(
  peak = gr,
  TxDb = txdb,
  tssRegion = c(-5000, 1000),
  verbose = FALSE
  )

  res <- tryCatch(
    paste(list(seq2gene(gr[i], tssRegion=c(-3000,1000), flankDistance=3000, txdb))[[1]], collapse = ","),
    error = function(e) return("Not available")
  )
  return(res)
}

y <- foreach(i = 1:length(gr), .combine = c) %dopar% {
  
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # dm6 genome annotation
  library(GenomicRanges)
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  peakAnno <- annotatePeak(
  peak = gr,
  TxDb = txdb,
  tssRegion = c(-5000, 1000),
  verbose = FALSE
  )

  res <- tryCatch(
     paste0(names(subsetByOverlaps(genes(txdb,single.strand.genes.only=FALSE),gr[i],maxgap=3000)),collapse=","),
    error = function(e) return("Not available")
  )
  return(res)
}

genes <- cbind(as.data.frame(paste(gr@seqnames,":",gr@ranges)),x)
colnames(genes)<-c("Regions","Genes")
allgenes <- cbind(as.data.frame(paste(gr@seqnames,":",gr@ranges)),y)
colnames(allgenes)<-c("Regions","Genes")

write.csv(genes, "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/Annotations/signac_peaks_annotated_3-1kb_seq2gene.csv",row.names=FALSE)
write.csv(allgenes, "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/Annotations/signac_peaks_annotated_3-1kb_all.csv",row.names=FALSE)

stopCluster(cluster)