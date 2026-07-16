# Supplemental Figure 2a
# By Nolan Origer

library(biomaRt)
library(Gviz)
library(GenomicRanges)

# Connect to the human genome dataset GRCh38.p14
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch SHROOM3 gene structure
s3.structure <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                              "ensembl_transcript_id", 
                                              "exon_chrom_start", 
                                              "exon_chrom_end",
                                              "start_position",
                                              "end_position"),
                           filters = "hgnc_symbol",
                           values = "SHROOM3",
                           mart = ensembl)
s3prot.structure <- dplyr::filter(s3.structure, ensembl_transcript_id %in% c("ENST00000296043","ENST00000646790","ENST00000644244"))
shrna.structure <- read.csv("~/shroom3-trc-shrna.csv", header = T)


# Create a GRanges object for exons
s3.gr <- GRanges(seqnames = "4", 
                    ranges = IRanges(start = s3.structure$exon_chrom_start, 
                                     end = s3.structure$exon_chrom_end),
                 transcript = s3.structure$ensembl_transcript_id)
genome(s3.gr) <- "GRCh38.p14"

s3prot.gr <- GRanges(seqnames = "4", 
                 ranges = IRanges(start = s3prot.structure$exon_chrom_start, 
                                  end = s3prot.structure$exon_chrom_end),
                 transcript = s3prot.structure$ensembl_transcript_id)
genome(s3prot.gr) <- "GRCh38.p14"

# Create a GRanges object for shRNA targets
shrna.gr <- GRanges(seqnames = "4", 
                    ranges = IRanges(start = shrna.structure$shrna_chrom_start, 
                                     end = shrna.structure$shrna_chrom_end),
                    transcript = shrna.structure$shrna_id)
genome(shrna.gr) <- "GRCh38.p14"

# Create track objects
gtrack <- GenomeAxisTrack()

s3.track <- GeneRegionTrack(s3.gr, 
                            chromosome = 4,
                            transcriptAnnotation = "transcript",
                            name = "SHROOM3",
                            rotation.title = 0,
                            cex.title = 0.8)
s3prot.track <- GeneRegionTrack(s3prot.gr, 
                            chromosome = 4,
                            transcriptAnnotation = "transcript",
                            name = "SHROOM3",
                            rotation.title = 0,
                            cex.title = 0.8)
shrna.track <- GeneRegionTrack(shrna.gr, 
                               chromosome = 4,
                               transcriptAnnotation = "transcript",
                               name = "shRNAs",
                               fill = "#FD5E53",
                               rotation.title = 0,
                               cex.title = 0.8)

# Plot the tracks
Gviz::plotTracks(list(s3prot.track, shrna.track, gtrack), 
                 chromosome = 4,
                 from = min(start(s3prot.track)), 
                 to = max(end(s3prot.track)),
                 extend.left = 50000,
                 extend.right = 1000,
                 col = NULL)