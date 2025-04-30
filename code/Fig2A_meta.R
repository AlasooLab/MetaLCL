library(ggrepel)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(readr)
library(rtracklayer)


edb <- EnsDb.Hsapiens.v86

sig_rsid <- "rs4819670"
missense_rsid <- "rs3180408"
target_gene <- "USP18"
flank_size <- 2.6e5  

# Load data
meta_locus <- read_tsv("Fig2A_data.tsv")

# Ensure only required columns are selected
meta_locus_lz <- meta_locus %>%
  dplyr::select(chrom, pos, snp, meta_nlog10p, ID, phenotype_id) %>%
  as.data.frame()

meta_locus_lz$p <- 10^(-meta_locus$meta_nlog10p)


loc_meta_data <- locus(data = meta_locus_lz, gene = target_gene, flank = flank_size, ens_db = edb,index_snp="rs4819670")
ld_data <- read.delim("GEUVADIS_USP18_range_LDs.tsv", header = TRUE, sep = "\t")
loc_meta_data[["data"]] <- loc_meta_data[["data"]] %>%
  left_join(ld_data, by = "ID")



#loc_meta_data_ld <- link_LD(loc_meta_data, token = "token",genome_build="grch38_high_coverage")

#recomb.hg38 <- import.bw("recomb1000GAvg.bw")
#loc_meta_data <- link_recomb(loc_meta_data, recomb = recomb.hg38)
pdf("fig2A_meta.pdf", width = 5.5, height = 5.5)
locus_plot(loc_meta_data,labels = c(sig_rsid, missense_rsid),label_x=c(4, -5),highlight = target_gene) 
dev.off()