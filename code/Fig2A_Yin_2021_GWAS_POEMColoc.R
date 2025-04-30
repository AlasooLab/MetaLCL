library(ggrepel)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(readr)


# Load Gene Data
edb <- EnsDb.Hsapiens.v86

sig_rsid <- "rs4819670"
missense_rsid <- "rs3180408"

target_gene <- "USP18"
flank_size <- 2.6e5 


sle_imputed_locus <- read_tsv("Fig2A_data.tsv")

sle_imputed_locus <- sle_imputed_locus %>%
  dplyr::select(chrom, pos, snp, log10p, ID, phenotype_id) %>%
  as.data.frame()

sle_imputed_locus$p <- 10^(-sle_imputed_locus$log10p)


loc_sle_data <- locus(data = sle_imputed_locus, gene = target_gene, flank = flank_size, ens_db = edb,index_snp="rs4819670")
ld_data <- read.delim("GEUVADIS_USP18_range_LDs.tsv", header = TRUE, sep = "\t")
loc_sle_data[["data"]] <- loc_sle_data[["data"]] %>%
  left_join(ld_data, by = "ID")

#recomb.hg38 <- import.bw("recomb1000GAvg.bw")
#loc_sle_data <- link_recomb(loc_sle_data, recomb = recomb.hg38)

#loc_sle_data_ld <- link_LD(loc_sle_data, token = "token",genome_build="grch38_high_coverage")
pdf("Fig2A_Yin_2021_GWAS_POEMColoc.pdf", width = 5.5, height = 5.5)
locus_plot(loc_sle_data,labels = c(sig_rsid, missense_rsid),label_x=c(4, -5),highlight = target_gene) 
dev.off()