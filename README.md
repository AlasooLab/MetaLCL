


We performed the largest _trans_-eQTL study in a single cell type where we profiled the expression of 18,792 genes in 3,737 individuals from nine cohorts. We then replicated these findings in an independent multi-ancestry MAGE cohort of 682 individuals.

## Datasets used in the analysis

| Cohort | Publication | Sample size | Assay | Raw expression data | Raw genotype data |
|---|---|---|---|---| ---|
ALSPAC | [Boyd 2013](https://doi.org/10.1093/ije/dys064), [Fraser 2013](https://doi.org/10.1093/ije/dys066), [Bryois 2014](https://doi.org/10.1371/journal.pgen.1004461) | 877 | microarray | [ALSPAC Data Access](http://www.bristol.ac.uk/alspac/researchers/access/) | [ALSPAC Data Access](http://www.bristol.ac.uk/alspac/researchers/access/) |
TwinsUK | [Buil 2015](http://dx.doi.org/10.1038/ng.3162) | 735 | RNA-seq | [EGAD00001001086](https://ega-archive.org/datasets/EGAD00001001086) | [TwinsUK Data Access](https://twinsuk.ac.uk/resources-for-researchers/access-our-data/) |
MAGE | [Taylor 2023](https://doi.org/10.1101/2023.11.04.565639) | 682 | RNA-seq | [PRJNA851328](https://www.ebi.ac.uk/ena/browser/view/PRJNA851328) | [1000 Genomes](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) |
CoLaus | [Firmann 2008](https://doi.org/10.1186/1471-2261-8-6), [Sönmez Flitman 2021](https://doi.org/10.1021/acs.jproteome.1c00585) | 553 | RNA-seq | [CoLaus Data Access](https://www.colaus-psycolaus.ch/professionals/how-to-collaborate/) | [CoLaus Data Access](https://www.colaus-psycolaus.ch/professionals/how-to-collaborate/) |
GEUVADIS | [Lappalainen 2013](https://doi.org/10.1038/nature12531) | 358 | RNA-seq | [PRJEB3366](https://www.ebi.ac.uk/ena/browser/view/PRJEB3366) | [1000 Genomes](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) |
MRCE | [Liang 2013](https://doi.org/10.1101%2Fgr.142521.112) | 484 | microarray | [E-MTAB-1428](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1428) | [EGAS00000000137](https://ega-archive.org/studies/EGAS00000000137) |
MRCA | [Liang 2013](https://doi.org/10.1101%2Fgr.142521.112) | 327 | microarray | [E-MTAB-1425](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1425) | [EGAS00000000137](https://ega-archive.org/studies/EGAS00000000137) |
GENCORD | [Gutierrez-Arcelus 2013](https://doi.org/10.7554/eLife.00523) | 187 | RNA-seq | [EGAD00001000425](https://ega-archive.org/datasets/EGAD00001000425) | [EGAD00001000428](https://ega-archive.org/datasets/EGAD00001000428) |
GTEx | [GTEx v8](https://www.science.org/doi/10.1126/science.aaz1776) | 113 | RNA-seq | [phs000424.v8.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2) | [phs000424.v8.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2) |
CAP | [CAP](https://doi.org/10.1186/s12864-020-06966-4) | 100 | RNA-seq | [phs000481.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000481.v3.p2) | [phs000481.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000481.v3.p2) |

## Accessing summary statistics

The full _trans_-eQTL meta-analysis summary statistics are publicly available on[ eQTL Catalogue FTP server](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/MetaLCL/).

Preferred URLs for downloading the summary statistics files can be found in the **summary_statistics_urls.tsv** file.

An example of how to download a summary statistics file:

```
curl -O ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/MetaLCL/Liang_2013_ALSPAC_CAP_eur_GENCORD_GEUVADIS_GTEx_TwinsUK_CoLausR2_naive_LCL_metaanalysis_ENSG00000138646_id.parquet

```


Description of the column names in summary statistics files can be found in **Columns.md**.



## Acknowledgements

If you find the results of the _trans_-eQTL meta-analysis useful in your research, please cite all nine participating studies **(TODO)** :

### ALSPAC
* Bryois, J. et al. Cis and trans effects of human genomic variants on gene expression. PLoS Genet. 10, e1004461 (2014).
* Boyd, A. et al. Cohort Profile: the ‘children of the 90s’--the index offspring of the Avon Longitudinal Study of Parents and Children. Int. J. Epidemiol. 42, 111–127 (2013).
* Fraser, A. et al. Cohort Profile: the Avon Longitudinal Study of Parents and Children: ALSPAC mothers cohort. Int. J. Epidemiol. 42, 97–110 (2013).

### TwinsUK
* Buil, A. et al. Gene-gene and gene-environment interactions detected by transcriptome sequence analysis in twins. Nat. Genet. 47, 88–91 (2015).

### CoLaus
* Firmann, M. et al. The CoLaus study: a population-based study to investigate the epidemiology and genetic determinants of cardiovascular risk factors and metabolic syndrome. BMC Cardiovasc. Disord. 8, 6 (2008).
* Sönmez Flitman, R. et al. Untargeted Metabolome- and Transcriptome-Wide Association Study Suggests Causal Genes Modulating Metabolite Concentrations in Urine. J. Proteome Res. 20, 5103–5114 (2021).

### GEUVADIS

* Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506–511 (2013).

### MRCA and MRCE

* Liang, L. et al. A cross-platform analysis of 14,177 expression quantitative trait loci derived from lymphoblastoid cell lines. Genome Res. 23, 716–726 (2013).

### GENCORD

* Gutierrez-Arcelus, M. et al. Passive and active DNA methylation and the interplay with genetic variation in gene regulation. Elife 2, e00523 (2013).

### GTEx v8

* Consortium, T. G. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318–1330 (2020).

### CAP
* Theusch, E., Chen, Y.-D. I., Rotter, J. I., Krauss, R. M. & Medina, M. W. Genetic variants modulate gene expression statin response in human lymphoblastoid cell lines. BMC Genomics 21, 555 (2020).

## Methods

More detailed descriptions of the methods used, along with citations, can be found in this [paper]() **(TODO)**.

### Genotype data quality control and imputation
#### Pre-imputation quality control
We lifted coordinates of the genotyped variants to the GRCh38 build with CrossMap v0.4.182. We aligned the strands of the genotyped variants to the 1000 Genomes 30x on GRCh38 reference panel using Genotype Harmonizer. We excluded genetic variants with Hardy-Weinberg p-value < 10-6, missingness > 0.05 and minor allele frequency < 0.01 from further analysis. We also excluded samples with more than 5% of their genotypes missing.

#### Genotype imputation and quality control
Most of the datasets were imputed using the 1000 Genomes reference panel based on the GRCh38 genome version. CoLaus has been imputed using the TOPMed reference panel, while still aligning with the same reference genome version. Additionally, GEUVADIS, GTEx and MAGE cohorts utilise whole genome sequencing data aligned to the GRCh38 reference genome.

We pre-phased and imputed the microarray genotypes to the 1000 Genomes 30x on GRCh38 reference panel using Eagle v2.4.185 and Minimac4. We used bcftools v1.9.0 to exclude variants with minor allele frequency (MAF) < 0.01 and imputation quality score R2 < 0.4 from downstream analysis. The genotype imputation and quality control steps are implemented in eQTL-Catalogue/genimpute (v22.01.1) workflow available from GitHub. We used QCTOOL v2.2.0 (https://www.chg.ox.ac.uk/~gav/qctool_v2/) to convert imputed genotypes to bgen format for analysis with regenie. 
### Gene expression data
#### Studies
We used gene expression data from seven RNA-seq studies (TwinsUK, CoLaus, GEUVADIS, GENCORD, GTEx v8, CAP, MAGE) and three microarray studies (ALSPAC, MRCA and MRCE).

#### RNA-seq quantification and normalisation
Quantification of the RNA-seq data was performed using version v22.05.1 of the eQTL-Catalogue/rnaseq workflow implemented in Nextflow. Before quantification, we used Trim Galore v0.5.0 to remove sequencing adapters from the fastq files. For gene expression quantification, we used HISAT2 v2.2.1 to align reads to the GRCh38 reference genome (Homo_sapiens.GRCh38.dna.primary_assembly.fa file downloaded from Ensembl). We counted the number of reads overlapping the genes in the GENCODE V30 reference transcriptome annotations with featureCounts v1.6.4. 

We excluded all samples that failed the quality control steps. We normalised the gene counts using the conditional quantile normalisation (cqn) R package v1.30.0 with gene GC nucleotide content as a covariate. We downloaded the gene GC content estimates from Ensembl biomaRt and calculated the exon-level GC content using bedtools v2.19.088. We also excluded lowly expressed genes, where 95 per cent of the samples within a dataset had transcripts per million (TPM)-normalised expression less than 1. Subsequently, we used the inverse normal transformation to standardise quantification estimates. Normalisation scripts together with containerised software are publicly available at https://github.com/eQTL-Catalogue/qcnorm.

#### Microarray data processing
Gene expression from 877 individuals in the ALSPAC cohort was profiled using Illumina Human HT-12 V3 BeadChips microarray. We used the normalised gene expression matrix from the original publication31. In the MRCA cohort, gene expression from 327 individuals was profiled using the Human Genome U133 Plus 2.0 microarray. We downloaded the raw CEL files from ArrayExpress (E-MTAB-1425) and normalised the data using the Robust Multi-Array Average (RMA) method from the affy Bioconductor package. In the MRCE cohort, gene expression from 484 individuals was profiled using the Illumina Human-6 v1 Expression BeadChip. As raw data was unavailable, we downloaded the processed gene expression matrix from ArrayExpress (E-MTAB-1428). In all three microarray datasets, we applied inverse normal transformation to each probe before performing trans-eQTL analysis. If there were multiple probes mapping to the same gene, the probe with the highest average expression was used.
### Trans-eQTL mapping and meta-analysis
We performed independent quality control and normalisation on all datasets and only included 18,792 protein coding genes in the analysis. Trans-eQTL analysis was conducted separately on each dataset with regenie. For studies containing related samples (TwinsUK, MRCA and MRCE) and ALSPAC, both step 1 and step 2 commands were employed, while for other datasets with a smaller number of unrelated samples, regenie was run in the linear regression mode (step 2 only). We used sex and six principal components of the phenotype matrix and six principal components of genotypic data as covariates in the trans-eQTL analysis. All scripts used to run trans-eQTL are publicly available at https://github.com/freimannk/regenie_analysis. Subsequently, we performed an inverse-variance weighted meta-analysis across studies. Meta-analysis workflow is available at https://github.com/freimannk/regenie_metaanalyse.

We used a cis window of ± 5Mb to assign identified eQTLs into cis and trans eQTLs. To determine significant loci, we excluded variants in close proximity (±1.5Mb) to the most highly associated variant per gene. This approach allowed us to identify distinct and robust signals while mitigating potential confounding effects from nearby variants. By applying these filters we found 79 trans-eQTLs loci at a suggestive p-value threshold of 1x10⁻¹¹.











