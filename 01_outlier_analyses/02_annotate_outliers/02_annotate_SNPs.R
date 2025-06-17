# This script will annotate SNPs of interest.

# Prior to running this code, make sure the SnpEff vcf file from the previous step and the SNPs of interest are all in a directory called "Annotation"


# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

#--------------------------------------------------------------------------------

# Load packages
install.packages(c('data.table', 'tidyverse', 'foreach'))
library(data.table)
library(tidyverse)
library(foreach)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

#--------------------------------------------------------------------------------

# Load the VCF file
vcf.fn <- "data/processed/outlier_analyses/snpeff/N.canaliculata_pops_SNPs_annotate.vcf"
# Parse the header
seqVCF_Header(vcf.fn)
# Convert VCF to GDS
seqVCF2GDS(vcf.fn, storage.option="ZIP_RA", "data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

#--------------------------------------------------------------------------------

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")
# Display the contents of the GDS file
genofile

# Structure of GDS file:
#+    [  ] *
#|--+ description   [  ] *
#|--+ sample.id   { Str8 19 ZIP_ra(111.6%), 84B } *
#|--+ variant.id   { Int32 14897468 ZIP_ra(34.6%), 19.7M } *
#|--+ position   { Int32 14897468 ZIP_ra(49.2%), 27.9M } *
#|--+ chromosome   { Str8 14897468 ZIP_ra(0.22%), 420.7K } *
#|--+ allele   { Str8 14897468 ZIP_ra(17.0%), 9.7M } *
#|--+ genotype   [  ] *
#|  |--+ data   { Bit2 2x19x14900846 ZIP_ra(27.4%), 37.0M } *
#|  |--+ extra.index   { Int32 3x0 ZIP_ra, 16B } *
#|  \--+ extra   { Int16 0 ZIP_ra, 16B }
#|--+ phase   [  ]
#|  |--+ data   { Bit1 19x14897468 ZIP_ra(0.10%), 33.6K } *
#|  |--+ extra.index   { Int32 3x0 ZIP_ra, 16B } *
#|  \--+ extra   { Bit1 0 ZIP_ra, 16B }
#|--+ annotation   [  ]
#|  |--+ id   { Str8 14897468 ZIP_ra(0.10%), 14.2K } *
#|  |--+ qual   { Float32 14897468 ZIP_ra(89.9%), 51.1M } *
#|  |--+ filter   { Int32,factor 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |--+ info   [  ]
#|  |  |--+ NS   { Int32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ DP   { Int32 14897468 ZIP_ra(44.8%), 25.4M } *
#|  |  |--+ DPB   { Float32 14897468 ZIP_ra(46.3%), 26.3M } *
#|  |  |--+ AC   { Int32 15106624 ZIP_ra(20.3%), 11.7M } *
#|  |  |--+ AN   { Int32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ AF   { Float32 15106624 ZIP_ra(21.2%), 12.2M } *
#|  |  |--+ RO   { Int32 14897468 ZIP_ra(45.6%), 25.9M } *
#|  |  |--+ AO   { Int32 15106624 ZIP_ra(39.5%), 22.8M } *
#|  |  |--+ PRO   { Float32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ PAO   { Float32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ QR   { Int32 14897468 ZIP_ra(55.0%), 31.3M } *
#|  |  |--+ QA   { Int32 15106624 ZIP_ra(54.9%), 31.6M } *
#|  |  |--+ PQR   { Float32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ PQA   { Float32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ SRF   { Int32 14897468 ZIP_ra(42.3%), 24.1M } *
#|  |  |--+ SRR   { Int32 14897468 ZIP_ra(42.3%), 24.1M } *
#|  |  |--+ SAF   { Int32 15106624 ZIP_ra(35.4%), 20.4M } *
#|  |  |--+ SAR   { Int32 15106624 ZIP_ra(35.4%), 20.4M } *
#|  |  |--+ SRP   { Float32 14897468 ZIP_ra(86.9%), 49.4M } *
#|  |  |--+ SAP   { Float32 15106624 ZIP_ra(63.1%), 36.4M } *
#|  |  |--+ AB   { Float32 15106624 ZIP_ra(67.8%), 39.0M } *
#|  |  |--+ ABP   { Float32 15106624 ZIP_ra(72.4%), 41.7M } *
#|  |  |--+ RUN   { Int32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ RPP   { Float32 15106624 ZIP_ra(62.8%), 36.2M } *
#|  |  |--+ RPPR   { Float32 14897468 ZIP_ra(85.7%), 48.7M } *
#|  |  |--+ RPL   { Float32 15106624 ZIP_ra(36.9%), 21.3M } *
#|  |  |--+ RPR   { Float32 15106624 ZIP_ra(36.8%), 21.2M } *
#|  |  |--+ EPP   { Float32 15106624 ZIP_ra(60.7%), 35.0M } *
#|  |  |--+ EPPR   { Float32 14897468 ZIP_ra(84.6%), 48.1M } *
#|  |  |--+ DPRA   { Float32 15106624 ZIP_ra(80.4%), 46.3M } *
#|  |  |--+ ODDS   { Float32 14897468 ZIP_ra(91.2%), 51.8M } *
#|  |  |--+ GTI   { Int32 14897468 ZIP_ra(10.6%), 6.0M } *
#|  |  |--+ TYPE   { Str8 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ CIGAR   { Str8 15106624 ZIP_ra(0.10%), 43.1K } *
#|  |  |--+ NUMALT   { Int32 14897468 ZIP_ra(0.83%), 483.7K } *
#|  |  |--+ MEANALT   { Float32 15106624 ZIP_ra(6.48%), 3.7M } *
#|  |  |--+ LEN   { Int32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ MQM   { Float32 15106624 ZIP_ra(45.1%), 26.0M } *
#|  |  |--+ MQMR   { Float32 14897468 ZIP_ra(56.6%), 32.1M } *
#|  |  |--+ PAIRED   { Float32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ PAIREDR   { Float32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ MIN_DP   { Int32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ END   { Int32 14897468 ZIP_ra(0.10%), 56.6K } *
#|  |  |--+ technology.NovaSeq   { Float32 15106624 ZIP_ra(0.10%), 57.4K } *
#|  |  |--+ ANN   { Str8 18564191 ZIP_ra(6.67%), 124.4M } *
#|  |  |--+ LOF   { Str8 13623 ZIP_ra(17.0%), 51.0K } *
#|  |  \--+ NMD   { Str8 6724 ZIP_ra(17.5%), 25.9K } *
#|  \--+ format   [  ]
#|     |--+ GQ   [  ] *
#|     |  \--+ data   { Float32 19x0 ZIP_ra, 16B } *
#|     |--+ GL   [  ] *
#|     |  \--+ data   { Float32 19x45323250 ZIP_ra(48.8%), 1.6G } *
#|     |--+ DP   [  ] *
#|     |  \--+ data   { VL_Int 19x14897468 ZIP_ra(72.8%), 222.4M } *
#|     |--+ AD   [  ] *
#|     |  \--+ data   { VL_Int 19x30004092 ZIP_ra(62.3%), 357.0M } *
#|     |--+ RO   [  ] *
#|     |  \--+ data   { VL_Int 19x14897468 ZIP_ra(73.8%), 219.0M } *
#|     |--+ QR   [  ] *
#|     |  \--+ data   { VL_Int 19x14897468 ZIP_ra(62.4%), 322.4M } *
#|     |--+ AO   [  ] *
#|     |  \--+ data   { VL_Int 19x15106624 ZIP_ra(44.9%), 124.4M } *
#|     |--+ QA   [  ] *
#|     |  \--+ data   { VL_Int 19x15106624 ZIP_ra(47.7%), 177.5M } *
#|     \--+ MIN_DP   [  ] *
#|        \--+ data   { VL_Int 19x0 ZIP_ra, 16B } *
#\--+ sample.annotation   [  ]

#--------------------------------------------------------------------------------

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        id=seqGetData(genofile, "variant.id")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

#--------------------------------------------------------------------------------

# Load SNPs of interest
bonf.sig.SNPs <- read.csv("data/processed/outlier_analyses/baypass/Outlier_SNPs/Nucella_outlier_SNPs_bonferroni", header=T)
SNPs.Interest.pval.001 <- read.csv("data/processed/outlier_analyses/baypass/Outlier_SNPs/Nucella_outlier_SNPs_pval0.001", header=T)

# Create SNP ID column
snps <- bonf.sig.SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
snps <- SNPs.Interest.pval.001 %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# Check structure
str(snps)

#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest

annotation = 
foreach(i=1:dim(snps)[1], 
.combine = "rbind",
.errorhandling = "remove")%do%{

message(i)
seqResetFilter(genofile)

tmp.i = snps[i,]$SNP_id

pos.tmp = snp.dt %>% filter(SNP_id %in% tmp.i) %>% .$id

seqSetFilter(genofile, variant.id = pos.tmp)

ann_data <- seqGetData(genofile, "annotation/info/ANN")$data

L = length(ann_data)

annotate.list =
foreach(k=1:L, 
.combine = "rbind")%do%{

  tmp = ann_data[k] 
  tmp2= str_split(tmp, "\\|")
  
  data.frame(
    id=pos.tmp,
    SNP_id = tmp.i,
    annotation.id=k,
    Allele = tmp2[[1]][1],
    Annotation = tmp2[[1]][2],
    Annotation_Impact = tmp2[[1]][3],
    Gene_Name = tmp2[[1]][4],
    Gene_ID = tmp2[[1]][5],
    Feature_Type = tmp2[[1]][6],
    Feature_ID = tmp2[[1]][7],
    Transcript_BioType = tmp2[[1]][8],
    Rank = tmp2[[1]][9],
    HGVS.c = tmp2[[1]][10],
    HGVS.p = tmp2[[1]][11],
    cDNA.pos.cDNA.length = tmp2[[1]][12],
    CDS.pos.CDS.length = tmp2[[1]][13],
    AA.pos.AA.length = tmp2[[1]][14],
    Distance = tmp2[[1]][15]

  )
}
return(annotate.list)
}



#--------------------------------------------------------------------------------

# Join annotation and SNP information
N.canaliculata_annotated_SNPs <- left_join(annotation, snp.dt, by = join_by(SNP_id))

# Write output (bonferroni)
write.csv(N.canaliculata_annotated_SNPs, "data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs_bonferroni_annotated.txt")
write.csv(N.canaliculata_annotated_SNPs, "data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs_bonferroni_annotated.csv")

# Write output (window rnp SNPs)
write.csv(N.canaliculata_annotated_SNPs, "data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs_rnp_annotated.txt")
write.csv(N.canaliculata_annotated_SNPs, "data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs_rnp_annotated.csv")
