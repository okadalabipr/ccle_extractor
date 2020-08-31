# Normalization

# get information of length

library(biomaRt)
library(dplyr)
library(stringr)

dataset <- "hsapiens_gene_ensembl"
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset)

listAttributes(mart)

gene.annotations <- biomaRt::getBM(
    mart=mart,
    attributes=c(
        "ensembl_gene_id", "start_position", "end_position"
    )
)
gene.annotations <- dplyr::transmute(
    gene.annotations, ensembl_gene_id, length=end_position-start_position
)   


# merge count file and length by ensembl gene id

CCLE_counts <- read.csv("CCLE_all_counts.csv")

TISSUE <- "BREAST"

CCLE_counts_selected <- data.frame(
    ensembl_gene_id=CCLE_counts$Name,
    Description=CCLE_counts$Description,
    CCLE_counts[,grep(TISSUE, colnames(CCLE_counts))]
)
CCLE_counts_selected[,1] <- str_sub(
    CCLE_counts_selected[,1], start=1, end=15
)
countsfile_length <- inner_join(
    gene.annotations,
    CCLE_counts_selected,
    by = "ensembl_gene_id"
)


# count TPM
library(edgeR)

counts_tpm <- (countsfile_length[,4:ncol(countsfile_length)]/countsfile_length$length)*1000
counts_tpm <- cbind(Description=countsfile_length$Description, counts_tpm) %>% as_tibble()

normfactor <- DGEList(
    counts=counts_tpm[,2:ncol(counts_tpm)],
    group=colnames(counts_tpm[,2:ncol(counts_tpm)])
)
normfactor <- calcNormFactors(normfactor, method="RLE")
normfactor_samples <- normfactor$samples
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors


for(i in 1:(dim(counts_tpm)[2]-1)){
    counts_tpm[,i+1] <- (counts_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
}
write.csv(counts_tpm, "tpm_RLE_normalized.csv")