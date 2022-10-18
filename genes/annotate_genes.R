if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
devtools::install_github("sinomem/QCtreeCNV")

library(data.table)

# get the correct mart
mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

# extract all the genes into a DT
genes <- biomaRt::getBM(attributes=c("chromosome_name", "start_position",
                                     "end_position", "ensembl_gene_id",
                                     "gene_biotype", "external_gene_name",
                                     "external_gene_source"), mart = mart)
setDT(genes)
# keep only those in the autosomes
setnames(genes, c("chromosome_name", "start_position", "end_position"),
                c("chr", "start", "end"))
genes <- genes[chr %in% 1:22,]
# list all available biotypes
unique(genes$gene_biotype)

# filter by biotype if wanted, e.g. only "protein_coding"
# genes <- genes[gene_biotype %in% biotypes, ]

# save genes table
fwrite(genes, "genes_table.tsv", sep = "\t")

# load the 30 loci table
DT <- fread("../pLI/data/30loci_pLI_LOEUF_annotated.tsv")
DT <- QCtreeCNV:::chr_uniform(DT)
setnames(DT, "stop", "end")

# annotate the table
DT[, ix := 1:.N]
for (cc in unique(DT$chr)) {
  DT_tmp <- DT[chr == cc, ]
  genes_tmp <- genes[chr == cc, ]
  for (i in 1:nrow(DT_tmp)) {
    st <- DT_tmp$start[i]
    en <- DT_tmp$end[i]
    index <- DT_tmp$ix[i]
    hits <- genes_tmp[between(start, st, en) | between(end, st, en) |
                        (start < st & end > en)]
    n_hits <- nrow(hits)
    if (n_hits > 0){
      DT[ix == index, `:=` (n_genes = n_hits,
                            genes = paste0(hits$ensembl_gene_id,
                                           collapse = "-"))]
    }
  }
}

DT[, ix := NULL]



fwrite(DT, "30loci_annotated_all_genes.tsv", sep = "\t")
