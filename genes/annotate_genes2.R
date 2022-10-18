library(data.table)

genes <- fread("./gencode19_genes.csv")
colnames(genes)
setnames(genes, c("seqname"), c("chr"))
genes <- QCtreeCNV:::chr_uniform(genes)
genes <- genes[chr %in% 1:22,]

# list all available biotypes
unique(genes$gene_type)

# filter by biotype if wanted, e.g. only "protein_coding"
# genes <- genes[gene_biotype %in% biotypes, ]

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
                            ens_genes = paste0(hits$gene_id, collapse = "-"),
                            genes = paste0(hits$gene_name, collapse = "-"))]
    }
  }
}

DT[, ix := NULL]

fwrite(DT, "30loci_annotated_all_genes.tsv", sep = "\t")

# load dosage data

dsg <- fread("../dosage_sens/Collins_rCNV_2022.dosage_sensitivity_scores.tsv")
setnames(dsg, "#gene", "gene")
dsg <- merge(genes, dsg, by.y = "gene", by.x = "gene_name")

DT <- fread("../pLI/data/30loci_pLI_LOEUF_annotated.tsv")
DT <- QCtreeCNV:::chr_uniform(DT)
setnames(DT, "stop", "end")
DT[, ix := 1:.N]
for (cc in unique(DT$chr)) {
  DT_tmp <- DT[chr == cc, ]
  genes_tmp <- dsg[chr == cc, ]
  for (i in 1:nrow(DT_tmp)) {
    st <- DT_tmp$start[i]
    en <- DT_tmp$end[i]
    index <- DT_tmp$ix[i]
    hits <- genes_tmp[between(start, st, en) | between(end, st, en) |
                        (start < st & end > en)]
    n_hits <- nrow(hits)
    if (n_hits > 0){
      DT[ix == index, `:=` (n_genes = n_hits,
                            pHaplosum = sum(hits$pHaplo),
                            pTriplosum = sum(hits$pTriplo),
                            ens_genes = paste0(hits$gene_id, collapse = "-"),
                            genes = paste0(hits$gene_name, collapse = "-"))]
    }
  }
}

fwrite(DT, "../dosage_sens/30loci_annotated_v2.tsv", sep = "\t")
