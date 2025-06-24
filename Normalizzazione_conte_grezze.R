# --- VERIFICA DELL'EFFICACIA DEL BATCH EFFECT ---
#BiocManager::install("DGEobj.utils")



# Caricamento dei pacchetti necessari
library(DESeq2)
library(DGEobj.utils)
library(biomaRt)
library(dplyr)
library(factoextra)
library(viridis)
library(vegan)
library(ComplexHeatmap)
# devtools::install_github("zhangyuqing/sva-devel")
library(sva)
library(SummarizedExperiment)
# --- CARICAMENTO DEI DATI ---

# Connettiti a Ensembl (esempio per umano)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Carica le informazioni sui batch
batch_info <- read.table("C:\\Users\\saran\\Downloads\\batch.txt", header = TRUE, sep = "\t")
# Carica i dati di espressione
rna <- read.delim("C:\\Users\\saran\\Downloads\\RawCounts123.txt", header = TRUE, sep = "\t")
dim(rna)#[1] 56185   124

# Assegna i nomi dei geni come rownames
rownames(rna) <- rna[, 1]
rna <- rna[, -1]

# Converte in matrice numerica
raw_counts <- as.matrix(rna)

# Rimuovi geni non espressi
raw_counts <- raw_counts[rowSums(raw_counts) > 0, ]

# Ora puoi usare rownames(raw_counts) per la mappatura con Ensembl
gene_symbols <- rownames(raw_counts)


dim(raw_counts)#[1] 39845   123
# --- PREPARAZIONE DEI DATI ---

# Converte in matrice numerica
Counts_matrix <- raw_counts

# Trasponi: campioni come righe, geni come colonne (per PCA, adonis2, kBET)
Counts_matrix_t <- t(Counts_matrix)
Counts_matrix_t <- Counts_matrix_t[, apply(Counts_matrix_t, 2, var) > 0]


batch <- as.factor(batch_info$batch)

# --- ANALISI PRIMA DELLA CORREZIONE ---

# PCA sui dati grezzi
pca_result_raw <- prcomp(Counts_matrix_t, scale. = TRUE)
pca_plot_raw <- fviz_pca_ind(pca_result_raw, label = "none", habillage = batch, addEllipses = TRUE) +
  labs(title = 'PCA - Raw Counts') +
  theme(plot.title = element_text(hjust = 0.5))
print(pca_plot_raw)

# Varianza spiegata dai batch (adonis2)
adonis_raw <- vegan::adonis2(Counts_matrix_t ~ batch, permutations = 999)
print("Varianza spiegata dai batch (RAW):")
print(adonis_raw$R2[1])
#[1] 0.04540064

# Distanza euclidea tra le medie dei batch 
batch_levels <- unique(batch)
batch_means_raw <- t(sapply(batch_levels, function(b) colMeans(Counts_matrix_t[batch == b, ])))
dist_euclidea_raw <- dist(batch_means_raw, method = "euclidean")
print("Distanza euclidea tra le medie dei batch (RAW):")
print(dist_euclidea_raw)
# 1071178

#Heatmap sui dati grezzi
Heatmap(cor(Counts_matrix), col = magma(100), name = "Correlation")

# --- CORREZIONE DEL BATCH EFFECT ---

# Prepara la matrice per ComBat_seq: geni come righe, campioni come colonne

dim(Counts_matrix) #[1] 39845   123

# Applica ComBat-seq
combat_rna <- sva::ComBat_seq(Counts_matrix, batch = batch, group = NULL)
dim(combat_rna)#[1] 39845   123
write.csv(combat_rna, file = "C:\\Users\\saran\\Downloads\\combat_rna_corrected_counts.csv", quote = FALSE)
# Trasponi di nuovo per PCA, adonis2, kBET, Heatmap
combat_rna_t <- t(combat_rna)

# --- ANALISI DOPO LA CORREZIONE ---
combat_rna_t <- combat_rna_t[, apply(combat_rna_t, 2, var) > 0]
# PCA sui dati corretti
pca_result_corr <- prcomp(combat_rna_t, scale. = TRUE)
pca_plot_corr <- fviz_pca_ind(pca_result_corr, label = "none", habillage = batch, addEllipses = TRUE) +
  labs(title = 'PCA - Batch Corrected') +
  theme(plot.title = element_text(hjust = 0.5))
print(pca_plot_corr)

# Varianza spiegata dai batch dopo la correzione
adonis_corr <- vegan::adonis2(combat_rna_t ~ batch, permutations = 999)
print("Varianza spiegata dai batch (CORRETTI):")
print(adonis_corr$R2[1])
#[1] 0.01696374

# Distanza euclidea tra le medie dei batch dopo la correzione
batch_means_corr <- t(sapply(batch_levels, function(b) colMeans(combat_rna_t[batch == b, ])))
dist_euclidea_corr <- dist(batch_means_corr, method = "euclidean")
print("Distanza euclidea tra le medie dei batch (CORRETTI):")
print(dist_euclidea_corr)
# 453790.3

#Heatmap sui dai batch dopo la correzione
Heatmap(cor(combat_rna), col = magma(100), name = "Correlation")

##----------------------------------------------------------------------------##
##                              Normalize counts                              ##
##----------------------------------------------------------------------------##

# do DESeq2 size factor normalization
sf <- estimateSizeFactorsForMatrix(combat_rna)
corces_rna_counts <- t( t(combat_rna) / sf )

# do +1 log2 transformation
rna_norm <- apply(corces_rna_counts + 1, 2, log2)

dim(rna_norm) #[1] 39845   123

##--------------------------SAVE NORMALIZED COUNTS----------------------------##
# Assicurati di NON avere oggetti chiamati "names" o "rownames"
rm(names, rownames)  # solo se li hai accidentalmente creati

# Converti la matrice in data frame
rna_norm_df <- as.data.frame(rna_norm)

# Recupera i nomi dei geni (da rownames) e aggiungili come prima colonna
rna_norm_df$gene <- rownames(rna_norm_df)

# Sposta la colonna 'gene' in prima posizione
rna_norm_df <- rna_norm_df[, c("gene", setdiff(colnames(rna_norm_df), "gene"))]

# Salva in CSV
write.csv(rna_norm_df,
          file = "C:\\Users\\saran\\Downloads\\rna_norm_with_genes.csv",
          quote = FALSE,
          row.names = FALSE)
