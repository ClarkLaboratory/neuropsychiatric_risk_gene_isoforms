# Script for combing all risk gene TPM counts and performing PCA (Supplementary Figures). Written by: J Gleeson (2023)

# Directory containing IsoLamp output files for each gene ending in '_TPM_values.csv'
setwd("<path to directory here>/PCA_data/")

library(gghalves)
library(ggplot2)
library(ggfortify)
library(data.table)
library(tidyr)
library(dplyr)
library(magrittr)

import_counts <- function(sample) {
  df <- fread(sample)
  dfm <- melt(df)
  dfm <- dfm %>% dplyr::mutate_at(c("variable"), .funs=toupper)
  dfm <- dfm %>% 
    separate(variable, sep="_", remove=F, c("individual", "region")) %>% 
    dplyr::rename("TPM" = "value")
  return(dfm)
}

# get sample paths from file list
samples <- list.files(pattern="_TPM_values.csv", recursive=T)

# create list of dfs for every sample
all_samples <- lapply(samples, import_counts)

# Merge into df
combined_all <- do.call("rbind", all_samples)

# need to create a unique id for each transcript and gene
combined_all$uniq_id <- paste0(combined_all$transcript_id, "_", combined_all$gene_id)

combined_all <- drop_na(combined_all)

# some brain regions weren't capitalised
combined_all <- combined_all %>% dplyr::mutate_at(c("region", "variable"), .funs=toupper)

# create df to plot with
df <- combined_all

dfmat <- dcast(df, uniq_id ~ variable, value.var = "TPM", fun.aggregate=sum)

# below to filter for novels only
#dfmat2 <- dfmat %>% dplyr::filter(startsWith(uniq_id, "tx"))

row.names(dfmat) <- dfmat$uniq_id
dfmat$uniq_id <- NULL

metadata <- read.csv("<path to directory here>/PCA_scripts/metadata_file_ISO.csv")
metadata <- metadata %>% dplyr::mutate_at(c("ids", "ind", "region"), .funs=toupper)

#PCA
pca_res <- prcomp(t(dfmat), center=T, scale=T)
pca_res_df <- as.data.frame(pca_res$x)
summary(pca_res)

metadata_ordered <- metadata[match(row.names(pca_res_df), metadata$ids),]

pca_res_metadata <- cbind(metadata_ordered,pca_res_df)

# plots
cpal <- c("#C82D76", "#615B9D", "#D49C2E","#BE501E", "#3A8C68", "#A7A7A7", "#8C6525")

#Plot PCA with RIN colour
pdf("PCs_pc1_pc3_RINe.pdf", width = 8, height = 6)
p <- ggplot(pca_res_metadata, aes(x=PC1, y=PC3, color=RIN, shape=ind)) + 
  geom_point(size=2) +
  scale_colour_gradient2(low = "#C82D76", mid = "#ffcc89", high = "#3A8C68", midpoint = 7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())
#xlab("PC1 27.8%") +
#ylab("PC2 11.9%") +
#scale_colour_manual(values=cpal)
p
dev.off()

#Plot PCA with age colour
pdf("PCs_pc4_pc5_age.pdf", width = 8, height = 6)
p <- ggplot(pca_res_metadata, aes(x=PC4, y=PC5, color=age, shape=ind)) + 
  geom_point(size=1.9) +
  #scale_colour_gradient2(low = "#C82D76", mid = "#ffcc89", high = "#3A8C68", midpoint = 7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())
#xlab("PC1 27.8%") +
#ylab("PC2 11.9%") +
#scale_colour_manual(values=cpal)
p
dev.off()

#Plot PCA with brain colour
pdf("PCs_pc1_pc2_brain.pdf", width = 8, height = 6)
p <- ggplot(pca_res_metadata, aes(x=PC1, y=PC2, color=region, shape=ind)) + 
  geom_point(size=2) +
  #scale_colour_gradient2(low = "#C82D76", mid = "#ffcc89", high = "#3A8C68", midpoint = 7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) 
#xlab("PC1 27.8%") +
#ylab("PC2 11.9%") +
#scale_colour_manual(values=cpal)
p
dev.off()

#Plot PCA with PMI colour
pdf("PCs_pc2_pc3_pmi.pdf", width = 8, height = 6)
p <- ggplot(pca_res_metadata, aes(x=PC2, y=PC3, color=pmi, shape=ind)) + 
  geom_point(size=1.9) +
  #scale_colour_gradient2(low = "#C82D76", mid = "#ffcc89", high = "#3A8C68", midpoint = 7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) 
#xlab("PC1 27.8%") +
#ylab("PC2 11.9%") +
#scale_colour_manual(values=cpal)
p
dev.off()
