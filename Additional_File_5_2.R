### Author: Dr Ricardo De Paoli-Iseppi. Creation date: 26/02/2021 ###
### Short script to count risk gene data from curated GWAS Catalog data ###
# Requires: GWAS Catalog of MHD of interest (.csv) and validation gene lists (.csv) e.g. MAGMA, TWAS, SMR (optional). 
# Validation evidence lists require minimum colnames 'PUBMEDID' and 'GENE'

library(tidyr)
library(dplyr)
library(stringr)

# Change csv file path below to most recent MHD filtered GWAS catalog.
setwd("<path to directory here>/Risk_gene_database_curation")
dt <- read.csv(file = 'BD_GWASCatalog_filtered_210226.csv', stringsAsFactors = FALSE)
#dt <- read.csv(file = 'ASD_GWASCatalog_filtered_210301.csv', stringsAsFactors = FALSE)
#dt <- read.csv(file = 'SZ_GWASCatalog_filtered_210304.csv', stringsAsFactors = FALSE)
#dt <- read.csv(file = 'MDD_GWASCatalog_filtered_210317.csv', stringsAsFactors = FALSE)
# Check.
head(dt)
str(dt)

# Separate rows on REPORTED.GENE and MAPPED_GENE where multiple genes are reported for each SNP. Multiple separators.
s_rg_dt <- dt %>% separate_rows(REPORTED.GENE, convert = TRUE)
s_mg_dt <- dt %>% separate_rows(MAPPED_GENE, convert = TRUE)

# Use group_by to summarise column of interest.
report_gene <- group_by(s_rg_dt, REPORTED.GENE)
mapped_gene <- group_by(s_mg_dt, MAPPED_GENE)
# Get gene count counts and associated appearance in studies, based on PUBMEDID.
rg_list <- summarize(report_gene, count_rg = n(), study_rg = n_distinct(PUBMEDID))
# Rename column where names is "REPORTED.GENE"
names(rg_list)[names(rg_list) == "REPORTED.GENE"] <- "GENE"
mg_list <- summarize(mapped_gene, count_mg = n(), study_mg = n_distinct(PUBMEDID))
# Rename column where names is "MAPPED_GENE"
names(mg_list)[names(mg_list) == "MAPPED_GENE"] <- "GENE"

# Combine tables for unique gene ID and reported/mapped gene counts. Keep all rows/non-matches.
comb_list <- full_join(rg_list, mg_list)
# Sort combined list by reported study count then by reported gene count.
sort_comb_list <- comb_list[order(-comb_list$study_rg, -comb_list$count_rg),]

# Write out to .csv
write.csv(sort_comb_list,"<path to directory here>/Risk_gene_database_curation/sorted_GWAS_counts.csv", row.names = FALSE)


## STEP 2: Inclusion of further evidence from MAGMA, TWAS and SMR + more ##
# Read in separate evidence files to dfs.
setwd("<path to directory here>/Risk_gene_database_curation/ASD")
#setwd("<path to directory here>/Risk_gene_database_curation/BD")
#setwd("<path to directory here>/Risk_gene_database_curation/SZ")
#setwd("<path to directory here>/Risk_gene_database_curation/MDD")

dt_MAGMA <- read.csv(file = 'MAGMA_evidence_MDD.csv', stringsAsFactors = FALSE)
dt_TWAS <- read.csv(file = 'TWAS_evidence_MDD.csv', stringsAsFactors = FALSE)
dt_SMR <- read.csv(file = 'SMR_evidence_MDD.csv', stringsAsFactors = FALSE)
#dt_DNAM <- read.csv(file = 'DNAm_evidence_SZ.csv', stringsAsFactors = FALSE)             #SZ specific
#dt_FINEMAP <- read.csv(file = 'FINEMAP_evidence_SZ.csv', stringsAsFactors = FALSE)       #SZ specific
dt_PPI <- read.csv(file = 'PPI_evidence_MDD.csv', stringsAsFactors = FALSE)               #SZ,MDD specific
dt_VALID <- read.csv(file = 'VALIDATION_evidence_MDD.csv', stringsAsFactors = FALSE)      #SZ,MDD specific
setwd("<path to directory here>/Risk_gene_database_curation")

# Get gene counts for each set of evidence. 
# MAGMA
gene_MAGMA <- group_by(dt_MAGMA, GENE)
MAGMA_list <- summarize(gene_MAGMA, count_magmag = n(), study_magmag = n_distinct(PUBMEDID))
# TWAS
gene_TWAS <- group_by(dt_TWAS, GENE)
TWAS_list <- summarize(gene_TWAS, count_twasg = n(), study_twasg = n_distinct(PUBMEDID))
# SMR
gene_SMR <- group_by(dt_SMR, GENE)
SMR_list <- summarize(gene_SMR, count_smrg = n(), study_smrg = n_distinct(PUBMEDID))
# DNAM (SZ specific)
gene_DNAM <- group_by(dt_DNAM, GENE)
DNAM_list <- summarize(gene_DNAM, count_dnamg = n(), study_dnamg = n_distinct(PUBMEDID))
# FINEMAP (SZ specific)
gene_FINEMAP <- group_by(dt_FINEMAP, GENE)
FINEMAP_list <- summarize(gene_FINEMAP, count_finemapg = n(), study_finemapg = n_distinct(PUBMEDID))
# PPI (SZ specific)
gene_PPI <- group_by(dt_PPI, GENE)
PPI_list <- summarize(gene_PPI, count_ppig = n(), study_ppig = n_distinct(PUBMEDID))
# VALID (SZ specific)
gene_VALID <- group_by(dt_VALID, GENE)
VALID_list <- summarize(gene_VALID, count_validg = n(), study_validg = n_distinct(PUBMEDID))


# Combine validation evidence lists with combined GWAS list created earlier. Clear data frame.
valid_list <- data.frame()
valid_list <- full_join(comb_list, MAGMA_list)
valid_list <- full_join(valid_list, TWAS_list)
valid_list <- full_join(valid_list, SMR_list)
#valid_list <- full_join(valid_list, DNAM_list)      #SZ specific
#valid_list <- full_join(valid_list, FINEMAP_list)   #SZ specific
valid_list <- full_join(valid_list, PPI_list)       #SZ,MDD specific
valid_list <- full_join(valid_list, VALID_list)     #SZ,MDD specific


# Count if gene has a value for each study/experiment type. In baseR.
rowLimit <- nrow(valid_list)
evidence_df <- data.frame(evidence=integer())
# Include specific columns to check (unique to each MHD depending on evidence).
for(i in 1:rowLimit){
  temp <- sum(!is.na(valid_list[i,c(2,4,6,8,10,12,14)]))
  evidence_df[i,] <- temp
}
valid_list <- cbind(valid_list,evidence_df)

# Sort by evidence across studies.
evidence_sort_list <- valid_list[order(-valid_list$evidence),]
# Write out to .csv
write.csv(evidence_sort_list,"<path to directory here>/Risk_gene_database_curation/sorted_evidence_counts.csv", row.names = FALSE)


## STEP 3: Combine evidence for risk genes across all disorders into a new data-set ##
# Read back-in generated evidence lists and prefix evidence column with MHD
SZ_evid_list <- read.csv(file = 'SZ/SZ_sorted_evidence_counts_210325.csv', stringsAsFactors = FALSE)
ASD_evid_list <- read.csv(file = 'ASD/ASD_sorted_evidence_counts_210304.csv', stringsAsFactors = FALSE)
BD_evid_list <- read.csv(file = 'BD/BD_sorted_evidence_counts_210303.csv', stringsAsFactors = FALSE)
MDD_evid_list <- read.csv(file = 'MDD/MDD_sorted_evidence_counts_210325.csv', stringsAsFactors = FALSE)

names(SZ_evid_list)[20] <- "SZ_evidence"
names(ASD_evid_list)[12] <- "ASD_evidence"
names(BD_evid_list)[12] <- "BD_evidence"
names(MDD_evid_list)[16] <- "MDD_evidence"

# Combine each MHD df on gene and then add up evidence for each to give total evidence for genes across each MHD. 
multi_list <- SZ_evid_list %>%
  select(GENE, SZ_evidence)
ASD_sub <- ASD_evid_list %>%
  select(GENE, ASD_evidence)
BD_sub <- BD_evid_list %>%
  select(GENE, BD_evidence)
MDD_sub <- MDD_evid_list %>%
  select(GENE, MDD_evidence)
multi_list <- full_join(multi_list, ASD_sub) # Already has SZ evidence column.
multi_list <- full_join(multi_list, BD_sub)
multi_list <- full_join(multi_list, MDD_sub)

rowLimit <- nrow(multi_list)
multi_MHD_evidence <- data.frame(multi_trait_evidence=integer())
for(i in 1:rowLimit){
  temp <- rowSums(multi_list[i,c(2,3,4,5)], na.rm=T)
  multi_MHD_evidence[i,] <- temp
}
multi_list <- cbind(multi_list,multi_MHD_evidence) 
# Sort by evidence across studies.
evidence_sort_multi_list <- multi_list[order(-multi_list$multi_trait_evidence),]
# Write out to .csv
write.csv(evidence_sort_multi_list,"<path to directory here>/Risk_gene_database_curation/sorted_evidence_multi_counts.csv", row.names = FALSE)

##### END #####