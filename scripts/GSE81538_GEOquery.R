# # installing GEOquery
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

# Loading packages
library('tidyverse')
library('GEOquery')

# Clearing Env
rm(list = ls())


####################################################################################
# Creating Clinical Annotation Table
####################################################################################

# # Downloading SOFT file and expression files for GSE81538 (Run only if needed)
# system('bash ./scripts/download_GSE81538.sh')

# loading SOFT file
GSE81538_SOFT <- getGEO(filename = 'data/GSE81538/SOFT/GSE81538_family.soft.gz' )
# Metadata
Meta(GSE81538_SOFT)
# Vector of GSM IDs
Sample_Names <- Meta(GSE81538_SOFT)["sample_id"][[1]]
# Initialising Dataframe
GSE81538_Clinical_DF <- NULL

# Running loop to extract clinical data from each sample
for (GSM_ID in Sample_Names) {
  # GEO sample ID
  Meta_GEO_Ac <- Meta(GSMList(GSE81538_SOFT)[[GSM_ID]])["geo_accession"]
  # Title ID
  Meta_Title_ID <- Meta(GSMList(GSE81538_SOFT)[[GSM_ID]])["title"]
  # Platform ID Info for Sample 1
  Meta_Plat_ID <- Meta(GSMList(GSE81538_SOFT)[[GSM_ID]])["platform_id"]
  # Sample Tissue
  Meta_Tissue <- Meta(GSMList(GSE81538_SOFT)[[GSM_ID]])["source_name_ch1"]
  
  # Clinical Info for Sample 1
  Meta_Clin <- Meta(GSMList(GSE81538_SOFT)[[GSM_ID]])["characteristics_ch1"] %>% as.data.frame()
  Meta_Clin$characteristics_ch1 <- Meta_Clin$characteristics_ch1 %>% as.character()
  Meta_Clin <- Meta_Clin %>% separate(characteristics_ch1, c("Characteristic","Value"), sep = ": ")
  Meta_Clin$Characteristic <- toupper(Meta_Clin$Characteristic) %>% gsub(" ", "_", .)
  # Building Vectors of Clinical Variables and Values
  Meta_Clin_Labels <- c("GSM_ID", "TITLE_ID", "PLATFORM", "TISSUE", Meta_Clin$Characteristic)
  Meta_Clin_Value <- c(Meta_GEO_Ac, Meta_Title_ID, Meta_Plat_ID, Meta_Tissue, Meta_Clin$Value) %>% as.character()
  # Combining vectors into DF
  Meta_Clin_DF <- data.frame(Meta_Clin_Value)
  
  # Adding to dataframe
  GSE81538_Clinical_DF <- rbind(GSE81538_Clinical_DF, Meta_Clin_Value)
  colnames(GSE81538_Clinical_DF) <- Meta_Clin_Labels
}

# Converting GSE81538_Clinical_DF Matrix to DF
GSE81538_Clinical_DF <- as.data.frame(GSE81538_Clinical_DF) %>% 
  # Converting Character NA to true NA
  na_if(., "NA")


# Remove NA from factors 
for (var in colnames(GSE81538_Clinical_DF)) {
  var_levels <- levels(GSE81538_Clinical_DF[, var])
  var_levels <- var_levels[!(var_levels %in% "NA")]
  GSE81538_Clinical_DF[, var] <- factor(GSE81538_Clinical_DF[, var], var_levels)
  
}

# Defining data type for variables
Integer_Vars  <- c('ER_IHC_ROUTINESTAIN_CLINICAL_READING',
                    'ER_IHC_ROUTINESTAIN_READING_2',
                    'ER_IHC_ROUTINESTAIN_READING_3',
                    'ER_IHC_RESTAIN_READING_1',
                    'ER_IHC_RESTAIN_READING_2',
                    'ER_IHC_RESTAIN_READING_3',
                    'ER_CONSENSUS',
                    'PGR_IHC_ROUTINESTAIN_CLINICAL_READING',
                    'PGR_ROUTINESTAIN_READING_2',
                    'PGR_ROUTINESTAIN_READING_3',
                    'PGR_RESTAIN_READING_1',
                    'PGR_RESTAIN_READING_2',
                    'PGR_RESTAIN_READING_3',
                    'PGR_CONSENSUS',
                    'HER2_IHC_CLINREADING',
                    'HER2_IHC_READING_2',
                    'HER2_IHC_READING_3',               
                    'HER2_CLINICAL_STATUS',
                    'HER2_SISH_READING_1',
                    'HER2_SISH_READING_2',
                    'HER2_SISH_READING_3',
                    'HER2_CONSENSUS',
                    'NHG_READING_1_TUBULAR',
                    'NHG_READING_1_PLEOMORPH',
                    'NHG_READING_1_MITOTIC',
                    'NHG_READING_2_TUBULAR',
                    'NHG_READING_2_PLEOMORPH',
                    'NHG_READING_2_MITOTIC',
                    'NHG_READING_3_TUBULAR',
                    'NHG_READING_3_PLEOMORPH',
                    'NHG_READING_3_MITOTIC',
                    'NHG_CONSENSUS')

for (Integer in Integer_Vars) {
  GSE81538_Clinical_DF[, Integer] <- as.integer(GSE81538_Clinical_DF[, Integer])
}



Numeric_Vars <-  c("KI67_READING_1", "KI67_READING_2", "KI67_READING_3", "KI67_CONSENSUS")
for (Numeric in Numeric_Vars) {
  GSE81538_Clinical_DF[, Numeric] <- as.numeric(GSE81538_Clinical_DF[, Numeric])
}

# Checking Final DF and saving
system('mkdir -p data/GSE81538/clean')
glimpse(GSE81538_Clinical_DF)
write.csv(GSE81538_Clinical_DF, file = "data/GSE81538/clean/GSE81538_Clinical_Annotation.csv", row.names = F)



####################################################################################
# Creating Gene Expression DF
####################################################################################
# Reading Gene expression file (values are FPKM)
GSE81538_GE <- read_csv("data/GSE81538/GSE81538_gene_expression_405_transformed.csv.gz")
# Checking that Tile ID from Clinical DF matches gene expression
colnames(GSE81538_GE[, 2:ncol(GSE81538_GE)]) %in% GSE81538_Clinical_DF$TITLE_ID %>% 
  table() %>% 
  print()

# Renaming First Column
colnames(GSE81538_GE)[1] <- "GENE"
glimpse(GSE81538_GE[, 1:10])

# Saving DF
write.csv(GSE81538_GE, file = "data/GSE81538/clean/GSE81538_Gene_Expression.csv", row.names = F)

####################################################################################
# Creating Transcript Expression DF
####################################################################################
# Reading Gene expression file (values are FPKM)
GSE81538_TE <- read_csv("data/GSE81538/GSE81538_transcript_expression_405.csv.gz")
glimpse(GSE81538_TE[, 1:10])
# Checking that Tile ID from Clinical DF matches gene expression
colnames(GSE81538_TE[, 2:ncol(GSE81538_TE)]) %in% GSE81538_Clinical_DF$TITLE_ID %>% 
  table() %>% 
  print()

# Renaming First Column
colnames(GSE81538_TE)[1] <- "TRANSCRIPT"
glimpse(GSE81538_TE[, 1:10])

# Saving DF
write.csv(GSE81538_TE, file = "data/GSE81538/clean/GSE81538_Transcript_Expression.csv", row.names = F)

#################
# Writing R.data 
#################
save(GSE81538_TE, GSE81538_GE, GSE81538_Clinical_DF, file = "data/GSE81538/clean/GSE81538_Cleaned.RData")
