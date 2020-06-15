# # installing GEOquery
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

# Loading packages
library('tidyverse')
library('GEOquery')


####################################################################################
# Creating Clinical Annotation Table
####################################################################################

# # Downloading SOFT file and expression files for GSE96058 (Run only if needed)
# system('bash ./scripts/download_GSE96058.sh')

# loading SOFT file
GSE96058_SOFT <- getGEO(filename = 'data/GSE96058/SOFT/GSE96058_family.soft.gz' )
# Metadata
Meta(GSE96058_SOFT)
# Vector of GSM IDs
Sample_Names <- Meta(GSE96058_SOFT)["sample_id"][[1]]
# Initialising Dataframe
GSE96058_Clinical_DF <- NULL

# Running loop to extract clinical data from each sample
for (GSM_ID in Sample_Names) {
  # GEO sample ID
  Meta_GEO_Ac <- Meta(GSMList(GSE96058_SOFT)[[GSM_ID]])["geo_accession"]
  # Title ID
  Meta_Title_ID <- Meta(GSMList(GSE96058_SOFT)[[GSM_ID]])["title"]
  # Platform ID Info for Sample 1
  Meta_Plat_ID <- Meta(GSMList(GSE96058_SOFT)[[GSM_ID]])["platform_id"]
  # Sample Tissue
  Meta_Tissue <- Meta(GSMList(GSE96058_SOFT)[[GSM_ID]])["source_name_ch1"]
  
  # Clinical Info for Sample 1
  Meta_Clin <- Meta(GSMList(GSE96058_SOFT)[[GSM_ID]])["characteristics_ch1"] %>% as.data.frame()
  Meta_Clin$characteristics_ch1 <- Meta_Clin$characteristics_ch1 %>% as.character()
  Meta_Clin <- Meta_Clin %>% separate(characteristics_ch1, c("Characteristic","Value"), sep = ": ")
  Meta_Clin$Characteristic <- toupper(Meta_Clin$Characteristic) %>% gsub(" ", "_", .)
  # Building Vectors of Clinical Variables and Values
  Meta_Clin_Labels <- c("GSM_ID", "TITLE_ID", "PLATFORM", "TISSUE", Meta_Clin$Characteristic)
  Meta_Clin_Value <- c(Meta_GEO_Ac, Meta_Title_ID, Meta_Plat_ID, Meta_Tissue, Meta_Clin$Value) %>% as.character()
  # Combining vectors into DF
  Meta_Clin_DF <- data.frame(Meta_Clin_Value)
  
  # Adding to dataframe
  GSE96058_Clinical_DF <- rbind(GSE96058_Clinical_DF, Meta_Clin_Value)
  colnames(GSE96058_Clinical_DF) <- Meta_Clin_Labels
}

# Converting GSE96058_Clinical_DF Matrix to DF
GSE96058_Clinical_DF <- as.data.frame(GSE96058_Clinical_DF) %>% 
  # Converting Character NA to true NA
  na_if(., "NA")


# Remove NA from factors 
for (var in colnames(GSE96058_Clinical_DF)) {
  var_levels <- levels(GSE96058_Clinical_DF[, var])
  var_levels <- var_levels[!(var_levels %in% "NA")]
  GSE96058_Clinical_DF[, var] <- factor(GSE96058_Clinical_DF[, var], var_levels)

}

# Defining data type for variables
Integer_Vars  <- c("ER_STATUS", "PGR_STATUS", "HER2_STATUS", "KI67_STATUS", "ER_PREDICTION_MGC", "PGR_PREDICTION_MGC", 
                   "HER2_PREDICTION_MGC", "KI67_PREDICTION_MGC", "ER_PREDICTION_SGC", "PGR_PREDICTION_SGC", "HER2_PREDICTION_SGC", 
                   "KI67_PREDICTION_SGC", "OVERALL_SURVIVAL_EVENT", "ENDOCRINE_TREATED", "CHEMO_TREATED")
for (Integer in Integer_Vars) {
  GSE96058_Clinical_DF[, Integer] <-  GSE96058_Clinical_DF[, Integer] %>% as.character() %>%  as.integer()
}



Numeric_Vars <-  c("AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "OVERALL_SURVIVAL_DAYS")
for (Numeric in Numeric_Vars) {
  GSE96058_Clinical_DF[, Numeric] <- GSE96058_Clinical_DF[, Numeric] %>% as.character() %>%  as.numeric()
}



# Converting receptor status to factored character
receptor_status_vars  <- c("ER_STATUS", "PGR_STATUS", "HER2_STATUS", "KI67_STATUS", "ER_PREDICTION_MGC", "PGR_PREDICTION_MGC", 
                           "HER2_PREDICTION_MGC", "KI67_PREDICTION_MGC", "ER_PREDICTION_SGC", "PGR_PREDICTION_SGC", "HER2_PREDICTION_SGC", 
                           "KI67_PREDICTION_SGC")
for (var in receptor_status_vars) {
  GSE96058_Clinical_DF[, var] <- ifelse(GSE96058_Clinical_DF[, var] == 0, "Negative", 
                                        ifelse(GSE96058_Clinical_DF[, var] == 1, "Positive", "NULL")) %>% 
                                 factor(levels = c("Negative", "Positive"))
}

# Converting overall survival in months variable
GSE96058_Clinical_DF$OVERALL_SURVIVAL_MONTHS <- (GSE96058_Clinical_DF$OVERALL_SURVIVAL_DAYS / 30.4167) %>% round(., digits = 1)

# Creating TNBC column in clinical DF
GSE96058_Clinical_DF$TNBC_STATUS <- ifelse(GSE96058_Clinical_DF$ER_STATUS == "Negative" & GSE96058_Clinical_DF$PGR_STATUS == "Negative" & GSE96058_Clinical_DF$HER2_STATUS == "Negative", "TNBC", "Non_TNBC")
GSE96058_Clinical_DF$TNBC_STATUS <- factor(GSE96058_Clinical_DF$TNBC_STATUS, c("TNBC", "Non_TNBC"))


# Creating folder of variables and levels
dir.create(path = "./data/GSE96058/clean/clinical_variable_factors", recursive = T, showWarnings = F)
for (clin_var in colnames(GSE96058_Clinical_DF)) {
  
  if (class(GSE96058_Clinical_DF[, clin_var]) == "factor") {
    var_levels <- levels(GSE96058_Clinical_DF[, clin_var])
    write(var_levels, file = paste("./data/GSE96058/clean/clinical_variable_factors/", clin_var, ".txt", sep = ""))
  }
}

# Checking Final DF and saving
system('mkdir -p data/GSE96058/clean')
glimpse(GSE96058_Clinical_DF)
write.csv(GSE96058_Clinical_DF, file = "data/GSE96058/clean/GSE96058_Clinical_Annotation.csv", row.names = F)



####################################################################################
# Creating Gene Expression DF
####################################################################################
# Reading Gene expression file (values are FPKM)
GSE96058_GE <- read_csv("data/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
# Checking that Tile ID from Clinical DF matches gene expression
colnames(GSE96058_GE[, 2:ncol(GSE96058_GE)]) %in% GSE96058_Clinical_DF$TITLE_ID %>% 
  table() %>% 
  print()

# Renaming First Column
colnames(GSE96058_GE)[1] <- "GENE"
glimpse(GSE96058_GE[, 1:10])

# Saving DF
write.csv(GSE96058_GE, file = "data/GSE96058/clean/GSE96058_Gene_Expression.csv", row.names = F)

####################################################################################
# Creating Transcript Expression DF
####################################################################################
# Reading Gene expression file (values are FPKM)
GSE96058_TE <- read_csv("data/GSE96058/GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz")
glimpse(GSE96058_TE[, 1:10])
# Checking that Tile ID from Clinical DF matches gene expression
colnames(GSE96058_TE[, 2:ncol(GSE96058_TE)]) %in% GSE96058_Clinical_DF$TITLE_ID %>% 
  table() %>% 
  print()

# Renaming First Column
colnames(GSE96058_TE)[1] <- "TRANSCRIPT"
glimpse(GSE96058_TE[, 1:10])

# Saving DF
write.csv(GSE96058_TE, file = "data/GSE96058/clean/GSE96058_Transcript_Expression.csv", row.names = F)


#################
# Writing R.data 
#################
save(GSE96058_TE, GSE96058_GE, GSE96058_Clinical_DF, file = "data/GSE96058/clean/GSE96058_Cleaned.RData")

