# reading in clinical matrix file
matrix_clin_df <- read.table(file = "data/GSE96058/GSE96058-GPL11154_series_matrix.txt", header = T)
matrix_clin_gse = GEOquery::getGEO(filename="data/GSE96058/GSE96058-GPL11154_series_matrix.txt")
View(matrix_clin_gse)


# Defining data type for variables
Integer_Vars  <- c("ER_STATUS", "PGR_STATUS", "HER2_STATUS", "KI67_STATUS", "ER_PREDICTION_MGC", "PGR_PREDICTION_MGC", 
                   "HER2_PREDICTION_MGC", "KI67_PREDICTION_MGC", "ER_PREDICTION_SGC", "PGR_PREDICTION_SGC", "HER2_PREDICTION_SGC", 
                   "KI67_PREDICTION_SGC", "OVERALL_SURVIVAL_EVENT", "ENDOCRINE_TREATED", "CHEMO_TREATED")
for (Integer in Integer_Vars) {
  GSE96058_Clinical_DF[, Integer] <- as.integer(GSE96058_Clinical_DF[, Integer])
}



Numeric_Vars <-  c("AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "OVERALL_SURVIVAL_DAYS")
for (Numeric in Numeric_Vars) {
  GSE96058_Clinical_DF[, Numeric] <- as.numeric(GSE96058_Clinical_DF[, Numeric])
}



# Converting clinical variables values (0 = neg 1 = pos)
for (Integer in Integer_Vars) {
  GSE96058_Clinical_DF[, Integer] <- (GSE96058_Clinical_DF[, Integer] - 1) %>% as.integer()
}

# Converting receptor status to factored character
matrix_clin_df  <- data.frame("GSM_ID" = matrix_clin_gse$geo_accession, 
                              "TITLE_ID" = matrix_clin_gse$title,
                              "ER_STATUS" = matrix_clin_gse$`er status:ch1`, 
                              "PGR_STATUS" = matrix_clin_gse$`pgr status:ch1`, 
                              "HER2_STATUS" = matrix_clin_gse$`her2 status:ch1`, 
                              "KI67_STATUS" = matrix_clin_gse$`ki67 status:ch1`, 
                              "ER_PREDICTION_MGC" = matrix_clin_gse$`er prediction mgc:ch1`, 
                              "PGR_PREDICTION_MGC" = matrix_clin_gse$`pgr prediction mgc:ch1`, 
                              "HER2_PREDICTION_MGC" = matrix_clin_gse$`her2 prediction mgc:ch1`, 
                              "KI67_PREDICTION_MGC" = matrix_clin_gse$`ki67 prediction mgc:ch1`, 
                              "ER_PREDICTION_SGC" = matrix_clin_gse$`er prediction sgc:ch1`, 
                              "PGR_PREDICTION_SGC" = matrix_clin_gse$`pgr prediction sgc:ch1`, 
                              "HER2_PREDICTION_SGC" = matrix_clin_gse$`her2 prediction sgc:ch1`, 
                              "KI67_PREDICTION_SGC" = matrix_clin_gse$`ki67 prediction sgc:ch1`, 
                              "OVERALL_SURVIVAL_DAYS" = matrix_clin_gse$`overall survival days:ch1`, 
                              "OVERALL_SURVIVAL_EVENT" = matrix_clin_gse$`overall survival event:ch1`, 
                              "ENDOCRINE_TREATED" = matrix_clin_gse$`endocrine treated:ch1`, 
                              "CHEMO_TREATED" = matrix_clin_gse$`chemo treated:ch1`, 
                              "PAM50_SUBTYPE" = matrix_clin_gse$`pam50 subtype:ch1`,
                              "AGE_AT_DIAGNOSIS" = matrix_clin_gse$`age at diagnosis:ch1`, 
                              "TUMOR_SIZE" = matrix_clin_gse$`tumor size:ch1`, 
                              "NODE_STATUS" = matrix_clin_gse$`lymph node status:ch1`, 
                              "NHG" = matrix_clin_gse$`nhg:ch1`)


Integer_Vars  <- c("ER_STATUS", "PGR_STATUS", "HER2_STATUS", "KI67_STATUS", "ER_PREDICTION_MGC", "PGR_PREDICTION_MGC", 
                   "HER2_PREDICTION_MGC", "KI67_PREDICTION_MGC", "ER_PREDICTION_SGC", "PGR_PREDICTION_SGC", "HER2_PREDICTION_SGC", 
                   "KI67_PREDICTION_SGC", "OVERALL_SURVIVAL_EVENT", "ENDOCRINE_TREATED", "CHEMO_TREATED")
for (Integer in Integer_Vars) {
  matrix_clin_df[, Integer] <- as.integer(matrix_clin_df[, Integer])
}


numeric_Vars  <- c("OVERALL_SURVIVAL_DAYS")
for (numeric in numeric_Vars) {
  matrix_clin_df[, numeric] <- as.numeric(matrix_clin_df[, numeric])
}


matrix_clin_df <- data.frame("GSM_ID" = matrix_clin_gse$geo_accession,
                              "OVERALL_SURVIVAL_DAYS" = matrix_clin_gse$`overall survival days:ch1`, 
                              "OVERALL_SURVIVAL_EVENT" = matrix_clin_gse$`overall survival event:ch1`,
                             "PAM")
