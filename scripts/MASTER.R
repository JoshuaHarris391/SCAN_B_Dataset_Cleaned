#################################################################
# Download Data
#################################################################
system('bash ./scripts/download_GSE96058.sh')
system('bash ./scripts/download_GSE81538.sh')


#################################################################
# Run Pre-processing Scripts
#################################################################
source(file = './scripts/GSE96058_GEOquery.R')
source(file = './scripts/GSE81538_GEOquery.R')

#################################################################
# Clear Environment and Load Cleaned Data
#################################################################
rm(list = ls())
load(file = "data/GSE96058/clean/GSE96058_Cleaned.RData")
load(file = "data/GSE81538/clean/GSE81538_Cleaned.RData")

