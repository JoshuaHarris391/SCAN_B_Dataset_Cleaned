# make destination dir
mkdir -p ./data/GSE81538/SOFT
# Download SOFT
wget -P data/GSE81538/SOFT 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/soft/GSE81538_family.soft.gz'
# Download supplementary files (contains FPKM results)
wget -P data/GSE81538 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_UCSC_Human_hg19_knownGenes_GTF_appended_10sep2012.gtf.gz'
wget -P data/GSE81538 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_gene_expression_405_transformed.csv.gz'
wget -P data/GSE81538 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_map_transcriptID_geneSymbol.csv.gz'
wget -P data/GSE81538 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_pathology_consensus_key.xlsx'
wget -P data/GSE81538 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_transcript_expression_405.csv.gz'
