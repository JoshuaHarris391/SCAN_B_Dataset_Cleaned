# make destination dir
mkdir -p ./data/GSE96058/SOFT
# Download SOFT
wget -P data/GSE96058/SOFT 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/soft/GSE96058_family.soft.gz'
# Download supplementary files (contains FPKM results)
wget -P data/GSE96058 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz'
wget -P data/GSE96058 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz'
wget -P data/GSE96058 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz'
