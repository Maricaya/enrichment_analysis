# load libraries
library("LOLA")
library("GenomicRanges")
library("data.table")

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript region_enrichment_analysis_LOLA.R <query_regions.bed> <database_path> <result_path>")
}

# Assign command-line arguments to variables
query_regions <- args[1]
database_path <- args[2]
result_path <- args[3]

# Define parameters (you can modify these based on your specific analysis requirements)
genome <- "hg38"  # or "mm10" depending on your use case
database_name <- basename(database_path)  # use the folder name as the database name
region_set <- basename(query_regions)  # use the file name as the region set name

### Load data

# Load query region sets
regionSet_query <- readBed(query_regions)

# Load background/universe region sets (e.g., consensus region set)
# If you have a specific background, add it here. For now, we'll assume the background is embedded in the database.
# Uncomment the following line if a specific background is required:
# regionSet_background <- readBed(background_regions)

# Load the database (requires resources downloaded from https://databio.org/regiondb)
database <- loadRegionDB(file.path(database_path))

###### LOLA

# Run LOLA
res <- runLOLA(regionSet_query, regionSet_query, database, cores=1)  # replace regionSet_background with your background if necessary

# Make description more descriptive
if (database_name == 'LOLACore') {
    res$description <- paste(res$description, res$cellType, res$antibody, sep='.')
} else {
    res$description <- paste(res$description, res$filename, sep='.')
}

# Ensure that the description values are unique
res$description <- make.names(res$description, unique=TRUE)

# Determine raw p-value
res$pValue <- 10^(-1 * res[['pValueLog']])

# Save results
fwrite(as.data.frame(res), file=file.path(result_path), row.names=FALSE)
