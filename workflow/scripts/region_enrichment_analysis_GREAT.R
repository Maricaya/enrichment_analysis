# load libraries
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
regions_file <- args[1]
background_file <- args[2]
database_path <- args[3]
result_path <- args[4]
genome <- args[5]
cores_n <- as.numeric(args[6])

# Validate genome argument
if (is.null(genome) || genome == "") {
    stop("Error: Genome is not specified.")
}

# set genome
if (genome == "hg19" | genome == "hg38") {
    orgdb <- "org.Hs.eg.db"
} else if (genome == "mm9" | genome == "mm10") {
    orgdb <- "org.Mm.eg.db"
} else {
    stop("Error: Unsupported genome version.")
}

# load query and background/universe region sets (e.g., consensus region set)
regionSet_query <- import(regions_file, format = "BED")
regionSet_background <- import(background_file, format = "BED")

# load database
database = read_gmt(database_path, from = "SYMBOL", to = "ENTREZ", orgdb = orgdb)

###### GREAT

# run GREAT
res <- great(
    gr = regionSet_query,
    gene_sets = database,
    tss_source = genome,
    biomart_dataset = NULL,
    min_gene_set_size = 5, #default: 5
    mode = "basalPlusExt",
    basal_upstream = 5000,
    basal_downstream = 1000,
    extension = 1000000,
    extended_tss = NULL,
    background = regionSet_background, #default: NULL
    exclude = "gap",
    cores = cores_n, #default: 1
    verbose = TRUE #default: great_opt$verbose
)

# get & save result table
tb <- getEnrichmentTable(res, min_region_hits = 0)
tb$description <- paste(tb$description, tb$id)
fwrite(as.data.frame(tb), file = result_path, row.names = FALSE)
