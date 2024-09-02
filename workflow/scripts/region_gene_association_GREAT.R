# load necessary libraries
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 解析参数
regions_file <- args[1]
database_path <- args[2]
gene_path <- args[3]
associations_table_path <- args[4]
associations_plot_path <- args[5]

# 你可以将 genome, great_params, 和 cores_n 设为默认值，或者继续从命令行中获取
genome <- "hg19"  # 默认值，或者从命令行参数中获取
great_params <- list(min_gene_set_size = 5, mode = "basalPlusExt", basal_upstream = 5, basal_downstream = 5, extension = 1000)
cores_n <- 1  # 默认值，或者从命令行参数中获取

# set genome
if (genome == "hg19" | genome == "hg38") {
    orgdb <- "org.Hs.eg.db"
} else if (genome == "mm9" | genome == "mm10") {
    orgdb <- "org.Mm.eg.db"
}

# load query region set
regionSet_query <- import(regions_file, format = "BED")

# load database
database = read_gmt(file.path(database_path), from = "SYMBOL", to = "ENTREZ", orgdb = orgdb)

###### GREAT

# run GREAT
res <- great(gr = regionSet_query,
      gene_sets = database,
      tss_source = genome,
      biomart_dataset = NULL,
      min_gene_set_size = great_params[["min_gene_set_size"]], # default: 5
      mode = great_params[["mode"]],
      basal_upstream = great_params[["basal_upstream"]],
      basal_downstream = great_params[["basal_downstream"]],
      extension = great_params[["extension"]],
      extended_tss = NULL,
      background = NULL,
      exclude = "gap",
      cores = cores_n, # default: 1
      verbose = TRUE # default: great_opt$verbose
     )

# plot gene-region association
pdf(file = file.path(associations_plot_path), width = 12, height = 4)
plotRegionGeneAssociations(res)
dev.off()

# get and save gene-region association
associations <- getRegionGeneAssociations(res)
fwrite(as.data.frame(associations), file = file.path(associations_table_path), row.names = TRUE)

# save unique associated genes by using mcols(), which returns a DataFrame object containing the metadata columns.
genes <- unique(unlist(mcols(associations)$annotated_genes))
write(genes, file.path(gene_path))
