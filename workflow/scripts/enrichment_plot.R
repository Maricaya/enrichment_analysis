# load libraries
library("ggplot2")
library("svglite")
library("data.table")

# Define evaltext function
evaltext <- function(x) {
    eval(parse(text = x))
}

# Define addline_format function (if needed)
addline_format <- function(x) gsub("\\s", "\n", x)

# Check if running interactively or within Snakemake
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Expected 3 arguments: feature_set, tool, db")
}

feature_set <- args[1]
tool <- args[2]
db <- args[3]

# Manually set the paths based on the provided arguments
result_path <- "test/results"  # Adjust this to your actual result path
config_path <- "test/config/example_enrichment_analysis_config.yaml"  # Adjust to your config file path

# Load the configuration file
config <- yaml::yaml.load_file(config_path)

# Set paths for enrichment result and plot based on input arguments
enrichment_result_path <- file.path(result_path, "enrichment_analysis", feature_set, tool, db, paste0(feature_set, "_", db, ".csv"))

enrichment_plot_path <- file.path(result_path, "enrichment_analysis", feature_set, tool, db, paste0(feature_set, "_", db, ".png"))


# Print paths for debugging
print(paste("Enrichment result path:", enrichment_result_path))
print(paste("Enrichment plot path:", enrichment_plot_path))

# Ensure enrichment result file exists
if (!file.exists(enrichment_result_path)) {
    stop("Enrichment result file not found:", enrichment_result_path)
}

# parameters
plot_cols <- config[["column_names"]][[tool]]

top_n <- plot_cols[["top_n"]]
pval_col <- plot_cols[["p_value"]]
adjp_col <- plot_cols[["adj_pvalue"]]
effect_col <- plot_cols[["effect_size"]]
overlap_col <- plot_cols[["overlap"]]
term_col <- plot_cols[["term"]]

# load enrichment result
if (file.size(enrichment_result_path) != 0L){
    enrichment_result <- data.frame(fread(enrichment_result_path, header=TRUE))
} else {
    file.create(enrichment_plot_path)
    quit(save = "no", status = 0)
}

top_n <- min(top_n, nrow(enrichment_result))

# evaluate overlap numerically if necessary
if(class(enrichment_result[[overlap_col]]) == "character") {
    enrichment_result[[overlap_col]] <- as.numeric(lapply(enrichment_result[[overlap_col]], evaltext))
}

# calculate comparable effect size either NES or odds-ratio/fold based
if (tool != "preranked_GSEApy" & tool != "pycisTarget" & tool != "RcisTarget") {
    # calculate log2(effect-size) and put in new column
    effect_col_new <- paste0("log2_", effect_col)
    enrichment_result[[effect_col_new]] <- log2(enrichment_result[[effect_col]])
    effect_col <- effect_col_new
}

# determine ranks
enrichment_result$PValue_Rnk <- if (tool != "pycisTarget" & tool != "RcisTarget") rank(enrichment_result[[pval_col]]) else rank(-enrichment_result[[pval_col]])
enrichment_result$Fold_Rnk <- rank(-abs(enrichment_result[[effect_col]]))
enrichment_result$Coverage_Rnk <- rank(-enrichment_result[[overlap_col]])
# calculate and sort by mean rank
enrichment_result$meanRnk <- rowMeans(enrichment_result[, c('PValue_Rnk', 'Fold_Rnk', 'Coverage_Rnk')])
enrichment_result <- enrichment_result[order(enrichment_result$meanRnk, decreasing = FALSE), ]

# format term column so that order is kept and values are unique
enrichment_result[[term_col]] <- make.unique(as.character(enrichment_result[[term_col]]), sep = "_")
enrichment_result[[term_col]] <- factor(enrichment_result[[term_col]], levels = enrichment_result[[term_col]])

# draw enrichment plot
do_enrichment_plot <- function(plot_data, title, x, y, size, colorBy, font.size, path, filename, top_n) {
    enr_p <- ggplot(plot_data, aes(x = !!sym(x), y = !!sym(y), size = !!sym(size), color = !!sym(colorBy))) +
        geom_point() +
        scale_color_continuous(low = "red", high = "blue", name = colorBy, guide = guide_colorbar(reverse = TRUE)) +
        ggtitle(title) +
        scale_size(range = c(3, 8)) +
        scale_y_discrete(limits = rev) +
        theme(axis.text.y = element_text(vjust = 0.6))

    ggsave(filename = file.path(path, paste0(filename, ".png")),
           plot = enr_p,
           width = 200,
           height = 10 * top_n,
           units = "mm")
}

# plot top_n terms by mean_rnk
do_enrichment_plot(plot_data = enrichment_result[1:top_n, ],
                   title = paste0(tool, ' results of ', feature_set, '\nin ', db),
                   x = effect_col,
                   y = term_col,
                   size = overlap_col,
                   colorBy = adjp_col,
                   font.size = 10,
                   path = dirname(enrichment_plot_path),
                   filename = paste0(feature_set, "_", db),
                   top_n = top_n
)
