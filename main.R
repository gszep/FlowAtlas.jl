# local environment variables
setwd("~/Documents/repos/tissue-immunology")
source("utils.R")
data_path = 'data'

################################################################################  loading cleaned batches and metadata
################################################################################
################################################################################

batches <- c()
batch_names <- c()

# load patients as batches one by one 
for ( patient in grep("^((?!(412C)).)*$", # patients to ignore
        Sys.glob( file.path(data_path,'*') ), perl=TRUE, value=TRUE)){

    file_paths <- grep("^((?!(Omentum|Oesophagus|Caecum|EDTA|Chest)).)*$", c( # tissues to ignore

        Sys.glob( file.path(patient,'*.cleaned.fcs') ),
        Sys.glob( file.path(patient,'blood','*.cleaned.fcs'))

    ), perl=TRUE, value=TRUE)

    dataSet <- read.flowSet( file_paths, transformation = FALSE, truncate_max_range = FALSE )
    dataSet <- Subset(dataSet,filter(dataSet, sampleFilter(size = 10000, filterId="dsFilter"))) # subsample files

    setInfo <- get_metadata(dataSet,file_paths)
    dataSet <- prepData(dataSet,setInfo$panel,setInfo$metadata)

    batch_names <- c(batch_names, unique(setInfo$metadata$patient_id))
    batches <- c(batches,dataSet); gc()

}

# merge batches into one object
dataSet <- sce_cbind( batches, method = "intersect", exprs = c("counts", "exprs"),
    colData_names = c("sample_id","tissue","patient_id"), batch_names = batch_names,
    cut_off_batch = 0.01, cut_off_overall = 0.01); gc()

# inherit metadata
rowData(dataSet) <- unique(lapply(batches,rowData))[1]

# separate out by tissue
tissueSets = list()
for ( tissue_name in unique(colData(dataSet)$tissue) ){
    tissueSets[[tissue_name]] <- filterSCE( dataSet, tissue == tissue_name)
}

################################################################################ diagnostic plots
################################################################################
################################################################################

pdf("figures/diagnostic.pdf")
print(plotCounts(dataSet, group_by = "tissue", color_by='patient_id'))
print(plotNRS(dataSet, color_by = "patient_id"))
print(pbMDS(dataSet, label_by = "tissue", shape_by = "patient_id", color_by = "tissue"))
dev.off()

################################################################################ clustering
################################################################################
################################################################################

set.seed(1234) # per tissue
for ( tissue_name in unique(colData(dataSet)$tissue) ){
    print(tissue_name)

    tissueSets[[tissue_name]] <- cluster(tissueSets[[tissue_name]],
        features = "type", xdim = 10, ydim = 10, maxK = 20); gc()
}

pdf("figures/clustering.pdf") # export heatmaps for expert annotation
for ( tissue_name in unique(colData(dataSet)$tissue) ){
    print(plotClusterExprs(tissueSets[[tissue_name]],features=NULL)+labs(title=tissue_name))
}
dev.off()

# import annotations from json
annotations <- fromJSON(file='figures/clustering.merge.json')
for ( tissue_name in unique(colData(dataSet)$tissue) ){

    cluster_ids <- c(1:20)
    for (cluster_name in names(annotations[[tissue_name]]) ){

        idx <- annotations[[tissue_name]][[cluster_name]]
        cluster_ids[idx] <- cluster_name
    }

    tissueSets[[tissue_name]] <- mergeClusters(tissueSets[[tissue_name]], k = "meta20", id = "label", overwrite = TRUE,
        table = data.frame( original_cluster = c(1:20), new_cluster = cluster_ids) )
}

pdf("figures/clustering.merge.pdf") # export merged clustermaps
for ( tissue_name in unique(colData(dataSet)$tissue) ){
    heatmap <- as.ggplot(plotExprHeatmap(tissueSets[[tissue_name]], hm_pal = rev(hcl.colors(10, "YlGnBu")),
        features = NULL, col_clust = FALSE, col_dend = FALSE, bar= TRUE,
        by = "cluster_id", m = "label"))

    print(heatmap+labs(title=tissue_name))
}
dev.off()

################################################################################ dimensionality reduction
################################################################################
################################################################################


dataSet <- runDR(dataSet, dr = "TSNE", cells = 500, features = "type")
dataSet <- runDR(dataSet, dr = "UMAP", cells = 1e3, features = "type")

tsne_plot <- plotDR(dataSet, "TSNE", color_by = "meta8")
umap_plot <- plotDR(dataSet, "UMAP", color_by = "meta8")

plot_grid(tsne_plot+theme(legend.position = "none")+xlab(expression('Nearest-Neighbour Embedding'))+ylab(""),
          umap_plot+theme(legend.position = "none")+xlab(expression('Manifold Approximation'))+ylab(""),
          get_legend(umap_plot+guides(color = guide_legend(title = "Cluster",override.aes=list(size=3)))),
          nrow = 1, rel_widths = c(5, 5, 2))

################################################################################ visualisation of results
################################################################################
################################################################################

dimensionality_reduction <- plotDR(dataSet, "UMAP", color_by = "meta8") + facet_wrap("tissue",ncol=3) + guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster")) + xlab(expression('Expression Manifold Projection X'["1"])) + ylab(expression('Expression Manifold Projection X'["2"]))
abundances <- plotAbundances(dataSet, k = "meta8", by = "cluster_id") + guides(color = guide_legend(nrow = 1, label.position = 'bottom', override.aes = list(size = 0.5), title = "Tissue")) + ylab("Relative Abundance Percentage") + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position="bottom")
plot_grid( dimensionality_reduction, abundances, ncol = 1, labels = c("A", "B"))