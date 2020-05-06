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
for (patient in grep("^((?!(412C)).)*$",
        Sys.glob( file.path(data_path,'*') ), perl=TRUE, value=TRUE)){

    file_paths <- grep("^((?!(Omentum|Oesophagus|Caecum)).)*$", c(

        Sys.glob( file.path(patient,'*.cleaned.fcs') ),
        Sys.glob( file.path(patient,'blood','*.cleaned.fcs'))

    ), perl=TRUE, value=TRUE)

    dataSet <- read.flowSet( file_paths, transformation = FALSE, truncate_max_range = FALSE )
    dataSet <- Subset(dataSet,filter(dataSet, sampleFilter(size = 10000, filterId="dsFilter")))

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

################################################################################ diagnostic plots
################################################################################
################################################################################

plotCounts(dataSet, group_by = "tissue")
plotExprs(dataSet, color_by = "tissue")

pbMDS(dataSet, label_by = "tissue", shape_by = "patient_id", color_by = "tissue")
plotNRS(dataSet, features = type_markers(dataSet), color_by = "tissue")

################################################################################ clustering
################################################################################
################################################################################

set.seed(1234)

plots = list()
for (tissue in unique(colData(dataSet)$tissue) ){
    tissueSet <- filterSCE( dataSet, tissue == tissue)

    tissueSet <- cluster(tissueSet, features = type_markers(tissueSet), xdim = 10, ydim = 10, maxK = 20); gc()
    heatmap <- as.ggplot(plotExprHeatmap(tissueSet, hm_pal = rev(hcl.colors(10, "YlGnBu")),
        features = NULL, row_anno = FALSE, col_dend = FALSE,
        by = "cluster_id" ))
    
    plots[[tissue]] = heatmap+labs(title=tissue)
}

pdf(paste("clustering","pdf",sep=".")) 
for (tissue in unique(colData(dataSet)$tissue) ){
    print(plots[[tissue]])
}
dev.off()

cluster_labels <- data.frame( original_cluster = c(1:20),
    new_cluster = c('Debris','CD8','Th17','Thf','CD4','CD4','Treg','CD4','Thf','Thf',
                    'Debris','Thf','Th22','Th22','CD8','Debris','Debris','Debris','Debris','Debris'))

tissueSet@metadata$cluster_codes$label <- NULL
tissueSet <- mergeClusters(tissueSet, k = "meta20", table = cluster_labels, id = "label")
plotExprHeatmap(tissueSet, hm_pal = rev(hcl.colors(10, "YlGnBu")),
    features = "type", by = "cluster_id", m = "label")

################################################################################ dimensionality reduction
################################################################################
################################################################################

dataSet <- runDR(dataSet, dr = "TSNE", cells = 500, features = "type")
dataSet <- runDR(dataSet, dr = "UMAP", cells = 1e3, features = "type")

tsne_plot <- plotDR(dataSet, "TSNE", color_by = "label")
umap_plot <- plotDR(dataSet, "UMAP", color_by = "label")

plot_grid(tsne_plot+theme(legend.position = "none")+xlab(expression('Nearest-Neighbour Embedding'))+ylab(""),
          umap_plot+theme(legend.position = "none")+xlab(expression('Manifold Approximation'))+ylab(""),
          get_legend(umap_plot+guides(color = guide_legend(title = "Cluster",override.aes=list(size=3)))),
          nrow = 1, rel_widths = c(5, 5, 2))

################################################################################ visualisation of results
################################################################################
################################################################################

dimensionality_reduction <- plotDR(dataSet, "UMAP", color_by = "label") + facet_wrap("tissue",ncol=7) + guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster")) + xlab(expression('Expression Manifold Projection X'["1"])) + ylab(expression('Expression Manifold Projection X'["2"]))
abundances <- plotAbundances(dataSet, k = "label", by = "cluster_id") + guides(color = guide_legend(nrow = 1, label.position = 'bottom', override.aes = list(size = 0.5), title = "Tissue")) + ylab("Relative Abundance Percentage") + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position="bottom")
plot_grid( dimensionality_reduction, abundances, ncol = 1, labels = c("A", "B"))


################################################################################ diffcyt
################################################################################
################################################################################

# design <- createDesignMatrix( metadata, cols_design = c("tissue", "patient_id") )
# contrast <- createContrast(c(0,1,0))

# differential_abundance <- diffcyt( dataSet,
#                                   design = design, contrast = contrast, analysis_type = "DA",
#                                   seed_clustering = 123 )

# display table of results for top cluster-marker combinations
# topTable(differential_abundance, format_vals = TRUE)

# calculate number of significant detected clusters at 10% false discovery
# table( topTable(differential_abundance, all = TRUE)$p_adj <=  0.1)
