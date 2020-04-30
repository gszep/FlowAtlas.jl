setwd("~/Documents/code/flow-cytometry")

require(flowCore)
require(flowViz)
require(flowAI)
source("utils.R")

################################################################################  anomaly filtering
################################################################################
################################################################################

data_path = 'data'
file_paths <- grep("^((?!filtered).)*$", c(

       Sys.glob( file.path(data_path,'*','*.fcs') ),
       Sys.glob( file.path(data_path,'*','blood','*.fcs') )),

    perl=TRUE, value=TRUE)

dataSet <- read.flowSet( file_paths, transformation = FALSE, truncate_max_range = FALSE ); gc()
setInfo <- get_metadata(dataSet,file_paths)

colnames(dataSet)
markernames(dataSet)

setInfo$panel

rgate <- rectangleGate("355 379/28-A"=c(1,4),"355 379/28-A"=c(1,4), filterId="Live T-Cells")
xyplot(`355 379/28-A`~`355 515/30-A`, xlab = 'Z-UV', ylab = 'CD3', smooth=FALSE, xbin=128, 
    data=transform("355 515/30-A"=asinh,"355 379/28-A"=asinh) %on% dataSet, filter=rgate )


xyplot( `SSC-A`~`FSC-A`, data=dataSet, ylab = 'Side Scatter', xlab = 'Forward Scatter',
        smooth=FALSE, xbin=128, filter=rgate)




plot(dataSet[[1]], c("FSC-A","SSC-A","355 515/30-A","405 585/15-A"), smooth = FALSE)

# for (file_path in c("data/423C/blood/Donor_Pre-Blood_023.fcs")) {

#     print(file_path)
#     flow_auto_qc(

#         file_path, second_fractionFR = 0.1, pen_valueFS = 200,  folder_results = dirname(file_path),
#         alphaFR = 10, fcs_QC = FALSE, fcs_highQ=".filtered", html_report = ".anomaly",
#         mini_report = FALSE, output = 0); gc()
# }

################################################################################  loading filtered and metadata
################################################################################
################################################################################

# file_paths <- c( Sys.glob(file.path(data_path,'*','*.filtered.fcs')), Sys.glob(file.path(data_path,'*','blood','*.filtered.fcs')) )
dataSet <- read.flowSet( file_paths, transformation = FALSE, truncate_max_range = FALSE); gc()
# dataSet <- Subset(dataSet, filter(dataSet, sampleFilter(size = 10000, filterId="dsFilter")))

require(CATALYST)
dataSet <- prepData(dataSet,panel,metadata)
gc()

################################################################################ diagnostic plots
################################################################################
################################################################################

plotCounts(dataSet, color_by = "condition")
plotExprs(dataSet, color_by = "condition")

CATALYST::plotMDS(dataSet, color_by = "condition")
plotNRS(dataSet, features = type_markers(dataSet), color_by = "condition")

################################################################################ clustering
################################################################################
################################################################################

set.seed(1234)
dataSet <- cluster(dataSet, features = type_markers(dataSet),
                   xdim = 10, ydim = 10, maxK = 20); gc()

cluster_labels <- data.frame( original_cluster = c(1:20),
    new_cluster = c('Debris','CD8','Th17','Thf','CD4','CD4','Treg','CD4','Thf','Thf',
                    'Debris','Thf','Th22','Th22','CD8','Debris','Debris','Debris','Debris','Debris'))

dataSet@metadata$cluster_codes$label <- NULL
dataSet <- mergeClusters(dataSet, k = "meta20", table = cluster_labels, id = "label")
plotClusterHeatmap(dataSet, hm2 = NULL, k = "meta20", m = "label", cluster_anno = FALSE)

require(ggplot2)
require(cowplot)
require(umap)

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

dimensionality_reduction <- plotDR(dataSet, "UMAP", color_by = "label") + facet_wrap("condition",ncol=7) + guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster")) + xlab(expression('Expression Manifold Projection X'["1"])) + ylab(expression('Expression Manifold Projection X'["2"]))
abundances <- plotAbundances(dataSet, k = "label", by = "cluster_id") + guides(color = guide_legend(nrow = 1, label.position = 'bottom', override.aes = list(size = 0.5), title = "Tissue")) + ylab("Relative Abundance Percentage") + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position="bottom")
plot_grid( dimensionality_reduction, abundances, ncol = 1, labels = c("A", "B"))


################################################################################ diffcyt
################################################################################
################################################################################

# design <- createDesignMatrix( metadata, cols_design = c("condition", "patient_id") )
# contrast <- createContrast(c(0,1,0))

# differential_abundance <- diffcyt( dataSet,
#                                   design = design, contrast = contrast, analysis_type = "DA",
#                                   seed_clustering = 123 )

# display table of results for top cluster-marker combinations
# topTable(differential_abundance, format_vals = TRUE)

# calculate number of significant detected clusters at 10% false discovery
# table( topTable(differential_abundance, all = TRUE)$p_adj <=  0.1)
