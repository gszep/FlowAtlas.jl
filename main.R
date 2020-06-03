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
for ( patient in grep("^((?!(390C)).)*$", # patients to ignore
        Sys.glob( file.path(data_path,'*') ), perl=TRUE, value=TRUE)){

    file_paths <- grep("^((?!(Omentum|Oesophagus|Caecum|EDTA|Chest)).)*$", c( # tissues to ignore

        Sys.glob( file.path(patient,'*.cleaned.fcs') ),
        Sys.glob( file.path(patient,'blood','*.cleaned.fcs'))

    ), perl=TRUE, value=TRUE)

    dataSet <- read.flowSet( file_paths, transformation = FALSE, truncate_max_range = FALSE ) # import as flowSet
    dataSet <- Subset(dataSet,filter(dataSet, sampleFilter(size = 10000, filterId="dsFilter"))) # subsample files

    # correction for systematic error in CD4 channel
    if ( patient==file.path(data_path,'390C') ) {
        for ( i in seq_along(dataSet) ) {
            dataSet[[i]]@exprs[,3] <- 0.15*dataSet[[i]]@exprs[,3]
        }
    }

    setInfo <- get_metadata(dataSet,file_paths) # extract metadata
    dataSet <- prepData(dataSet,setInfo$panel,setInfo$metadata,cofactor=250.0) # apply logicale transform

    batch_names <- c(batch_names, unique(setInfo$metadata$patient_id))
    batches <- c(batches,dataSet); gc()

}

# merge batches into one object
dataSet <- sce_cbind( batches, method = "intersect", exprs = c("counts", "exprs"),
    colData_names = c("sample_id","tissue","patient_id"), batch_names = batch_names,
    cut_off_batch = 0, cut_off_overall = 0); gc()

# inherit metadata
rowData(dataSet) <- unique(lapply(batches,rowData))[1]

# manual batch corrections
expression <- assay(dataSet,"exprs")

patient <- dataSet$patient_id=='390C'

channel <- rownames(dataSet)=='CD3'
assay(dataSet,"exprs")[channel,patient] <- expression[channel,patient] + 0.3

channel <- rownames(dataSet)=='CD8'
assay(dataSet,"exprs")[channel,patient] <- 0.9*expression[channel,patient]+0.5

channel <- rownames(dataSet)=='CD45RA'
assay(dataSet,"exprs")[channel,patient] <- 2*expression[channel,patient]-1

channel <- rownames(dataSet)=='CD69'
assay(dataSet,"exprs")[channel,patient] <- 1.2*expression[channel,patient]-0.5

channel <- rownames(dataSet)=='CXCR3'
assay(dataSet,"exprs")[channel,patient] <- 0.6*expression[channel,patient]-0.5

patient <- dataSet$patient_id=='403C'

channel <- rownames(dataSet)=='CD69'
assay(dataSet,"exprs")[channel,patient] <- 1.2*expression[channel,patient]

patient <- dataSet$patient_id=='423C'

channel <- rownames(dataSet)=='HLA'
assay(dataSet,"exprs")[channel,patient] <- 0.6*expression[channel,patient]

channel <- rownames(dataSet)=='CCR4'
assay(dataSet,"exprs")[channel,patient] <- 0.6*expression[channel,patient]

channel <- rownames(dataSet)=='Helios'
assay(dataSet,"exprs")[channel,patient] <- 0.8*expression[channel,patient]

patient <- dataSet$patient_id=='428C'

channel <- rownames(dataSet)=='CD8'
assay(dataSet,"exprs")[channel,patient] <- 0.8*expression[channel,patient]+1

channel <- rownames(dataSet)=='CD45RA'
assay(dataSet,"exprs")[channel,patient] <- 1.5*expression[channel,patient]

channel <- rownames(dataSet)=='Helios'
assay(dataSet,"exprs")[channel,patient] <- expression[channel,patient]+0.25

plotExprs( dataSet, color_by='tissue')
plotExprs( filterSCE( dataSet, tissue == 'Blood'), color_by='patient_id')

# failed attempts to use batchnorm algorithm
# scMerge(dataSet,batch_name="patient_id", cell_type=dataSet$tissue,
#     ctl=c('CD3','HLA'), marker=rownames(dataSet), marker_list=rownames(dataSet),
#     exprs="exprs",assay_name="exprs", verbose = TRUE)

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
        features = NULL, xdim = 10, ydim = 10, maxK = 20); gc()
}

tissue_name <- 'BM'
tissueSets[[tissue_name]] <- cluster(tissueSets[[tissue_name]],
        features = NULL, xdim = 10, ydim = 10, maxK = 20, seed=10)
plotClusterExprs(tissueSets[[tissue_name]],
    features=NULL, color_by = NULL)

pdf("figures/clustering.pdf") # export heatmaps for expert annotation
for ( tissue_name in unique(colData(dataSet)$tissue) ){
    print(plotClusterExprs(tissueSets[[tissue_name]],features=NULL,color_by=NULL)+labs(title=tissue_name))
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

dataSet <- cluster(dataSet,
    features = NULL, xdim = 10, ydim = 10, maxK = 20); gc()
plotClusterExprs( filterSCE( dataSet, tissue == 'Spleen'), k='meta20',
    features=NULL, color_by = NULL)


dataSet <- runDR(dataSet, dr = "TSNE", cells = 500, features = NULL)
dataSet <- runDR(dataSet, dr = "UMAP", cells = 1e3, features = NULL)

tsne_plot <- plotDR(dataSet, "TSNE", color_by = "meta20")
umap_plot <- plotDR(dataSet, "UMAP", color_by = "meta20")

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


xy <- reducedDim(dataSet,'UMAP')
colnames(xy) <- c("x", "y")
df <- data.frame(colData(dataSet), xy)
df <- df[!(is.na(df$x) | is.na(df$y)), ]
 
dataSet
sce <- SingleCellExperiment(df,assay=c('exprs'))
sce

cluster(sce,features = NULL, xdim = 10, ydim = 10, maxK = 20)