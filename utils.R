library(SingleCellExperiment)
library(S4Vectors)

library(rjson)
library(reshape2)
library(magrittr) 
library(dplyr)

require(flowCore)
require(CATALYST)

require(BiocSingular)
require(BiocParallel)
require(scMerge)

require(ggridges)
require(ggplot2)
require(ggplotify)
require(cowplot)

require(umap)


.check_sce <- CATALYST:::.check_sce
.check_cd_factor <- CATALYST:::.check_cd_factor
.check_k <- CATALYST:::.check_k

.get_features <- CATALYST:::.get_features
.agg <- CATALYST:::.agg

# override prepData from catalyst package
prepData <- function(fs, panel, md, features = NULL, transform = TRUE, cofactor = 150,
    panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
    md_cols = list(file = "file_name", id = "sample_id", factors = c("tissue", "patient_id")) ) {

    # assure panel & metadata are data.frames
    for (u in c("panel", "md"))
        if (!is.null(get(u)))
            assign(u, data.frame(get(u), 
                check.names = FALSE, 
                stringsAsFactors = FALSE))
    
    # fill up missing values with function defaults
    stopifnot(is.list(panel_cols), is.list(md_cols))
    args <- as.list(environment())
    for (i in c("md_cols", "panel_cols")) {
        defs <- as.list(formals("prepData")[[i]][-1])
        miss <- !names(defs) %in% names(args[[i]])    
        if (any(miss)) assign(i, c(args[[i]], defs[miss])[names(defs)])
    }
    stopifnot(
        c("channel", "antigen") %in% names(panel_cols),
        c("file", "id", "factors") %in% names(md_cols))
    
    # check channels listed in panel 
    stopifnot(panel[[panel_cols$channel]] %in% colnames(fs))
    
    if (is.null(features)) {
        # default to channels listed in panel
        features <- as.character(panel[[panel_cols$channel]])
    } else {
        # check validity of 'features'
        chs <- colnames(fs)
        check1 <- is.logical(features) && length(features) == length(chs)
        check2 <- is.integer(features) && all(features %in% seq_along(chs))
        check3 <- all(features %in% chs)
        if (!any(check1, check2, check3))
            stop("Invalid argument 'features'. Should be either", 
                " a logial vector,\n  a numeric vector of indices, or",
                " a character vector of column names.")
        # subset panel & reorder features accordingly
        m <- match(panel[[panel_cols$channel]], features, nomatch = 0)
        panel <- panel[m != 0, , drop = FALSE]
        features <- features[m]
    }
    
    # check that filenames or identifiers 
    # match b/w 'flowSet' & metadata
    ids0 <- md[[md_cols$file]]
    ids1 <- fsApply(fs, identifier)
    ids2 <- keyword(fs, "FILENAME")
    if (length(unlist(ids2)) == length(fs))
        ids2 <- basename(ids2)
    check1 <- all(ids1 %in% ids0)
    check2 <- all(ids2 %in% ids0)
    ids_use <- which(c(check1, check2))[1]
    ids <- list(ids1, ids2)[[ids_use]]
    if (is.null(ids)) {
        stop("Couldn't match 'flowSet'/FCS filenames\n", 
            "with those listed in 'md[[md_cols$file]]'.")
    } else {
        # reorder 'flowSet' frames according to metadata table
        fs <- fs[match(md[[md_cols$file]], ids)]
    }
    
    # assure correctness of formats
    k <- c(md_cols$id, md_cols$factors)
    md <- md[, k, drop = FALSE] %>% 
        mutate_all(factor) %>% 
        rename("sample_id" = md_cols$id)
    
    # replace problematic characters
    as <- panel[[panel_cols$antigen]]
    as[is.na(as)] <- panel$fcs_colname[is.na(as)]
    
    # column & panel subsetting
    fs <- fs[, features]
    chs0 <- colnames(fs)
    
    # replace channel w/ antigen names
    m1 <- match(panel[[panel_cols$channel]], chs0, nomatch=0)
    m2 <- match(chs0, panel[[panel_cols$channel]], nomatch=0)
    as <- as[m2]
    ns <- table(as)
    for (a in names(ns)) if (ns[a] > 1)
        as[as == a] <- paste(a, seq_len(ns[a]), sep = ".")
    flowCore::colnames(fs)[m1] <- as
    chs <- colnames(fs)
    
    # get exprs.
    es <- matrix(fsApply(fs, exprs), byrow = TRUE,
        nrow = length(chs), dimnames = list(chs, NULL))
    
    # fix event times
    t <- grep("time", colnames(fs), ignore.case = TRUE)
    if (length(t) != 0) {
        ns <- fsApply(fs, nrow)
        t0 <- c(1, cumsum(ns) + 1)
        tx <- t0[-1] - 1
        for (i in seq_along(fs)[-1]) {
            idx <- seq(t0[i], tx[i])
            es[t, idx] <- es[t, idx] + es[t, tx[i - 1]]
        }
    }

    # get & check marker classes if provided
    mcs <- c("type", "state", "none")
    if (is.null(panel_cols$class) | is.null(panel[[panel_cols$class]])) {
        mcs <- factor("none", levels = mcs)
    } else {
        mcs <- factor(panel[[panel_cols$class]], levels = mcs)
        if (any(is.na(mcs)))
            stop("Invalid marker classes detected;",
                " valid classes are 'type', 'state', and 'none'.")
    }
    
    # construct row/colData & int_metadata
    rd <- DataFrame(
        row.names = chs, channel_name = chs0, 
        marker_name = chs, marker_class = mcs)
    m <- match(chs0, panel[[panel_cols$channel]], nomatch = 0)
    rd$use_channel <- panel$use_channel

    md$n_cells <- as.numeric(fsApply(fs, nrow))
    k <- setdiff(names(md), "n_cells")
    cd <- DataFrame(lapply(md[k], function(u) {
        v <- as.character(rep(u, md$n_cells))
        factor(v, levels = levels(u))
    }), row.names = NULL) 
    
    # construct SCE
    sce <- SingleCellExperiment(
        assays = list(counts = es), 
        rowData = rd, colData = cd,
        metadata = list(experiment_info = md))
    
    ds <- keyword(fs[[1]])
    l <- list(cyt = "\\$CYT$", sn = "\\$CYTSN$")
    keep <- lapply(l, grep, names(ds))
    int_metadata(sce)$description <- ds[unlist(keep)]
    
    # (optionally) do arcsinh-transformation & return SCE
    if (transform) .transform(sce, cofactor) else sce
}

# logicale transofrmation of signals
.transform <- function(x, cf, ain = "counts", aout = "exprs", dir = c("forwards", "backwards")) {
    dir <- match.arg(dir)
    chs <- channels(x)
    stopifnot(is.numeric(cf), cf > 0)
    if (length(cf) == 1) {
        int_metadata(x)$cofactor <- cf
        cf <- rep(cf, nrow(x))
    } else {
        stopifnot(!is.null(names(cf)), chs %in% names(cf))
        cf <- cf[match(chs, names(cf))]
        int_metadata(x)$cofactor <- cf
    }
    if (dir == "forwards") {
        fun <- asinh
        op <- "/"
    } else {
        fun <- sinh
        op <- "*"
    }
    y <- assay(x, ain)
    y <- fun(sweep(y, 1, cf, op))
    assay(x, aout, FALSE) <- y
    return(x)
}

# extract panel, tissue and patient ids
get_metadata <- function(dataSet,file_paths,functional=c('CD69','CD103','HLA','PD') ){

    file_name <- pData(dataSet)$name

    # doner and tissue grouping
    sample_id <- regmatches( file_name, regexpr("(.*)(?=.fcs)", file_name, perl=TRUE))
    patient_id <- regmatches( file_paths, regexpr("(?<=/)([0-9]+C)(?=/)", file_paths, perl=TRUE))

    tissue <- regmatches( file_name, regexpr("(?<=_)(.*)(?=_)", file_name, perl=TRUE))
    tissue <- gsub( '.+-Blood','Blood', tissue ) # all blood samples have one label

    # panel definition
    fcs_colname <- colnames(dataSet)
    markers <- c(1:length(fcs_colname))

    marker_names <- Filter( function(x) length(x)>0, markernames(dataSet))
    antigen <- c(unlist(unique(lapply( marker_names, function(x) gsub("[^[:alnum:]]\\w+","",x) ))))
    functional <- match(functional,antigen)

    lineage <- setdiff(markers,functional)
    marker_class <- rep("none",length(fcs_colname))
    marker_class[lineage] <- "type"; marker_class[functional] <- "state"

    marker_class <- factor(marker_class, levels = c("type", "state", "none"))
    panel <- data.frame( fcs_colname, antigen, marker_class, stringsAsFactors = FALSE )
    print(panel)
    rownames(panel) <- c()

    # organising data into one object
    metadata <- data.frame( file_name = file_name, sample_id, tissue, patient_id, stringsAsFactors = FALSE )
    panel <- panel[panel$marker_class!="none",]
    return(list("metadata" = metadata, "panel" = panel))
}

# override plotClusterExprs from catalyst package
plotClusterExprs <- function(x, k = "meta20", features = "type", color_by = "condition") {
    
    # check validity of input arguments
    .check_sce(x, TRUE)
    .check_cd_factor(x, color_by)

    k <- .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    features <- .get_features(x, features)

    # order clusters according to hierarchical 
    # clustering on median feature expressions 
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    d <- dist(ms, method="euclidean")
    o <- hclust(d, method="average")$order
    
    # construct data.frame of expression matrix include cell metadata
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, 
        id.vars = names(cd),
        variable.name = "antigen", 
        value.name = "expression")
    # add average across all clusters as referebce
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    
    # compute cluster frequencies
    fq <- tabulate(x$cluster_id) / ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    
    # reorder clusters
    df$cluster_id <- factor(df$cluster_id, 
        levels = rev(c("avg", levels(x$cluster_id)[o])),
        labels = rev(c("average", paste0(names(fq), " (", fq, "%)")[o])))
    
    ggplot(df, aes_string(
        x = "expression", y = "cluster_id", 
        col = color_by)) + 
        facet_wrap(~antigen, scales = "free_x", nrow = 2) + 
        geom_density_ridges(alpha = 0.2) + 
        theme_ridges() + theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"))
}

# override batch normalisation algorithm
scMerge <- function(sce_combine, ctl = NULL, kmeansK = NULL, 
    exprs = "logcounts", hvg_exprs = "counts", batch_name = "batch", marker = NULL, 
    marker_list = NULL, ruvK = 20, replicate_prop = 1, cell_type = NULL, 
    cell_type_match = FALSE, cell_type_inc = NULL, BSPARAM = ExactParam(), 
    svd_k = 50, dist = "cor", WV = NULL, WV_marker = NULL, 
    BPPARAM = SerialParam(), return_all_RUV = FALSE, BACKEND = NULL,
    assay_name = NULL, plot_igraph = TRUE, verbose = FALSE) {
    
    ## Checking input expression
    if (is.null(exprs)) {
        stop("exprs is NULL.")
    }
    
    
    ## Checking input expression assay name in SCE object
    if (!exprs %in% SummarizedExperiment::assayNames(sce_combine)) {
        stop(paste("No assay named", exprs))
    }
    
    colsum_exprs = DelayedMatrixStats::colSums2(SummarizedExperiment::assay(sce_combine, exprs))
    rowsum_exprs = DelayedMatrixStats::rowSums2(SummarizedExperiment::assay(sce_combine, exprs))
    if(any(colsum_exprs == 0) | any(rowsum_exprs == 0)){
        message("Automatically removing ", sum(colsum_exprs == 0), " cells and ",
                sum(rowsum_exprs == 0), " genes that are all zeroes in the data")
        sce_combine = sce_combine[rowsum_exprs != 0, colsum_exprs != 0]
    }
    
    ## Checking if the cell names are non-unique
    cellNames = colnames(sce_combine)
    
    if (length(cellNames) != length(unique(cellNames))) {
        stop("Please make sure column names are unique.")
    }
    
    if (is.null(assay_name)) {
        stop("assay_name is NULL, please provide a name to store the results under")
    }
    
    if (length(ruvK) > 1){
        message("You chose more than one ruvK. The argument return_all_RUV is forced to be TRUE.")
        return_all_RUV = TRUE
    }
    
    if (return_all_RUV) {
        message("You chose return_all_RUV = TRUE. The result will contain all RUV computations. This could be a very large object.")
        ## We need an assay_name for every ruvK, if return_all_RUV is
        ## TRUE
        if (length(assay_name) != length(ruvK)) {
            stop("You chose return_all_RUV = TRUE. In this case, the length of assay_name must be equal to the length of ruvK")
        }
    }
    
    
    ## Extracting data matrix from SCE object
    exprs_mat <- SummarizedExperiment::assay(sce_combine, exprs)
    if (!is.matrix(exprs_mat)) {
        # stop(paste0("The assay named '", exprs, "' must be of class 'matrix', please convert this."))
    }
    
    
    sce_rownames <- rownames(sce_combine)
    
    # if (is.null(colnames(exprs_mat)) | 
    #     length(colnames(exprs_mat)) != length(unique(colnames(exprs_mat)))) {
    #     stop("colnames of exprs is NULL or contains duplicates")
    # }
    
    
    hvg_exprs_mat <- SummarizedExperiment::assay(sce_combine, hvg_exprs)
    if (!is.matrix(hvg_exprs_mat)) {
        # stop(paste0("The assay named '", hvg_exprs, "' must be of class 'matrix', please convert this."))
    }
    
    ## Checking negative controls input
    if (is.null(ctl)) {
        stop("Negative control genes are needed. \n")
    } else {
        if (is.character(ctl)) {
            ctl <- which(sce_rownames %in% ctl)
        }
        if (length(ctl) == 0) {
            stop("Could not find any negative control genes in the row names of the expression matrix", 
                call. = FALSE)
        }
    }
    
    ## Checking the batch info
    if (!(batch_name %in% colnames(colData(sce_combine)))) {
        stop("Could not find a 'batch' column in colData(sce_combine)", 
            call. = FALSE)
    }
    
    sce_combine$batch = sce_combine[[batch_name]]
    
    if (is.factor(sce_combine$batch)) {
        batch <- droplevels(sce_combine$batch)
    } else {
        batch <- sce_combine$batch
    }
    
    
    
    # ## If the user supplied a parallelParam class, then regardless
    # ## of parallel = TRUE or FALSE, we will use that class Hence
    # ## no if statement for this case.
    # if (!is.null(parallelParam)) {
    #     message("Step 1: Computation will run in parallel using supplied parameters")
    # }
    # 
    # ## If parallel is TRUE, but user did not supplied a
    # ## parallelParam class, then we set it to bpparam()
    # if (parallel & is.null(parallelParam)) {
    #     message("Step 1: Computation will run in parallel using BiocParallel::bpparam()")
    #     parallelParam = BiocParallel::bpparam()
    # }
    # 
    # ## If parallel is FALSE, or the user did not supplied a
    # ## parallelParam class, we will use SerialParam()
    # if (!parallel | is.null(parallelParam)) {
    #     message("Step 1: Computation will run in serial")
    #     parallelParam = BiocParallel::SerialParam()
    # }
    
    
    ## Finding pseudo-replicates
    t1 <- Sys.time()
    repMat <- scReplicate(sce_combine = sce_combine, batch = batch, 
        kmeansK = kmeansK, exprs = exprs, hvg_exprs = hvg_exprs, 
        marker = marker, marker_list = marker_list, replicate_prop = replicate_prop, 
        cell_type = cell_type, cell_type_match = cell_type_match, 
        cell_type_inc = cell_type_inc, dist = dist, WV = WV, 
        WV_marker = WV_marker, BPPARAM = BPPARAM, 
        BSPARAM = BSPARAM, plot_igraph = plot_igraph, verbose = verbose)
    t2 <- Sys.time()
    
    timeReplicates <- t2 - t1
    
    cat("Dimension of the replicates mapping matrix: \n")
    print(dim(repMat))
    
    ## Performing RUV normalisation
    
    message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")
    
    ruv3res <- scRUVIII(Y = exprs_mat, M = repMat, ctl = ctl, 
        k = ruvK, batch = batch, fullalpha = NULL, cell_type = cell_type, 
        return_all_RUV = return_all_RUV, BSPARAM = BSPARAM, BPPARAM = BPPARAM,
        svd_k = svd_k)
    t3 <- Sys.time()
    
    timeRuv <- t3 - t2
    
    sce_combine <- sce_combine
    
    if (return_all_RUV) {
        ## if return_all_RUV is TRUE, then the previous check ensured
        ## assay_name is not NULL and matches the length of ruvK And
        ## the scRUVIII should've just returned with a single result
        ## (ruv3res_optimal)
        listNewY <- lapply(ruv3res[names(ruv3res) != "optimal_ruvK"], 
                           function(x) {t(x$newY)})
        
        for (i in seq_len(length(listNewY))) {
            SummarizedExperiment::assay(sce_combine, assay_name[i]) <- listNewY[[i]]
        }
    } else {
        ## If return_all_RUV is FALSE, then scRUVIII should've just
        ## returned with a single result (ruv3res_optimal)
        SummarizedExperiment::assay(sce_combine, assay_name) <- t(ruv3res$newY)
    }
    
    if(is.null(BACKEND)){
        
    } else {
        for(i in SummarizedExperiment::assayNames(sce_combine)){
            SummarizedExperiment::assay(sce_combine, i) = DelayedArray::realize(SummarizedExperiment::assay(sce_combine, i), BACKEND)
        }
    }
    
    S4Vectors::metadata(sce_combine) <- c(S4Vectors::metadata(sce_combine), 
        list(ruvK = ruvK, ruvK_optimal = ruv3res$optimal_ruvK, 
            scRep_res = repMat, timeReplicates = timeReplicates, 
            timeRuv = timeRuv))
    
    message("scMerge complete!")
    
    return(sce_combine)
}