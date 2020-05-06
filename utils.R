library(SingleCellExperiment)
library(S4Vectors)

library(magrittr) 
library(dplyr)

require(flowCore)
require(CATALYST)
require(scMerge)

require(ggplot2)
require(ggplotify)
require(cowplot)

require(umap)

# override prepData from catalyst package
prepData <- function(fs, panel, md, features = NULL, transform = TRUE, cofactor = 5,
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

.transform <- function(x, cf, 
    ain = "counts", aout = "exprs",
    dir = c("forwards", "backwards")) {
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
    rownames(panel) <- c()

    # organising data into one object
    metadata <- data.frame( file_name = file_name, sample_id, tissue, patient_id, stringsAsFactors = FALSE )
    panel <- panel[panel$marker_class!="none",]
    return(list("metadata" = metadata, "panel" = panel))
}