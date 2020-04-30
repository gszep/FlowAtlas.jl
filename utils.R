
get_metadata <- function(dataSet,file_paths,
    markers=c(5:24),lineage=c(5,6,9,10,13,15,17,18,19,21,24) ){

    file_name <- pData(dataSet)$name

    # doner and tissue grouping
    sample_id <- regmatches( file_name, regexpr("(.*)(?=.fcs)", file_name, perl=TRUE))
    patient_id <- regmatches( file_paths, regexpr("(?<=/)([0-9]+C)(?=/)", file_paths, perl=TRUE))

    condition <- regmatches( file_name, regexpr("(?<=_)(.*)(?=_)", file_name, perl=TRUE))
    condition <- gsub( '.+-Blood','Blood', condition ) # all blood samples have one label

    # panel definition
    fcs_colname <- colnames(dataSet)

    marker_names <- Filter( function(x) length(x)>0, markernames(dataSet))
    marker_names <- unlist(unique(lapply( marker_names, function(x) gsub("[^[:alnum:]]\\w+","",x) )))
    antigen <- c("FSC_A","FSC_W","SSC_A","SSC_W",marker_names,"Time")

    functional <- setdiff(markers, lineage)
    marker_class <- rep("none",length(fcs_colname))
    marker_class[lineage] <- "type"; marker_class[functional] <- "state"

    marker_class <- factor(marker_class, levels = c("type", "state", "none"))
    panel <- data.frame( fcs_colname, antigen, marker_class, stringsAsFactors = FALSE )

    # organising data into one object
    metadata <- data.frame( file_name = file_paths, sample_id, condition, patient_id, stringsAsFactors = FALSE )
    panel <- panel[panel$marker_class!="none",]
    return(list("metadata" = metadata, "panel" = panel))
}