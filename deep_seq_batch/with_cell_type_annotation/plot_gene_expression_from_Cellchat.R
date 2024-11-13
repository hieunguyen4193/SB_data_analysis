plotGeneExpression <- function (object, features = NULL, signaling = NULL, enriched.only = TRUE, 
          type = c("violin", "dot", "bar"), color.use = NULL, group.by = NULL, 
          ...) 
{
  type <- match.arg(type)
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  if (!identical(rownames(meta), colnames(object@data.signaling))) {
    cat("The cell barcodes in 'meta' is ", head(rownames(meta)), 
        "\n")
    warning("The cell barcodes in 'meta' is different from those in the used data matrix.\n              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
    rownames(meta) <- colnames(object@data.signaling)
  }
  w10x <- Seurat::CreateSeuratObject(counts = object@data.signaling, 
                                     meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  }
  else if (!is.null(signaling)) {
    res <- extractEnrichedLR(object, signaling = signaling, 
                             geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- VlnPlot(w10x, features = feature.use)
  }
  else if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, color.use = color.use, 
                  ...)
  }
  else if (type == "bar") {
    gg <- barPlot(w10x, features = feature.use, color.use = color.use, 
                  ...)
  }
  return(gg)
}
