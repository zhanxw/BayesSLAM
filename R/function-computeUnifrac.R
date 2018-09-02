#' Compute Unifrac distances using GUnifrac for MetaPhlAn count table
#'
#' The default tree structure are provided by curatedMetagenomicData package
#' 
#' @param otu.tab.rff a matrix, rownames are MetaPhlAn tip species
#' @param tree an R phylo object. If NULL is provided, the default phylo object for MetaPhlAn data is used.
#' @param alpha weight for the sharing lineages used for calculation of UniFrac distances
#'
#' @return a list of UniFrac distance matrix using alpha 0, 0.5 and 1
#' @export
#'
#' @examples
#' prefix = system.file("extdata", package = "BayesSLAM")
#' load(file.path(prefix, "Castro-NallarE_2015.Rdata"))
#' countsTable <- otuTable
#' otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
#' otu.tab.rff <- countsTable
#' unifracs = computeUnifrac(otu.tab.rff)
#'    
computeUnifrac <- function(otu.tab.rff, tree = NULL, alpha = c(0, 0.5, 1)) {
  if (is.null(tree)) {
    tree <- computeMetaphlanTree(rownames(otu.tab.rff))
  }

  otu.tab.rff <- otu.tab.rff [rownames(otu.tab.rff) %in% tree$tip.label, ]
  
  unifracs <-
    GUniFrac(t(otu.tab.rff), tree, alpha)$unifracs
  unifracs
}
