
#' Compute a subset of tree for a count table
#' 
#' Based on a MetaPhlAn tree, only keep the tip nodes that appear in the count table row names.
#'
#' @param otu.tab.rff counts table (taxa by sample), provided by MetaPhlAn
#' @param tree R phylo object
#'
#' @return a phylo tree
#' @export
#'
#' @examples
#' library(ape)
#' prefix = system.file("extdata", package = "BayesSLAM")
#' load(file.path(prefix, "Castro-NallarE_2015.Rdata"))
#' countsTable <- otuTable
#' tree <- computeMetaphlanTree(countsTable)
computeMetaphlanTree <- function(otu.tab.rff, tree = NULL) {
  if (is.null(tree)) {
  # ot <- otu.tab.rff
  # Calculate the UniFracs
  # this preprocessed file was in curateMetagenomicData, we save and reuse its copy
  tree.file <- system.file("extdata/metaphlan2_selected.tree.reroot.nwk.bz2", package = "BayesSLAM")
  stopifnot(file.exists(tree.file))
  tree <- ape::read.tree(tree.file)
  }
  stopifnot(!is.null(tree))
  
  last_level <- 
    lapply(tree$tip.label, function(x){tmp <- strsplit(x, "\\|")[[1]]; 
    if(length(tmp) == 1) {NA} else{tmp[length(tmp)-1]}})
  last_level <- do.call(c, last_level)
  tree$tip.label <- last_level
  
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(otu.tab.rff)))
  tree
}