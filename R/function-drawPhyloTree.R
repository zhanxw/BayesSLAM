## Draw phylogenetic tree from taxonomy table
## Xiaowei Zhan 
## August 2018

if (FALSE) {
  options(stringsAsFactors = F)
  suppressPackageStartupMessages(library(phyloseq))
  suppressPackageStartupMessages(library(curatedMetagenomicData))
  suppressPackageStartupMessages(library(data.tree))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggtree))
}
#' Draw phylogenetic tree in circular layout
#'
#' @param taxa phyloseq taxonomyTable class or a data frame
#' @param highlight a vector of node names to be colored, e.g. "k__Bacteria", 
#'                  or a data frame with optionally columns (col, size, alpha)
#' @param title Title for the plot
#'
#' @return an ggplot2 object
#' @export
#'
#' @examples
#' if (FALSE) {
#'   res <- curatedMetagenomicData("Castro-*metaphlan*", dryrun=FALSE) #one dataset
#'   taxa<-tax_table(ExpressionSet2phyloseq(res[[1]]))
#' } else {
#'   prefix = system.file("extdata", package = "BayesSLAM")
#'   load(file.path(prefix, "Castro-NallarE_2015.Rdata"))
#'   taxa <- taxaTable
#' }
#' drawPhyloTree(taxa)
#' drawPhyloTree(taxa, highlight = "p__Bacteroidetes")
#' drawPhyloTree(taxa, highlight = data.frame(taxa = c("p__Bacteroidetes", "s__Prevotella_salivae"), size = c(8, 5), col = c("red", "blue"), alpha = c(.6, .7)))
drawPhyloTree <- function(taxa, highlight = NULL, title = NULL) {
  stopifnot(class(taxa) %in% c("taxonomyTable", "data.frame"))
  if (class(taxa) == "taxonomyTable") {
    taxa <- data.frame(taxa@.Data)
  }
  stopifnot(all(c("Kingdom", "Phylum", "Class", 
                  "Order", "Family", "Genus", "Species") %in% colnames(taxa) ))
  # Convert taxonomyTable to data.tree
  taxa$pathString <- sprintf("Root/k__%s/p__%s/c__%s/o__%s/f__%s/g__%s/s__%s",
                             taxa$Kingdom,
                             taxa$Phylum,
                             taxa$Class,
                             taxa$Order,
                             taxa$Family, 
                             taxa$Genus,
                             taxa$Species)
  taxa$pathString <- gsub(pattern = "/.__NA", replacement = "", taxa$pathString)
  taxa <- subset(taxa, !grepl('^t__', rownames(taxa)))
  root = data.tree::as.Node(taxa)
  
  # Convert to Newick format
  root.newick <- data.tree::ToNewick(root,heightAttribute = NULL)
  # Codes below are for customizing tree branches lengths
  # root.newick <- ToNewick(root, heightAttribute = function(x){
  #   #print(str(x$name)); 
  #   name = substr(x$name, 1, 1); 
  #   #print (x);
  #   ret = switch(name, T = 9, # Tree is the root
  #                k = 7, 
  #                p = 6, c = 5, o = 4, f = 3, g = 2, s = 1  );
  #   cat(x$name, " ", name, " ", ret, "\n");
  #   ret})
  tmp.fn <- tempfile()
  cat(root.newick, file = tmp.fn)
  
  # Read the read into the phylo format
  tree <- ape::read.tree(tmp.fn)
  # Plot the tree
  p <- ggtree(tree, layout = "circular")
  if (!is.null(highlight)) {
    if (class(highlight) == "character") {
      # Specify node attributes to highlight
      dd <-  data.frame(taxa = highlight)
    } else if (class(highlight) == "data.frame") {
      stopifnot(ncol(highlight) >= 1)
      stopifnot(!colnames(highlight)[1] %in% c("col", "size", "alpha"))
      dd <- highlight
      dd$taxa = highlight[,1]
    }
    # specify color, size and alpha
    cname = colnames(dd)
    if (! "col" %in% cname) {
      if (length(highlight) > 0) {
        dd$col = "red" 
      } else {
        dd$col = highlight
      }
    }
    if (! "size" %in% cname) {
      dd$size = 10
    }
    if (! "alpha" %in% cname) {
      dd$alpha = .8
    }
    p <- p %<+% dd
    groups <- unique(dd[, c("size", "col", "alpha")])
    ## note, below is tedious, as we cannot add geom_point2 programmatically
    extractAes <- function(dd, groups, p, i) {
      . <- subset(dd, 
                  size == groups[i, "size"] &
                    col == groups[i, "col"] &
                    alpha == groups[i, "alpha"])
      print(.)
      ret <- p$data
      ret$selected = FALSE
      ret$col = ret$size = ret$alpha = rep(NA, nrow(ret))
      nodeIdx = p$data$node[p$data$label %in% .$taxa ]
      ret$selected[nodeIdx] = TRUE
      ret$col[nodeIdx] = groups[i, "col"]
      ret$size[nodeIdx] = groups[i, "size"]
      ret$alpha[nodeIdx] = groups[i, "alpha"]
      return (ret)
    }
    cmd <- c('p <- p ', rep('', nrow(groups)))
    for (i in 1:nrow(groups)) {
      assign(sprintf("dat%d", i), extractAes(dd, groups, p, i))
      cmd[i+1] <- (sprintf(" + geom_point2(data = dat%d,
                           aes(subset = selected,
                               col = col, 
                               size = size,
                               alpha = alpha))", i))
    }
    eval(parse(text = (paste(cmd, collapse = ''))))
    p <- p + scale_size_identity() + scale_alpha_identity() + scale_color_identity()
  }
  p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5));
  return(p)
}

if (FALSE) {
  library(cowplot)
  load("~/Downloads/pvals.RData")
  str(as.data.frame(all_taxa_df))
  all_taxa_df <- all_taxa_df[!duplicated(all_taxa_df),]
  plots <- lapply(seq_along(pvals), function(i) {
    #plots <- lapply(seq(2), function(i) {
    print(i)
    p <- pvals[[i]]
    highlight = data.frame(taxa = names(p [ p < 0.05 & !is.na(p)]),
                           col = "red", alpha = .6, size = 5)
    highlight$col[highlight$taxa %in%
                    c("s__Faecalibacterium_prausnitzii")] = "blue"
    highlight$col[highlight$taxa %in%
                    c("s__Holdemania_filiformis")] = "green"
    highlight$col[highlight$taxa %in%
                    c("s__Bacteroides_thetaiotaomicron")] = "yellow"
    
    drawPhyloTree(taxa = as.data.frame(all_taxa_df),
                  highlight = highlight,
                  title = fn[i])
  })  
  plots2 <- plot_grid(plotlist = plots, ncol = 2)
  ggsave("plot.findings.scale.t5.pdf", plots2, width = 8, height = 16)
}