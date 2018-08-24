## Draw phylogenetic tree from taxonomy table
## Xiaowei Zhan 
## August 2018
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(curatedMetagenomicData))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))

#' Draw phylogenetic tree in circular layout
#'
#' @param taxa phyloseq taxonomyTable class or a data frame
#' @param highlight a vector of node names to be colored, e.g. "k__Bacteria", 
#'                  or a data frame with optionally columns (col, size, alpha)
#'
#' @return an ggplot2 object
#' @export
#'
#' @examples
#' res <- curatedMetagenomicData("Castro-*metaphlan*", dryrun=FALSE) #one dataset
#' taxa<-tax_table(ExpressionSet2phyloseq(res[[1]]))
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
  root.newick <- ToNewick(root,heightAttribute = NULL)
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
  tree <- read.tree(tmp.fn)
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
    aesList <- list()
    extractAes <- function(dd, groups, p, i) {
      . <- subset(dd, 
                  size == groups[i, "size"] &
                    col == groups[i, "col"] &
                    alpha == groups[i, "alpha"])
      print(.)
      return (list(nodeIdx = p$data$node[p$data$label %in% .$taxa ],
                   col = groups[i, "col"],
                   size = groups[i, "size"],
                   alpha = groups[i, "alpha"]))
    }
    if (nrow(groups) == 1) {
      aesList[[1]] <- extractAes(dd, groups, p, 1)
      p <- p + 
        geom_point2(aes(subset = (node %in% aesList[[1]]$nodeIdx), 
                        col = aesList[[1]]$col, 
                        size = aesList[[1]]$size,
                        alpha = aesList[[1]]$alpha))
    } else if (nrow(groups) == 2) {
      aesList[[1]] <- extractAes(dd, groups, p, 1)
      aesList[[2]] <- extractAes(dd, groups, p, 2)
      p <- p + 
        geom_point2(aes(subset = (node %in% aesList[[1]]$nodeIdx), 
                        col = aesList[[1]]$col, 
                        size = aesList[[1]]$size,
                        alpha = aesList[[1]]$alpha)) +
        geom_point2(aes(subset = (node %in% aesList[[2]]$nodeIdx), 
                        col = aesList[[2]]$col, 
                        size = aesList[[2]]$size,
                        alpha = aesList[[2]]$alpha))
    } else if (nrow(groups) == 3) {
      aesList[[1]] <- extractAes(dd, groups, p, 1)
      aesList[[2]] <- extractAes(dd, groups, p, 2)
      aesList[[3]] <- extractAes(dd, groups, p, 3)
      p <- p + 
        geom_point2(aes(subset = (node %in% aesList[[1]]$nodeIdx), 
                        col = aesList[[1]]$col, 
                        size = aesList[[1]]$size,
                        alpha = aesList[[1]]$alpha)) +
        geom_point2(aes(subset = (node %in% aesList[[2]]$nodeIdx), 
                        col = aesList[[2]]$col, 
                        size = aesList[[2]]$size,
                        alpha = aesList[[2]]$alpha)) +
        geom_point2(aes(subset = (node %in% aesList[[3]]$nodeIdx), 
                        col = aesList[[3]]$col, 
                        size = aesList[[3]]$size,
                        alpha = aesList[[3]]$alpha))    
    } else {
      warning("Only the fist 3 aes themes are implemented")
      aesList[[1]] <- extractAes(dd, groups, p, 1)
      aesList[[2]] <- extractAes(dd, groups, p, 2)
      aesList[[3]] <- extractAes(dd, groups, p, 3)
      p <- p + 
        geom_point2(aes(subset = (node %in% aesList[[1]]$nodeIdx), 
                        col = aesList[[1]]$col, 
                        size = aesList[[1]]$size,
                        alpha = aesList[[1]]$alpha)) +
        geom_point2(aes(subset = (node %in% aesList[[2]]$nodeIdx), 
                        col = aesList[[2]]$col, 
                        size = aesList[[2]]$size,
                        alpha = aesList[[2]]$alpha)) +
        geom_point2(aes(subset = (node %in% aesList[[3]]$nodeIdx), 
                        col = aesList[[3]]$col, 
                        size = aesList[[3]]$size,
                        alpha = aesList[[3]]$alpha))  
    }
    print(aesList)
    print(p)
    # dd %>% group_by(size, col, alpha) %>% do({
    #   nodeIdx <- p$data$node[p$data$label %in% .$taxa ]
    #   col_ = .$col[1]
    #   size_ = .$size[1]
    #   alpha_ = .$alpha[1]
    #   p <<-       p + 
    #            geom_point2(aes(subset = (node %in% nodeIdx), 
    #                     col = col_, 
    #                     size = size_,
    #                     alpha = alpha_))
    # 
    #       cat("highlight", nodeIdx, " with ", col_, size_, alpha_, "\n")
    #       print(p)
    #       data.frame()
    #     })
    p <- p + scale_size_identity() + scale_alpha_identity() + scale_color_identity()
  }
  p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5));
  return(p)
}

if (FALSE) {
  load("~/Downloads/pvals.RData")
  str(as.data.frame(all_taxa_df))
  all_taxa_df <- all_taxa_df[!duplicated(all_taxa_df),]
  plots <- lapply(seq_along(pvals), function(i) {
    # plots <- lapply(seq(2), function(i) {
    print(i)
    p <- pvals[[i]]
    # if (i != ) {
    #   highlight = data.frame(taxa = names(p [ p < 0.05 & !is.na(p)]),
    #                          col = "red", alpha = .6, size = 5)
    #   drawPhyloTree(taxa = as.data.frame(all_taxa_df),
    #                 highlight = highlight,
    #                 title = fn[i, 1])
    # } else {
    highlight = data.frame(taxa = names(p [ p < 0.05 & !is.na(p)]),
                           col = "red", alpha = .6, size = 5)
    highlight$col[highlight$taxa %in%
                    c("s__Faecalibacterium_prausnitzii",
                      "s__Holdemania_filiformis",
                      "s__Bacteroides_thetaiotaomicron")] = "blue"
    highlight$size[highlight$taxa %in%
                     c("s__Faecalibacterium_prausnitzii",
                       "s__Holdemania_filiformis",
                       "s__Bacteroides_thetaiotaomicron")] = 6
    
    drawPhyloTree(taxa = as.data.frame(all_taxa_df),
                  highlight = highlight,
                  title = fn[i, 1])
    # }
  })  
  plots2 <- plot_grid(plotlist = plots, ncol = 2)
  ggsave("plot.findings.t5.pdf", plots2, width = 8, height = 16)
}
