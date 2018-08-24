## countsTable: sample by leafOTU
## @return list pvalue, orig
#' Unified statistical tests for microbiome data
#'
#' @param testName specify the name of the statistical test
#' @param countsTable OTU feature by sample matrix 
#' @param pheno vector, factor, 2 or 3 groups
#' @param cov matrix, can be NULL
#' @param taxaTable phyloseq taxonomyTable (created by tax_table) 
#'
#' @return a list of pvalue and detail. The detail element stores original results
#' @export
#'
#' @examples
#' load("~/Dropbox/xiaowei/data/real_data/curatedMetagenomicData.Rdata", verbose = T)
#'  ## two groups
#' otuTable <- otu_list[[1]]
#' taxaTable <- tax_list[[1]]
#' pheno <- z_list[[1]]
#' countsTable <- otuTable
#' res <- mbTest("permanova", otuTable, pheno, NULL, taxaTable)
#' res <- mbTest("permanovaG", otuTable, pheno, NULL, taxaTable)
#' res <- mbTest("mircat", otuTable, pheno, NULL, taxaTable)
#' res <- mbTest("qcat", otuTable, pheno, NULL, taxaTable)
#' res <- mbTest("mispu", otuTable, pheno, NULL, taxaTable)
#' ## three groups
#' otuTable <- otu_list[[15]]
#' taxaTable <- tax_list[[15]]
#' pheno <- z_list[[15]]
#' countsTable <- otuTable
#' res <- mbTest("anova", otuTable, pheno, NULL, NULL)
#' res <- mbTest("kruskal", otuTable, pheno, NULL, NULL)
#' res <- mbTest("deseq2", otuTable, pheno, NULL, NULL)
#' res <- mbTest("metagenomeSeq", otuTable, pheno, NULL, NULL)
#' res <- mbTest("metagenomeSeq.zig", otuTable, pheno, NULL, NULL)
mbTest <- function(testName = c(
  # test all OTUs
  "t", "anova",
  "wilcox", "kruskal",
  "deseq", "deseq2", "edger",
  "metagenomeSeq",
  "metagenomeSeq.zig", # old method for metagneomeseq zero-inflated gaussian mixture model
  # test composition difference
  "permanova", 
  "permanovaG", 
  "mircat", 
  "qcat",
  "zidgm",
  "mispu"), 
  countsTable, pheno, cov=NULL, 
  taxaTable = NULL) {
  # check data validity
  stopifnot(ncol(countsTable) == length(pheno))
  stopifnot(length(unique(pheno)) >= 2)
  
  # set up the return value
  pheno.level <- levels(factor(pheno))
  startTime <- Sys.time()
  ret <- list(p.value = NA, detail = NA, timing = NA)
  
  # main body
  if (testName == "t") {
    countsTable <- apply(countsTable, 2, function(x) { 
      s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
    ret$detail <- apply(countsTable, 1, function(x) {
      t.test(x[pheno == pheno.level[1]], x[pheno == pheno.level[2]])
    })
    ret$p.value <- sapply(ret$detail, function(x) {x$p.value})
  } else if (testName == "anova") {
    countsTable <- apply(countsTable, 2, function(x) { 
      s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
    ret$detail <- apply(countsTable, 1, function(x) {
      aov( x ~ pheno )
    })
    ret$p.value <- lapply(ret$detail, function(x) {summary(x)[[1]]$`Pr(>F)`[1]})
  } else if (testName == "wilcox") {
    countsTable <- apply(countsTable, 2, function(x) { 
      s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
    ret$detail <- apply(countsTable, 1, function(x) {
      wilcox.test(x[pheno == pheno.level[1]], x[pheno == pheno.level[2]])
    })
    ret$p.value <- sapply(ret$detail, function(x) {x$p.value})
  } else if (testName == "kruskal") {
    countsTable <- apply(countsTable, 2, function(x) { 
      s = sum(x); if (s == 0) {warning("Sample has no counts")}; x / sum(x)} )
    ret$detail <- apply(countsTable, 1, function(x) {
      kruskal.test(x, pheno)
    })
    ret$p.value <- sapply(ret$detail, function(x) {x$p.value})
  } else if (testName == "deseq") {
    warning("DESeq may not work properly")
    library(DESeq)
    #designfac = factor(gsub("[[:print:]]+\\;", "", sample_names(physeq)))
    designfac = pheno
    # Enforce Orientation
    #if (!taxa_are_rows(physeq)) {
    #  physeq <- t(physeq)
    #}
    # Convert to matrix, round up to nearest integer
    #x = ceiling(as(otu_table(physeq), "matrix")) + 1
    x = ceiling(as((countsTable), "matrix")) + 1
    dim(x)
    # Add taxonomy data
    #taxADF = as(data.frame(as(tax_table(physeq), "matrix")), "AnnotatedDataFrame")
    #taxaTable = taxaTable[rownames(taxaTable)%in%rowanames(countsTable),]
    #dim(taxaTable)
    #taxADF = as(data.frame(as(taxaTable, "matrix")), "AnnotatedDataFrame")
    #dim(taxADF)
    # Initalize the count data sets.
    cds = newCountDataSet(x, conditions = designfac) 
    # First, estimate size factors, then estimate dispersion.  Size factors
    cds = estimateSizeFactors(cds)  #, locfunc=genefilter::shorth)
    sizeFactors(cds)
    # Now dispersions Variance estimation, passing along additional options
    cds = estimateDispersions(cds, fitType = c("local"))
    res = nbinomTest(cds, levels(designfac)[1], levels(designfac)[2])
    ret$detail <- (res)
    ret$p.value = res$pval
  } else if (testName == "deseq2") {
    library(DESeq2)
    x = ceiling(as((countsTable), "matrix")) + 1
    dim(x)
    # Add taxonomy data
    #taxADF = as(data.frame(as(tax_table(physeq), "matrix")), "AnnotatedDataFrame")
    #taxaTable = taxaTable[rownames(taxaTable)%in%rownames(countsTable),]
    #dim(taxaTable)
    # Initalize the count data sets.
    dds <- DESeqDataSetFromMatrix(x, DataFrame(pheno), ~pheno)
    dds <- DESeq(dds)
    res <- results(dds)
    ret$detail <- res
    ret$p.value <- res[, "pvalue"]
    names(ret$p.value) <- rownames(res)
  } else if (testName == "edger") {
    library(edgeR)
    y = DGEList(counts = countsTable, group = pheno)  #, remove.zeros=TRUE)
    z = edgeR:::calcNormFactors(y, method = "RLE")
    print(z[[2]]$norm.factors);
    z[[2]]$norm.factors <- size_factor_estimator(t(countsTable), method = "RLE");
    print(size_factor_estimator(t(countsTable), method = "RLE"));
    if (!all(is.finite(z$samples$norm.factors))) {
      stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
    }
    # Estimate dispersions
    z1 = estimateCommonDisp(z)
    z2 = estimateTagwiseDisp(z1)
    et = exactTest(z2)
    ret$detail <- et
    ret$p.value <- et$table$PValue
    names(ret$p.value) <- rownames(et$table)
  } else if (testName == "metagenomeSeq") {
    require("metagenomeSeq")
    
    if(length(unique(pheno)) != 2) {
      stop("metagenomeSeq does not support multiple groups in fitFeatureModel()")
    }
    OTU = as(countsTable, "matrix")
    if (is.null(rownames(countsTable))) {
      rownames(countsTable) <-sprintf("OTU%03d", seq(nrow(countsTable)))
    }
    # Convert sample_data to AnnotatedDataFrame
    ADF = AnnotatedDataFrame(data.frame(pheno, row.names = colnames(countsTable)))
    # define dummy 'feature' data for OTUs, using their name Helps with
    # extraction and relating to taxonomy later on.
    TDF = AnnotatedDataFrame(data.frame(OTUname = rownames(countsTable), row.names =  rownames(countsTable)))
    # Create the metagenomeSeq object
    MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = cumNorm(MGS)
    # fit = fitZig(MGS, model.matrix(~pheno))
    # # You need to specify all OTUs to get the full table from MRfulltable.
    # x = MRfulltable(fit, number = nrow(assayData(MGS)$counts))
    # 
    fit = fitFeatureModel(MGS, model.matrix(~pheno))
    ret$detail <- fit
    ret$p.value <- fit$pvalues
  } else if (testName == "metagenomeSeq.zig") {
    require("metagenomeSeq")
    
    OTU = as(countsTable, "matrix") + 1
    if (is.null(rownames(countsTable))) {
      rownames(countsTable) <-sprintf("OTU%03d", seq(nrow(countsTable)))
    }
    # Convert sample_data to AnnotatedDataFrame
    ADF = AnnotatedDataFrame(data.frame(pheno, row.names = colnames(countsTable)))
    # define dummy 'feature' data for OTUs, using their name Helps with
    # extraction and relating to taxonomy later on.
    TDF = AnnotatedDataFrame(data.frame(OTUname = rownames(countsTable), row.names =  rownames(countsTable)))
    # Create the metagenomeSeq object
    MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = cumNorm(MGS)
    # fit = fitZig(MGS, model.matrix(~pheno))
    # # You need to specify all OTUs to get the full table from MRfulltable.
    # x = MRfulltable(fit, number = nrow(assayData(MGS)$counts))
    # 
    if (length(levels(pheno)) == 2) {
      fit = fitZig(MGS, model.matrix(~pheno))
      ret$detail <- fit
      ret$p.value <- topTable(eBayes(fit$fit), number=nrow(countsTable), sort.by = "none")$P.Value
      names(ret$p.value) <- rownames(countsTable)
    } else {
      levels(pheno) <- sprintf("L%s", levels(pheno))
      mod = model.matrix(~pheno)
      colnames(mod) = levels(pheno)
      fit = fitZig(obj = MGS, mod = mod)
      #fit = fitZig(MGS, model.matrix(~pheno))
      zigFit = fit$fit
      finalMod = fit$fit$design
      
      nLevel <- levels(pheno)
      m.contrast <- outer(nLevel, nLevel, paste)
      contrast <- sub(" ", "-", m.contrast[upper.tri(m.contrast)])
      contrast.matrix = makeContrasts(contrasts = contrast, levels = finalMod) 
      fit2 = contrasts.fit(zigFit, contrast.matrix)
      fit2 = eBayes(fit2)
      res <- topTable(fit2, number=nrow(countsTable), sort.by = "none")
      ret$detail <- res
      ret$p.value <- res$P.Value
      names(ret$p.value) <- rownames(res)
    }
  } else if (testName == "permanova") {
    library(GUniFrac)
    #countsTable <- otuTable
    # Rarefaction
    otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
    source("function-computeUnifrac.R")
    unifracs = computeUnifrac(otu.tab.rff)
    dw <- unifracs[, , "d_1"] # Weighted UniFrac
    du <- unifracs[, , "d_UW"] # Unweighted UniFrac
    dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
    d0 <- unifracs[, , "d_0"]      # GUniFrac with alpha 0
    d5 <- unifracs[, , "d_0.5"]    # GUniFrac with alpha 0.5
    # Permanova - Distance based multivariate analysis of variance
    ret$detail <- adonis(as.dist(d5) ~ pheno)
    ret$p.value <- ((ret$detail)$aov.tab)$`Pr(>F)` [1]
  } else if (testName == "permanovaG") {
    library(GUniFrac)
    #countsTable <- otuTable
    # Rarefaction
    otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
    source("function-computeUnifrac.R")
    unifracs = computeUnifrac(otu.tab.rff)
    unifracs <-
      GUniFrac(t(otu.tab.rff), tree, alpha = c(0, 0.5, 1))$unifracs
    dw <- unifracs[, , "d_1"] # Weighted UniFrac
    du <- unifracs[, , "d_UW"] # Unweighted UniFrac
    dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
    d0 <- unifracs[, , "d_0"]      # GUniFrac with alpha 0
    d5 <- unifracs[, , "d_0.5"]    # GUniFrac with alpha 0.5
    
    ret$detail <- PermanovaG(unifracs[, , c("d_1", "d_UW")] ~ pheno)
    ret$p.value <- ret$detail$p.tab$omni.p.value
  } else if (testName == "mircat") {
    library(MiRKAT)
    # Rarefaction
    otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
    source("function-computeUnifrac.R")
    unifracs = computeUnifrac(otu.tab.rff)
    D.weighted = unifracs[, , "d_1"]
    D.unweighted = unifracs[, , "d_UW"]
    D.BC = as.matrix(vegdist(t(otu.tab.rff) , method = "bray"))
    K.weighted = D2K(D.weighted)
    K.unweighted = D2K(D.unweighted)
    K.BC = D2K(D.BC)
    
    # Y = matrix(rnorm(n * 3, 0, 1), n , 3)
    #MiRKAT(Y = t(countsTable),
    #        K = K.weighted,
    #        X = cov))
    Ks <- list(K.unweighted, K.BC)
    ret$detail <- MiRKAT(as.numeric(pheno) - 1, X = NULL, Ks = Ks, out_type = "D")
    ret$p.value <- ret$detail$omnibus_p
  } else if (testName == "qcat") {
    library(miLineage)
    ## TODO: QCAT support continuous outcomes only?
    ret$detail <- QCAT(t(countsTable), matrix(as.numeric(pheno), nc = 1), 1, NULL, fdr.alpha = 0.05)
    ret$p.value <-  ret$detail$pval
  } else if (testName == "zidgm") {
    library(miLineage)
    ret$detail <- ZIGDM(
      t(countsTable),
      NULL, # cov related to presence-absence
      matrix(as.numeric(pheno), nc = 1), # cov related to mean abundance
      NULL, # cov related to dispersion 
      test.type = "Mean",
      1,
      ZI.LB = 10,
      NULL,
      fdr.alpha = 0.05
    )
    ret$p.value = ret$detail$global.pval
  } else if (testName == "mispu") {
    library(MiSPU)
    otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
    
    tree.file <- system.file("extdata/metaphlan2_selected.tree.reroot.nwk.bz2", package = "curatedMetagenomicData")
    stopifnot(file.exists(tree.file))
    tree <- read.tree(tree.file)
    last_level <- 
      lapply(tree$tip.label, function(x){tmp <- strsplit(x, "\\|")[[1]]; 
      if(length(tmp) == 1) {NA} else{tmp[length(tmp)-1]}})
    last_level <- do.call(c, last_level)
    tree$tip.label <- last_level
    
    tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(ot)))
    otu.tab.rff <- otu.tab.rff [rownames(otu.tab.rff) %in% tree$tip.label, ]
    
    ret$detail = MiSPU(
      as.numeric(pheno) - 1,
      t(otu.tab.rff),
      tree,
      NULL,
      model =  "binomial",
      pow = c(2:8, Inf),
      n.perm = 1000
    )
    ret$p.value <- ret$detail$aMiSPU$pvalue
  }
  endTime <- Sys.time()
  ret$timing <- difftime(endTime, startTime, units = "secs")
  ret
}

if (FALSE) {
  ## append taxonomy info
  setwd("~/Dropbox/xiaowei/data/real_data/")
  load("curatedMetagenomicData.Rdata", verbose = T)
  library(phyloseq)
  str(tax_list[[14]])
  head(tax_list[[14]]@.Data)
  str(Y_list[[15]])
  tax_list[[15]] <- tax_list[[14]] ## change this
  str(tax_list[[15]])
}

# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
