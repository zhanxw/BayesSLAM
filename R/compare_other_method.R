## countsTable: sample by leafOTU
## @return list pvalue, orig
#' Unified statistical tests for microbiome data
#'
#' @param testName specify the name of a statistical test
#' @param countsTable OTU feature by sample matrix 
#' @param pheno vector, factor, 2 or 3 groups
#' @param cov matrix, can be NULL
#' @param taxaTable phyloseq taxonomyTable (created by tax_table) 
#'
#' @return a list of pvalue, detail and timing. The detail element stores original results
#' @export
#'
#' @examples
#'  ## two groups
#' if (FALSE) {
#'   load("~/Dropbox/xiaowei/data/real_data/curatedMetagenomicData.Rdata", verbose = TRUE)
#'   otuTable <- otu_list[[1]]
#'   taxaTable <- tax_list[[1]]
#'   pheno <- z_list[[1]]
#' } else {
#'   prefix = system.file("extdata", package = "BayesSLAM")
#'   load(file.path(prefix, "Castro-NallarE_2015.Rdata"))
#' }
#' countsTable <- otuTable
#' res <- slam("anova", otuTable, pheno, NULL, NULL)
#' res <- slam("kruskal", otuTable, pheno, NULL, NULL)
#' res <- slam("deseq2", otuTable, pheno, NULL, NULL)
#' res <- slam("edger", otuTable, pheno, NULL, NULL)
#' res <- slam("metagenomeSeq", otuTable, pheno, NULL, NULL)
#' res <- slam("metagenomeSeq.zig", otuTable, pheno, NULL, NULL)
#'
#' res <- slam("permanova", otuTable, pheno, NULL, taxaTable)
#' res <- slam("permanovaG", otuTable, pheno, NULL, taxaTable)
#' res <- slam("mircat", otuTable, pheno, NULL, taxaTable)
#' res <- slam("qcat", otuTable, pheno, NULL, taxaTable)
#' res <- slam("mispu", otuTable, pheno, NULL, taxaTable)
#' 
#' ## three groups
#' if (FALSE) {
#'   load("~/Dropbox/xiaowei/data/real_data/curatedMetagenomicData.Rdata", verbose = TRUE)
#'   otuTable <- otu_list[[15]]
#'   taxaTable <- tax_list[[15]]
#'   pheno <- z_list[[15]]
#' } else {
#'   prefix = system.file("extdata", package = "BayesSLAM")
#'   load(file.path(prefix, "ZellerG_2014.Rdata"))
#' }
#' countsTable <- otuTable
#' res <- slam("anova", otuTable, pheno, NULL, NULL)
#' res <- slam("kruskal", otuTable, pheno, NULL, NULL)

slam <- function(testName = c(
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
  if (!checkPrereq(testName)) {
    warning("Please install required packages for ", testName)    
    return (NULL)
  }
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
    cds = DESeq::newCountDataSet(x, conditions = designfac) 
    # First, estimate size factors, then estimate dispersion.  Size factors
    cds = DESeq::estimateSizeFactors(cds)  #, locfunc=genefilter::shorth)
    DESeq::sizeFactors(cds)
    # Now dispersions Variance estimation, passing along additional options
    cds = DESeq::estimateDispersions(cds, fitType = c("local"))
    res = DESeq::nbinomTest(cds, levels(designfac)[1], levels(designfac)[2])
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
    dds <- DESeq2::DESeqDataSetFromMatrix(x, DataFrame(pheno), ~pheno)
    dds <- DESeq2::DESeq(dds)
    res <- results(dds)
    ret$detail <- res
    ret$p.value <- res[, "pvalue"]
    names(ret$p.value) <- rownames(res)
  } else if (testName == "edger") {
    library(edgeR)
    y = edgeR::DGEList(counts = countsTable, group = pheno)  #, remove.zeros=TRUE)
    z = edgeR::calcNormFactors(y, method = "RLE")
    # We make sure edgeR normalization makes sense
    print(z[[2]]$norm.factors);
    z[[2]]$norm.factors <- size_factor_estimator(t(countsTable), method = "RLE");
    print(size_factor_estimator(t(countsTable), method = "RLE"));
    if (!all(is.finite(z$samples$norm.factors))) {
      stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
    }
    # Estimate dispersions
    z1 = edgeR::estimateCommonDisp(z)
    z2 = edgeR::estimateTagwiseDisp(z1)
    et = edgeR::exactTest(z2)
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
    MGS = metagenomeSeq::newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = metagenomeSeq::cumNorm(MGS)
    # fit = fitZig(MGS, model.matrix(~pheno))
    # # You need to specify all OTUs to get the full table from MRfulltable.
    # x = MRfulltable(fit, number = nrow(assayData(MGS)$counts))
    # 
    fit = metagenomeSeq::fitFeatureModel(MGS, model.matrix(~pheno))
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
    MGS = metagenomeSeq::newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = metagenomeSeq::cumNorm(MGS)
    # fit = fitZig(MGS, model.matrix(~pheno))
    # # You need to specify all OTUs to get the full table from MRfulltable.
    # x = MRfulltable(fit, number = nrow(assayData(MGS)$counts))
    # 
    if (length(levels(pheno)) == 2) {
      fit = metagenomeSeq::fitZig(MGS, model.matrix(~pheno))
      ret$detail <- fit
      ret$p.value <- topTable(eBayes(fit$fit), number=nrow(countsTable), sort.by = "none")$P.Value
      names(ret$p.value) <- rownames(countsTable)
    } else {
      levels(pheno) <- sprintf("L%s", levels(pheno))
      mod = model.matrix(~pheno)
      colnames(mod) = levels(pheno)
      fit = metagenomeSeq::fitZig(obj = MGS, mod = mod)
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
    # source("function-computeUnifrac.R")
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
    # source("function-computeUnifrac.R")
    unifracs = computeUnifrac(otu.tab.rff)
    # unifracs <-
    #   GUniFrac(t(otu.tab.rff), tree, alpha = c(0, 0.5, 1))$unifracs
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
    # source("function-computeUnifrac.R")
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
    ret$detail <- QCAT(t(countsTable), matrix(as.numeric(pheno), ncol = 1), 1, NULL, fdr.alpha = 0.05)
    ret$p.value <-  ret$detail$pval
  } else if (testName == "zidgm") {
    library(miLineage)
    ret$detail <- ZIGDM(
      t(countsTable),
      NULL, # cov related to presence-absence
      matrix(as.numeric(pheno), ncol = 1), # cov related to mean abundance
      NULL, # cov related to dispersion 
      test.type = "Mean",
      1,
      ZI.LB = 10,
      NULL,
      fdr.alpha = 0.05
    )
    ret$p.value = ret$detail$global.pval
  } else if (testName == "mispu") {
    library(GUniFrac)
    library(MiSPU)
    otu.tab.rff <- Rarefy(countsTable)$otu.tab.rff
    
    tree <- computeMetaphlanTree(rownames(otu.tab.rff))
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


#' Check if the required packages are installed
#'
#' @param testName a name of a statistical test
#'
#' @return TRUE only all required packages are installed
#' @export
#'
#' @examples
#' checkPrereq("deseq")
#' checkPrereq("mispu")
checkPrereq <- function(testName) {
  d <- "test package
t stats
  anova stats
  wilcox stats
  kruskal stats
  deseq DESeq
  deseq2 DESeq2
  edger edgeR
  metagenomeseq metagenomeSeq
  metagenomeSeq.zig metagenomeSeq
  permanova GUniFrac
  permanovaG GUniFrac
  mircat MiRKAT
  qcat miLineage
  zidgm miLineage
  mispu MiSPU,GUniFrac"
  d <- read.table(text = d, header = TRUE, stringsAsFactors = F)
  if (!testName %in% d$test) {
    warning("Specified test [", testName, "] is not yet supported.")
    return(FALSE)
  }
  
  needed.pkg <- strsplit(x = d$package[d$test == testName], ",")[[1]]
  installed <- installed.packages()
  isSucc = TRUE
  for (pkg in needed.pkg) {
    if (pkg %in% installed[, "Package"]) {next}
    warning("Please install ", pkg, " in order to use ", testName)
    isSucc = FALSE
  }
  return(isSucc)
}
# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
