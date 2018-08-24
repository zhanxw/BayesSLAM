real_data_preprocessor = function(otu, tax) {
  # Find out the basic taxonomic level
  n <- dim(tax)[1];
  p <- dim(tax)[2];
  index_base <- rep(TRUE, n);
  names(index_base) <- rownames(tax);
  for (i in 1:n) {
    for (j in 1:(p - 1)) {
      if (is.na(tax[i, j + 1])) {
        break;
      } else {
        index_base[paste0(tolower(substr(colnames(tax)[j], 1, 1)), "__", tax[i, j])] <- FALSE;
      }
    }
  }
  
  # Generate the aggregate structure
  S <- matrix(0L, ncol = sum(index_base), nrow = sum(1 - index_base));
  colnames(S) <- rownames(tax)[which(index_base)];
  rownames(S) <- rownames(tax)[which(!index_base)];
  for (i in which(index_base)) {
    for (j in 1:p) {
      if (!is.na(tax[i, j]) && rownames(tax)[i] != paste0(tolower(substr(colnames(tax)[j], 1, 1)), "__", tax[i, j])) {
        S[paste0(tolower(substr(colnames(tax)[j], 1, 1)), "__", tax[i, j]), rownames(tax)[i]] <- 1L;
      }
    }
  }
  S_x <- matrix(0L, nrow = dim(S)[1], ncol = max(rowSums(S)));
  for (ii in 1:dim(S)[1]) {
    S_x[ii, 1:sum(S[ii,] == 1)] <- which(S[ii,] == 1)
  }
  rownames(S_x) <- rownames(S);
  
  # Return the count matrix in bottom taxonomic level
  Y <- t(otu[which(index_base),]);
  
  return(list(Y = Y, S = S_x));
}

b_estimator = function(Y, z, S, aggregate, s_hat, a = 2, a_0 = 2) {
  n <- dim(Y)[1];
  p <- dim(Y)[2];
  K <- length(unique(z));
  if (aggregate) {
    p_ext <-dim(S)[1];
    pp <- dim(S)[2];
  } else {
    p <- 0;
  }
  Y_temp <- matrix(0, nrow = n, ncol = p + p_ext);
  Y_temp[, 1:p] <- Y;
  if (aggregate) {
    for (j in 1:p_ext) {
      index <- S[j, S[j,] != 0];
      if (length(index) == 1) {
        Y_temp[, p + j] <- Y[, index];
      } else {
        Y_temp[, p + j] <- rowSums(Y[, index]);
      }
    }
  }
  A <- Y_temp/s_hat;
  A[which(A == 0)] <- NA;
  logA <- log(A);
  b_0 <- (a_0 - 1)*apply(logA, 2, var, na.rm = TRUE);
  b_0[which(is.na(b_0))] <- 1;
  B <- matrix(1, nrow = K, ncol = p + p_ext);
  for (k in 1:K) {
    B[k,] <- apply(logA[which(z == k - 1),], 2, var, na.rm = TRUE);
  }
  B[which(is.na(B))] <- 1;
  
  return(list(B = B, b_0 = b_0));
}

BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
  }
  return(PPI_sorted[k])
}


# =============================================================================================
# The function of estimating size factors from count data
# ---------------------------------------------------------------------------------------------
# Input:  Y, a n-by-p count matrix, where n is the number of samples and p is the number of 
#         features
# Input:  method, a categorical variable chosen from
#         TC:   total counts
#         CSS:  cumulative-sum scaling by Paulson et al., 2013
#         TMM:  trimmed mean by M-values by Robinson and Oshlack, 2010
#         RLE:  relative log expression by Anders and Huber, 2010
#         Q75:  the 75th percentile by Bullard et al., 2010
#         none: all size factors are set to 1
# Input:  constraint, a boolean variable, with TRUE rescaling size factors to their sum equal 1
# ---------------------------------------------------------------------------------------------
# Output: s, a n-dimensional vector, where each element is the estimated size factor for the 
#         corresponding sample
# =============================================================================================
size_factor_estimator = function(Y, method = c("TSS", "CSS", "TMM", "RLE", "Q75")) {
  n <- dim(Y)[1];
  s <- rep(NA, n);
  if (method == "TSS") {
    s <- rowSums(Y);
  } else if (method == "CSS") {
    p <- cumNormStatFast_copy(t(Y));
    for (i in 1:n) {
      temp <- Y[i, which(Y[i,] != 0)];
      s[i] <- sum(Y[i, which(Y[i,] <= quantile(temp, p))]);
    }
  } else if (method == "TMM") {
    s <- calcNormFactors_copy(t(Y), method = method, refColumn = NULL, logratioTrim = 0.3, 
                              sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10);
  } else if (method == "RLE") {
    s <- calcNormFactors_copy(t(Y), method = method);
    # if (sum(is.na(s)) > 0) {
    #   Y_temp <- Y;
    #   Y_temp[which(Y_temp == 0)] <- NA;
    #   geomeans <- exp(colMeans(log(Y_temp), na.rm = TRUE));
    #   Y_temp <- t(t(Y_temp)/geomeans);
    #   s <- apply(Y_temp[,which(geomeans > 0)], 1, median, na.rm = TRUE);
    # }
  } else if (method == "Q75") {
    s <- apply(Y, 1, quantile, 0.75);
    s <- pmax(s, 1);
    # s <- calcNormFactors_copy(t(Y), method = "upperquartile");
  }
  # if (constraint) {
  #   s <- s/sum(s);
  # }
  names(s) <- rownames(Y);
  return(s);
}

# =============================================================================================
# The functions attached to size_factor_estimator(...)
# =============================================================================================
rescale <- function(x, constraint = c("sum", "product")) {
  if (constraint == "sum") {
    return(x/sum(x));
  } else if (constraint == "product") {
    return(x/exp(mean(log(x))))
  }
} 

# Download from https://github.com/HCBravoLab/metagenomeSeq/blob/master/R/cumNormStatFast.R on 
# April 12, 2018
cumNormStatFast_copy <- function(mat, pFlag = FALSE, rel = 0.1, ...){
  # mat = returnAppropriateObj(obj, FALSE, FALSE)
  smat = lapply(1:ncol(mat), function(i) {
    sort(mat[which(mat[, i] > 0), i], decreasing = TRUE)
  })
  leng = max(sapply(smat, length))
  if(any(sapply(smat, length) == 1)) {
    stop("Warning sample with one or zero features")
  }
  smat2 = array(NA, dim = c(leng, ncol(mat)))
  for(i in 1:ncol(mat)){
    smat2[leng:(leng - length(smat[[i]]) + 1), i] = smat[[i]]
  }
  rmat2 = sapply(1:ncol(smat2), function(i){
    quantile(smat2[, i], p = seq(0, 1, length.out = nrow(smat2)), na.rm = TRUE)
  })
  smat2[is.na(smat2)] = 0
  ref1 = rowMeans(smat2)
  ncols = ncol(rmat2)
  diffr = sapply(1:ncols, function(i) {
    ref1 - rmat2[, i]
  })
  diffr1 = matrixStats::rowMedians(abs(diffr))
  if(pFlag == TRUE){
    plot(abs(diff(diffr1))/diffr1[-1], type = "h", ...)
    abline(h = rel)
    axis(1, at = seq(0, length(diffr1), length.out = 5), labels = seq(0, 1, length.out = 5))
  }
  x= which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
  if(x <= 0.50){
    message("Default value being used")
    x = 0.50
  }
  # if(class(obj) == "MRexperiment"){
  #   obj@expSummary$cumNormStat = x		
  # }
  return(x)
}

# Download from https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R on April 12, 2018
calcNormFactors_copy <- function(object, lib.size = NULL, method = c("TMM", "RLE", "upperquartile",
                                                                     "none"), 
                                 refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, 
                                 doWeighting = TRUE, Acutoff = -1e10, p = 0.75, ...) {
  x <- as.matrix(object)
  if (any(is.na(x))) {
    stop("NA counts not permitted")
  }
  if (is.null(lib.size)) {
    lib.size <- colSums(x)
  }
  if (any(is.na(lib.size))) {
    stop("NA lib.sizes not permitted")
  }
  method <- match.arg(method)
  allzero <- .rowSums(x>0, nrow(x), ncol(x)) == 0
  if(any(allzero)) {
    x <- x[!allzero, , drop = FALSE]
  }
  if(nrow(x) == 0 || ncol(x) == 1) {
    method = "none"
  }
  f <- switch(method,
              TMM = {
                f75 <- .calcFactorQuantile_copy(data = x, lib.size = lib.size, p = 0.75)
                if(is.null(refColumn)) {
                  refColumn <- which.min(abs(f75 - mean(f75)))
                }
                if(length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x)) {
                  refColumn <- 1
                }
                f <- rep(NA, ncol(x))
                for(i in 1:ncol(x)) {
                  f[i] <- .calcFactorWeighted_copy(obs = x[, i], ref = x[, refColumn], 
                                                   libsize.obs = lib.size[i], 
                                                   libsize.ref = lib.size[refColumn], 
                                                   logratioTrim = logratioTrim, sumTrim = sumTrim, 
                                                   doWeighting = doWeighting, Acutoff = Acutoff)
                }
                # f
                f*lib.size
              },
              # RLE = .calcFactorRLE(x)/lib.size,
              RLE = .calcFactorRLE_copy(x),
              upperquartile = .calcFactorQuantile_copy(x, lib.size, p = p),
              none = rep(1, ncol(x))
  )
  # f <- f/exp(mean(log(f)))
  f
}

.calcFactorRLE_copy <- function (data) {
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile_copy <- function (data, lib.size, p = 0.75)
{
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
}

.calcFactorWeighted_copy <- function(obs, ref, libsize.obs = NULL, libsize.ref = NULL, 
                                     logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
                                     Acutoff = -1e10) {
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs)) {
    nO <- sum(obs) 
  } else {
    nO <- libsize.obs
  }
  if (is.null(libsize.ref)) {
    nR <- sum(ref)
  } else {
    nR <- libsize.ref
  }
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  if(max(abs(logR)) < 1e-6) {
    return(1)
  }
  n <- length(logR)
  loL <- floor(n*logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n*sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)
  if (doWeighting) {
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep], na.rm = TRUE)
  } else {
    f <- mean(logR[keep], na.rm=TRUE)
  }
  if(is.na(f)) {
    f <- 0
  }
  2^f
}