library(zoo);
library(lattice);

####### gamma PPI #######
gamma.PPI.draw = function(gamma.PPI, tax.name.full, alpha.sig = 0.05){
  library(ggplot2)
  p.size = length(tax.name.full)
  dim(gamma.PPI) = c( 1,length(gamma.PPI))
  if(is.null(colnames(gamma.PPI))){
    colnames(gamma.PPI) = tax.name.full
    taxa.nam.ori = tax.name.full
  }else{
    taxa.nam.ori = colnames(gamma.PPI)
  }
  # (1)order the taxa by taxonomic level
  rank.od = c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  ind.by = match(substring(taxa.nam.ori, 1, 3), rank.od)
  ord = order(ind.by)
  ord.nam = taxa.nam.ori[ord]
  PPI.order = gamma.PPI[ord]; dim(PPI.order) = c( 1,length(gamma.PPI)); colnames(PPI.order) = ord.nam
  #taxon.name.sel = ord.nam[which(gamma.PPI > cutoff.value)]
  # 
  rank.group.size = rep(NA, length(rank.od))
  for(i in 1:length(rank.od)){
    rank.group.size[i] = length(grep(rank.od[i], ord.nam, value=TRUE))
    x.str = c(0, cumsum(rank.group.size)+1)[1:length(rank.od)] 
    x.end = cumsum(rank.group.size) + 1
  }
  rect.df = data.frame(x1 = x.str, x2 = x.end, y1 = rep(0, length(rank.od)), y2 = rep(1, length(rank.od)),
                       rank = c("kingdom", "phylum", "class", "order", 
                                "family", "genus", "species"))
  # # (2)order each taxonomic level by Gamma PPI
  # level.pnt = data.frame(str = rect.df$x1+1, ed = rect.df$x2)
  # # modification:
  # level.pnt$str[2] = level.pnt$str[2]-1; level.pnt$ed[c(1, dim(level.pnt)[1])] = level.pnt$ed[c(1, dim(level.pnt)[1])] - 1
  # 
  # for(i in 1:length(rank.od)){
  #   PPI.tmp = PPI.order[,level.pnt$str[i]:level.pnt$ed[i]]
  #   PPI.order[,level.pnt$str[i]:level.pnt$ed[i]] = PPI.tmp[order(PPI.tmp,decreasing = T)]
  # }
  # #
  cutoff.value = BayFDR(gamma.PPI, alpha.sig)
  ggplot() +
    geom_rect(data=rect.df,
              aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=rank), 
              alpha=0.35, size = 0) +
    geom_segment( color='gray60',aes(x=seq_along(PPI.order),
                                     xend=seq_along(PPI.order), y=0, yend = as.vector(PPI.order))) +
    geom_point(aes(x=seq_along(PPI.order), y= as.vector(PPI.order))) +
    xlab("Taxa index") + ylab("Posterior probabilities of inclusion") +
    geom_hline(aes(yintercept = cutoff.value), colour="#990000", linetype="dashed") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer( limits=c("kingdom","phylum","class","order", "family", "genus","species"), 
                       palette="Spectral") 
}

gmm_ppi = function(alpha, z, resolution = 100) {
  # Set hyperparameters
  a = 2;
  a_0 = 2;
  
  # Estimate b and b_0
  b_0_hat <- (a_0 - 1)*var(alpha);
  temp <- rep(NA, length(unique(z)));
  count <- 1;
  for (k in unique(z)) {
    temp[count] <- var(alpha[z == k]);
    count <- count + 1;
  }
  b_hat <- (a - 1)*max(temp);
  
  
  B <- seq(0.1, max(1, b_hat*2), length.out = resolution);
  B_0 <- seq(0.1, max(1, b_0_hat*2), length.out = resolution);
  PPI <- matrix(0, nrow = length(B), ncol = length(B_0));
  par(mfrow = c(2, 2), oma = c(2, 0, 2, 0));
  h <- 1;
  h_0 <- h;
  for (i in 1:length(B)) {
    b <- B[i];
    for (j in 1:length(B_0)) {
      b_0 <- B_0[j];
      PPI[i, j] <- 1/(1 + exp(-lklh(alpha, z, a, b, h, a_0, b_0, h_0)));
    }
  }
  data = data.frame(x = rep(B, each = length(B_0)), y = rep(B_0, length(B)), z = c(t(PPI)));
  myPanel <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...);
    panel.points(b_hat, b_0_hat, cex = 1.5, col = 2, pch = 16);
    panel.points(1, 1, cex = 1.5, col = 3, pch = 16);
  }
  print(levelplot(z~x*y, data = data, col.regions = gray(0:100/100)[100:1], aspect = "iso", xlab = expression(b), ylab = expression(b[0]), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.05), labels = list(cex = 1.5)), panel = myPanel, main = expression(h==1)), split=c(1, 1, 2, 2));
  h <- 10;
  h_0 <- h;
  for (i in 1:length(B)) {
    b <- B[i];
    for (j in 1:length(B_0)) {
      b_0 <- B_0[j];
      PPI[i, j] <- 1/(1 + exp(-lklh(alpha, z, a, b, h, a_0, b_0, h_0)));
    }
  }
  data = data.frame(x = rep(B, each = length(B_0)), y = rep(B_0, length(B)), z = c(t(PPI)));
  print(levelplot(z~x*y, data = data, col.regions = gray(0:100/100)[100:1], aspect = "iso", xlab = expression(b), ylab = expression(b[0]), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.05), labels = list(cex = 1.5)), panel = myPanel, main = expression(h==10)), split=c(2, 1, 2, 2), newpage=FALSE);
  h <- 100;
  h_0 <- h;
  for (i in 1:length(B)) {
    b <- B[i];
    for (j in 1:length(B_0)) {
      b_0 <- B_0[j];
      PPI[i, j] <- 1/(1 + exp(-lklh(alpha, z, a, b, h, a_0, b_0, h_0)));
    }
  }
  data = data.frame(x = rep(B, each = length(B_0)), y = rep(B_0, length(B)), z = c(t(PPI)));
  print(levelplot(z~x*y, data = data, col.regions = gray(0:100/100)[100:1], aspect = "iso", xlab = expression(b), ylab = expression(b[0]), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.05), labels = list(cex = 1.5)), panel = myPanel, main = expression(h==100)), split=c(1, 2, 2, 2), newpage=FALSE);
  h <- 1000;
  h_0 <- h;
  for (i in 1:length(B)) {
    b <- B[i];
    for (j in 1:length(B_0)) {
      b_0 <- B_0[j];
      PPI[i, j] <- 1/(1 + exp(-lklh(alpha, z, a, b, h, a_0, b_0, h_0)));
    }
  }
  data = data.frame(x = rep(B, each = length(B_0)), y = rep(B_0, length(B)), z = c(t(PPI)));
  print(levelplot(z~x*y, data = data, col.regions = gray(0:100/100)[100:1], aspect = "iso", xlab = expression(b), ylab = expression(b[0]), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.05), labels = list(cex = 1.5)), panel = myPanel, main = expression(h==1000)), split=c(2, 2, 2, 2), newpage=FALSE);
  ppi <- rep(NA, 4);
  names(ppi) <- c(1, 10, 100, 1000);
  count <- 1;
  for(h in c(1, 10, 100, 1000)) {
    h_0 <- h;
    ppi[count] <- 1/(1 + exp(-lklh(alpha, z, a, b_hat, h, a_0, b_0_hat, h_0)));
    ppi[count] <- 1/(1 + exp(-lklh(alpha, z, a, 1, h, a_0, 1, h_0)));
    count <- count + 1;
  }
  par(mfrow = c(1, 1));
  return(ppi);
}

lklh = function(alpha, z, a, b, h, a_0, b_0, h_0) {
  logr = 0;
  for (k in unique(z)) {
    logr = logr - 0.5*log(sum(z == k)*h + 1) + lgamma(a + sum(z == k)/2) - lgamma(a) + a*log(b) - (a + sum(z == k)/2)*log(b + 0.5*(sum((alpha[z == k])^2) - sum(alpha[z == k])^2/(sum(z == k) + 1/h)));
  }
  logr = logr - (-0.5*log(length(z)*h_0 + 1) + lgamma(a_0 + length(z)/2) - lgamma(a_0) + a_0*log(b_0) - (a_0 + length(z)/2)*log(b_0 + 0.5*(sum(alpha^2) - sum(alpha)^2/(length(z) + 1/h_0))));
  return(logr);
}


rdirichlet = function(alpha) {
  p <- length(alpha);
  temp <- rep(NA, p);
  for (j in 1:p) {
    temp[j] <- rgamma(1, alpha[j], 1);
  }
  temp <- temp/sum(temp);
  return (temp);
}

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
  
  # Return the count matrix in basic taxonomic level
  Y <- t(otu[which(index_base),]);
  
  return(list(Y = Y, S = S));
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



data_generator_dm = function(zt, gammat, sigma_b, sigma_w, seed) {
  n <- length(zt);
  p <- length(gammat);
  K <- max(zt) + 1;
  set.seed(seed);
  N <- floor(runif(n, 5000, 10000));
  # N <- floor(runif(n, 20000000, 60000000));
  Mt <- matrix(NA, nrow = K, ncol = p);
  Mt[, which(gammat == 0)] <- rep(runif(sum(gammat == 0), 0, 2), each = K);
  # Mt[, which(gammat == 0)] <- rep(runif(sum(gammat == 0), -10, 25), each = K);
  temp <- seq(1 - (K - 1)/2*sigma_b, 1 - (K - 1)/2*sigma_b + (K - 1)*sigma_b, by = sigma_b);
  for (j in which(gammat == 1)) {
    Mt[sample(1:K, K), j] <- temp;
    # Mt[sample(1:K, K), j] <- temp + runif(1, -10, 25);
  }
  At <- matrix(NA, nrow = n, ncol = p);
  Phit <- matrix(NA, nrow = n, ncol = p);
  Y <- matrix(NA, nrow = n, ncol = p);
  colnames(Y) <- 1:p;
  rownames(Y) <- paste0(1:n, "-grp", zt + 1);
  for (j in 1:p) {
    if (gammat[j] == 1) {
      colnames(Y)[j] <- paste0(colnames(Y)[j], "-TP");
    }
  }
  for (i in 1:n) {
    At[i,] <- exp(rnorm(p, Mt[zt[i] + 1,], sigma_w));
    Phit[i,] <- rdirichlet(At[i,]);
    Y[i,] <- rmultinom(1, N[i], Phit[i,]);
  }
  index <- which(colSums(Y) == 0);
  if (length(index) > 0) {
    Y <- Y[, -index];
    gammat <- gammat[-index];
    At <- At[, -index];
    Phit <- Phit[, -index];
  }
  return(list(Y = Y, gammat = gammat, zt = zt, At = At, Phit = Phit));
}



data_generator_zinb = function(zt, gammat, sigma_b, sigma_w, seed, model = c("zinb", "zip", "nb", "poisson")) {
  n <- length(zt);
  p <- length(gammat);
  K <- max(zt) + 1;
  set.seed(seed);
  s <- runif(n, 0.5, 4);
  g <- rep(1, p);
  if (model == "zinb" || model == "nb") {
    psi <- rexp(p, 1/10);
  } else if (model == "zip" || model == "poisson") {
    psi <- rep(NA, p);
  }
  Mt <- matrix(NA, nrow = K, ncol = p);
  Mt[, which(gammat == 0)] <- rep(runif(sum(gammat == 0), 0, 2), each = K);
  # Mt[, which(gammat == 0)] <- rep(runif(sum(gammat == 0), 8, 13), each = K);
  temp <- seq(1 - (K - 1)/2*sigma_b, 1 - (K - 1)/2*sigma_b + (K - 1)*sigma_b, by = sigma_b);
  for (j in which(gammat == 1)) {
    Mt[sample(1:K, K), j] <- temp;
    # Mt[sample(1:K, K), j] <- temp + runif(1, 8, 13);
  }
  A <- matrix(NA, nrow = n, ncol = p);
  for (i in 1:n) {
    A[i,] <- exp(rnorm(p, Mt[zt[i] + 1,], sigma_w));
  }
  Y <- matrix(NA, nrow = n, ncol = p);
  colnames(Y) <- 1:p;
  rownames(Y) <- paste0(1:n, "-grp", zt + 1);
  for (j in 1:p) {
    if (gammat[j] == 1) {
      colnames(Y)[j] <- paste0(colnames(Y)[j], "-TP");
    }
  }
  for (j in 1:p) {
    for (i in 1:n) {
      if (model == "zinb" || model == "nb") {
        Y[i, j] <- rnbinom(1, mu = s[i]*g[j]*A[i, j], size = psi[j]);
      } else if (model == "zip" || model == "poisson") {
        Y[i, j] <- rpois(1, s[i]*g[j]*A[i, j]);
      }
    }
  }
  H <- matrix(0L, nrow = n, ncol = p);
  if (model == "zinb" || model == "zip") {
    H[sample(1:(n*p), floor(n*p/2))] <- 1;
    Y[which(H == 1)] <- 0;
  }
  index <- which(colSums(Y) == 0);
  if (length(index) > 0) {
    Y <- Y[, -index];
    gammat <- gammat[-index];
    g <- g[-index];
    psi <- psi[-index];
    H <- H[, -index];
    A <- A[, -index];
  }
  return(list(Y = Y, gammat = gammat, zt = zt, st = s, gt = g, phit = psi, Ht = H, At = A));
}



# =============================================================================================
# The function of evaluating classification performance
# ---------------------------------------------------------------------------------------------
# Input:  t, a p-dimensional binary vector, which denotes the truth
# Input:  p, a p-dimensional vector, which gives probabilities or p-values
# Input:  cutoff, a vector, which defines the resolution to plot ROC curve
# Input:  alpha, a positive value between 0 and 1, which defines the significance level
# Input:  probability, a boolean variable, with FALSE meaning p is a vector of p-values
# ---------------------------------------------------------------------------------------------
# Output: a list of classification performance metrics
#         roc:                 the points to draw a receiver operating characteristic curve
#         auc:                 the area under the roc curve
#         specificity
#         sensitivity (recall)
#         precision
#         f1:                  F-1 score
#         mcc:                 Matthews correlation coefficient
#         acc:                 accuracy
# =============================================================================================
classification_performance_evaluator = function(t, p, cutoff = seq(0, 1.01, by = 0.01), 
                                                 alpha = 0.05, probability = TRUE) {
  if (!probability) {
    p <- p.adjust(p, method = "fdr");
  }
  result <- matrix(NA, nrow = length(cutoff), ncol = 4);
  colnames(result) <- c("tp", "fp", "fn", "tn");
  for (i in 1:length(cutoff)) {
    if (probability) {
      tab <- tabulate_error(t, p >= cutoff[i]);
    } else {
      tab <- tabulate_error(t, p <= cutoff[i]);
    }
    result[i, "tn"] <- tab[1, 1];
    result[i, "fp"] <- tab[2, 1];
    result[i, "fn"] <- tab[1, 2];
    result[i, "tp"] <- tab[2, 2];
  }
  roc <- matrix(NA, nrow = length(cutoff), ncol = 2);
  colnames(roc) <- c("fpr", "tpr");
  roc[, "fpr"] <- result[, "fp"]/(result[, "fp"] + result[, "tn"]);
  roc[, "tpr"] <-result[, "tp"]/(result[, "tp"] + result[, "fn"])
  auc <- abs(sum(diff(rev(roc[, "fpr"]))*rollmean(rev(roc[, "tpr"]), 2)));
  pauc <- abs(sum(diff(rev(roc[which(roc[, "fpr"] <= 0.1), "fpr"]))*rollmean(rev(roc[which(roc[, "fpr"] <= 0.1), "tpr"]), 2)));
  pauc <- 0.5*(1 + (pauc-0.01/2)/(0.1 - 0.01/2));
  
  if (probability) {
    c <- BayFDR(p, alpha);
    tab <- tabulate_error(t, p >= c);
  } else {
    c <- alpha;
    # tab <- tabulate_error(t, p.adjust(p, method = "fdr") <= c);
    tab <- tabulate_error(t, p <= c);
  }
  tp <- tab[2, 2];
  fp <- tab[2, 1];
  fn <- tab[1, 2];
  tn <- tab[1, 1];
  sp <- tn/(fp + tn);
  se <- tp/(tp + fn);
  pr <- tp/(tp + fp);
  f1 <- 2*(pr*se)/(pr + se);
  mcc <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
  acc <- (tp + tn)/(tp + fp + fn + tn);
  return (list(roc = roc, auc = auc, pauc = pauc, specificity = sp, sensitivity = se, precision = pr, f1 = f1, mcc = mcc, acc = acc));
}

# =============================================================================================
# The functions attached to classification_performance_evaluator(...)
# =============================================================================================
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

tabulate_error = function(gamma_true, gamma) {
  table = matrix(0L, 2, 2);
  p <- length(gamma_true);
  for (i in 1:p) {
    table[gamma[i] + 1, gamma_true[i] + 1] <- table[gamma[i] + 1, gamma_true[i] + 1] + 1;
  }
  return (table);
}

roc = function(gamma_true, mpv, cutoff, increasing) {
  result <- matrix(NA, nrow = length(cutoff), ncol = 2);
  for (i in 1:length(cutoff)) {
    if (increasing) {
      tab <- tabulate_error(gamma_true, mpv >= cutoff[i]);
    } else {
      tab <- tabulate_error(gamma_true, mpv < cutoff[i]);
    }
    result[i, 1] <- tab[2, 1]/sum(tab[, 1]);
    result[i, 2] <- tab[2, 2]/sum(tab[, 2]);
  }
  if (increasing) {
    result <- rbind(c(1, 1), result);
    result <- rbind(result, c(0, 0));
  } else {
    result <- rbind(c(0, 0), result);
    result <- rbind(result, c(1, 1));
  }
  auc <- sum(diff(rev(result[, 1]))*rollmean(rev(result[, 2]), 2));
  return (list(roc = result, auc = auc));
}

# Download from https://github.com/cran/zoo/blob/master/R/rollmean.R
# rollmean <- function(x, k, fill = if (na.pad) NA, na.pad = FALSE, 
#                      align = c("center", "left", "right"), ...) {
#   UseMethod("rollmean")
# }



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