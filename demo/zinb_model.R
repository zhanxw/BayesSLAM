# ============================================================================================
# Set working directory
# setwd("~/Dropbox/xiaowei");
prefix = system.file("extdata", package = "BayesSLAM")

# Load library
library(ggplot2);

# Load functions
if (FALSE) {
  Rcpp::sourceCpp('code/core_zinb_x4.cpp');   # The MCMC algorithm of fitting a zero-inflated negative binomial (ZINB) model written in C++
  Rcpp::sourceCpp('code/core_dm_x4.cpp');   # The MCMC algorithm of fitting a Dirichlet-multinomial (DM) model written in C++
  source('code/functions.R');   # Other related functions written in R
}

# ============================================================================================



# ============================================================================================
# Load data
if (FALSE) {
  load("data/real_data/curatedMetagenomicData.Rdata");
  # Let's analyze TettAJ_2016 dataset, which contains 97 samples in 2 groups (control vs. psoriasis)
  otu <- otu_list[["TettAJ_2016"]];   # Abundance table at all taxonomic levels (682 speices, 93 samples in this case)
  tax <- tax_list[["TettAJ_2016"]];   # Taxonomy table
  z <- as.integer(z_list[["TettAJ_2016"]]) - 1;   # Group indicator (49 control and 44 experimental samples in this case)
  save(list = c("otu", "tax", "z"), file = "TettAJ_2016.Rdata")
} else {
  print(prefix)
  print
  load(file.path(prefix, "data/TettAJ_2016.Rdata"))
}

# Preprocess data
data <- real_data_preprocessor(otu, tax);   # Generate data passed to our algorithm
Y <- data$Y;   # Abundance table at bottom taxonomic level (379 species and 9 genus in this case) 
S <- data$S;   # Phylogenetic tree structure that can be read by our core functions
dim(otu)
dim(Y)
dim(S) ## why S has 294 rows ( 294 != 682 - 379) ??
S_to_Sx <- function(S) {
  if (any(S > 2)) {
    warning("Is your S a boolean matrix?")
  }
  S_x <- matrix(0L, nrow = dim(S)[1], ncol = max(rowSums(S)));
  for (ii in 1:dim(S)[1]) {
    S_x[ii, 1:sum(S[ii,] == 1)] <- which(S[ii,] == 1)
  }
  rownames(S_x) <- rownames(S);
  S_x
}
S <- S_to_Sx(S)
# ============================================================================================



# ============================================================================================
# Fit ZINB model
iter <- 1000;
aggregate <- TRUE;
n <- dim(Y)[1];   # Number of samples
p <- dim(Y)[2];   # Number of features
K <- length(unique(z));   # Number of groups

# Estimate plug-in size factors (optional)
s_tss <- rescale(size_factor_estimator(Y, method = "TSS"), constraint = "product"); 
# s_q75 <- rescale(size_factor_estimator(Y, method = "Q75"), constraint = "product");  
# s_tmm <- rescale(size_factor_estimator(Y, method = "TMM"), constraint = "product");  
# s_css <- rescale(size_factor_estimator(Y, method = "CSS"), constraint = "product");  
s_1 <- rep(1, n);   # Starting points of size factors if choosing "DPP"

# # Estimate hyperparameter b and b_0
# # b_0 <- rep(1, p + dim(S)[1]);
# # B <- matrix(1, nrow = K, ncol = p + dim(S)[1]);
# temp <- b_estimator(Y, z, S, aggregate, s_tss);
# b_0 <- temp$b_0;
# B <- temp$B;

# Run MCMC
start_time <- proc.time();
R <- zinb_model_estimator(Y, z, s_1, iter, TRUE, S, aggregate, TRUE, 2, 1000);
end_time <- proc.time();
time <- end_time - start_time;
# save(R, time, file = paste0("result/real_data_local/id=", id, "_chain=", chain, ".Rdata"));
# ============================================================================================



# ============================================================================================
# Fit DM model
iter <- 1000;
aggregate <- TRUE;
n <- dim(Y)[1];   # Number of samples
p <- dim(Y)[2];   # Number of features
K <- length(unique(z));   # Number of groups

# # Estimate hyperparameter b and b_0
# b_0 <- rep(1, p + dim(S)[1]);
# B <- matrix(1, nrow = K, ncol = p + dim(S)[1]);

# Run MCMC
start_time <- proc.time();
R <- dm_model_estimator(Y, z, iter, S, aggregate, TRUE);
end_time <- proc.time();
time <- end_time - start_time;
# ============================================================================================



# ============================================================================================
# Visualize results
PPIs <- R$gamma_ppi;
names(PPIs) <- c(colnames(Y), row.names(S));
qplot(1:length(PPIs), PPIs) + geom_bar(stat = "identity") + ylim(0, 1) + ylab("Posterior probabilities of inclusion") + xlab("Taxa index") + geom_hline(yintercept = BayFDR(PPIs, 0.05), colour = "red", linetype="dashed")

index <- which(PPIs >= BayFDR(PPIs, 0.05));
if (K == 2) {
  data_temp <- data.frame(cbind(apply(R$logafc[, index], 2, quantile, 0.025), apply(R$logafc[, index], 2, quantile, 0.5), apply(R$logafc[, index], 2, quantile, 0.975)));
  names(data_temp) <- c("lower", "median", "upper");
  data_temp$taxa <- names(index);
  data_temp$taxa = factor(data_temp$taxa, levels = data_temp$taxa[order(data_temp$median)], ordered = TRUE);
  ggplot(data = data_temp, aes(x = taxa)) + ylab(expression(log(alpha[j1]/alpha[j0]))) + geom_point(aes(y = median), colour="red") + geom_errorbar(mapping = aes(ymin = lower, ymax = upper, width = 0.25)) + coord_flip() + geom_hline(yintercept = 0, linetype="dashed");
}
# ============================================================================================