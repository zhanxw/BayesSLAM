## strip "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter" to ""g__Methanobrevibacter"
lastTaxonomy <- function(x) {
  sapply(x, function(v) {
    res <- strsplit(v, '\\|')[[1]]
    res[length(res)]
  })
}

## convert full_taxa_names (e.g.k__Archaea|p__Euryarchaeota|c__Methanobacteria)
## to a data frame
createTaxaTable <- function(full_taxa) {
  last_taxa <- lastTaxonomy(full_taxa)
  tmp <- lapply(last_taxa, function(x) {
    idx = grep(pattern = sprintf("%s$", x), x = full_taxa)
    ret = strsplit(full_taxa[idx],"\\|")[[1]]
    # print(ret)
    ret <- c(ret, rep(NA, 8 - length(ret)))
    (ret)
  })
  # print(str(tmp))
  tmp <- do.call(rbind, tmp)
  colnames(tmp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  abbrev <- c("k", "p", "c", "o", "f", "g", "s", "t")
  for (i in 1:8) {
    pattern = sprintf("^%s__", abbrev[i])
    tmp[,i] <- sub(pattern = pattern, replacement = "", tmp[,i])
  }
  rownames(tmp) <- last_taxa
  tmp
}