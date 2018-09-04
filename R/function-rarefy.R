#' Efficient rarefy count matrix
#'
#' @param otu a taxa-by-sample count matrix
#'
#' @return rarefy results
#' @export
#'
#' @examples
#'   prefix = system.file("extdata", package = "BayesSLAM")
#'   load(file.path(prefix, "Castro-NallarE_2015.Rdata"))
#'   otuTable <- rarefy(otuTable)
rarefy <- function(otu) {
  depth = min(colSums(otu))
  ret <- apply(otu, 2, rarefy.vec, depth)
  rownames(ret) <- rownames(otu)
  ret
}

rarefy.vec <- function(x, depth){
  input <- (x)
  ord <- (order(-input))
  out <- rep(0, length(input))
  k = depth # 288710134 ## total depth after rarefy
  s = sum(input) ## total depth before rarefy
  for (i in 1:length(input)) { ## length(input)) { 
    count <- input[ord[i]]
    stopifnot(s == sum(input[ord[i:length(input)]]))
    r <- rhyper(1, m = count, n = s - count, k = k)
    out[ord[i]] <- r
    k <- k - r
    s <- s - count
    # cat(i, 'count=',count, 'r=',r, 'k=', k, 's=', s, "\n")
    if (count == 0 || k == 0 || s == 0) {
      break
    }
  }
  out
}

if (FALSE) {
  system.time(tmp2 <- rarefy(otu))
  rownames(tmp2) <- rownames(otu)
  head(tmp2[,1])
  
  library(GUniFrac)
  system.time(tmp3 <- GUniFrac::Rarefy(t(otu), depth = 288710134)) # takes 69s
  head(tmp3$otu.tab.rff[1,])
  library(rtk)
  help(package = "rtk")
  
  library(rtk)
  system.time(tmp4 <- rtk(otu, repeats = 1, depth = 288710134, ReturnMatrix = 1, verbose= TRUE, threads = 1))
  head(tmp4$raremat[[1]][,1])
  
  ## multithread does not improve running time...
  system.time(tmp5 <- rtk(otu, repeats = 1, depth = 288710134, ReturnMatrix = 1, verbose= TRUE, threads = 4))
  head(tmp5$raremat[[1]][,1])
  
  library(extraDistr)
  summary(otu[,1])
  summary(otu[,2])
  colSums(otu)
  length(table(otu[,1]))
  system.time(res.h <- rmvhyper(1, n = table(otu[,1]), k = 288710134)) ## this will not work, return NA
  
  
}