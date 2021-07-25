# Selbach, Nature 2011, suppl. Eq. S4
calcKd <- function(log2ri, ti, tcc) {
  # log2ri : log2 H/L
  # ti : timepoint (hours)
  # tcc : duration of cell cycle (hours)
  ok <- is.finite(log2ri)
  ri <- 2^log2ri
  (sum(log((ri[ok]+1)*ti[ok]))/sum(ti[ok]^2)) - (log(2)/tcc)
}
# Selbach, Nature 2011, suppl. Eq. S5
calcThalf <- function(log2ri, ti, tcc) {
  kdp <- calcKd(log2ri,ti,tcc)
  log(2)/kdp
}