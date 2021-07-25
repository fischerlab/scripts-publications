### helper functions from Michael Stadler
map2col <- function(x, xlim=NULL, do.log2p1=FALSE, cols=rev(brewer.pal(11,'Spectral'))) {
  if(do.log2p1)
    x <- log2(x+1)
  if(!is.null(xlim)) {
    x[x < xlim[1]] <- xlim[1]
    x[x > xlim[2]] <- xlim[2]
  }
  tmp <- colorRamp(cols)((x-min(x,na.rm=TRUE))/diff(range(x,na.rm=TRUE)))
  rgb(tmp[,1], tmp[,2], tmp[,3], max=255)
}

####