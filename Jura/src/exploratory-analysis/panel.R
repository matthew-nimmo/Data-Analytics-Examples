panel.cor <- function(x, y)
{
  usr <- par("usr")
  par(usr=c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits=3)[1]
  text(0.5, 0.5, txt, cex=1.5, col=ifelse(r>0.3, "blue", ifelse(r<(-0.3), "red", "black")))
}

panel.hist <- function(x)
{
  usr <- par("usr")
  par(usr=c(usr[1:2], 0, 1.5))
  h <- hist(x, plot=FALSE, breaks="FD")
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="darkgrey")
}

panel.scatter <- function (x, y) 
{
  dens <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(dens)[1,] + 1L
  col <- viridis::viridis(256)[dens]
  points(x, y, pch=16, col=col, bg=NA, cex=0.8)
}
