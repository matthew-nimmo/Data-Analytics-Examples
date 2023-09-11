#' cumulative probability plot.
#'
#' Generates a log-probability plot of numeric variable.
#'
#' For simple cumprob plots, \code{\link{cumprob.default}} will be used. However, there are
#' \code{cumprob} methods for formula and list R objects.
#' @param x numeric vector.
#' @param ... other parameters.
#'
#' @return None
#'
#' @keywords cumprob
#' @export
cumprob <- function(x, ...) UseMethod('cumprob')

#' cumulative probability plot.
#'
#' Generates a log-probability plot of numeric variable.
#'
#' @param x numeric vector.
#' @param horizontal a logical value indicating if the plot is rotated.
#' @param add a logical indicating if the plot should be added to the current cumprob plot.
#' @param main a main title for the plot.
#' @param sub a sub title for the plot.
#' @param xlab a label for the x axis, defaults to description of x.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param col the colour of the series.
#' @param pch plotting symbol.
#' @param cex symbol expansion.
#' @param type of the plot. see \code{\link{plot.default}}.
#' @param ... other parameters.
#'
#' @return None
#'
#' @keywords cumprob.default
#' @export
cumprob.default <- function(x, add=FALSE, main=NULL, sub=NULL, xlab=NULL,
                            xlim=NULL, col="black", pch=20, cex=1, type="p", ...)
{
	if (!is.numeric(x))
		stop(paste(sQuote("x"))," Needs to be numeric")

	if (is.null(xlab))
		xlab <- ToLabel(deparse(substitute(x)))

	i <- !is.na(x)
	x <- sort(x[i])
	n <- length(x)

	if (n <= 1)
		return(NA)

	if ((max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) == 0)
		return(NA)

	if (is.null(xlim))
		xlim <- range(x, na.rm=TRUE)

	y <- qnorm((1:n) / (n+1))
	ylab <- "Cumulative Probability"
	ylim <- range(y[is.finite(y)])
	if (any(!is.finite(ylim)))
		ylim <- qnorm(0.001, 0.999)

	if (is.null(main))
		main <- "Log probability plot"

	o <- par("mar")
	if (main == "") {
	  mar <- o
	  mar[3] <- 1.1
	  par(mar=mar)
	}

  if (!add)
	{
		plot(x, y, type="n", main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		     log="x", frame.plot=FALSE, axes=FALSE)

    axis(1)

	  labels <- c(0.01, 0.1, 0.5, 0.9, 0.99)
	  aty <- qnorm(labels)
    abline(h=aty, col="grey")
    mtext(labels, side=2, line=1, at=aty, cex=par("cex.axis")*par("cex"), las=1)
  }

	if (type=="l" | type=="b") {
		lines(x, y, col=col, ...)
	}

	if (type=="p" | type=="b") {
		points(x, y, pch=pch, cex=cex, col=col, ...)
	}

	par(mar=o)

	return(y)
}

#' cumulative probability plot.
#'
#' Generates a log-probability plot of numeric variable.
#'
#' @param formula a formula, such as y ~ grp, where y is a numeric vector of data values to be
#' split into groups according to the grouping variable \code{grp} (usually a factor).
#' @param data a data.frame from which the variables in \code{formula} should be taken.
#' @param xlab a label for the x axis, defaults to description of x.
#' @param col the colour of the series.
#' @param plot.legend plot the legend. Defaults to TRUE.
#' @param ... other parameters.
#'
#' @return None
#'
#' @keywords cumprob.formula
#' @export
cumprob.formula <- function(formula, data=NULL, xlab=NULL, col=NULL, plot.legend=FALSE, title=NULL, ...)
{
	mf <- model.frame(formula=formula, data=data)
  xlabels <- attr(terms(mf), "term.labels")
	ylabels <- names(mf)[attr(terms(mf), "response")]

	for (i in xlabels)
	{
		if (length(ylabels) == 0)
		{
			if (is.null(col))
				col <- "black"
			cumprob(mf[,i], xlab=i, col=col, ...)
		}
		else
		{
			for (j in ylabels)
			{
				x <- mf[,i]
				y <- mf[,j]
				if (!is.factor(x))
					x <- factor(x)

				s <- split(y, x, drop=TRUE)
				xlim <- range(y[y>0], na.rm=TRUE)

				labs <- levels(x)
				olabs <- vapply(s, length, FUN.VALUE=0)
				ord <- order(olabs, decreasing=TRUE)
				olabs <- names(olabs[ord])

				n <- length(labs)
				if (is.null(col)) {
					col <- 1:n
				}
				if (length(col) == 1) {
				  col <- rep(col, n)
				}
				ocol <- col[ord]

				indx <- 1
				first <- TRUE
				v <- data.frame(x=rep(NA,length(olabs)),y=rep(NA,length(olabs)), col=ocol)
				row.names(v) <- olabs
				for (k in olabs)
				{
					z <- cumprob(s[[k]], add=!first, xlim=xlim, xlab=ifelse(is.null(xlab), j, xlab), col=ocol[indx], ...)
			    v[k,1] <- max(s[[k]])
			    v[k,2] <- max(z)
					first <- FALSE
					indx <- indx + 1
				}

				basicPlotteR::addTextLabels(v$x, v$y, labels=row.names(v), cex.label=0.9, col.label=v$col, keepLabelsInside=FALSE)

				if (plot.legend) {
				  title <- ifelse(is.null(title), xlabels, title)
				  legend("topleft", labs, title=title, col=col, pch=20, bty="o", bg="white", cex=0.7*par("cex.axis"))
				}
				#par(o)
			}
		}
	}

	return(invisible())
}

#' cumulative probability plot.
#'
#' Generates a log-probability plot of numeric variable.
#'
#' @param list a list of variable names to plot.
#' @param data a data.frame from which the variables in \code{list} should be taken.
#' @param col the colour (can be a vector of colours) of the series.
#' @param ... other parameters.
#'
#' @return None
#'
#' @keywords cumprob.list
#' @export
cumprob.list <- function(list, data, col=NULL, ...)
{
	list <- list[list %in% names(data)]
	for (i in list)
	{
		if (is.null(col))
			col <- "black"
		cumprob(data[,i], xlab=i, col=col, ...)
	}

	return(invisible())
}

ToLabel <- function(x)
{
  if (missing(x))
    return("")

  lab <- deparse(substitute(x))

  start <- regexpr("\\$", lab)[1] + 1
  if(start==0)
    start=1

  end <- regexpr("\\[", lab)[1] - 1
  if(end<0)
    end <- nchar(lab)

  s <- gsub(".\\$","",deparse(substitute(x)))
  s
}
