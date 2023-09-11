declus.for <- function(grade, x, y, z=NULL, ncells=1, cell=c(1,10), anis=c(1,1), offsets=5, minmax=0)
{
	if (is.null(z))
		z <- rep(1,length(x))

	i <- !is.na(grade)
	x <- x[i]
	y <- y[i]
	z <- z[i]
	grade <- grade[i]
	n <- length(x)
	wt <- vector(mode="numeric", length=n)

	# define a "lower" origin to use for the cell sizes:
	xo1 <- min(x) - 0.01
	yo1 <- min(y) - 0.01
	zo1 <- min(z) - 0.01
	xcs = cell[1]
	ycs = (cell[1] * anis[1])
	zcs = (cell[1] * anis[2])
	ncellx <- as.integer((max(x)-(xo1-xcs))/xcs) + 1
	ncelly <- as.integer((max(y)-(yo1-ycs))/ycs) + 1
	ncellz <- as.integer((max(z)-(zo1-zcs))/zcs) + 1
	ncellt <- (ncellx * ncelly * ncellz)

	ret <- .Fortran("declus", nd=as.integer(n), nc=as.integer(ncellt), vr=as.double(grade), x=as.double(x), y=as.double(y), z=as.double(z),
           wtopt=as.double(wt), yanis=as.double(anis[1]), zanis=as.double(anis[2]), as.integer(ncells),
		   as.double(cell[1]), as.double(cell[2]), as.integer(offsets), as.integer(minmax))

	return(ret$wtopt)
}

declus <- function(vr, x, y, z=NULL, yanis=1, zanis=1, ncell=24, cmin=1, cmax=25, noff=5, minmax=0)
{
	nd <- length(vr)
	if(is.null(z)) {
		z <- rep(0, nd)
	}

	wt <- vector(mode="numeric", length=nd)
	wtopt <- vector(mode="numeric", length=nd)
	ijk <- vector(mode="numeric", length=nd)

	imin <- TRUE
	if (minmax == 1) {
	  imin <- FALSE
	}

	roff <- noff

	# compute min, max, and average:
	vrav <- mean(vr, na.rm=TRUE)
	vrmin <- min(vr, na.rm=TRUE)
	vrmax <- max(vr, na.rm=TRUE)
	xmin <- min(x, na.rm=TRUE)
	xmax <- max(x, na.rm=TRUE)
	ymin <- min(y, na.rm=TRUE)
	ymax <- max(y, na.rm=TRUE)
	zmin <- min(z, na.rm=TRUE)
	zmax <- max(x, na.rm=TRUE)

	# initialize the "best" weight values:
	vrop <- vrav
	best <- 0.0

	# define a "lower" origin to use for the cell sizes:
	xo1 <- xmin - 0.01
	yo1 <- ymin - 0.01
	zo1 <- zmin - 0.01

	# define the increment for the cell size:
	xinc <- (cmax-cmin) / ncell
	yinc <- yanis * xinc
	zinc <- zanis * xinc

	# loop over "ncell+1" cell sizes in the grid network:
	xcs <- cmin - xinc
	ycs <- (cmin*yanis) - yinc
	zcs <- (cmin*zanis) - zinc

	# MAIN LOOP over cell sizes:
	for (lp in 1:ncell)
	{
		xcs <- xcs + xinc
		ycs <- ycs + yinc
		zcs <- zcs + zinc

		# initialize the weights to zero:
		wt <- rep(0.0, nd)

		# determine the maximum number of grid cells in the network:
		ncellx <- as.integer((xmax - (xo1 - xcs)) / xcs) + 1
		ncelly <- as.integer((ymax - (yo1 - ycs)) / ycs) + 1
		ncellz <- as.integer((zmax - (zo1 - zcs)) / zcs) + 1
		ncellt <- (ncellx * ncelly * ncellz)
		cellwt <- vector(mode="numeric", length=ncellt)

		# loop over all the origin offsets selected:
		xfac <- min((xcs/roff), (0.5 * (xmax - xmin)))
		yfac <- min((ycs/roff), (0.5 * (ymax - ymin)))
		zfac <- min((zcs/roff), (0.5 * (zmax - zmin)))

		for (kp in 1:noff)
		{
			xo <- xo1 - (kp - 1.0) * xfac
			yo <- yo1 - (kp - 1.0) * yfac
			zo <- zo1 - (kp - 1.0) * zfac

			# initialize the cumulative weight indicators:
      cellwt <- rep(0.0, ncellt)

			# determine which cell each datum is in:
			icellx <- as.integer((x - xo) / xcs) + 1
			icelly <- as.integer((y - yo) / ycs) + 1
			icellz <- as.integer((z - zo) / zcs) + 1
			icell <- icellx + (icelly - 1) * ncellx  + (icellz -1 ) * ncelly * ncellx
			ijk <- icell

      k <- table(icell)
			cellwt[as.numeric(names(k))] <- as.numeric(k)

			# The weight assigned to each datum is inversely proportional to the
			# number of data in the cell.  We first need to get the sum of weights
			# so that we can normalize the weights to sum to one:
			sumw <- 1 / sum(ifelse(cellwt>0, 1/cellwt, 0))

			# Accumulate the array of weights (that now sum to one):
			for (i in 1:nd)
			{
				ipoint <- ijk[i]
				wt[i] <- wt[i] + (1.0 / cellwt[ipoint]) * sumw
			}
		}

		# compute the weighted average for this cell size:
		sumw  <- sum(wt, na.rm=TRUE)
		sumwg <- sum(wt*vr, na.rm=TRUE)
		vrcr <- sumwg / sumw

		# see if this weighting is optimal:
		if ((imin & vrcr < vrop) | (!imin & vrcr > vrop) | (lp == 1))
		{
			best <- xcs
			vrop <- vrcr
			wtopt <- wt
		}
	}

	# Get the optimal weights:
	sumw <- sum(wtopt, na.rm=TRUE)

	if (sumw != 0) {
	  facto <- nd / sumw
	  wtopt <- wtopt * facto
	}

	return(wtopt)
}
