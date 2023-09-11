build_m_vario <- function() {
  data.frame(
    Cd = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 140),
    Co = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 40),
    Cr = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 140),
    Cu = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 40),
    Ni = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 40),
    Pb = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 40),
    Zn = c(0.1, 0.6, 0.5, 0.5, 0.5, 0.3, 1.5, 1.5, 1.5, 140)
  )
}

build_m_grid <- function(c_geo) {
  xinc <- 0.25
  yinc <- 0.25

  xr <- range(c_geo$x, na.rm=TRUE)
  xr <- c(xinc*floor(xr[1]/xinc), xinc*ceiling(xr[2]/xinc))
  yr <- range(c_geo$y, na.rm=TRUE)
  yr <- c(yinc*floor(yr[1]/yinc), yinc*ceiling(yr[2]/yinc))

  gr <- GridTopology(cellcentre.offset = c(xr[1],yr[1]),
                     cellsize = c(xinc,yinc),
                     cells.dim = c(diff(xr)/xinc, diff(yr)/yinc)+1)
  gr <- SpatialGrid(gr)
  gr <- raster(gr)

  return(gr)
}

build_m_ind <- function(c_geo, b_bnd) {
  flds_rock <- c("SEQ","KIM","QUA","ARG","POR")
  gr <- raster(ncols=3*54, nrows=3*66, xmn=0, xmx=5.4, ymn=0, ymx=6.6)

  x <- c_geo
  coordinates(x) <- ~ x + y
  z <- as.data.frame(coordinates(gr))
  names(z) <- c("x","y")
  coordinates(z) <- ~ x + y

  set.seed(9)

  bm <- lapply(flds_rock, function(i) {
    z <- krige(as.formula(paste(i, "~ 1")), x, newdata=z, nsim=0,
               model=vgm(0.3, "Sph", 1.5, add.to=vgm(0.6, "Sph", 0.5, 0.1)),
               maxdist=3, nmin=1, nmax=10, omax=3)
    z@data$var1.pred[is.na(z@data$var1.pred)] <- 0
    z@data$var1.pred <- pmax(0, pmin(1, z@data$var1.pred))
    raster(SpatialPixelsDataFrame(z@coords, z@data))
  })
  bm <- brick(bm)
  names(bm) <- flds_rock
  bm <- rasterize(b_bnd, bm, mask=TRUE)

  flds <- names(bm)
  jura_rock <- bm
  jura_rock <- which.max(jura_rock)
  jura_rock <- ratify(jura_rock)
  levels(jura_rock) <- list(data.frame(ID=1:length(flds), rock=flds))
  bm$rock <- jura_rock

  return(bm)
}

build_m_id <- function(ppmt, m_ind, c_flds) {
  # Inverse Distance
  flds <- c_flds$c_geo$independ
  df <- ppmt$Y
  z <- sapply(flds, function(i) {
    x <- df[, c("x","y",i)]
    x <- x[complete.cases(x), ]
    z <- as.data.frame(coordinates(m_ind))
    names(z) <- c("x","y")

    y <- idw(as.formula(paste(i, " ~ 1")), ~ x + y, data=x, newdata=z,
             idp=2, maxdist=6, nmin=1, nmax=10, omax=5,
             block=c(res(m_ind)[1]/3, res(m_ind)[2]/3))
    y$var1.pred
  })
  z <- as.data.frame(z)
  z <- as.data.frame(ppmt_inv(z, ppmt))
  names(z) <- flds
  z <- SpatialPixelsDataFrame(coordinates(m_ind), z)
  z <- brick(z)

  # Aggregate to parent cell.
  x <- raster::aggregate(z, fact=3)

  # Add back to sub-cell.
  bm <- raster::extract(x, coordinates(m_ind))
  bm <- as.data.frame(bm)
  bm <- SpatialPixelsDataFrame(coordinates(m_ind), bm)
  bm <- brick(bm)
  bm <- mask(bm, m_ind$rock)
  bm$rock <- m_ind$rock

  # Export model
  x <- as.data.frame(rasterToPoints(bm))
  x$rock <- (levels(bm$rock)[[1]]$rock)[x$rock]
  write.csv(x, "./data/O_outputs/jura_model_id.csv", na="", row.names=FALSE)
  saveRDS(bm, "./data/O_outputs/jura_model_id.rds")
  writeRaster(bm, "./data/O_outputs/jura_model_id.grd", format="raster", overwrite=TRUE)
  writeRaster(bm, "./data/O_outputs/jura_model_id.nc", format="CDF", overwrite=TRUE)

  return(bm)
}

build_m_sk <- function(ppmt, m_ind, m_vario, c_flds) {
  browser()
  flds <- c_flds$c_geo$independ
  x <- rbind(m_vario, round(sapply(m_vario, function(x) x[9] / x[8]), 2))
  x <- t(x)
  vario <- lapply(flds, function(i) {
    vgm(x[i,6], "Sph", x[i,8], anis=c(x[i,10], x[i,11]),
        add.to=vgm(x[i,2], "Sph", x[i,4], x[i,1], anis=c(x[i,10], x[i,11])))
  })
  names(vario) <- flds
  rm(x)

  # Simple Kriging
  df <- ppmt$Y
  z <- sapply(flds, function(i) {
    x <- df[, c("x","y",i)]
    x <- x[complete.cases(x), ]
    coordinates(x) <- ~ x + y
    z <- as.data.frame(coordinates(m_ind))
    names(z) <- c("x","y")
    coordinates(z) <- ~ x + y

    y <- krige(as.formula(paste(i, " ~ 1")), x, z, beta=0,
               model=vario[[i]], maxdist=6, nmin=1, nmax=10, omax=5,
               block=c(res(m_ind)[1]/3, res(m_ind)[2]/3))
    y$var1.pred
  })
  z <- as.data.frame(z)
  z <- as.data.frame(ppmt_inv(z, ppmt))
  names(z) <- flds
  z <- SpatialPixelsDataFrame(coordinates(m_ind), z)
  z <- brick(z)

  # Aggregate to parent cell.
  x <- raster::aggregate(z, fact=3)

  # Add back to sub-cell.
  bm <- raster::extract(x, coordinates(m_ind))
  bm <- as.data.frame(bm)
  bm <- SpatialPixelsDataFrame(coordinates(m_ind), bm)
  bm <- brick(bm)
  bm <- mask(bm, m_ind$rock)
  bm$rock <- m_ind$rock

  # Export model
  x <- as.data.frame(rasterToPoints(bm))
  x$rock <- (levels(bm$rock)[[1]]$rock)[x$rock]
  write.csv(x, "./data/O_outputs/jura_model_sk.csv", na="", row.names=FALSE)
  saveRDS(bm, "./data/O_outputs/jura_model_sk.rds")
  writeRaster(bm, "./data/O_outputs/jura_model_sk.grd", format="raster", overwrite=TRUE)
  writeRaster(bm, "./data/O_outputs/jura_model_sk.nc", format="CDF", overwrite=TRUE)

  return(bm)
}

build_m_ok <- function(ppmt, m_ind, m_vario, c_flds) {
  browser()
  flds <- c_flds$c_geo$independ
  x <- rbind(m_vario, round(sapply(m_vario, function(x) x[9] / x[8]), 2))
  x <- t(x)
  vario <- lapply(flds, function(i) {
    vgm(x[i,6], "Sph", x[i,8], anis=c(x[i,10], x[i,11]),
        add.to=vgm(x[i,2], "Sph", x[i,4], x[i,1], anis=c(x[i,10], x[i,11])))
  })
  names(vario) <- flds
  rm(x)

  # Simple Kriging
  df <- ppmt$Y
  z <- sapply(flds, function(i) {
    x <- df[, c("x","y",i)]
    x <- x[complete.cases(x), ]
    coordinates(x) <- ~ x + y
    z <- as.data.frame(coordinates(m_ind))
    names(z) <- c("x","y")
    coordinates(z) <- ~ x + y

    y <- krige(as.formula(paste(i, " ~ 1")), x, z,
               model=vario[[i]], maxdist=6, nmin=1, nmax=10, omax=5,
               block=c(res(m_ind)[1]/3, res(m_ind)[2]/3))
    y$var1.pred
  })
  z <- as.data.frame(z)
  z <- as.data.frame(ppmt_inv(z, ppmt))
  names(z) <- flds
  z <- SpatialPixelsDataFrame(coordinates(m_ind), z)
  z <- brick(z)

  # Aggregate to parent cell.
  x <- raster::aggregate(z, fact=3)

  # Add back to sub-cell.
  bm <- raster::extract(x, coordinates(m_ind))
  bm <- as.data.frame(bm)
  bm <- SpatialPixelsDataFrame(coordinates(m_ind), bm)
  bm <- brick(bm)
  bm <- mask(bm, m_ind$rock)
  bm$rock <- m_ind$rock

  # Export model
  x <- as.data.frame(rasterToPoints(bm))
  x$rock <- (levels(bm$rock)[[1]]$rock)[x$rock]
  write.csv(x, "./data/O_outputs/jura_model_ok.csv", na="", row.names=FALSE)
  saveRDS(bm, "./data/O_outputs/jura_model_ok.rds")
  writeRaster(bm, "./data/O_outputs/jura_model_ok.grd", format="raster", overwrite=TRUE)
  writeRaster(bm, "./data/O_outputs/jura_model_ok.nc", format="CDF", overwrite=TRUE)

  return(bm)
}
