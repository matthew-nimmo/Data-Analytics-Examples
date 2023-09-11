v3_factory <- function() {
  list(
    tar_target(v3_umap, build_v3_umap(c_geo, c_flds)),
    tar_target(v3_grps, build_v3_grps(v3_umap))
  )
}

build_v3_umap <- function(c_geo, c_flds) {
  x <- c_geo[, c_flds$c_geo$independ]
  x <- x[complete.cases(x), ]
  x <- x[!duplicated(x), ]
  x <- scale(x)

  set.seed(21)
  um <- umap::umap(x, n_components=2, n_neighbors=15, min_dist=1e-10)
  saveRDS(um, "./models/v3/v3_umap.Rds")

  return(um)
}

build_v3_gmm <- function(c_all, c_flds) {
  fit_gmm <- function(flds) {
    df <- c_all[, flds]
    df <- df[complete.cases(df), ]
    bic <- mclustBIC(df, G=1:20)
    m <- Mclust(df, x=bic)
    return(m)
  }

  set.seed(17)
  m <- list(
    gmm_ass = fit_gmm(c_flds$c_geo$independ),
    gmm_com = fit_gmm(c_flds$c_com$depend),
    gmm_met = fit_gmm(c_flds$c_met$depend)
  )

  saveRDS(m, "./models/v3/v3_gmm.Rds")

  return(m)
}

build_v3_grps <- function(v3_umap) {
  k_clusters <- 7
  set.seed(123)

  x <- v2_umap$layout
  ms_h <- round(quantile(dist(x), 0.01), 2)
  m <- list(
    km = kmeans(x, centers=k_clusters, iter.max=100, nstart=3),
    pam = pam(x, k=k_clusters, metric="euclidean", stand=FALSE),
    clara = clara(x, k=k_clusters, metric="euclidean",
                  stand=FALSE, samples=50, pamLike=FALSE),
    gmm = GMM(x, k_clusters, "eucl_dist", "random_subset", 10, 10),
    ms = meanShift(x, x, kernelType="EPANECHNIKOV",
                   iterations=30, epsilonCluster=ms_h),
    ap0 = apcluster(apcluster::negDistMat(r=2), x),
    apK = apclusterK(apcluster::negDistMat(r=2), x, K=k_clusters)
  )

  saveRDS(m, "./models/v3/v3_groups.Rds")

  return(m)
}

build_v3_class <- function(c_geo, c_flds, v2_grps) {
  x <- c_geo[, c_flds$c_geo$independ]
  x <- x[complete.cases(x), ]
  x$group <- factor(v2_grps$km$cluster)

  m <- rpart(group ~ ., data=x, model=TRUE)
  saveRDS(m, "./models/v3/v3_groups-dt.Rds")

  return(m)
}

build_v3_bn <- function(c_all, c_flds, v2_gmm) {
  return(NULL)
}

build_v3_reg <- function(d_geo, d_com, d_met, d_flds) {
  set.seed(307)
  ctrl <- trainControl(method="repeatedcv", number=10, repeats=10)
  cb_grid <- expand.grid(committees=1:20, neighbors=1:9)

  #library(doParallel)
  #cl <- makeCluster(detectCores() - 1)
  #registerDoParallel(cl)

  m <- lapply(1:length(c_flds), function(i) {
    u <- get(names(c_flds)[i])
    x <- u[, c_flds[[i]]$independ, drop=FALSE]
    y <- u[, c_flds[[i]]$depend, drop=FALSE]

    apply(y, 2, function(y) {
      j <- complete.cases(cbind(x,y))
      train(x[j,], y[j], method="cubist", trControl=ctrl)
    })
  })
  names(m) <- names(c_flds)

  #stopCluster(cl)
  #registerDoSEQ()

  saveRDS(m, "./models/v3/v3_regression.Rds")

  return(m)
}
