v2_factory <- function() {
  list(
    tar_target(v2_data, prep_v2_data(b_geo, b_com, b_met, b_flds)),
    tar_target(v2_umap, build_v2_umap(v2_data, b_flds)),
    tar_target(v2_gmm, build_v2_gmm(v2_data, b_flds)),
    tar_target(v2_scov, build_v2_scov(b_bm, v2_data, v2_umap, b_flds)),
    tar_target(v2_ucov, build_v2_ucov(b_bm, v2_data, v2_umap, b_flds)),
    tar_target(v2_grps, build_v2_grps(v2_umap, b_flds)),
    tar_target(v2_class, build_v2_class(v2_data, b_flds, v2_grps)),
    tar_target(v2_bn, build_v2_bn(v2_data, b_flds, v2_gmm)),
    tar_target(v2_reg, build_v2_reg(v2_data, b_flds))
  )
}

prep_v2_data <- function(b_geo, b_com, b_met, b_flds) {
  flds <- unlist(b_flds$b_geo)
  df1 <- b_geo[, flds]
  df1 <- df1[complete.cases(df1), ]
  df1$source <- rep("b_geo", nrow(df1))
  df1$x <- b_geo$x
  df1$y <- b_geo$y

  flds <- unlist(b_flds$b_com)
  df2 <- b_com[, flds]
  df2 <- df2[complete.cases(df2), ]
  df2$source <- rep("b_com", nrow(df2))
  df2$x <- b_com$x
  df2$y <- b_com$y

  flds <- unlist(b_flds$b_met)
  df3 <- b_met[, flds]
  df3 <- df3[complete.cases(df3), ]
  df3$source <- rep("b_met", nrow(df3))
  df3$x <- b_met$x
  df3$y <- b_met$y

  df <- dplyr::bind_rows(df1, df2, df3)
  flds <- Reduce(union, sapply(b_flds, "[[", "factors"))
  for (i in flds) {
    df[[i]] <- factor(df[[i]])
  }

  #write_csv(df, "./models/v2/jura_v2_data.csv", na="")
  #saveRDS(df, "./models/v2/jura_v2_data.Rds")

  return(df)
}

build_v2_umap <- function(v2_data, b_flds) {
  x <- v2_data[, b_flds$b_geo$independ]
  x <- scale(x)

  set.seed(21)
  um <- umap::umap(x, n_components=2, n_neighbors=15, min_dist=1e-10)
  saveRDS(um, "./models/v2/v2_umap.Rds")

  return(um)
}

build_v2_scov <- function(b_bm, v2_data, v2_umap, b_flds) {
  set.seed(1)

  flds <- colnames(v2_umap$data)
  z <- b_bm[[flds]]
  z <- as.data.frame(z)
  z <- cbind(coordinates(b_bm), z)
  indx <- complete.cases(z) & !duplicated(z)
  z <- z[indx, ]

  y <- z[, c("x", "y")]

  x_com <- v2_data[v2_data$source=="b_com", c("x", "y", b_flds$b_com$depend[1])]
  m_com <- gausspr(x_com[, 1:2], x_com[,3], kernel="rbfdot", scaled=FALSE, variance.model=TRUE)
  nn_com <- get.knnx(x_com[, 1:2], y, k=1, algorithm="kd_tree")
  m_com.cor <- sapply(1:nrow(y), function(i) {
    m_com@kernelf(as.numeric(y[i,]), as.numeric(x_com[nn_com$nn.index[i],1:2]))
  })
  o_com <- order(m_com.cor, decreasing=TRUE)
  z_com <- predict(m_com, y,type="sdeviation")

  x_met <- v2_data[v2_data$source=="b_met", c("x", "y", b_flds$b_met$depend[1])]
  m_met <- gausspr(x_met[, 1:2], x_met[,3], kernel="rbfdot", scaled=FALSE, variance.model=TRUE)
  nn_met <- get.knnx(x_met[, 1:2], y, k=1, algorithm="kd_tree")
  m_met.cor <- sapply(1:nrow(y), function(i) {
    m_com@kernelf(as.numeric(y[i,]), as.numeric(x_met[nn_met$nn.index[i],1:2]))
  })
  o_met <- order(m_met.cor, decreasing=TRUE)
  z_met <- predict(m_met, y,type="sdeviation")

  return(list(y=y, indx=indx,
              b_com=list(nn=nn_com, x=x_com, o=o_com, z=z_com, m.cor=m_com.cor, m=m_com),
              b_met=list(nn=nn_met, x=x_met, o=o_met, z=z_met, m.cor=m_met.cor, m=m_met)))
}

build_v2_ucov <- function(b_bm, v2_data, v2_umap, b_flds) {
  set.seed(1)

  flds <- colnames(v2_umap$data)
  z <- b_bm[[flds]]
  z <- as.data.frame(z)
  z <- cbind(coordinates(b_bm), z)
  indx <- complete.cases(z) & !duplicated(z)
  z <- z[indx, ]

  y <- z[, flds]
  y <- scale(y, center=attr(v2_umap$data, "scaled:center"), scale=attr(v2_umap$data, "scaled:scale"))
  y <- as.data.frame(predict(v2_umap, y))

  x_com <- v2_data[v2_data$source=="b_com", flds]
  x_com <- scale(x_com, center=attr(v2_umap$data, "scaled:center"), scale=attr(v2_umap$data, "scaled:scale"))
  x_com <- as.data.frame(predict(v2_umap, x_com))
  x_com <- cbind(x_com, v2_data[v2_data$source=="b_com", b_flds$b_com$depend[1]])
  names(x_com) <- c("x", "y", b_flds$b_com$depend[1])
  m_com <- gausspr(x_com[, 1:2], x_com[,3], kernel="rbfdot", scaled=FALSE, variance.model=TRUE)
  nn_com <- get.knnx(x_com[, 1:2], y, k=1, algorithm="kd_tree")
  m_com.cor <- sapply(1:nrow(y), function(i) {
    m_com@kernelf(as.numeric(y[i,]), as.numeric(x_com[nn_com$nn.index[i],1:2]))
  })
  o_com <- order(m_com.cor, decreasing=TRUE)
  z_com <- predict(m_com, y, type="sdeviation")

  x_met <- v2_data[v2_data$source=="b_met", flds]
  x_met <- scale(x_met, center=attr(v2_umap$data, "scaled:center"), scale=attr(v2_umap$data, "scaled:scale"))
  x_met <- as.data.frame(predict(v2_umap, x_met))
  x_met <- cbind(x_met, v2_data[v2_data$source=="b_met", b_flds$b_met$depend[1]])
  names(x_met) <- c("x", "y", b_flds$b_met$depend[1])
  m_met <- gausspr(x_met[, 1:2], x_met[,3], kernel="rbfdot", scaled=FALSE, variance.model=TRUE)
  nn_met <- get.knnx(x_met[, 1:2], y, k=1, algorithm="kd_tree")
  m_met.cor <- sapply(1:nrow(y), function(i) {
    m_met@kernelf(as.numeric(y[i,]), as.numeric(x_met[nn_met$nn.index[i],1:2]))
  })
  o_met <- order(m_met.cor, decreasing=TRUE)
  z_met <- predict(m_met, y, type="sdeviation")

  return(list(y=y, indx=indx,
              b_com=list(nn=nn_com, x=x_com, o=o_com, z=z_com, m.cor=m_com.cor, m=m_com),
              b_met=list(nn=nn_met, x=x_met, o=o_met, z=z_met, m.cor=m_met.cor, m=m_met)))
}

build_v2_gmm <- function(v2_data, b_flds) {
  fit_gmm <- function(flds) {
    df <- v2_data[, flds]
    df <- df[complete.cases(df), ]
    bic <- mclustBIC(df, G=1:20)
    m <- Mclust(df, x=bic)
    return(m)
  }

  set.seed(17)
  m <- list(
    gmm_ass = fit_gmm(b_flds$b_geo$independ),
    gmm_com = fit_gmm(b_flds$b_com$depend),
    gmm_met = fit_gmm(b_flds$b_met$depend)
  )

  saveRDS(m, "./models/v2/v2_gmm.Rds")

  return(m)
}

build_v2_grps <- function(v2_umap, b_flds) {
  x <- as.data.frame(v2_umap$layout)

  #k_clusters <- parameters::n_clusters_silhouette(x)
  #k_clusters <- as.numeric(k_clusters$n_Clusters)[which.max(k_clusters$Silhouette)]
  k_clusters <- 5

  ms_h <- round(quantile(dist(x), 0.01), 2)
  set.seed(123)

  m <- list(
    km = kmeans(x, centers=k_clusters, iter.max=100, nstart=3),
    pam = pam(x, k=k_clusters, metric="euclidean", stand=FALSE),
    clara = clara(x, k=k_clusters, metric="euclidean",
                  stand=FALSE, samples=50, pamLike=FALSE),
    gmm = GMM(x, k_clusters, "eucl_dist", "random_subset", 10, 10),
    ms = meanShift(as.matrix(x), kernelType="EPANECHNIKOV",
                   iterations=30, epsilonCluster=ms_h),
    ap0 = apcluster(apcluster::negDistMat(r=2), x),
    apK = apclusterK(apcluster::negDistMat(r=2), x, K=k_clusters)
  )

  saveRDS(m, "./models/v2/v2_groups.Rds")

  return(m)
}

build_v2_class <- function(v2_data, b_flds, v2_grps) {
  x <- v2_data[, b_flds$b_geo$independ]
  x <- x[complete.cases(x), ]
  x$group <- factor(v2_grps$km$cluster)

  m <- rpart(group ~ ., data=x, model=TRUE)
  saveRDS(m, "./models/v2/v2_groups-dt.Rds")

  return(m)
}

build_v2_bn <- function(v2_data, b_flds, v2_gmm) {
  flds1 <- unlist(b_flds$b_geo)
  flds2 <- b_flds$b_com$depend
  flds3 <- b_flds$b_met$depend
  df <- v2_data[, c(flds1, flds2, flds3)]

  # Learn BN for c_geo
  flds <- unlist(b_flds$b_geo)
  x <- df[, flds]
  x <- x[complete.cases(x), ]
  bn <- bnlearn::hc(x)
  bn <- bn.fit(bn, x)
  df.imputed <- impute(bn, df[, flds])
  z <- predict(v2_gmm$gmm_ass, newdata=df.imputed[, colnames(v2_gmm$gmm_ass$data)])
  df.imputed$gmm_ass <- factor(z$classification)

  # Learn BN for c_com
  flds <- unlist(b_flds$b_com)
  x <- df[, flds]
  x <- x[complete.cases(x), ]
  bn <- bnlearn::hc(x)
  bn <- bn.fit(bn, x)
  y <- cbind(df.imputed, df[, setdiff(names(x), names(df.imputed))])
  y <- impute(bn, y[, names(x)])
  df.imputed <- cbind(df.imputed, y[, setdiff(names(y), names(df.imputed))])

  # Learn BN for c_met
  flds <- unlist(b_flds$b_met)
  flds <- flds[!grepl("calc", names(flds))]
  x <- df[, flds]
  x <- x[complete.cases(x), ]
  bn <- bnlearn::hc(x)
  bn <- bn.fit(bn, x)
  y <- cbind(df.imputed, df[, setdiff(names(x), names(df.imputed))])
  y <- impute(bn, y[, names(x)])
  df.imputed <- cbind(df.imputed, y[, setdiff(names(y), names(df.imputed))])

  # Blacklist
  bl <- tiers2blacklist(list(flds1, flds2))
  bl <- rbind(bl, tiers2blacklist(list(c(flds1,flds2), flds3)))

  # Whitelist
  wl <- tiers2blacklist(list(flds1, "gmm_ass"))

  # Learn DAG
  set.seed(17)
  bn <- bnlearn::hc(df.imputed, whitelist=wl, blacklist=bl)

  saveRDS(bn, "./models/v2/v2_bn.Rds")

  return(bn)
}

build_v2_reg <- function(v2_data, b_flds) {
  set.seed(761)
  ctrl <- trainControl(method="repeatedcv", number=10, repeats=10, returnData=FALSE)

  #library(doParallel)
  #cl <- makeCluster(detectCores() - 1)
  #registerDoParallel(cl)

  flds <- sapply(b_flds, "[[", "depend")
  flds <- sapply(flds, length)
  flds <- b_flds[flds>0]
  m <- lapply(flds, function(i) {
    x <- v2_data[, c(i$factors, i$independ), drop=FALSE]
    y <- v2_data[, i$depend, drop=FALSE]

    apply(y, 2, function(yy) {
      k <- complete.cases(cbind(x,yy))
      train(x[k,], yy[k], method="cubist", trControl=ctrl)
    })
  })
  names(m) <- names(flds)

  #stopCluster(cl)
  #registerDoSEQ()

  saveRDS(m, "./models/v2/v2_jura-regression.Rds")

  return(m)
}
