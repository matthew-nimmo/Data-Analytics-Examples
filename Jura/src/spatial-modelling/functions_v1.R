v1_factory <- function() {
  list(
    tar_target(v1_data, prep_v1_data(v0_data)),
    tar_target(v1_grps, build_v1_grps(v1_data, a_flds)),
    tar_target(v1_class, build_v1_class(v1_data, a_flds, v1_grps)),
    tar_target(v1_reg, build_v1_reg(v1_data, a_flds)),
    tar_target(v1_bm, deploy_v1_bm(a_bm, v1_reg, a_calc, "cubist")),
    tar_target(v1_model, deploy_v1_model(v1_data, v0_umap, v0_bn, v0_gmm, v1_grps, v1_class, v1_reg, a_flds, a_calc))
  )
}

prep_v1_data <- function(v0_data) {
  return(v0_data)
}

build_v1_gmm <- function(v1_data, a_flds) {
  fit_gmm <- function(flds) {
    df <- v1_data[v1_data$KEEP, flds]
    bic <- mclustBIC(df, G=1:20)
    m <- Mclust(df, x=bic)
    return(m)
  }

  set.seed(17)
  m <- list(
    gmm_ass = fit_gmm(a_flds$a_geo$independ),
    gmm_com = fit_gmm(a_flds$a_com$depend),
    gmm_met = fit_gmm(a_flds$a_met$depend)
  )

  saveRDS(m, "./models/v1/v1_jura-gmm.Rds")

  return(m)
}

build_v1_grps <- function(v1_data, a_flds) {
  x <- v1_data[v1_data$KEEP, a_flds$a_geo$independ]
  x <- scale(x)
  x <- as.data.frame(x)

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

  saveRDS(m, "./models/v1/v1_jura-groups.Rds")

  return(m)
}

build_v1_class <- function(v1_data, a_flds, v1_grps) {
  set.seed(761)
  ctrl1 <- trainControl(method="repeatedcv", number=10, repeats=10, returnData=FALSE)
  ctrl2 <- trainControl(method="LGOCV", number=10, p=0.75, returnData=FALSE)

  x <- v1_data[v1_data$KEEP, a_flds$a_geo$independ]
  x <- as.matrix(x)
  y <- factor(v1_grps$km$cluster)

  models <- c("null", "naive_bayes", "pls", "CSimca", "lda",
              "kknn", "lvq", "gaussprPoly", "svmRadial", "rpart", "xgbTree")
  m <- lapply(models, \(m) {
    print(m)
    tryCatch(train(x, y, method=m, trControl=ctrl1),
             silent=TRUE,
             error = \(e) {
               train(x, y, method=m, trControl=ctrl2)
             })
  })
  names(m) <- models

  saveRDS(m, "./models/v1/v1_jura-groups-class.Rds")

  return(m)
}

build_v1_reg <- function(v1_data, a_flds) {
  set.seed(761)
  ctrl1 <- trainControl(method="repeatedcv", number=10, repeats=10, returnData=FALSE)
  ctrl2 <- trainControl(method="LGOCV", number=10, p=0.75, returnData=FALSE)

  #library(doParallel)
  #cl <- makeCluster(detectCores() - 1)
  #registerDoParallel(cl)

  flds <- sapply(a_flds, "[[", "depend")
  flds <- sapply(flds, length)
  flds <- a_flds[flds>0]
  models <- c("null", "lm", "enet", "rqnc", "earth", "bayesglm", "gamSpline", "brnn",
              "gaussprRadial", "rvmRadial", "BstLm", "cubist", "bagEarth", "ranger", "gamboost")
  m <- lapply(flds, function(i) {
    x <- v1_data[, i$independ, drop=FALSE]
    x <- as.matrix(x)
    y <- v1_data[, i$depend, drop=FALSE]

    apply(y, 2, function(y) {
      j <- complete.cases(cbind(x,y))
      mm <- lapply(models, \(m) {
        print(m)
        tryCatch(train(x[j,], y[j], method=m, trControl=ctrl1),
                 silent=TRUE,
                 error = \(e) {
                   train(x[j,], y[j], method=m, trControl=ctrl2)
                 })
      })
      names(mm) <- models
      mm
    })
  })
  names(m) <- names(flds)

  #stopCluster(cl)
  #registerDoSEQ()

  saveRDS(m, "./models/v1/v1_jura-regression.Rds")

  return(m)
}

deploy_v1_bm <- function(a_bm, v1_reg, a_calc, model.name) {
  m <- lapply(names(v1_reg), \(m) lapply(v1_reg[[m]], "[[", model.name))
  m <- Reduce("c", m)

  y <- as.data.frame(a_bm)
  indx <- complete.cases(y)
  names(y) <- gsub("rock_rock", "Rock", names(y))
  y$Rock <- factor(y$Rock)

  # Add predictions.
  for (i in names(m)) {
    z <- as.numeric(rep(NA, nrow(y)))
    z[indx] <- predict(m[[i]], y[indx, ])
    y[[i]] <- z
  }

  # Apply calculated variables.
  mc <- Reduce("c", a_calc)
  for (i in names(mc)) {
    y[[i]] <- eval(parse(text=mc[[i]]), y)
  }

  bm <- stack(SpatialPixelsDataFrame(coordinates(a_bm), as.data.frame(y)))
  saveRDS(bm, "./models/v1/v1_jura-block-model.Rds")

  y <- cbind(as.data.frame(coordinates(a_bm)), as.data.frame(y))
  write_csv(y[indx,], "./models/v1/v0_jura-block-model.csv", na="")

  return(bm)
}

deploy_v1_model <- function(v1_data, v0_umap, v0_bn, v0_gmm, v1_grps, v1_class, v1_reg, a_flds, a_calc) {
  x <- v1_data[v1_data$KEEP, colnames(v0_umap$a_geo$data)]
  x <- scale(x, center=attr(v0_umap$a_geo$data, "scaled:center"), scale=attr(v0_umap$a_geo$data, "scaled:scale"))
  um <- as.data.frame(predict(v0_umap$a_geo, x))
  names(um) <- c("x_umap","y_umap")
  x <- cbind(v1_data[v1_data$KEEP, ], um)
  x$GMM <- as.character(v0_gmm$classification)
  x$GRP <- as.character(v1_grps$km$cluster)
  flds <- c("x", "y", "x_umap", "y_umap", "source",
            a_flds$a_geo$factors, "GMM", "GRP",
            a_flds$a_geo$independ,
            Reduce(union, sapply(a_flds, "[[", "depend")),
            names(Reduce("c", a_calc)))
  x <- x[, flds]

  m <- list(
    flds = list(factors = c("GMM", "GRP", a_flds$a_geo$factors),
                independent = a_flds$a_geo$independ,
                dependent = Reduce(union, sapply(a_flds, "[[", "depend")),
                calc = Reduce("c", a_calc)),
    data = x,
    umap = v0_umap$a_geo,
    bn = v0_bn,
    models = list(class = v1_class$xgbTree,
                  reg = lapply(names(v1_reg), \(m) lapply(v1_reg[[m]], "[[", "cubist")))
  )

  predict_reg <- function(x) {
    y <- as.data.frame(x)
    indx <- complete.cases(y[, m$flds$independ])
    y$Rock <- factor(y$Rock)

    # Add predictions.
    for (i in names(m$models$reg)) {
      z <- as.numeric(rep(NA, nrow(y)))
      z[indx] <- predict(m$models$reg[[i]], y[indx, ])
      y[[i]] <- z
    }

    # Apply calculated variables.
    mc <- Reduce("c", a_calc)
    for (i in names(mc)) {
      y[[i]] <- eval(parse(text=mc[[i]]), y)
    }

    return(y)
  }

  predict_class <- function(x) {
    y <- x
    indx <- complete.cases(y[, m$flds$independ])
    y$GRP <- rep(as.numeric(NA), nrow(y))
    y$GRP[indx] <- predict(m$models$class$xgbTree, y[indx, m$flds$independ])

    return(y)
  }

  predict_bn <- function(x) {
    y <- x
    flds <- setdiff(names(m$bn), names(x))
    flds <- union(flds, m$flds$dependent)
    for (i in flds) {
      y[[i]] <- rep(as.numeric(NA), nrow(y))
    }
    flds <- names(m$bn)[sapply(m$bn, "class") == "bn.fit.dnode"]
    for (i in flds) {
      cf <- coef(m$bn[[i]])
      if (length(dim(cf)) > 1) {
        y[[i]] <- factor(y[[i]], levels=row.names(cf))
      } else {
        y[[i]] <- factor(y[[i]], levels=names(cf))
      }
    }

    z <- bnlearn::impute(m$bn, y[, names(m$bn)])
    y <- cbind(x[, setdiff(names(x), names(z))], z)

    # Apply calculated variables.
    mc <- Reduce("c", a_calc)
    for (i in names(mc)) {
      y[[i]] <- eval(parse(text=mc[[i]]), y)
    }

    return(y)
  }

  m$predict <- function(x=m$data, type=c("reg", "class", "bn")) {
    type <- match.arg(type)
    y <- switch(type,
      reg = predict_bn(x),
      class = predict_class(x),
      bn = predict_bn(x)
    )
    return(y)
  }

  saveRDS(m, "./models/v1/v1_geomet_model.Rds")
  return(m)
}
