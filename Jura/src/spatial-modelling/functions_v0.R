source("./src/digits.R")
source("./src/sum_stats.R")
source("./src/plot_missing.R")

my.pal <- colorRampPalette(
  rgb(c(251,169,0,205,245),
      c(234,36,0,205,245),
      c(124,24,0,205,245), maxColorValue=255), alpha=1)

prep_a_geo <- function(a_geo_file) {
  df <- read.csv(a_geo_file, stringsAsFactors=FALSE, na=c("","NA","-9","-98","-99")) %>%
    dplyr::rename(x=Xloc, y=Yloc) %>%
    dplyr::select("x", "y", "Rock", "Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn", "UCS") %>%
    dplyr::mutate(source = "a_geo")
  df <- as.data.frame(df)

  return(df)
}

prep_a_com <- function(a_com_file) {
  df <- read.csv(a_com_file, stringsAsFactors=FALSE, na=c("","NA")) %>%
    dplyr::select("x", "y", "Rock", "Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn",
                  "UCS", "Ab", "BBMWi") %>%
    dplyr::mutate(source = "a_com")
  df <- as.data.frame(df)

  return(df)
}

prep_a_met <- function(a_met_file) {
  df <- read.csv(a_met_file, stringsAsFactors=FALSE, na=c("","NA")) %>%
    dplyr::select("x", "y", "Rock", "Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn",
                  "Cu_recov", "Ni_recov", "Pb_recov", "Zn_recov", "Cd_recov", "Co_recov", "Cr_recov") %>%
    dplyr::mutate(source = "a_met")
  df <- as.data.frame(df)

  return(df)
}

prep_a_bm <- function(a_bm_file) {
  df <- readRDS(a_bm_file)
  return(df)
}

prep_v0_data <- function(a_geo, a_com, a_met) {
  df <- dplyr::bind_rows(a_geo, a_com, a_met)
  write_csv(df, "./data/A_raw/jura_a_dat.csv", na="")
  saveRDS(df, "./data/A_raw/jura_a_dat.Rds")

  return(df)
}

build_v0_umap <- function(v0_data, flds, n_neighbors=15) {
  set.seed(21)

  x <- v0_data[, flds]
  x <- x[complete.cases(x), ]
  x <- scale(x)
  um <- umap::umap(x, n_components=2, n_neighbors=n_neighbors, min_dist=1e-10)

  return(um)
}

build_v0_ucov <- function(a_bm, v0_data, v0_umap) {
  set.seed(1)

  flds <- colnames(v0_umap$data)
  z <- a_bm[[flds]]
  z <- as.data.frame(z)
  z <- cbind(coordinates(a_bm), z)
  indx <- complete.cases(z) & !duplicated(z)
  z <- z[indx, ]

  y <- z[, flds]
  y <- scale(y, center=attr(v0_umap$data, "scaled:center"), scale=attr(v0_umap$data, "scaled:scale"))
  y <- as.data.frame(predict(v0_umap, y))

  x_met <- v0_data[v0_data$source=="a_met", flds]
  x_met <- scale(x_met, center=attr(v0_umap$data, "scaled:center"), scale=attr(v0_umap$data, "scaled:scale"))
  x_met <- as.data.frame(predict(v0_umap, x_met))
  x_met <- cbind(x_met, v0_data[v0_data$source=="a_met", colnames(v0_umap$data)[1]])
  names(x_met) <- c("x", "y", colnames(v0_umap$data)[1])

  m_met <- gausspr(x_met[, 1:2], x_met[,3], kernel="rbfdot", scaled=FALSE, variance.model=TRUE)
  nn_met <- get.knnx(x_met[, 1:2], y, k=1, algorithm="kd_tree")
  m_met.cor <- sapply(1:nrow(y), function(i) {
    m_met@kernelf(as.numeric(y[i,]), as.numeric(x_met[nn_met$nn.index[i],1:2]))
  })
  o_met <- order(m_met.cor, decreasing=TRUE)
  z_met <- predict(m_met, y, type="sdeviation")

  return(data.frame(id=(1:nrow(z))[indx], x=y[,1], y=y[,2], z=m_met.cor))
}

build_v0_grps <- function(v0_umap, k=3) {
  set.seed(21)
  m <- kmeans(v0_umap$layout, k)

  return(m)
}

build_v0_fa <- function(v0_data, v0_flds) {
  x <- v0_data[, v0_flds$assays]
  m_ass <- list(
    twofactor = fa(x, nfactors=2, rotate="oblimin", fm="minres"),
    threefactor = fa(x, nfactors=3, rotate="oblimin", fm="minres")
  )

  x <- v0_data[v0_data$source=="a_met", v0_flds$met]
  m_met <- list(
    twofactor = fa(x, nfactors=2, rotate="oblimin", fm="minres"),
    threefactor = fa(x, nfactors=3, rotate="oblimin", fm="minres")
  )

  return(list(assays=m_ass, met=m_met))
}

build_v0_bn <- function(v0_data, v0_umap1, v0_grps1, v0_flds) {
  flds0 <- v0_flds$factors
  flds1 <- v0_flds$assays
  flds2 <- v0_flds$met
  flds <- c(flds0, flds1, flds2)

  x <- v0_data[, flds]
  for (i in flds0) {
    x[[i]] <- factor(x[[i]])
  }
  x$GMM1 <- as.character(rep(NA, nrow(x)))
  x$x1 <- as.numeric(rep(NA, nrow(x)))
  x$y1 <- as.numeric(rep(NA, nrow(x)))

  i <- as.numeric(names(v0_grps1$cluster))
  x$GMM1[i] <- as.character(v0_grps1$cluster)
  x$GMM1 <- factor(x$GMM1)

  i <- as.numeric(row.names(v0_umap1$data))
  x$x1[i] <- v0_umap1$layout[,1]
  x$y1[i] <- v0_umap1$layout[,2]

  z <- x[v0_data$source=="a_met", ]

  wl <- tiers2blacklist(list(c(flds1,flds2), c("GMM1","x1","y1")))
  wl <- rbind(wl, data.frame(from=flds1, to=paste0(flds1, "_recov")))
  bl <- tiers2blacklist(list(flds1,flds2))

  set.seed(1711)
  dag <- empty.graph(nodes=names(z))
  bn <- bn.fit(dag, z, replace.unidentifiable=TRUE)
  bn0 <- structural.em(z, start=bn, return.all=TRUE, impute="bayes-lw", impute.args=list(n=500),
                       maximize.args=list(whitelist=wl, blacklist=bl))

  bn <- structural.em(x, start=bn0$fitted, return.all=TRUE, impute="bayes-lw", impute.args=list(n=500),
                      maximize.args=list(whitelist=wl, blacklist=bl))

  bn <- structural.em(z, start=bn$fitted, return.all=TRUE, impute="bayes-lw", impute.args=list(n=500),
                      maximize.args=list(whitelist=wl, blacklist=bl))

  return(bn)
}

compare_v0_bn <- function(v0_bn, v0_data, v0_flds) {
  flds <- names(v0_bn$fitted)
  y <- v0_bn$imputed

  i <- as.numeric(row.names(y))
  x <- v0_data[i, v0_flds$met]
  for (i in setdiff(flds, names(x))) {
    x[[i]] <- as.numeric(rep(NA, nrow(x)))
  }
  for (i in names(y)[sapply(y, is.factor)]) {
    x[[i]] <- factor(x[[i]], levels=levels(y[[i]]))
  }
  x <- x[, flds]
  x <- impute(v0_bn$fitted, data=x)

  return(list(x=x, y=y))
}

impute_v0_bn <- function(v0_bn, v0_data) {
  flds <- names(v0_bn$fitted)
  x <- v0_data
  for (i in setdiff(flds, names(x))) {
    x[[i]] <- as.numeric(rep(NA, nrow(x)))
  }
  for (i in flds[sapply(v0_bn$imputed, is.factor)]) {
    x[[i]] <- factor(x[[i]], levels=levels(v0_bn$imputed[[i]]))
  }
  y <- impute(v0_bn$fitted, data=x[, flds])
  flds <- setdiff(names(x), flds)
  y <- cbind(x[, flds], y)

  return(y)
}

create_v0_fig_miss <- function(v0_data) {
  plot_missing(v0_data, seg=20, horiz=TRUE) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.background = element_rect(color=1, size=1))
}

create_v0_fig_umap1 <- function(a_bm, v0_umap, v0_data) {
  x_bm <- a_bm[[colnames(v0_umap$data)]]
  x_bm <- as.data.frame(x_bm)
  x_bm <- x_bm[complete.cases(x_bm), ]
  x_bm <- x_bm[!duplicated(x_bm), ]
  x_bm <- scale(x_bm, center=attr(v0_umap$data, "scaled:center"), scale=attr(v0_umap$data, "scaled:scale"))
  x_bm <- predict(v0_umap, x_bm)
  x_bm <- as.data.frame(x_bm)
  x_bm$source <- "bm"

  #x_geo <- as.data.frame(v0_umap$layout)
  #x_geo$source <- v0_data$source[as.numeric(row.names(v0_umap$data))]
  #x_geo$source <- gsub("a_", "", x_geo$source)

  x_geo <- v0_data[, colnames(v0_umap$data)]
  x_geo <- scale(x_geo, center=attr(v0_umap$data, "scaled:center"), scale=attr(v0_umap$data, "scaled:scale"))
  x_geo <- as.data.frame(x_geo)
  i <- complete.cases(x_geo)
  x_geo <- x_geo[i, ]
  x_geo <- predict(v0_umap, x_geo)
  x_geo <- as.data.frame(x_geo)
  x_geo$source <- gsub("a_", "", v0_data$source[i])

  x <- rbind(x_bm, x_geo)
  x <- x[!duplicated(x), ]

  ggplot(data=x, mapping=aes(x=V1, y=V2, color=source, shape=source, size=source)) +
    geom_point() +
    scale_colour_manual(values=c("bm"="grey90","geo"="grey50","com"="blue","met"="red")) +
    scale_shape_manual(values=c("bm"=15, "geo"=15,"com"=3,"met"=17)) +
    scale_size_manual(values=c("bm"=3, "geo"=2,"com"=2,"met"=2)) +
    labs(x=NULL, y=NULL) +
    coord_fixed() +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_tufte()
}

create_v0_fig_dag <- function() {
  geomet_oa_dag <- dagify(
    g ~ Na + Mg + K + Ca + C + Fe + Al,
    p ~ Sc + Ti + V + Cr + Y + Zr + Nb + Hf + Ta + La + Ce + Th + Al,
    z ~ Cu + Ni + Pb + Zn + Cd + Co,
    pr ~ C + Ca,
    ac ~ Ca,
    cl ~ K + Na,
    log ~ lith + altn,
    mt ~ log + g + p,
    dom ~ mt + z,
    dens ~ g + z,
    h ~ g + z + dens,
    ab ~ dom + h + dens,
    bmwi ~ dom + h + dens,
    ai ~ dom,
    gs ~ ab + bmwi,
    r ~ z + dom + gs,
    labels = c(g = "gangue",
               p = "protolith",
               mt = "material",
               z = "mineralization",
               r = "recovery",
               pr = "preg-robbing",
               ac = "acid consuming",
               cl = "clays"),
    exposure = "g",
    latent = "mt",
    outcome = "r"
  )

  geomet_oa_dag %>%
    tidy_dagitty(layout="fr", seed=14) %>%
    ggplot(aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_dag_edges_link(aes(start_cap=ggraph::circle(3,"mm"),
                            end_cap=ggraph::circle(3,"mm"))) +
    geom_dag_point(size=8) +
    geom_dag_text(size=2) +
    geom_dag_label_repel(
      aes(label=label, fill=label),
      segment.alpha = 1,
      size = 2.5,
      fill = list("lightblue","lightblue","lightblue","lightblue","white","white","white","salmon"),
      show.legend = FALSE,
      box.padding = grid::unit(0.5, "lines")
    ) +
    theme_dag() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.background = element_rect(color=1, size=1))

}

create_v0_fig_ucov1 <- function(v0_data, v0_umap, v0_ucov) {
  x <- v0_ucov[, -1]
  k <- round(100 * sum(x$z>0.9) / nrow(x), 0)
  x$z <- cut(x$z, breaks=c(0, 0.8, 0.9, 1))
  levels(x$z) <- c("gap", "partial gap", "no gap")
  x$z <- as.character(x$z)

  y <- v0_data[v0_data$source=="a_met", colnames(v0_umap$data)]
  y <- scale(y, center=attr(v0_umap$data, "scaled:center"), scale=attr(v0_umap$data, "scaled:scale"))
  y <- predict(v0_umap, y)
  y <- as.data.frame(y)
  names(y) <- c("x", "y")
  y$z <- "sample"

  z <- rbind(x, y)
  names(z) <- c("x","y","coverage")

  ggplot(data=z, mapping=aes(x=x, y=y, color=coverage, shape=coverage)) +
    geom_point(size=2) +
    scale_colour_manual(values=c("gap"="grey20","partial gap"="grey50","no gap"="grey90","sample"="red")) +
    scale_shape_manual(values=c("gap"=15,"partial gap"=15,"no gap"=15,"sample"=17)) +
    ggtitle(paste0("UMAP Coverage = ", k, "%")) +
    labs(x=NULL, y=NULL) +
    coord_fixed() +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_tufte()
}

create_v0_fig_eda1 <- function(v0_data, v0_flds) {
  p <- lapply(v0_flds$met, \(i) {
    y <- v0_data[, c(i, "Rock")]
    y <- y[complete.cases(y), ]
    ggplot(y, aes_string(x=i, y="Rock")) +
      geom_boxplot(notch=FALSE, orientation="y", fill="grey70") +
      labs(y=NULL) +
      theme_tufte()
  })
  p <- arrangeGrob(grobs=p, ncol=5, nrow=2)

  return(p)
}

create_v0_fig_eda2 <- function(v0_data, v0_flds, v0_grps2){
  x <- v0_data[v0_data$source=="a_met",]
  x$grps <- v0_grps2$cluster
  p <- lapply(v0_flds$assays, \(i) {
    y <- x[, c(i, paste0(i, "_recov"), "grps")]
    y <- y[complete.cases(y), ]
    ggplot(y, aes_string(x=i, y=paste0(i, "_recov"), color="grps")) +
      geom_point(show.legend=FALSE) +
      theme_tufte()
  })
  p <- arrangeGrob(grobs=p, ncol=4, nrow=2)

  return(p)
}

create_v0_fig_eda3 <- function(v0_data, v0_umap2, v0_grps2) {
  x <- as.data.frame(v0_umap2$layout)
  names(x) <- c("x_umap", "y_umap")
  x <- cbind(x, v0_data[as.numeric(row.names(v0_umap2$layout)),])
  x$grps <- v0_grps2$cluster

  p <- lapply(colnames(v0_umap2$data), \(i) {
    ggplot(data=x, mapping=aes_string(x="x_umap", y="y_umap")) +
      geom_voronoi_tile(aes(fill=grps), show.legend=FALSE) +
      scale_fill_gradientn(colors=grey.colors(4, start=0.6)) +
      geom_voronoi_segment(colour="grey80") +
      geom_point(shape=15, size=2, mapping=aes_string(color=i)) +
      scale_color_gradientn(colors=rev(my.pal(20))) +
      labs(x=NULL, y=NULL) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      theme_tufte()+
      theme(legend.key.size = unit(3, 'mm'))
  })
  p <- arrangeGrob(grobs=p, ncol=3, nrow=3)

  return(p)
}

create_v0_fig_eda4 <- function(v0_data, v0_flds, v0_grps2) {
  x <- v0_data[as.numeric(names(v0_grps2$cluster)), ]
  x$GRP <- factor(v0_grps2$cluster)
  p <- lapply(v0_flds$met, \(i) {
    y <- x[, c(i, "GRP")]
    y <- y[complete.cases(y), ]
    ggplot(y, aes_string(x=i, y="GRP")) +
      geom_boxplot(notch=FALSE, orientation="y", fill="grey70") +
      labs(x=i, y="Group") +
      theme_tufte()
  })
  p <- arrangeGrob(grobs=p, ncol=4, nrow=2)

  return(p)
}
