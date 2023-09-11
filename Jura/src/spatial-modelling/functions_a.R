source("./src/sum_stats.R")
source("./src/digits.R")
source("./src/model2coefs.R")

a_data_factory <- function() {
  list(
    tar_target(a_bnd_file, "../../data/boundaries/bound.shp", format="file"),
    tar_target(a_geo_file, "../../data/samples/jura_modified.csv", format="file"),
    tar_target(a_com_file, "../../data/metallurgy/jura_comminution.csv", format="file"),
    tar_target(a_met_file, "../../data/metallurgy/jura_flotation.csv", format="file"),
    tar_target(a_bm_file, "../../data/resource/jura_model_ok.rds", format="file"),

    tar_target(a_bnd, readOGR(a_bnd_file)),
    tar_target(a_geo, prep_a_geo(a_geo_file)),
    tar_target(a_com, prep_a_com(a_com_file)),
    tar_target(a_met, prep_a_met(a_met_file)),
    tar_target(a_bm, readRDS(a_bm_file)),
    tar_target(a_meta, a_profile(a_geo, a_com, a_met)),
    tar_target(a_map, map_a_map()),
    tar_target(a_flds, map_a_flds()),
    tar_target(a_calc, map_a_calc())
  )
}

prep_a_geo <- function(a_geo_file) {
  # Load data.
  df <- read.csv(a_geo_file, stringsAsFactors=FALSE, na=c("","NA","-9","-98","-99"))
  df <- as.data.frame(df)
  #write_csv(df, "./data/A_raw/jura_a_geo.csv", na="")
  #saveRDS(df, "./data/A_raw/jura_a_geo.Rds")

  return(df)
}

prep_a_com <- function(a_com_file) {
  # Load data.
  df <- read.csv(a_com_file, stringsAsFactors=FALSE, na=c("","NA"))
  df <- as.data.frame(df)
  #write_csv(df, "./data/A_raw/jura_a_com.csv")
  #saveRDS(df, "./data/A_raw/jura_a_com.Rds")

  return(df)
}

prep_a_met <- function(a_met_file) {
  df <- read.csv(a_met_file, stringsAsFactors=FALSE, na=c("","NA"))
  df <- as.data.frame(df)
  #write_csv(df, "./data/A_raw/jura_a_met.csv", na="")
  #saveRDS(df, "./data/A_raw/jura_a_met.Rds")

  return(df)
}

a_profile <- function(a_geo, a_com, a_met) {
  profile_data <- function(x) {
    # Variables.
    flds_num <- names(x)[sapply(x, is.numeric)]
    flds_char <- setdiff(names(x), flds_num)

    # Summary stats.
    if (length(flds_num) > 0) {
      x.stats <- sapply(x[, flds_num], sum_stats)
    } else {
      x.stats <- NULL
    }

    # Outliers.
    row.names(x) <- 1:nrow(x)
    x <- x[, sapply(x, is.numeric)]
    x <- x[complete.cases(x), ]
    k <- max(3, floor(0.05*nrow(x)))
    score.lof <- lofactor(x, k=k)
    names(score.lof) <- row.names(x)
    cut_out <- median(score.lof, na.rm=TRUE) + 3 * IQR(score.lof, na.rm=TRUE)
    outliers.lof <- names(score.lof)[score.lof > cut_out]
    set.seed(123)
    n <- min(64, ceiling(0.6*nrow(x)))
    tr <- isolationForest(X=x, n_trees=100, Phi=n, parallel=TRUE,
                          future_plan="multicore", extension_level=1, vanilla=TRUE)
    score.if <- invisible(predict.isolationForest(tr, newdata=x))
    names(score.if) <- row.names(x)
    cut_out <- median(score.if, na.rm=TRUE) + 3 * IQR(score.if, na.rm=TRUE)
    outliers.if <- names(score.if)[score.if > cut_out]
    outliers.all <- union(outliers.lof, outliers.if)

    y <- list(
      dim = list(cols=ncol(x), rows=nrow(x)),
      flds = list(numeric=length(flds_num), char=length(flds_char)),
      stats = x.stats,
      missing_rate = round(100 * sum(is.na(x)) / (ncol(x)*nrow(x)), 2),
      outliers = outliers.all
    )

    return(y)
  }

  x <- list(
    a_geo = profile_data(a_geo),
    a_com = profile_data(a_com),
    a_met = profile_data(a_met)
  )

  return(x)
}

map_a_flds <- function() {
  list(
    a_geo = list(
      factors = c("Rock"),
      independ = c("Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn"),
      depend = c("UCS"),
      calc = NULL),

    a_com = list(
      factors = c("Rock"),
      independ = c("Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn", "UCS"),
      depend = c("Ab", "BBMWi"),
      calc = NULL),

    a_met = list(
      factors = c("Rock"),
      independ = c("Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn"),
      depend = c("Cu_tail", "Ni_tail", "Pb_tail", "Zn_tail", "Cd_tail", "Co_tail", "Cr_tail"),
      calc = c("Cu_recov", "Ni_recov", "Pb_recov", "Zn_recov", "Cd_recov", "Co_recov", "Cr_recov")),

    a_bm = list(
      factors = c("rock"),
      independ = c("Cd", "Co", "Cr", "Cu", "Ni", "Pb", "Zn"),
      depend = NULL,
      calc = NULL)
  )
}

map_a_calc <- function(a_flds) {
  m <- list(
    a_geo = NULL,
    a_com = NULL,
    a_met = list(
      Cu_recov="100 * (Cu - Cu_tail) / Cu",
      Ni_recov="100 * (Ni - Ni_tail) / Ni",
      Pb_recov="100 * (Pb - Pb_tail) / Pb",
      Zn_recov="100 * (Zn - Zn_tail) / Zn",
      Cd_recov="100 * (Cd - Cd_tail) / Cd",
      Co_recov="100 * (Co - Co_tail) / Co",
      Cr_recov="100 * (Cr - Cr_tail) / Cr"
    )
  )

  return(m)
}

map_a_map <- function() {
  list(
    a_geo = c(x="Xloc", y="Yloc", rock="Rock", cd_pct="Cd",
              co_pct="Co", cr_pct="Cr", cu_pct="Cu",
              ni_pct="Ni", pb_pct="Pb", zn_pct="Zn", ucs="UCS"),

    a_com = c(x="x", y="y", cd_pct="Cd", co_pct="Co", cr_pct="Cr",
              cu_pct="Cu", ni_pct="Ni", pb_pct="Pb", zn_pct="Zn",
              rock="Rock", ucs="UCS", axb="Ab", bbmwi="BBMWi"),

    a_met = c(x="x", y="y", cd_pct="Cd", co_pct="Co", cr_pct="Cr", cu_pct="Cu",
              ni_pct="Ni", pb_pct="Pb", zn_pct="Zn", rock="Rock",
              cu_tail="Cu_tail", ni_tail="Ni_tail", pb_tail="Pb_tail", zn_tail="Zn_tail",
              cd_tail="Cd_tail", co_tail="Co_tail", cr_tail="Cr_tail",
              cu_recov="Cu_recov", ni_recov="Ni_recov", pb_recov="Pb_recov",
              zn_recov="Zn_recov", cd_recov="Cd_recov", co_recov="Co_recov", cr_recov="Cr_recov"),

    a_bm = c(cd_pct="Cd", co_pct="Co", cr_pct="Cr", cu_pct="Cu",
             ni_pct="Ni", pb_pct="Pb", zn_pct="Zn", rock="rock"),

    a_pore = NULL,
    a_pmill = NULL
  )
}
