b_data_factory <- function() {
  list(
    tar_target(b_fldsnz, find_b_fldsnz(a_geo, a_com, a_met, a_flds)),
    tar_target(b_fldscor, find_b_fldscor(a_geo, a_com, a_met, a_flds)),
    tar_target(b_fldsln, find_b_fldsln(a_geo, a_com, a_met, a_flds)),
    tar_target(b_geo, prep_b_geo(a_geo, a_map)),
    tar_target(b_com, prep_b_com(a_com, a_map)),
    tar_target(b_met, prep_b_met(a_met, a_map)),
    tar_target(b_bm, prep_b_bm(a_bm, a_map)),
    tar_target(b_flds, map_b_flds())
  )
}

# Find near-zero variance variables.
find_b_fldsnz <- function(a_geo, a_com, a_met, a_flds) {
  flds <- c(a_flds$a_geo$factors, a_flds$a_geo$independ)
  flds_geo <- nearZeroVar(a_geo[, flds], freqCut=85/15, uniqueCut=10)
  flds_geo <- ifelse(length(flds_geo)>0, flds[flds_geo], NA)

  flds <- c(a_flds$a_com$factors, a_flds$a_com$independ)
  flds_com <- nearZeroVar(a_com[, flds], freqCut=85/15, uniqueCut=10)
  flds_com <- ifelse(length(flds_com)>0, flds[flds_com], NA)

  flds <- c(a_flds$a_met$factors, a_flds$a_met$independ)
  flds_met <- nearZeroVar(a_met[, flds], freqCut=85/15, uniqueCut=10)
  flds_met <- ifelse(length(flds_met)>0, flds[flds_met], NA)

  flds <- list(a_geo=flds_geo, a_com=flds_com, a_met=flds_met)

  return(flds)
}

# Find variables with high correlation (>70%).
find_b_fldscor <- function(a_geo, a_com, a_met, a_flds) {
  flds <- a_flds$a_geo$independ
  z <- cor(a_geo[, flds], use="pairwise.complete.obs")
  z[is.na(z)] <- 0
  flds_geo <- findCorrelation(z, cutoff=0.7)
  flds_geo <- ifelse(length(flds_geo)>0, flds[flds_geo], NA)

  flds <- a_flds$a_com$independ
  z <- cor(a_com[, flds], use="pairwise.complete.obs")
  z[is.na(z)] <- 0
  flds_com <- findCorrelation(z, cutoff=0.7)
  flds_com <- ifelse(length(flds_com)>0, flds[flds_com], NA)

  flds <- a_flds$a_met$independ
  z <- cor(a_met[, flds], use="pairwise.complete.obs")
  z[is.na(z)] <- 0
  flds_met <- findCorrelation(z, cutoff=0.7)
  flds_met <- ifelse(length(flds_met)>0, flds[flds_met], NA)

  flds <- list(a_geo=flds_geo, a_com=flds_com, a_met=flds_met)

  return(flds)
}

# Find variables that are linear combinations of other variables.
find_b_fldsln <- function(a_geo, a_com, a_met, a_flds) {
  flds <- a_flds$a_geo$independ
  z <- a_geo[, flds]
  z <- z[complete.cases(z), ]
  flds_geo <- findLinearCombos(z)$remove
  flds_geo <- ifelse(length(flds_geo)>0, flds[flds_geo], NA)

  flds <- a_flds$a_com$independ
  z <- a_com[, flds]
  z <- z[complete.cases(z), ]
  flds_com <- findLinearCombos(z)$remove
  flds_com <- ifelse(length(flds_com)>0, flds[flds_com], NA)

  flds <- a_flds$a_met$independ
  z <- a_met[, flds]
  z <- z[complete.cases(z), ]
  flds_met <- findLinearCombos(z)$remove
  flds_met <- ifelse(length(flds_met)>0, flds[flds_met], NA)

  flds <- list(a_geo=flds_geo, a_com=flds_com, a_met=flds_met)

  return(flds)
}

prep_b_geo <- function(a_geo, a_map) {
  # Map variables.
  flds <- unlist(a_map$a_geo)
  df <- a_geo[, flds]
  names(df) <- names(flds)

  # Add extra code here.
  df$rock <- toupper(df$rock)
  flds <- names(df)
  flds <- flds[sapply(df, is.numeric)]
  if (length(flds) > 0) {
    x <- df[, flds]
    x[x < 0] <- abs(x[x < 0])
    df[, flds] <- x
  }
  df$ucs <- ifelse(df$ucs == 0, NA, df$ucs)

  # Fill missing values
  i <- is.na(df)
  if (any(i)) {
    set.seed(17)
    x <- mice(df, m=1, maxit=30)
    df <- complete(x)
  }

  #write_csv(df, "./data/B_cleaned/jura_b_geo.csv")
  saveRDS(df, "./data/B_cleaned/jura_b_geo.Rds")

  return(df)
}

prep_b_com <- function(a_com, a_map) {
  # Map variables.
  flds <- unlist(a_map$a_com)
  df <- a_com[, flds]
  names(df) <- names(flds)

  # Fill missing values
  i <- is.na(df)
  if (any(i)) {
    set.seed(17)
    x <- mice(df, m=1, maxit=30)
    df <- complete(x)
  }

  #write_csv(df, "./data/B_cleaned/jura_b_com.csv")
  saveRDS(df, "./data/B_cleaned/jura_b_com.Rds")

  return(df)
}

prep_b_met <- function(a_met, a_map) {
  # Map variables.
  flds <- unlist(a_map$a_met)
  df <- a_met[, flds]
  names(df) <- names(flds)

  # Fill missing values
  i <- is.na(df)
  if (any(i)) {
    set.seed(17)
    x <- mice(df, m=1, maxit=30)
    df <- complete(x)
  }

  #write_csv(df, "./data/B_cleaned/jura_b_met.csv")
  saveRDS(df, "./data/B_cleaned/jura_b_met.Rds")

  return(df)
}

prep_b_bm <- function(a_bm, a_map) {
  # Map variables.
  flds <- unlist(a_map$a_bm)
  df <- a_bm[[flds]]
  names(df) <- names(flds)

  #write_csv(df, "./data/B_cleaned/jura_b_block-model.csv")
  #saveRDS(df, "./data/B_cleaned/jura_b_block-model.Rds")

  return(df)
}

map_b_flds <- function() {
  list(
    b_geo = list(
      factors = c("rock"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct", "ni_pct", "pb_pct", "zn_pct"),
      depend = c("ucs"),
      calc = NULL),

    b_com = list(
      factors = c("rock"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct", "ni_pct", "pb_pct", "zn_pct", "ucs"),
      depend = c("axb", "bbmwi"),
      calc = NULL),

    b_met = list(
      factors = c("rock"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct", "ni_pct", "pb_pct", "zn_pct"),
      depend = c("cu_tail", "ni_tail", "pb_tail", "zn_tail", "cd_tail", "co_tail", "cr_tail"),
      calc = c("cu_recov", "ni_recov", "pb_recov", "zn_recov", "cd_recov", "co_recov", "cr_recov")),

    b_bm = list(
      factors = c("rock"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct", "ni_pct", "pb_pct", "zn_pct"),
      depend = NULL,
      calc = NULL)
  )
}

map_b_calc <- function(b_flds) {
  m <- list(
    b_geo = NULL,
    b_com = NULL,
    b_met = list(
      cu_recov="100 * (cu_pct - cu_tail) / cu_pct",
      ni_recov="100 * (ni_pct - ni_tail) / ni_pct",
      pb_recov="100 * (pb_pct - pb_tail) / pb_pct",
      zn_recov="100 * (zn_pct - zn_tail) / zn_pct",
      cd_recov="100 * (cd_pct - cd_tail) / cd_pct",
      co_recov="100 * (co_pct - co_tail) / co_pct",
      cr_recov="100 * (cr_pct - cr_tail) / cr_pct"
    )
  )

  return(m)
}
