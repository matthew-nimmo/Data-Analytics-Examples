c_data_factory <- function() {
  list(
    tar_target(c_geo, prep_c_geo(b_geo, b_flds)),
    tar_target(c_com, prep_c_com(b_com, b_flds)),
    tar_target(c_met, prep_c_met(b_met, b_flds)),
    tar_target(c_all, prep_c_all(c_geo, c_com, c_met)),
    tar_target(c_flds, map_c_flds())
  )
}

prep_c_geo <- function(b_geo, b_flds) {
  flds <- unlist(b_flds$b_geo)
  df <- b_geo[, flds]

  #write_csv(df, "./data/C_cleaned/jura_c_geo.csv")
  saveRDS(df, "./data/C_cleaned/jura_c_geo.Rds")

  return(df)
}

prep_c_com <- function(b_com, b_flds) {
  flds <- unlist(b_flds$b_com)
  df <- b_com[, flds]

  #write_csv(df, "./data/C_cleaned/jura_c_com.csv")
  saveRDS(df, "./data/C_cleaned/jura_c_com.Rds")

  return(df)
}

prep_c_met <- function(b_met, b_flds) {
  flds <- unlist(b_flds$b_met)
  df <- b_met[, flds]

  #write_csv(df, "./data/C_cleaned/jura_c_met.csv")
  saveRDS(df, "./data/C_cleaned/jura_c_met.Rds")

  return(df)
}

prep_c_bm <- function(b_bm, b_flds) {
  return(NULL)
}

prep_c_pore <- function(b_pore, b_flds) {
  return(NULL)
}

prep_c_pmill <- function(b_mill, b_flds) {
  return(NULL)
}

prep_c_all <- function(c_geo, c_com, c_met) {
  df1 <- c_geo
  df1$source <- rep("GEO", nrow(df1))

  df2 <- c_com
  df2$source <- rep("COM", nrow(df2))

  df3 <- c_met
  df3$source <- rep("MET", nrow(df3))

  df <- dplyr::bind_rows(df1, df2, df3)

  #write_csv(df, "./data/C_cleaned/jura_c_all.csv")
  saveRDS(df, "./data/C_cleaned/jura_c_all.Rds")

  return(df)
}

map_c_flds <- function() {
  list(
    c_geo = list(
      factors = c("lith"),
      independ = c("cd_pct", "co_pct", "cr_pct",
                   "cu_pct", "ni_pct", "pb_pct", "zn_pct"),
      depend = c("ucs")),

    c_com = list(
      factors = c("lith"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct",
                   "ni_pct", "pb_pct", "zn_pct"),
      depend = c("axb", "bbmwi")),

    c_met = list(
      factors = c("lith"),
      independ = c("cd_pct", "co_pct", "cr_pct", "cu_pct",
                   "ni_pct", "pb_pct", "zn_pct"),
      depend = c("cu_tail", "ni_tail", "pb_tail", "zn_tail", "cd_tail", "co_tail", "cr_tail",
                 "cu_recov", "ni_recov", "pb_recov", "zn_recov", "cd_recov", "co_recov", "cr_recov"))
  )
}
