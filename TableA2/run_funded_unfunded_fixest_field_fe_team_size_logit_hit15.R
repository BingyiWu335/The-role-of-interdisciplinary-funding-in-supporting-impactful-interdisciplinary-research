#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

setFixest_nthreads(0)

input_dir <- "./data"
out_dir <- "./TableA2/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

controls <- paste(
  "ln_ref",
  "ln_team_size",
  "ln_institution_count",
  "ln_A_HI_mean_reg",
  "ln_J_HI_reg",
  sep = " + "
)

specs <- list(
  MDP = file.path(input_dir, "mdp_reg_input.parquet"),
  IDP = file.path(input_dir, "idp_reg_input.parquet")
)

fit_one <- function(data, outcome) {
  fml <- as.formula(paste0(
    outcome,
    " ~ funded * RS + ",
    controls,
    " | field + paper_year"
  ))
  if (outcome == "Hit15") {
    feglm(fml, data = data, family = binomial(link = "logit"), vcov = ~ year_source_cluster)
  } else {
    feols(fml, data = data, vcov = ~ year_source_cluster)
  }
}

safe_r2 <- function(fit, type) {
  out <- tryCatch(as.numeric(r2(fit, type)), error = function(e) NA_real_)
  if (length(out) == 0) NA_real_ else out[1]
}

rows <- list()
samples <- list()
stats <- list()

for (group in names(specs)) {
  path <- specs[[group]]
  message("[", Sys.time(), "] Reading ", group, ": ", path)
  dt <- as.data.table(nanoparquet::read_parquet(path))
  dt[, funded := as.integer(funded)]
  dt[, field := factor(field)]
  dt[, paper_year := factor(paper_year)]
  dt[, year_source_cluster := factor(year_source_cluster)]

  samples[[group]] <- data.table(
    group = group,
    n = nrow(dt),
    funded_n = dt[funded == 1L, .N],
    control_n = dt[funded == 0L, .N],
    field_n = uniqueN(dt$field),
    year_n = uniqueN(dt$paper_year),
    source_n = uniqueN(dt$source_id),
    cluster_n = uniqueN(dt$year_source_cluster)
  )

  for (outcome in c("Hit15", "Log1p_NC1")) {
    message("[", Sys.time(), "] Fitting ", group, " ", outcome)
    fit <- fit_one(dt, outcome)
    ct <- as.data.table(coeftable(fit), keep.rownames = "term")
    setnames(ct, old = c("Estimate", "Std. Error"), new = c("estimate", "std_error"), skip_absent = TRUE)
    stat_col <- intersect(names(ct), c("t value", "z value"))
    p_col <- intersect(names(ct), c("Pr(>|t|)", "Pr(>|z|)"))
    if (length(stat_col) == 1) setnames(ct, stat_col, "statistic")
    if (length(p_col) == 1) setnames(ct, p_col, "p_value")

    keep_terms <- c("funded", "RS", "funded:RS", "RS:funded")
    ct <- ct[term %in% keep_terms]
    ct[, group := group]
    ct[, outcome := outcome]
    ct[, model_family := ifelse(outcome == "Hit15", "binomial_logit", "ols")]
    ct[, nobs := nobs(fit)]
    model_r2 <- if (outcome == "Hit15") safe_r2(fit, "pr2") else safe_r2(fit, "ar2")
    ct[, r2 := model_r2]
    ct[, odds_ratio := fifelse(outcome == "Hit15", exp(estimate), NA_real_)]
    ct[, odds_ratio_ci95_low := fifelse(outcome == "Hit15", exp(estimate - 1.96 * std_error), NA_real_)]
    ct[, odds_ratio_ci95_high := fifelse(outcome == "Hit15", exp(estimate + 1.96 * std_error), NA_real_)]
    rows[[paste(group, outcome, sep = "_")]] <- ct[, .(
      group, outcome, model_family, term, estimate, std_error, statistic, p_value,
      odds_ratio, odds_ratio_ci95_low, odds_ratio_ci95_high, nobs, r2
    )]
    stats[[paste(group, outcome, sep = "_")]] <- data.table(
      group = group,
      outcome = outcome,
      model_family = ifelse(outcome == "Hit15", "binomial_logit", "ols"),
      nobs = nobs(fit),
      r2 = model_r2,
      r2_type = ifelse(outcome == "Hit15", "pseudo_r2", "adjusted_r2")
    )
  }
}

core <- rbindlist(rows, fill = TRUE)
core[, stars := fifelse(p_value < 0.001, "***",
                 fifelse(p_value < 0.01, "**",
                 fifelse(p_value < 0.05, "*",
                 fifelse(p_value < 0.1, "+", ""))))]

fwrite(core, file.path(out_dir, "fixest_core_terms.csv"))
fwrite(rbindlist(samples), file.path(out_dir, "fixest_samples.csv"))
fwrite(rbindlist(stats), file.path(out_dir, "fixest_model_stats.csv"))

print(core)
print(rbindlist(samples))
print(rbindlist(stats))
