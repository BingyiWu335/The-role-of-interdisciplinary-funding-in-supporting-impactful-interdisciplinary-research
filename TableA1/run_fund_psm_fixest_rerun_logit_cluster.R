#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

setFixest_nthreads(0)

args <- commandArgs(trailingOnly = TRUE)
base_name <- if (length(args) >= 1) args[[1]] else "psm_team_rs_1v1_norepl_opt"
out_suffix <- if (length(args) >= 2) args[[2]] else base_name

out_dir <- file.path("./TableA1/output", out_suffix)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

matched_csv <- file.path("./data", paste0(base_name, "_matched.parquet"))

message("[", Sys.time(), "] Reading matched sample: ", matched_csv)
dt <- as.data.table(nanoparquet::read_parquet(matched_csv))
dt <- dt[ps_weight > 0]

dt[, ID015F := as.integer(ID015F)]
dt[, Funder_id := as.integer(Funder)]
dt[, Funder_fe := factor(Funder_id)]
dt[, field := factor(field_id)]
dt[, years := factor(year)]
dt[, match_cluster := interaction(Funder_id, field_id, year, drop = TRUE)]

dt[, RSC := RS - mean(RS, na.rm = TRUE)]
dt[, reference_count_ln := log1p(reference_count)]
dt[, team_size_ln := log1p(team_size)]
dt[, institution_count_ln := log1p(institution_count)]
dt[, A_HI_mean_ln := log1p(A_HI_mean)]
dt[, J_HI_ln := log1p(J_HI)]
dt[, team_prior_pub_mean_ln := log1p(team_prior_pub_mean)]
dt[, NC1_ln := log1p(NC1)]

controls <- paste(
  "reference_count_ln",
  "team_size_ln",
  "institution_count_ln",
  "A_HI_mean_ln",
  "J_HI_ln",
  "team_prior_pub_mean_ln",
  sep = " + "
)

model_specs <- list(
  pooled_hit15 = list(funder = NA_integer_, outcome = "Hit15", type = "hit"),
  pooled_nc1 = list(funder = NA_integer_, outcome = "NC1_ln", type = "nc"),
  nsf_hit15 = list(funder = 1L, outcome = "Hit15", type = "hit"),
  nsf_nc1 = list(funder = 1L, outcome = "NC1_ln", type = "nc"),
  nsfc_hit15 = list(funder = 2L, outcome = "Hit15", type = "hit"),
  nsfc_nc1 = list(funder = 2L, outcome = "NC1_ln", type = "nc"),
  nih_hit15 = list(funder = 3L, outcome = "Hit15", type = "hit"),
  nih_nc1 = list(funder = 3L, outcome = "NC1_ln", type = "nc")
)

build_formula <- function(outcome, rhs, pooled) {
  fe_part <- if (pooled) "years + field + Funder_fe" else "years + field"
  as.formula(paste0(outcome, " ~ ", rhs, " + ", controls, " | ", fe_part))
}

fit_model <- function(formula, data_i, type) {
  if (type == "hit") {
    feglm(
      formula,
      data = data_i,
      weights = ~ ps_weight,
      family = binomial(link = "logit"),
      vcov = ~ match_cluster
    )
  } else {
    feols(
      formula,
      data = data_i,
      weights = ~ ps_weight,
      vcov = ~ match_cluster
    )
  }
}

extract_terms <- function(fit, terms, family_label, spec_name, model_name, data_i) {
  ct <- as.data.table(coeftable(fit), keep.rownames = "term")
  setnames(ct, old = c("Estimate", "Std. Error"), new = c("estimate", "std_error"), skip_absent = TRUE)
  stat_col <- intersect(names(ct), c("z value", "t value"))
  p_col <- intersect(names(ct), c("Pr(>|z|)", "Pr(>|t|)"))
  if (length(stat_col) == 1) setnames(ct, stat_col, "statistic")
  if (length(p_col) == 1) setnames(ct, p_col, "p_value")
  ct <- ct[term %in% terms]
  ct[, specification := spec_name]
  ct[, model := model_name]
  ct[, family := family_label]
  ct[, nobs := nobs(fit)]
  ct[, weighted_n := sum(data_i$ps_weight, na.rm = TRUE)]
  ct[, cluster_n := uniqueN(data_i$match_cluster)]
  if (family_label == "binomial_logit") {
    ct[, odds_ratio := exp(estimate)]
    ct[, odds_ratio_ci95_low := exp(estimate - 1.96 * std_error)]
    ct[, odds_ratio_ci95_high := exp(estimate + 1.96 * std_error)]
  } else {
    ct[, odds_ratio := NA_real_]
    ct[, odds_ratio_ci95_low := NA_real_]
    ct[, odds_ratio_ci95_high := NA_real_]
  }
  ct[, .(
    specification, model, family, term, estimate, std_error, statistic, p_value,
    odds_ratio, odds_ratio_ci95_low, odds_ratio_ci95_high, nobs, weighted_n, cluster_n
  )]
}

sample_rows <- list()
term_rows <- list()

for (nm in names(model_specs)) {
  spec <- model_specs[[nm]]
  data_i <- if (is.na(spec$funder)) dt else dt[Funder_id == spec$funder]
  pooled <- is.na(spec$funder)
  family_label <- if (spec$type == "hit") "binomial_logit" else "ols"

  sample_rows[[nm]] <- data.table(
    model = nm,
    n = nrow(data_i),
    weighted_n = sum(data_i$ps_weight, na.rm = TRUE),
    treated_idp_n = data_i[ID015F == 1L, .N],
    control_mdp_n = data_i[ID015F == 0L, .N],
    cluster_n = uniqueN(data_i$match_cluster),
    funder_n = uniqueN(data_i$Funder_id),
    field_n = uniqueN(data_i$field),
    year_n = uniqueN(data_i$years)
  )

  message("[", Sys.time(), "] Fitting ", nm, " RS model")
  f_rs <- build_formula(spec$outcome, "RSC * ID015F", pooled)
  fit_rs <- fit_model(f_rs, data_i, spec$type)
  term_rows[[paste0(nm, "_rs")]] <- extract_terms(
    fit_rs,
    c("RSC", "ID015F", "RSC:ID015F", "ID015F:RSC"),
    family_label,
    "rs_interaction",
    nm,
    data_i
  )
}

core <- rbindlist(term_rows, fill = TRUE)
core[, stars := fifelse(p_value < 0.001, "***",
                 fifelse(p_value < 0.01, "**",
                 fifelse(p_value < 0.05, "*",
                 fifelse(p_value < 0.1, "+", ""))))]

samples <- rbindlist(sample_rows, fill = TRUE)

fwrite(core, file.path(out_dir, paste0(base_name, "_rs_only_fixest_core_terms.csv")))
fwrite(samples, file.path(out_dir, paste0(base_name, "_rs_only_fixest_samples.csv")))

message("[", Sys.time(), "] Wrote outputs to ", out_dir)
print(core)
print(samples)
