#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(fixest)
})

input_csv <- "./data/paper2_main_fund_lag_input.parquet"
out_dir <- "./Sec5_4/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[", Sys.time(), "] Reading input: ", input_csv)
df <- nanoparquet::read_parquet(input_csv)

df$Funder <- as.factor(df$Funder)
df$field_id <- as.factor(df$field_id)
df$year <- as.factor(df$year)

controls_base <- paste(
  "log1p(reference_count)",
  "log1p(team_size)",
  "log1p(institution_count)",
  "log1p(A_HI_mean)",
  "log1p(J_HI)",
  sep = " + "
)

rhs_base <- paste("RS * ID015F", controls_base, sep = " + ")
rhs_lag <- paste(rhs_base, "fund_lag", sep = " + ")

specs <- list(
  list(sample = "All", funder = NA, fe = "year + field_id + Funder"),
  list(sample = "NSF", funder = 1, fe = "year + field_id"),
  list(sample = "NSFC", funder = 2, fe = "year + field_id"),
  list(sample = "NIH", funder = 3, fe = "year + field_id")
)

outcomes <- list(
  list(outcome = "NC", dep = "log1p(NC1)", family = "ols"),
  list(outcome = "Hit15", dep = "Hit15", family = "logit")
)

vcov_type <- "hetero"

lag_specs <- list(
  list(label = "baseline_full_sample", rhs = rhs_base, needs_lag = FALSE),
  list(label = "baseline_lag_sample", rhs = rhs_base, needs_lag = TRUE),
  list(label = "plus_fund_lag", rhs = rhs_lag, needs_lag = TRUE)
)

extract_fit <- function(fit, model_id, sample_name, outcome, family, lag_control) {
  ct <- as.data.frame(coeftable(fit))
  ct$term <- rownames(ct)
  names(ct) <- sub("Estimate", "estimate", names(ct), fixed = TRUE)
  names(ct) <- sub("Std. Error", "std_error", names(ct), fixed = TRUE)
  stat_col <- intersect(names(ct), c("z value", "t value"))
  p_col <- intersect(names(ct), c("Pr(>|z|)", "Pr(>|t|)"))
  if (length(stat_col) == 1) names(ct)[names(ct) == stat_col] <- "statistic"
  if (length(p_col) == 1) names(ct)[names(ct) == p_col] <- "p_value"
  ct$model_id <- model_id
  ct$sample <- sample_name
  ct$outcome <- outcome
  ct$family <- family
  ct$lag_control <- lag_control
  ct$nobs <- nobs(fit)
  ct$odds_ratio <- if (family == "logit") exp(ct$estimate) else NA_real_
  ct$stars <- ifelse(ct$p_value < 0.001, "***",
    ifelse(ct$p_value < 0.01, "**",
      ifelse(ct$p_value < 0.05, "*",
        ifelse(ct$p_value < 0.1, "+", "")
      )
    )
  )
  ct[, c(
    "model_id", "sample", "outcome", "family", "lag_control",
    "term", "estimate", "std_error", "statistic", "p_value",
    "odds_ratio", "stars", "nobs"
  )]
}

safe_r2 <- function(fit, type) {
  out <- tryCatch(as.numeric(r2(fit, type)), error = function(e) NA_real_)
  if (length(out) == 0) NA_real_ else out[1]
}

fits <- list()
coef_rows <- list()
stat_rows <- list()
sample_rows <- list()

model_i <- 0
for (sp in specs) {
  d0 <- df
  if (!is.na(sp$funder)) {
    d0 <- d0[d0$Funder == as.character(sp$funder), ]
  }

  for (oc in outcomes) {
    for (ls in lag_specs) {
      model_i <- model_i + 1
      model_id <- paste(sp$sample, oc$outcome, ls$label, sep = "_")

      needed <- c(
        "RS", "ID015F", "reference_count", "team_size", "institution_count",
        "A_HI_mean", "J_HI", "year", "field_id", "Funder"
      )
      if (oc$outcome == "NC") needed <- c(needed, "NC1")
      if (oc$outcome == "Hit15") needed <- c(needed, "Hit15")
      if (ls$needs_lag) needed <- c(needed, "fund_lag")
      needed <- unique(needed)
      if (!is.na(sp$funder)) needed <- setdiff(needed, "Funder")

      d <- d0[complete.cases(d0[, needed]), ]
      sample_rows[[model_id]] <- data.frame(
        model_id = model_id,
        sample = sp$sample,
        outcome = oc$outcome,
        lag_control = ls$label,
        rows_before_complete_cases = nrow(d0),
        nobs_input_complete_cases = nrow(d),
        stringsAsFactors = FALSE
      )

      fml <- as.formula(paste0(oc$dep, " ~ ", ls$rhs, " | ", sp$fe))
      message("[", Sys.time(), "] Fitting ", model_id, " n=", nrow(d), " formula=", deparse(fml))

      fit <- if (oc$family == "ols") {
        feols(fml, data = d, vcov = vcov_type, notes = FALSE)
      } else {
        feglm(fml, data = d, family = binomial(link = "logit"), vcov = vcov_type, notes = FALSE)
      }
      fits[[model_id]] <- fit
      coef_rows[[model_id]] <- extract_fit(fit, model_id, sp$sample, oc$outcome, oc$family, ls$label)
      stat_rows[[model_id]] <- data.frame(
        model_id = model_id,
        sample = sp$sample,
        outcome = oc$outcome,
        family = oc$family,
        lag_control = ls$label,
        nobs = nobs(fit),
        adj_r2 = if (oc$family == "ols") safe_r2(fit, "ar2") else NA_real_,
        pseudo_r2 = if (oc$family == "logit") safe_r2(fit, "pr2") else NA_real_,
        year_fe = "Yes",
        field_fe = "Yes",
        funder_fe = ifelse(is.na(sp$funder), "Yes", "No"),
        vcov = vcov_type,
        stringsAsFactors = FALSE
      )
    }
  }
}

coef_out <- do.call(rbind, coef_rows)
core_terms <- coef_out[coef_out$term %in% c("ID015F", "RS", "ID015F:RS", "RS:ID015F", "fund_lag"), ]
stats_out <- do.call(rbind, stat_rows)
samples_out <- do.call(rbind, sample_rows)

write.csv(coef_out, file.path(out_dir, "main_regression_fund_lag_coefficients.csv"), row.names = FALSE)
write.csv(core_terms, file.path(out_dir, "main_regression_fund_lag_core_terms.csv"), row.names = FALSE)
write.csv(stats_out, file.path(out_dir, "main_regression_fund_lag_model_stats.csv"), row.names = FALSE)
write.csv(samples_out, file.path(out_dir, "main_regression_fund_lag_samples.csv"), row.names = FALSE)

sink(file.path(out_dir, "main_regression_fund_lag_summaries.txt"))
cat("Paper2 main regressions with fund_lag robustness\n")
cat("Input:", input_csv, "\n")
cat("RS variable: RS\n")
cat("IDP variable: ID015F\n")
cat("fund_lag: publication year - earliest exact work-funder fund_year from Paper6/ALL_3vec_cos.par\n")
cat("SE:", vcov_type, "\n")
cat("baseline_lag_sample uses the same non-missing fund_lag sample as plus_fund_lag.\n\n")
for (nm in names(fits)) {
  cat("\n====================", nm, "====================\n")
  print(summary(fits[[nm]]))
}
sink()

print(core_terms)
message("[", Sys.time(), "] Done. Outputs in ", out_dir)
