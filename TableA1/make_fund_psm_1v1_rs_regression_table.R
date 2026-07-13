#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

setFixest_nthreads(0)

base_name <- "psm_team_rs_1v1_norepl_opt"
data_dir <- "./data"
matched_csv <- file.path(data_dir, paste0(base_name, "_matched.parquet"))
out_dir <- file.path(data_dir, "fund_psm_1v1_rs_only_table")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

controls <- paste(
  "reference_count_ln",
  "team_size_ln",
  "institution_count_ln",
  "A_HI_mean_ln",
  "J_HI_ln",
  "team_prior_pub_mean_ln",
  sep = " + "
)

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

model_specs <- list(
  list(col = "(1)", key = "pooled_nc", funder = NA_integer_, funder_label = "Pooled", outcome = "NC1_ln", dv = "NC"),
  list(col = "(2)", key = "pooled_hit15", funder = NA_integer_, funder_label = "Pooled", outcome = "Hit15", dv = "Hit15"),
  list(col = "(3)", key = "nsf_nc", funder = 1L, funder_label = "NSF", outcome = "NC1_ln", dv = "NC"),
  list(col = "(4)", key = "nsf_hit15", funder = 1L, funder_label = "NSF", outcome = "Hit15", dv = "Hit15"),
  list(col = "(5)", key = "nsfc_nc", funder = 2L, funder_label = "NSFC", outcome = "NC1_ln", dv = "NC"),
  list(col = "(6)", key = "nsfc_hit15", funder = 2L, funder_label = "NSFC", outcome = "Hit15", dv = "Hit15"),
  list(col = "(7)", key = "nih_nc", funder = 3L, funder_label = "NIH", outcome = "NC1_ln", dv = "NC"),
  list(col = "(8)", key = "nih_hit15", funder = 3L, funder_label = "NIH", outcome = "Hit15", dv = "Hit15")
)

stars <- function(p) {
  fifelse(p < 0.001, "***",
    fifelse(p < 0.01, "**",
      fifelse(p < 0.05, "*",
        fifelse(p < 0.1, "+", ""))))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_int <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

safe_r2 <- function(fit, type) {
  out <- tryCatch(as.numeric(r2(fit, type = type)), error = function(e) NA_real_)
  if (length(out) == 0) NA_real_ else out[[1]]
}

build_formula <- function(outcome, pooled) {
  fe_part <- if (pooled) "years + field + Funder_fe" else "years + field"
  as.formula(paste0(outcome, " ~ RSC * ID015F + ", controls, " | ", fe_part))
}

fit_model <- function(formula, data_i, outcome) {
  if (outcome == "Hit15") {
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

extract_coef <- function(fit, spec) {
  ct <- as.data.table(coeftable(fit), keep.rownames = "term")
  setnames(ct, old = c("Estimate", "Std. Error"), new = c("estimate", "std_error"), skip_absent = TRUE)
  stat_col <- intersect(names(ct), c("z value", "t value"))
  p_col <- intersect(names(ct), c("Pr(>|z|)", "Pr(>|t|)"))
  if (length(stat_col) == 1) setnames(ct, stat_col, "statistic")
  if (length(p_col) == 1) setnames(ct, p_col, "p_value")
  ct <- ct[term %in% c("ID015F", "RSC", "RSC:ID015F", "ID015F:RSC")]
  ct[term %in% c("RSC:ID015F", "ID015F:RSC"), term := "ID015F:RSC"]
  ct[, `:=`(
    column = spec$col,
    model_key = spec$key,
    funder = spec$funder_label,
    dv = spec$dv,
    outcome = spec$outcome,
    model_family = ifelse(spec$outcome == "Hit15", "binomial_logit", "ols"),
    nobs = nobs(fit),
    weighted_n = sum(fit$weights, na.rm = TRUE),
    cluster_n = uniqueN(fit$fixef_vars) # overwritten below in stats table
  )]
  ct
}

coef_rows <- list()
stats_rows <- list()

for (spec in model_specs) {
  data_i <- if (is.na(spec$funder)) dt else dt[Funder_id == spec$funder]
  pooled <- is.na(spec$funder)
  message("[", Sys.time(), "] Fitting ", spec$key)
  fit <- fit_model(build_formula(spec$outcome, pooled), data_i, spec$outcome)
  coef_rows[[spec$key]] <- extract_coef(fit, spec)

  adj_r2 <- if (spec$outcome == "Hit15") NA_real_ else safe_r2(fit, "ar2")
  pseudo_r2 <- if (spec$outcome == "Hit15") safe_r2(fit, "pr2") else NA_real_
  stats_rows[[spec$key]] <- data.table(
    column = spec$col,
    model_key = spec$key,
    funder = spec$funder_label,
    dv = spec$dv,
    outcome = spec$outcome,
    model_family = ifelse(spec$outcome == "Hit15", "binomial_logit", "ols"),
    n = nrow(data_i),
    nobs = nobs(fit),
    weighted_n = sum(data_i$ps_weight, na.rm = TRUE),
    treated_idp_n = data_i[ID015F == 1L, .N],
    control_mdp_n = data_i[ID015F == 0L, .N],
    cluster_n = uniqueN(data_i$match_cluster),
    field_n = uniqueN(data_i$field),
    year_n = uniqueN(data_i$years),
    controls = "Yes",
    year_fe = "Yes",
    field_fe = "Yes",
    funder_fe = ifelse(pooled, "Yes", "No"),
    se_type = "Clustered by Funder x field x year",
    adj_r2 = adj_r2,
    pseudo_r2 = pseudo_r2
  )
}

coef_dt <- rbindlist(coef_rows, fill = TRUE)
coef_dt[, stars := stars(p_value)]
stats_dt <- rbindlist(stats_rows, fill = TRUE)

fwrite(coef_dt, file.path(out_dir, "fund_psm_1v1_rs_only_coefficients.csv"))
fwrite(stats_dt, file.path(out_dir, "fund_psm_1v1_rs_only_model_stats.csv"))

value_for <- function(term_name, col_name) {
  row <- coef_dt[column == col_name & term == term_name]
  if (nrow(row) == 0) return("")
  paste0(fmt_num(row$estimate[1]), row$stars[1])
}

se_for <- function(term_name, col_name) {
  row <- coef_dt[column == col_name & term == term_name]
  if (nrow(row) == 0) return("")
  paste0("(", fmt_num(row$std_error[1]), ")")
}

stat_for <- function(field, col_name) {
  row <- stats_dt[column == col_name]
  if (nrow(row) == 0) return("")
  row[[field]][1]
}

model_cols <- vapply(model_specs, function(x) x$col, character(1))
table_rows <- list(
  c("Funder", vapply(model_cols, function(col) stat_for("funder", col), character(1))),
  c("DV", vapply(model_cols, function(col) stat_for("dv", col), character(1))),
  c("IDP", vapply(model_cols, function(col) value_for("ID015F", col), character(1))),
  c("", vapply(model_cols, function(col) se_for("ID015F", col), character(1))),
  c("RS", vapply(model_cols, function(col) value_for("RSC", col), character(1))),
  c("", vapply(model_cols, function(col) se_for("RSC", col), character(1))),
  c("IDP * RS", vapply(model_cols, function(col) value_for("ID015F:RSC", col), character(1))),
  c("", vapply(model_cols, function(col) se_for("ID015F:RSC", col), character(1))),
  c("Controls", vapply(model_cols, function(col) stat_for("controls", col), character(1))),
  c("Year FE", vapply(model_cols, function(col) stat_for("year_fe", col), character(1))),
  c("Field FE", vapply(model_cols, function(col) stat_for("field_fe", col), character(1))),
  c("Funder FE", vapply(model_cols, function(col) stat_for("funder_fe", col), character(1))),
  c("SE", rep("Match cluster", length(model_cols))),
  c("Nb. Obs.", vapply(model_cols, function(col) fmt_int(stat_for("nobs", col)), character(1))),
  c("Adj. R2", vapply(model_cols, function(col) fmt_num(stat_for("adj_r2", col)), character(1))),
  c("Pseudo R2", vapply(model_cols, function(col) fmt_num(stat_for("pseudo_r2", col)), character(1)))
)

wide <- as.data.table(do.call(rbind, table_rows))
setnames(wide, c("", model_cols))
fwrite(wide, file.path(out_dir, "fund_psm_1v1_rs_only_regression_table.csv"))

md_lines <- c(
  "# Fund PSM 1v1 Regression Table (RS-only)",
  "",
  "|  | (1) | (2) | (3) | (4) | (5) | (6) | (7) | (8) |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(wide))) {
  vals <- as.character(wide[i])
  vals <- gsub("\\|", "\\\\|", vals)
  md_lines <- c(md_lines, paste0("| ", paste(vals, collapse = " | "), " |"))
}
md_lines <- c(
  md_lines,
  "",
  "Notes: NC is log1p(NC1). Hit15 models are fixed-effects logit models. IDP equals 1 for IDP-funded papers and 0 for matched MDP-funded papers. RS is the Bwu RS measure centered at the matched-sample mean before interacting with IDP. Controls include log1p(reference_count), log1p(team_size), log1p(institution_count), log1p(A_HI_mean), log1p(J_HI), and log1p(team_prior_pub_mean). Standard errors in parentheses are clustered by the matching layer Funder x field x year. + p<0.10, * p<0.05, ** p<0.01, *** p<0.001."
)
writeLines(md_lines, file.path(out_dir, "fund_psm_1v1_rs_only_regression_table.md"))

process_lines <- c(
  "# Fund 数据集 1v1 PSM 回归流程说明（RS-only）",
  "",
  "## 1. 样本口径",
  "",
  "- 使用 Paper2 fund 数据集中的 funded 论文，处理组为 IDP funded papers（ID015F = 1），对照组为 MDP funded papers（ID015F = 0）。",
  "- 本轮最终口径只使用 Bwu 表中的 RS 作为跨学科性指标；不使用 ref_rao_stirling_subfield，也不因为该变量缺失而筛样本。",
  "- 使用的 1v1 PSM 样本为 psm_team_rs_1v1_norepl_opt，这是基于 RS-only base 构造的 1:1 不放回匹配样本。",
  "",
  "## 2. PSM 匹配流程",
  "",
  "- 精确匹配层：Funder x field x year。",
  "- 倾向得分模型的协变量：ln_team_size、ln_institution_count、ln_A_HI_mean、ln_team_prior_pub_mean。",
  "- 倾向得分 caliper：logit propensity score 的 0.2 SD。",
  "- 在每个精确层内做 global-optimal 1:1 no-replacement matching。",
  "- 最终匹配到 1,202,196 个 IDP treated 和 1,202,196 个 MDP controls，共 2,404,392 行。",
  "",
  "## 3. 回归模型",
  "",
  "NC 模型：",
  "",
  "```r",
  "feols(log1p(NC1) ~ RSC * ID015F + controls | year + field, weights = ps_weight, vcov = ~ match_cluster)",
  "```",
  "",
  "Hit15 模型：",
  "",
  "```r",
  "feglm(Hit15 ~ RSC * ID015F + controls | year + field, weights = ps_weight, family = binomial(link = \"logit\"), vcov = ~ match_cluster)",
  "```",
  "",
  "- RSC = RS - mean(RS)，即以匹配样本均值中心化后的 RS。",
  "- Pooled 模型额外加入 Funder fixed effects；分 funder 模型不加入 Funder fixed effects。",
  "- Controls 包括 log1p(reference_count)、log1p(team_size)、log1p(institution_count)、log1p(A_HI_mean)、log1p(J_HI)、log1p(team_prior_pub_mean)。",
  "- 标准误聚类在匹配层：Funder x field x year。",
  "",
  "## 4. 输出文件",
  "",
  paste0("- 回归表：", file.path(out_dir, "fund_psm_1v1_rs_only_regression_table.md")),
  paste0("- 表格 CSV：", file.path(out_dir, "fund_psm_1v1_rs_only_regression_table.csv")),
  paste0("- 系数明细：", file.path(out_dir, "fund_psm_1v1_rs_only_coefficients.csv")),
  paste0("- 模型统计量：", file.path(out_dir, "fund_psm_1v1_rs_only_model_stats.csv"))
)
writeLines(process_lines, file.path(out_dir, "fund_psm_1v1_rs_only_process_zh.md"))

message("[", Sys.time(), "] Wrote outputs to ", out_dir)
print(wide)
