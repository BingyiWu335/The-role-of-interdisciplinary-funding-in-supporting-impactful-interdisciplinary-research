#!/usr/bin/env bash
set -u
P2=/data/coding/Paper2
R=/data/coding/Paper2_replication
cd "$P2/derived_psm" || exit 1

echo "=== [1/3] convert matched parquet -> csv (estimator reads *_matched.csv) ==="
python3 - <<'PY'
import duckdb, os
os.chdir("/data/coding/Paper2/derived_psm")
for b in ["psm_team_rs_1v1_norepl_opt","psm_team_rsnn_norepl_opt","psm_team_rsnn_repl","psm_team_refsubfield_cosine_norepl"]:
    src, dst = b+"_matched.parquet", b+"_matched.csv"
    if os.path.exists(dst):
        print("exists", dst); continue
    if os.path.exists(src):
        duckdb.sql(f"COPY (SELECT * FROM '{src}') TO '{dst}' (HEADER, DELIMITER ',')")
        print("converted", dst, os.path.getsize(dst))
    else:
        print("MISSING", src)
PY

echo "=== [2/3] run A.1 estimator per panel base-name ==="
for b in psm_team_rs_1v1_norepl_opt psm_team_rsnn_norepl_opt psm_team_rsnn_repl psm_team_refsubfield_cosine_norepl; do
  echo "----- estimator: $b -----"
  Rscript "$R/code/s2_tableA1_psm/run_fund_psm_fixest_rerun_logit_cluster.R" "$b" 2>&1 \
    | grep -viE 'NOTE|singleton|observations removed|fixed-effect|^\s*$' | tail -3
done

echo "=== [3/3] render Panel A table ==="
Rscript "$R/code/s2_tableA1_psm/make_fund_psm_1v1_rs_regression_table.R" 2>&1 | tail -3

echo "=== STAGE2 DONE ==="
