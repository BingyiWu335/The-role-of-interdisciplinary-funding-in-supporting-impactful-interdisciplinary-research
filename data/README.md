# Anonymized replication data

Column-selected copies of the reproduction inputs: only `work_id` + the variables used in the
table models are kept. Removed: grant/funding identifiers and metadata (`award_id`, `fund_num`,
`fund_year*`, `grant_row_count`, `funder_label`, `*_project_n`, ...) plus unused indicator columns.
`Funder` (NSF/NSFC/NIH, already public in the paper) is kept — it is needed for the pooled Funder FE.

The Parquet files in this directory are version-controlled with Git Large File Storage (Git LFS).

Files: `paper2_{nsf,nsfc,nih,all}.parquet` (Tables 2–5 and A.3–A.6);
`psm_team_*_matched.parquet` (Table A.1, three panels); `{mdp,idp}_reg_input.parquet`
(Table A.2); and `paper2_main_{fund_lag,pure_idp_mdp}_input.parquet` (Section 5.4).
