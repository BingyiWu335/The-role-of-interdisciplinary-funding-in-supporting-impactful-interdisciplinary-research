# Analysis code for TableA5 of Wu et al. (2026).

library(readr); library(dplyr); library(stringr); library(lmtest); library(sandwich)
library(car); library(interactions); library(MASS); library(fixest); library(ggplot2)
library(marginaleffects); library(modelsummary); library(pandoc); library(bestNormalize)
library(officer); library(flextable); library(svglite); library(scales)

dir.create("./TableA5", recursive = TRUE, showWarnings = FALSE)

### ----- New Hit -----

df1 <- nanoparquet::read_parquet("./data/paper2_nsf.parquet")
df2 <- nanoparquet::read_parquet("./data/paper2_nsfc.parquet")
df3 <- nanoparquet::read_parquet("./data/paper2_nih.parquet")
df4 <- nanoparquet::read_parquet("./data/paper2_all.parquet")

# Pre-compute log-transformed control variables.
.add_ln <- function(d) dplyr::mutate(d,
  ln_reference_count   = log1p(reference_count),
  ln_team_size         = log1p(team_size),
  ln_institution_count = log1p(institution_count),
  ln_A_HI_mean         = log1p(A_HI_mean),
  ln_J_HI              = log1p(J_HI),
  ln_NC1               = log1p(NC1))
df1 <- .add_ln(df1); df2 <- .add_ln(df2); df3 <- .add_ln(df3); df4 <- .add_ln(df4)

df1$RSC <- scale(df1$RS, center = TRUE, scale = FALSE)
df2$RSC <- scale(df2$RS, center = TRUE, scale = FALSE)
df3$RSC <- scale(df3$RS, center = TRUE, scale = FALSE)
df4$RSC <- scale(df4$RS, center = TRUE, scale = FALSE)

df1$D015C <- scale(df1$D015F, center = TRUE, scale = FALSE)
df2$D015C <- scale(df2$D015F, center = TRUE, scale = FALSE)
df3$D015C <- scale(df3$D015F, center = TRUE, scale = FALSE)
df4$D015C <- scale(df4$D015F, center = TRUE, scale = FALSE)

df1$D010C <- scale(df1$D010F, center = TRUE, scale = FALSE)
df2$D010C <- scale(df2$D010F, center = TRUE, scale = FALSE)
df3$D010C <- scale(df3$D010F, center = TRUE, scale = FALSE)
df4$D010C <- scale(df4$D010F, center = TRUE, scale = FALSE)

df1$D020C <- scale(df1$D020F, center = TRUE, scale = FALSE)
df2$D020C <- scale(df2$D020F, center = TRUE, scale = FALSE)
df3$D020C <- scale(df3$D020F, center = TRUE, scale = FALSE)
df4$D020C <- scale(df4$D020F, center = TRUE, scale = FALSE)





nsf_rs_hit  <- feglm(Hit15 ~ RS + D010F + ln_reference_count + ln_team_size + ln_institution_count + 
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df1, vcov = 'HC3')

nsfc_rs_hit  <- feglm(Hit15 ~ RS + D010F + ln_reference_count + ln_team_size + ln_institution_count +
                        ln_A_HI_mean + ln_J_HI | year + field_id,
                      family = binomial(link = "logit"), data = df2, vcov = 'HC3')

nih_rs_hit  <- feglm(Hit15 ~ RS + D010F + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df3, vcov = 'HC3')

all_rs_hit  <- feglm(Hit15 ~ RS + D010F + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id + Funder,
                     family = binomial(link = "logit"), data = df4, vcov = 'HC3')

tab <- msummary(list("nsf_rs_hit" = nsf_rs_hit, "nsfc_rs_hit" = nsfc_rs_hit,
                     "nih_rs_hit" = nih_rs_hit, "all_rs_hit" = all_rs_hit),
                fmt = 3, stars = TRUE, gof_omit = "AIC|FE|Within|RMSE|^R2$", output = "flextable") %>% autofit() %>% fontsize(size = 7)
flextable::save_as_docx(tab, path = "./TableA5/Hit_RS_and_D010F.docx")






nsf_rs_hit  <- feglm(Hit15 ~ RSC*D010C + ln_reference_count + ln_team_size + ln_institution_count + 
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df1, vcov = 'HC3')

nsfc_rs_hit  <- feglm(Hit15 ~ RSC*D010C + ln_reference_count + ln_team_size + ln_institution_count +
                        ln_A_HI_mean + ln_J_HI | year + field_id,
                      family = binomial(link = "logit"), data = df2, vcov = 'HC3')

nih_rs_hit  <- feglm(Hit15 ~ RSC*D010C + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df3, vcov = 'HC3')

all_rs_hit  <- feglm(Hit15 ~ RSC*D010C + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id + Funder,
                     family = binomial(link = "logit"), data = df4, vcov = 'HC3')

tab <- msummary(list("nsf_rs_hit" = nsf_rs_hit, "nsfc_rs_hit" = nsfc_rs_hit,
                     "nih_rs_hit" = nih_rs_hit, "all_rs_hit" = all_rs_hit),
                fmt = 3, stars = TRUE, gof_omit = "AIC|FE|Within|RMSE|^R2$", output = "flextable") %>% autofit() %>% fontsize(size = 7)
flextable::save_as_docx(tab, path = "./TableA5/Hit_RS_D010F.docx")





nsf_rs_hit  <- feglm(Hit15 ~ RS + D020F + ln_reference_count + ln_team_size + ln_institution_count + 
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df1, vcov = 'HC3')

nsfc_rs_hit  <- feglm(Hit15 ~ RS + D020F + ln_reference_count + ln_team_size + ln_institution_count +
                        ln_A_HI_mean + ln_J_HI | year + field_id,
                      family = binomial(link = "logit"), data = df2, vcov = 'HC3')

nih_rs_hit  <- feglm(Hit15 ~ RS + D020F + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df3, vcov = 'HC3')

all_rs_hit  <- feglm(Hit15 ~ RS + D020F + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id + Funder,
                     family = binomial(link = "logit"), data = df4, vcov = 'HC3')

tab <- msummary(list("nsf_rs_hit" = nsf_rs_hit, "nsfc_rs_hit" = nsfc_rs_hit,
                     "nih_rs_hit" = nih_rs_hit, "all_rs_hit" = all_rs_hit),
                fmt = 3, stars = TRUE, gof_omit = "AIC|FE|Within|RMSE|^R2$", output = "flextable") %>% autofit() %>% fontsize(size = 7)
flextable::save_as_docx(tab, path = "./TableA5/Hit_RS_and_D020F.docx")





nsf_rs_hit  <- feglm(Hit15 ~ RSC*D020C + ln_reference_count + ln_team_size + ln_institution_count + 
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df1, vcov = 'HC3')

nsfc_rs_hit  <- feglm(Hit15 ~ RSC*D020C + ln_reference_count + ln_team_size + ln_institution_count +
                        ln_A_HI_mean + ln_J_HI | year + field_id,
                      family = binomial(link = "logit"), data = df2, vcov = 'HC3')

nih_rs_hit  <- feglm(Hit15 ~ RSC*D020C + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id,
                     family = binomial(link = "logit"), data = df3, vcov = 'HC3')

all_rs_hit  <- feglm(Hit15 ~ RSC*D020C + ln_reference_count + ln_team_size + ln_institution_count +
                       ln_A_HI_mean + ln_J_HI | year + field_id + Funder,
                     family = binomial(link = "logit"), data = df4, vcov = 'HC3')

tab <- msummary(list("nsf_rs_hit" = nsf_rs_hit, "nsfc_rs_hit" = nsfc_rs_hit,
                     "nih_rs_hit" = nih_rs_hit, "all_rs_hit" = all_rs_hit),
                fmt = 3, stars = TRUE, gof_omit = "AIC|FE|Within|RMSE|^R2$", output = "flextable") %>% autofit() %>% fontsize(size = 7)
flextable::save_as_docx(tab, path = "./TableA5/Hit_RS_D020F.docx")
