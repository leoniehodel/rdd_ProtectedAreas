# ===========================================================================
# RDD ANALYSIS
# Structure:
#   1. Setup: outcome labels, donut widths (computed ONCE for all tables)
#   2. Table S2: Summary Stats
#   2. Table S3: RDD Overall
#.  3. Table S4: RDD established pasture
#   4. Table S5: RDD by pasture cohort
#   5. Table S7: Cell-level vs PA-level means comparison table
#   6. Table S8: RDD Indigenous Lands
#   7. Table S9: RDD Indigenous Lands, excluding buffers around settlements
#   8. Table S10: RDD Strictly Protected Areas
#   Leonie Hodel
#  usethis::use_mit_license("Leonie Hodel")
# ===========================================================================
library(dplyr)      # filter, mutate, summarise, case_when, recode, bind_rows
library(tidyr)      # pivot_wider, replace_na
library(purrr)      # map, map_dfr, imap_dfr, safely
library(tibble)     # tibble()
library(readr)      # write_csv
library(rdrobust)   # rdrobust
library(gt)
library(arrow)
library(MatchIt)
source('rdd_helpers.R')
mem.maxVSize(50000)
# ---------------------------------------------------------------------------
# 1. SETUP
# ---------------------------------------------------------------------------

##### Load analysis data ------------------------
analysis_data <- read_parquet("data/rdd_analysis_data.parquet")
##### Define Output directories------------------------

dir_fig_rdd <- "results/figures/rdd_figures"  # main RDD / cattle figures
dir_fig_si  <- "results/figures/SI"           # supplementary figures
dir_tables  <- "results/tables"               # summary tables + raw model csv
# create results folder
invisible(lapply(c(dir_fig_rdd, dir_fig_si, dir_tables),
                 function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

var_labels_rdd <- c(
  pasture_area_perc                = "(1) Pasture (%)",
  deforestation_area_perc          = "(2) Deforestation (%)",
  stocking_rate_ha                 = "(3) Cattle density (all pasture)",
  stocking_rate_ha_deforested_cell = "(4) Cattle density (recently cleared)"
  #cattle_on_established_pasture    = "Cattle density (old pasture)"
)

vars_cattle <- c("stocking_rate_ha", "stocking_rate_ha_deforested_cell",
                 "cattle_on_established_pasture")

## 1.1 Donut widths -------------------------------------------------------
# Donut half-width per outcome = turning point of a polynomial fitted to the
# donut-sensitivity curve (quadratic vertex for pasture, cubic minimum
# otherwise). Computed once and reused by every table below.

# donut widths full sample
donuts_m <- seq(0, 1800, by = 120)
run_donut <- function(gap_m, dat, yvar, xkm = "signed_distance_external_border_km") {
  d <- dat %>% filter(!(get(xkm) >= -gap_m/1000 & get(xkm) < 0))
  if (nrow(d) < 200 || var(d[[xkm]], na.rm = TRUE) == 0 || var(d[[yvar]], na.rm = TRUE) == 0)
    return(tibble(gap_m, est = NA_real_))
  res <- purrr::safely(rdrobust, otherwise = NULL)(y = d[[yvar]], x = d[[xkm]], p = 1, h = 2, vce = "hc0")
  if (is.null(res$result)) return(tibble(gap_m, est = NA_real_))
  tibble(gap_m = gap_m, est = res$result$Estimate[1, "tau.us"])
}
donut_range_km <- map_dfr(names(var_labels_rdd), function(v) {
  message("donut sensitivity: ", v)
  map_dfr(donuts_m, run_donut, dat = analysis_data, yvar = v) %>% mutate(variable = v)
}) %>% filter(!is.na(est), gap_m < 1500) %>%
  group_by(variable) %>%
  group_modify(~{
    qc <- coef(lm(est ~ poly(gap_m, 2, raw = TRUE), data = .x)); quad <- -qc[2] / (2 * qc[3])
    cc <- coef(lm(est ~ poly(gap_m, 3, raw = TRUE), data = .x))
    b1 <- cc[2]; b2 <- cc[3]; b3 <- cc[4]; D <- 4 * b2^2 - 12 * b3 * b1
    cmin <- if (is.na(D) || D < 0) NA_real_ else (-2 * b2 + sqrt(D)) / (6 * b3)
    tibble(donut_m = if (.y$variable == "pasture_area_perc") quad else cmin)
  }) %>% ungroup()

donut_range_km<- donut_range_km  %>% mutate(donut_m = pmax(replace_na(donut_m, 0), 0)) %>% { setNames(.$donut_m / 1000, .$variable) }

# a shared donut keeps this additinal variable comparable and avoids selecting a second width on a smaller subsample
donut_range_km["cattle_on_established_pasture"] <- donut_range_km["stocking_rate_ha"]

cat(" ---------- donut width per outcome (m) ----------------")
print(round(donut_range_km * 1000))


# ---------------------------------------------------------------------------
# 2. TABLE S2 - Summary statistics
# ---------------------------------------------------------------------------

# label, column, decimals...
var_spec <- tribble(
  ~label,                                                                ~col,                               ~digits,
  "Forest (km²)",                                                       "forest_area_km2",                        4,
  "Pasture (km²)",                                                      "pasture_area_km2",                       4,
  "Pasture (%)",                                                        "pasture_area_perc",                      1,
  "Deforestation (km²)",                                                "deforestation_area_km2",                 4,
  "Deforestation (%, if forested in 2017)",                             "deforestation_area_perc",                1,
  "Cattle density (animals/ha pasture)",                                "stocking_rate_ha",                       2,
  "Cattle density on deforested cell (animals/ha deforested pasture)",  "stocking_rate_ha_deforested_cell",       2
)

# category rows: share of cells + number of PAs
cat_rows <- analysis_data %>%
  group_by(category) %>%
  summarise(N = n(), n_pa = n_distinct(pa_id), .groups = "drop") %>%
  mutate(
    Variable = sprintf("... %s (%d PAs)", category, n_pa),
    Mean     = sprintf("%.0f%%", 100 * N / sum(N)),
    SD = "", Min = "", P25 = "", P75 = "", Max = ""
  ) %>%
  select(Variable, N, Mean, SD, Min, P25, P75, Max) %>%
  mutate(N = format(N, big.mark = ","))

# continuous rows: one pass per column, NAs dropped per variable
num_rows <- pmap_dfr(var_spec, function(label, col, digits) {
  x <- analysis_data[[col]]
  x <- x[!is.na(x)]
  q <- unname(quantile(x, c(0.25, 0.75)))
  fmt <- function(v) formatC(v, format = "f", digits = digits, big.mark = ",")
  tibble(
    Variable = label,
    N    = format(length(x), big.mark = ","),
    Mean = fmt(mean(x)), SD  = fmt(sd(x)),
    Min  = fmt(min(x)),  P25 = fmt(q[1]),
    P75  = fmt(q[2]),    Max = fmt(max(x))
  )
})

tab_s2 <- bind_rows(cat_rows, num_rows)
n_cat  <- nrow(cat_rows)

gt_s2 <- tab_s2 %>%
  gt() %>%
  tab_row_group(label = md("*Protected area type*"), rows = 1:n_cat) %>%
  tab_row_group(label = md("*Grid-cell characteristics*"),
                rows = (n_cat + 1):nrow(tab_s2)) %>%
  row_group_order(groups = c("*Protected area type*", "*Grid-cell characteristics*")) %>%
  cols_label(N = "N (grid cells)", SD = md("Std. Dev."),
             P25 = "Pctl. 25", P75 = "Pctl. 75") %>%
  cols_align("left",  columns = Variable) %>%
  cols_align("right", columns = N:Max) %>%
  tab_options(table.font.size = px(13), data_row.padding = px(4)) %>%
  opt_row_striping()
gt_s2
gtsave(gt_s2, file.path(dir_tables, "2_summary_stats.html"))

# ---------------------------------------------------------------------------
# 3. TABLE S3 - Overall RDD, full sample
#    Four outcomes; both the MSE-optimal bandwidth (h) and twice it (2h).
# ---------------------------------------------------------------------------

vars_to_test <- c("pasture_area_perc", "deforestation_area_perc",
                  "stocking_rate_ha", "stocking_rate_ha_deforested_cell")

rdd_donut_all <- rdd_grid(list(All = analysis_data), vars_to_test,
                          donut_range_km)

rdd_clean <- save_rdd(rdd_donut_all, "3_rdd_summary_table_overall",
                      raw_file = "3_rdd_raw_overall.csv",
                      drop = "subset")
rdd_clean

# ---------------------------------------------------------------------------
# 4. TABLE S4 - Cattle on established pasture (no recent clearing)
#    Same specification as Table 3, restricted to one outcome.
#    `cattle_on_established_pasture` is defined at the top of the master
#    script: stocking rate in cells with pasture but no recent deforestation.
# ---------------------------------------------------------------------------

rdd_donut_oldpasture <- rdd_grid(list(All = analysis_data),
                                 "cattle_on_established_pasture",
                                 donut_range_km)

rdd_clean_oldpasture <- save_rdd(
  rdd_donut_oldpasture, "4_rdd_summary_table_overall_established_pasture",
  raw_file = "4_rdd_raw_overall_established_pasture.csv",
  drop = "subset"
)
rdd_clean_oldpasture


# ---------------------------------------------------------------------------
# 5. TABLE 5 - RDD by pasture cohort
#    One RDD per (pasture cohort x outcome); MSE-optimal bandwidth only.
# ---------------------------------------------------------------------------

pasture_cohorts <- analysis_data %>%
  filter(!is.na(pasture_group)) %>%
  split(.$pasture_group)

rdd_donut_results <- rdd_grid(
  pasture_cohorts,
  c("pasture_area_perc", "deforestation_area_perc",
    "stocking_rate_ha", "stocking_rate_ha_deforested_cell"),
  donut_range_km,
  subset_col = "pasture_type"
)

rdd_clean_pasture <- save_rdd(
  rdd_donut_results, "5_rdd_summary_table_pasture_group",
  raw_file = "5_rdd_raw_pasture_group.csv",
  rowname_col = "pasture_type", keep_h_only = TRUE, drop = "bandwidth_type"
)
rdd_clean_pasture

# ---------------------------------------------------------------------------
# 6. TABLE S6  -
# ---------------------------------------------------------------------------

pa_group_category <- analysis_data %>%
  filter(!is.na(pasture_group), !is.na(category)) %>%
  distinct(pa_id, category, pasture_group)

il_sp_group_table <- table(pa_group_category$category,
                           pa_group_category$pasture_group)
il_sp_group_table

il_sp_group_df <- as.data.frame.matrix(addmargins(il_sp_group_table)) %>%
  tibble::rownames_to_column("category")
write_csv(il_sp_group_df, file.path(dir_tables, "6_il_sp_pasture_group_counts.csv"))



# -----------------------------------------------------------------------------
# 7. TABLE S7 -
# Cell-level vs PA-level means comparison table
# -----------------------------------------------------------------------------
bin_width_f =2500
table_inside_only <- TRUE
match_wide <- analysis_data %>%
  mutate(distance_bin = bin_width_f * round(signed_distance_external_border / bin_width_f)) %>%
  filter(distance_bin < 0) %>%
  group_by(pa_id, category, distance_bin) %>%
  summarise(mean_pasture       = mean(pasture_area_perc, na.rm = TRUE),
            mean_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
            mean_forest        = mean(forest_area_km2, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(category), !is.na(mean_pasture)) %>%
  mutate(category = factor(category)) %>%
  pivot_wider(names_from = distance_bin,
              values_from = c(mean_pasture, mean_deforestation, mean_forest)) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))
names(match_wide) <- make.names(names(match_wide))
match_wide <- match_wide[!is.na(match_wide$category) & !is.na(match_wide$mean_forest_.2500), ]

m_out <- matchit(category ~ mean_forest_.2500, data = match_wide,
                 method = "nearest", distance = "glm")
matched_il_ids <- match_data(m_out) %>%
  filter(category == "Indigenous Lands") %>% pull(pa_id)

# duplicate the matched-IL cells as an extra display category
add_matched <- function(df) {
  bind_rows(df,
            df %>%
              filter(category == "Indigenous Lands", pa_id %in% matched_il_ids) %>%
              mutate(category = "Indigenous Lands (matched)"))
}

tab_data <- if (table_inside_only) {
  analysis_data %>% filter(within_pa == 1)   # inside-PA cells only
} else {
  analysis_data
}
tab_data <- add_matched(tab_data)

# cell-level: every grid cell weighted equally (mean and SD across cells);
# n counts only cells with a non-NA value for that variable

var_labels <- c(
  pasture_area_perc         = "Pasture (%)",
  deforestation_area_perc   = "Deforestation (%)",
  stocking_rate_ha          = "Cattle density (animals/pasture ha)",
  stocking_rate_ha_deforested_cell = "Cattle density (animals/recently cleared pasture ha)"
)
vars <- names(var_labels)

cell_stats <- tab_data %>%
  group_by(category) %>%
  summarise(across(all_of(vars), list(cell_mean = ~ mean(.x, na.rm = TRUE),
                                      cell_sd   = ~ sd(.x, na.rm = TRUE),
                                      cell_n    = ~ sum(!is.na(.x))),
                   .names = "{.col}__{.fn}"),
            .groups = "drop")

# PA-level: one value per protected area, then summarise (mean, SD, median across PAs)
pa_level_tab <- tab_data %>%
  group_by(pa_id, category) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# n counts only PAs with a non-NA (non-NaN) mean for that variable
pa_stats <- pa_level_tab %>%
  group_by(category) %>%
  summarise(across(all_of(vars),
                   list(pa_mean   = ~ mean(.x, na.rm = TRUE),
                        pa_sd     = ~ sd(.x, na.rm = TRUE),
                        pa_median = ~ median(.x, na.rm = TRUE),
                        pa_n      = ~ sum(!is.na(.x))),
                   .names = "{.col}__{.fn}"),
            .groups = "drop")

# bin-weighted: mean per 120 m distance bin, then average the inside bins
# (distance_bin > 0) with each bin weighted equally -- matches the SP/IL figure.
bin_width_fig <- 120

bin_stats <- add_matched(analysis_data) %>%
  mutate(distance_bin = bin_width_fig * round(signed_distance_external_border / bin_width_fig)) %>%
  group_by(category, distance_bin) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(distance_bin > 0) %>%
  group_by(category) %>%
  summarise(across(all_of(vars), list(bin_mean = ~ mean(.x, na.rm = TRUE),
                                      bin_sd   = ~ sd(.x, na.rm = TRUE)),
                   .names = "{.col}__{.fn}"), .groups = "drop")

to_long <- function(df) {
  df %>%
    pivot_longer(matches("__"),
                 names_to = c("variable", "stat"), names_sep = "__",
                 values_to = "value") %>%
    select(category, variable, stat, value)
}

comparison <- bind_rows(to_long(cell_stats), to_long(bin_stats), to_long(pa_stats)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(variable = factor(var_labels[variable], levels = var_labels)) %>%
  arrange(variable, category) %>%
  select(variable, category, n_PAs = pa_n, n_cells = cell_n,
         bin_mean, bin_sd, cell_mean, cell_sd, pa_mean, pa_sd, pa_median)

write_csv(comparison, file.path(dir_tables, "7_cell_vs_pa_comparison.csv"))

# combine mean +/- SD for display (2 dp for pasture/deforestation, 3 dp for cattle)
fmt_n <- function(x, d) mapply(function(v, dd) formatC(v, format = "f", digits = dd), x, d)
comparison_disp <- comparison %>%
  mutate(
    dec  = ifelse(grepl("Cattle", as.character(variable)), 3, 2),
    bin  = paste0(fmt_n(bin_mean, dec),  " ± ", fmt_n(bin_sd, dec)),
    cell = paste0(fmt_n(cell_mean, dec), " ± ", fmt_n(cell_sd, dec)),
    pa   = paste0(fmt_n(pa_mean, dec),   " ± ", fmt_n(pa_sd, dec)),
    pamed = fmt_n(pa_median, dec)
  ) %>%
  select(variable, category, n_PAs, n_cells, bin, cell, pa, pamed)

comparison_gt <- comparison_disp %>%
  gt(rowname_col = "category", groupname_col = "variable") %>%
  fmt_number(columns = c(n_PAs, n_cells), decimals = 0, use_seps = TRUE) %>%
  cols_label(
    n_PAs   = md("**N PAs**"),
    n_cells = md("**N cells**"),
    bin     = md("**Bin-weighted mean ± SD**"),
    cell    = md("**Cell-weighted mean ± SD**"),
    pa      = md("**PA-level mean ± SD**"),
    pamed   = md("**PA-level median**")
  ) %>%
  cols_align(align = "center", columns = c(n_PAs, n_cells, bin, cell, pa, pamed)) %>%
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = 12,
    heading.align = "center",
    data_row.padding = px(3),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    column_labels.background.color = "#f2f2f2",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

gtsave(comparison_gt, file.path(dir_tables, "7_cell_vs_pa_comparison.html"))
print(comparison)
# ---------------------------------------------------------------------------
# 7. TABLE S8 - Indigenous Lands
# ---------------------------------------------------------------------------
var_labels_rdd <- c(
  pasture_area_perc                = "(1) Pasture (%)",
  deforestation_area_perc          = "(2) Deforestation (%)",
  stocking_rate_ha                 = "(3) Cattle density (all pasture)",
  stocking_rate_ha_deforested_cell = "(4) Cattle density (recently cleared)",
  cattle_on_established_pasture    = "(5) Cattle density (established pasture)"
)
il  <- filter(analysis_data, category == "Indigenous Lands")
spa <- filter(analysis_data, category == "Strictly Protected Areas")
stopifnot(nrow(il) > 0, nrow(spa) > 0)

s8_tab <- rdd_grid(list(`Indigenous Lands` = il),
                 c("pasture_area_perc", "deforestation_area_perc", vars_cattle),
                 donut_range_km)

s8 <- save_rdd(s8_tab, "8_rdd_IL_pasture_defo",
                rowname_col = "variable_label", groupname_col = NULL,
                keep_h_only = TRUE, drop = c("bandwidth_type", "subset"))
s8

# -----------------------------------------------------------------------------
# 8. TABLE S9 - Indigenous Lands, excluding buffers around settlements
#    Cattle outcomes only, at two exclusion radii.
# ----------------------------------------------------------------------------

excl <- c(`Excl. <1000 m of settlement` = 1000,
          `Excl. <2000 m of settlement` = 2000)

# exclude the datapoints outside of the PA for this
il_excl <- map(excl, function(thr) {
  keep <- il$within_pa == 0 |
    is.na(il$dist_to_next_aldeia) |
    il$dist_to_next_aldeia > thr
  message(sprintf("thr = %s m: dropped %d of %d inside cells (%.1f%%)",
                  thr,
                  sum(!keep), sum(il$within_pa == 1),
                  100 * sum(!keep) / sum(il$within_pa == 1)))
  il[keep, ]
})

s9_tab <- rdd_grid(il_excl, vars_cattle, donut_range_km,
                 subset_col = "exclusion")

s9 <- save_rdd(s9_tab, "9_rdd_IL_cattle_settlement",
                rowname_col = "exclusion", keep_h_only = TRUE,
                drop = "bandwidth_type")
s9

# -----------------------------------------------------------------------------
# 9. TABLE S10 - Strictly Protected Areas
#    Same specification as Table S8, different subset.
# -----------------------------------------------------------------------------

s10_tab <- rdd_grid(list(`Strictly Protected Areas` = spa),
                 c("pasture_area_perc", "deforestation_area_perc", vars_cattle),
                 donut_range_km)

s10 <- save_rdd(s10_tab, "10_rdd_SPA_pasture_defo_cattle",
                rowname_col = "variable_label", groupname_col = NULL,
                keep_h_only = TRUE, drop = c("bandwidth_type", "subset"))
s10

