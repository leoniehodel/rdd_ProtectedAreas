# ===========================================================================
# RDD FIGURES
# Structure:
#   1. Define vars to test
#.  2. Figure 2: Overall RDD plot
#   3. Figure 3: Pasture Cohorts RDD plot
#   4. Figure 4: IL vs SPA RDD plot (+ Figure S6: SP_IL_comparison_full.jpg")
#   5. Figure S4: Donut range estimation plot
#   6. Figure S5: Pasture Cohorts histogram
#   Leonie Hodel
# ===========================================================================

library(dplyr)
library(tidyr)
library(readr)
library(arrow)
library(purrr)
library(MatchIt)
library(ggplot2)
library(patchwork)
library(stringr)
library(gt)
library(rdrobust)

# ---------------------------------------------------------------------------
# 1. Define vars to test
# ---------------------------------------------------------------------------
##### Load analysis data ----------------------------------------------------
analysis_data <- read_parquet("data/rdd_analysis_data.parquet")

##### Define Output directories----------------------------------------------

dir_fig_rdd <- "results/figures/rdd_figures"  # main RDD / cattle figures
dir_fig_si  <- "results/figures/SI"           # supplementary figures
dir_tables  <- "results/tables"               # summary tables + raw model csv
# create results folder
invisible(lapply(c(dir_fig_rdd, dir_fig_si, dir_tables),
                 function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# how the variables were originally defined

# forest_threshold <- grid_size / 10
#stocking_rate_km2 = n_cattle / grid_size
#stocking_rate_ha  = stocking_rate_km2 / 100
#stocking_rate_ha_deforested_cell <- ifelse(
#   deforestation_area_perc > 0,
#   stocking_rate_ha,
#   NA
# )
#deforestation_area_perc = ifelse(
#         forest_area_km2 > forest_threshold,
#         pmin(100 * deforestation_area_km2 / forest_area_km2, 100),
#         NA
#       )

#analysis_data$cattle_on_established_pasture <- ifelse(
#  analysis_data$pasture_area_km2 > 0 &
#    analysis_data$deforestation_area_km2 <= 0,
#  analysis_data$stocking_rate_ha, NA_real_)

#analysis_data$distance_external_border_km <-analysis_data$distance_to_external_border/1000
#analysis_data$signed_distance_external_border_km<- analysis_data$signed_distance_external_border/1000

vars_to_test <- c(
  "deforestation_area_perc",
  "pasture_area_perc" ,
  "stocking_rate_ha" ,
  "stocking_rate_ha_deforested_cell"
)

# ---------------------------------------------------------------------------
# 2. Figure 2: Overall RDD plot
# ---------------------------------------------------------------------------

# general width of the bins
bin_width <- 120

binned_overall <- analysis_data %>%
  dplyr::mutate(distance_bin = bin_width * round(signed_distance_external_border / bin_width)) %>%
  dplyr::group_by(distance_bin) %>%
  dplyr::summarise(
    perc_pasture = mean(pasture_area_perc, na.rm = TRUE),
    perc_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    mean_cattle_occurence = mean(stocking_rate_ha, na.rm = TRUE),
    mean_cattle_deforestation = mean(stocking_rate_ha_deforested_cell, na.rm = TRUE),
    .groups = "drop"
  )

# labes and colors for the different variables
overall_config <- list(
  list(var = "perc_pasture",              label = "Pasture (%)",                            color = "yellowgreen"),
  list(var = "perc_deforestation",        label = "Deforestation (%)",                      color = "orangered"),
  list(var = "mean_cattle_occurence",     label = "Cattle density\n(animals/ha)",            color = "blue"),
  list(var = "mean_cattle_deforestation", label = "Cattle on deforested\narea (animals/ha)", color = "#9467BD")
)

make_overall_col <- function(cfg, df, show_x = TRUE, show_arrows = FALSE) {
  p <- ggplot(df, aes(x = distance_bin, y = .data[[cfg$var]])) +
    geom_point(size = 0.15, color = "black", alpha = 0.5) +
    geom_smooth(data = df %>% filter(distance_bin <= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2, linewidth = 0.5) +
    geom_smooth(data = df %>% filter(distance_bin >= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    labs(y = cfg$label) +
    theme_minimal(base_size = 7) +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      axis.title.y = element_text(size = 6),
      axis.text = element_text(size = 5),
      plot.margin = margin(2, 2, 2, 2)
    )
  if (show_arrows) {
    y_top <- max(df[[cfg$var]], na.rm = TRUE)
    arrow_y <- y_top * 0.55; label_y <- y_top * 0.45
    p <- p +
      annotate("segment", x = 3000, xend = 7000, y = arrow_y, yend = arrow_y,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.3, color = "grey30") +
      annotate("text", x = 5000, y = label_y, label = "Inside PA", size = 1.8, color = "grey30") +
      annotate("segment", x = -3000, xend = -7000, y = arrow_y, yend = arrow_y,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.3, color = "grey30") +
      annotate("text", x = -5000, y = label_y, label = "Outside PA", size = 1.8, color = "grey30")
  }
  if (show_x) {
    p <- p + labs(x = "Distance to boundary (m)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
  }
  p
}

overall_panels <- lapply(seq_along(overall_config), function(i) {
  make_overall_col(overall_config[[i]], binned_overall, show_arrows = (i == 1))
})

combined_overall <- overall_panels[[1]] | overall_panels[[2]] | overall_panels[[3]] | overall_panels[[4]]
combined_overall <- combined_overall +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 8, face = "bold")))
combined_overall
ggsave(file.path(dir_fig_rdd, "Figure2_overall_combined.png"),
       combined_overall, width = 170, height = 50, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------
# 4. Figure 3: Pasture Cohort RDD plot
# ---------------------------------------------------------------------------

# print out the number of PAs in each cohort
analysis_data %>%
  dplyr::group_by(pasture_group) %>%
  summarise(
    mean_pasture_in_group = mean(pasture_area_perc, na.rm = TRUE),
    n_PAs = n_distinct(pa_id)
  )
pasture_labels <- unique(analysis_data$pasture_group)
binned <- analysis_data %>%
  mutate(
    distance_bin = bin_width * round(signed_distance_external_border / bin_width),
    pasture_group = factor(pasture_group)) %>%
  group_by(pasture_group, distance_bin) %>%
  summarise(
    perc_pasture = mean(pasture_area_perc, na.rm = TRUE),
    perc_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    mean_cattle_occurence = mean(stocking_rate_ha, na.rm = TRUE),
    mean_cattle_deforestation = mean(stocking_rate_ha_deforested_cell, na.rm = TRUE),
    .groups = "drop"
  )

var_config <- list(
  list(var = "perc_pasture",              label = "Pasture (%)",                          color = "yellowgreen", ylim = NULL),
  list(var = "perc_deforestation",        label = "Deforestation (%)",                    color = "orangered",   ylim = NULL),
  list(var = "mean_cattle_occurence",     label = "Cattle density (animals/ha)",           color = "blue",        ylim = NULL),
  list(var = "mean_cattle_deforestation", label = "Cattle on deforested area (animals/ha)", color = "#9467BD",     ylim = c(0, 0.3))
)
make_row <- function(cfg, show_x = FALSE, show_strip = FALSE) {
  p <- ggplot(binned, aes(x = distance_bin, y = .data[[cfg$var]])) +
    geom_point(size = 0.3, color = "black", alpha = 0.7) +
    geom_smooth(data = binned %>% filter(distance_bin <= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2) +
    geom_smooth(data = binned %>% filter(distance_bin >= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~pasture_group, ncol = 6) +
    labs(y = cfg$label) +
    theme_minimal(base_size = 10) +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
  if (show_strip) {
    p <- p + theme(strip.text = element_text(size = 7, lineheight = 0.9))
    # annotate only the first facet panel
    first_level <- levels(binned$pasture_group)[1]
    y_top <- max(binned[[cfg$var]], na.rm = TRUE)
    arrow_y <- y_top * 0.85
    label_y <- y_top * 0.75
    p <- p + geom_segment(
      data = data.frame(pasture_group = factor(first_level, levels = levels(binned$pasture_group)),
                         x = 3000, xend = 7000, y = arrow_y, yend = arrow_y),
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
      linewidth = 0.3, color = "grey30", inherit.aes = FALSE
    ) +
    geom_text(
      data = data.frame(pasture_group = factor(first_level, levels = levels(binned$pasture_group)),
                         x = 5000, y = label_y),
      aes(x = x, y = y, label = "Inside PA"),
      size = 1.8, color = "grey30", inherit.aes = FALSE
    ) +
    geom_segment(
      data = data.frame(pasture_group = factor(first_level, levels = levels(binned$pasture_group)),
                         x = -3000, xend = -7000, y = arrow_y, yend = arrow_y),
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
      linewidth = 0.3, color = "grey30", inherit.aes = FALSE
    ) +
    geom_text(
      data = data.frame(pasture_group = factor(first_level, levels = levels(binned$pasture_group)),
                         x = -5000, y = label_y),
      aes(x = x, y = y, label = "Outside PA"),
      size = 1.8, color = "grey30", inherit.aes = FALSE
    )
  } else {
    p <- p + theme(strip.text = element_blank())
  }
  if (!is.null(cfg$ylim)) p <- p + coord_cartesian(ylim = cfg$ylim)
  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  } else {
    p <- p + labs(x = "Distance to boundary (m)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  p
}

# generate and save the plot
panels <- lapply(seq_along(var_config), function(i) {
  make_row(var_config[[i]], show_x = (i == length(var_config)), show_strip = (i == 1))
})

combined <- panels[[1]] / panels[[2]] / panels[[3]] / panels[[4]] + plot_annotation(tag_levels = "A")

combined
# the 22 rows removed are from panel 4, some NAs bc some bins dont have data.
ggsave(file.path(dir_fig_rdd, "Figure3_pasture_groups_combined.png"), combined, width = 12, height = 10, dpi = 300)

# ---------------------------------------------------------------------------
# 5. Figure 4: Matching and IL vs SPA RDD plot
# ---------------------------------------------------------------------------

##### Match Indigenous Lands and Strictly Protected Areas ------------------------


# definition of forest bin width around PAs
bin_width_f =2500
binned_policy_match <- analysis_data %>%
  mutate(distance_bin = bin_width_f * round(signed_distance_external_border / bin_width_f)) %>%
  group_by(category, distance_bin, pa_id) %>%
  summarise(
    mean_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    mean_forest = mean(forest_area_km2, na.rm = TRUE),
    mean_pasture = mean(pasture_area_perc, na.rm = TRUE),
    perc_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    perc_pasture = mean(pasture_area_perc, na.rm = TRUE),
    mean_cattle_occurence = mean(stocking_rate_ha, na.rm = TRUE),
    # cattle number on
    mean_cattle_deforestation = mean(stocking_rate_ha_deforested_cell, na.rm = TRUE),
    n = n(),
    n_cattle_occurence = sum(!is.na(stocking_rate_ha)),
    n_cattle_deforestation = sum(!is.na(stocking_rate_ha_deforested_cell)),
    .groups = "drop"
  )

pasture_summary_pas <- binned_policy_match %>%
  group_by(pa_id, category,distance_bin) %>% filter(distance_bin<0) %>%
  summarise(mean_deforestation = mean(perc_deforestation, na.rm = TRUE),
            mean_pasture = mean(perc_pasture, na.rm = TRUE),
            mean_forest = mean(mean_forest, na.rm = TRUE) ) %>%
  ungroup()  %>%  filter(!is.na(category),!is.na(mean_pasture)) %>%
  mutate(category = factor(category), mean_pasture= mean_pasture)

wide_df <- pasture_summary_pas %>%
  pivot_wider(
    names_from = distance_bin,
    values_from = c(mean_pasture, mean_deforestation, mean_forest)
  ) %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

names(wide_df) <- make.names(names(wide_df))
clean_df <- wide_df[!is.na(wide_df$category) & !is.na(wide_df$mean_forest_.2500), ]

m.out1 <- matchit(category ~ mean_forest_.2500, data = clean_df, method = "nearest", distance = "glm")
summary(m.out1, un = FALSE)
m.data <- match_data(m.out1)

# select the data that includes the matched samples
analysis_data_matched <- analysis_data %>% filter(pa_id %in% m.data$pa_id)

binned_policy <- analysis_data_matched %>%
  mutate(distance_bin = bin_width * round(signed_distance_external_border / bin_width)) %>%
  group_by(category, distance_bin) %>%
  summarise(
    mean_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    mean_pasture = mean(pasture_area_perc, na.rm = TRUE),
    perc_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    perc_pasture = mean(pasture_area_perc, na.rm = TRUE),
    mean_cattle_occurence = mean(stocking_rate_ha, na.rm = TRUE),
    # cattle number on
    mean_cattle_deforestation = mean(stocking_rate_ha_deforested_cell, na.rm = TRUE),
    n = n(),
    n_cattle_occurence = sum(!is.na(stocking_rate_ha)),
    n_cattle_deforestation = sum(!is.na(stocking_rate_ha_deforested_cell)),

    .groups = "drop"
  )

##### Indigenous Lands vs Strictly Protected Areas RDD plots  ------------------------

il_sp_config <- list(
  list(var = "mean_pasture",              label = "Pasture (%)",                            color = "yellowgreen"),
  list(var = "mean_deforestation",        label = "Deforestation (%)",                      color = "orangered"),
  list(var = "mean_cattle_occurence",     label = "Cattle density\n(animals/ha)",            color = "blue"),
  list(var = "mean_cattle_deforestation", label = "Cattle on deforested\narea (animals/ha)", color = "#9467BD")
)

# function to generate the Figure 4 annotated subplots
make_il_sp_row <- function(cfg, df, show_x = FALSE, show_strip = FALSE) {
  label_data <- df %>%
    group_by(category) %>%
    summarise(
      global_max_y = max(df[[cfg$var]], na.rm = TRUE),
      min_x = min(distance_bin), max_x = max(distance_bin),
      mean_left = mean(.data[[cfg$var]][distance_bin < 0], na.rm = TRUE),
      mean_right = mean(.data[[cfg$var]][distance_bin > 0], na.rm = TRUE),
      sd_left = sd(.data[[cfg$var]][distance_bin < 0], na.rm = TRUE),
      sd_right = sd(.data[[cfg$var]][distance_bin > 0], na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(df, aes(x = distance_bin, y = .data[[cfg$var]])) +
    geom_point(size = 0.15, color = "black", alpha = 0.5) +
    geom_smooth(data = df %>% filter(distance_bin <= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2, linewidth = 0.5) +
    geom_smooth(data = df %>% filter(distance_bin >= 0),
                method = "lm", se = TRUE, color = cfg$color, fill = cfg$color, alpha = 0.2, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_text(data = label_data,
              aes(x = min_x + 100, y = global_max_y * 1.05,
                  label = paste0(round(mean_left, 2), " ± ", round(sd_left, 2))),
              hjust = 0, size = 1.6, fontface = "italic", inherit.aes = FALSE) +
    geom_text(data = label_data,
              aes(x = max_x - 100, y = global_max_y * 1.05,
                  label = paste0(round(mean_right, 2), " ± ", round(sd_right, 2))),
              hjust = 1, size = 1.6, fontface = "italic", inherit.aes = FALSE) +
    expand_limits(y = max(df[[cfg$var]], na.rm = TRUE) * 1.15) +
    facet_wrap(~category, ncol = 2) +
    labs(y = cfg$label) +
    theme_minimal(base_size = 7) +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      axis.title.y = element_text(size = 6),
      axis.text = element_text(size = 5),
      plot.margin = margin(1, 2, 1, 2)
    )
  if (show_strip) {
    p <- p + theme(strip.text = element_text(size = 7, face = "bold"))
    first_cat <- levels(factor(df$category))[1]
    y_top <- max(df[[cfg$var]], na.rm = TRUE)
    arrow_y <- y_top * 0.85; label_y <- y_top * 0.75
    p <- p +
      geom_segment(data = data.frame(category = first_cat, x = 3000, xend = 7000, y = arrow_y, yend = arrow_y),
        aes(x=x,xend=xend,y=y,yend=yend), arrow=arrow(length=unit(1.5,"mm"),type="closed"),
        linewidth=0.3, color="grey30", inherit.aes=FALSE) +
      geom_text(data = data.frame(category = first_cat, x = 5000, y = label_y),
        aes(x=x,y=y,label="Inside PA"), size=1.8, color="grey30", inherit.aes=FALSE) +
      geom_segment(data = data.frame(category = first_cat, x = -3000, xend = -7000, y = arrow_y, yend = arrow_y),
        aes(x=x,xend=xend,y=y,yend=yend), arrow=arrow(length=unit(1.5,"mm"),type="closed"),
        linewidth=0.3, color="grey30", inherit.aes=FALSE) +
      geom_text(data = data.frame(category = first_cat, x = -5000, y = label_y),
        aes(x=x,y=y,label="Outside PA"), size=1.8, color="grey30", inherit.aes=FALSE)
  } else {
    p <- p + theme(strip.text = element_blank())
  }
  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  } else {
    p <- p + labs(x = "Distance to boundary (m)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
  }
  p
}

il_sp_panels <- lapply(seq_along(il_sp_config), function(i) {
  make_il_sp_row(il_sp_config[[i]], binned_policy, show_x = (i == length(il_sp_config)), show_strip = (i == 1))
})

combined_figure <- il_sp_panels[[1]] / il_sp_panels[[2]] / il_sp_panels[[3]] / il_sp_panels[[4]] +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 8, face = "bold")))
combined_figure
ggsave(file.path(dir_fig_rdd, "Figure4_SP_IL_comparison_38each.jpg"),
       combined_figure, width = 85, height = 145, units = "mm", dpi = 600)

## IL SPA plots — FULL dataset (all PAs, not the matched 38-each sample)
binned_policy_full <- analysis_data %>%
  mutate(distance_bin = bin_width * round(signed_distance_external_border / bin_width)) %>%
  group_by(category, distance_bin) %>%
  summarise(
    mean_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    mean_pasture = mean(pasture_area_perc, na.rm = TRUE),
    perc_deforestation = mean(deforestation_area_perc, na.rm = TRUE),
    perc_pasture = mean(pasture_area_perc, na.rm = TRUE),
    mean_cattle_occurence = mean(stocking_rate_ha, na.rm = TRUE),
    mean_cattle_deforestation = mean(stocking_rate_ha_deforested_cell, na.rm = TRUE),
    n = n(),
    n_cattle_occurence = sum(!is.na(stocking_rate_ha)),
    n_cattle_deforestation = sum(!is.na(stocking_rate_ha_deforested_cell)),
    .groups = "drop"
  )

il_sp_panels_full <- lapply(seq_along(il_sp_config), function(i) {
  make_il_sp_row(il_sp_config[[i]], binned_policy_full,
                 show_x = (i == length(il_sp_config)), show_strip = (i == 1))
})

combined_figure_full <- il_sp_panels_full[[1]] / il_sp_panels_full[[2]] / il_sp_panels_full[[3]] / il_sp_panels_full[[4]] +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 8, face = "bold")))

ggsave(file.path(dir_fig_si, "6_SP_IL_comparison_full.jpg"),
       combined_figure_full, width = 85, height = 145, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------
# 3. Figure S4: Donut range estimation plot
# ---------------------------------------------------------------------------
donuts_m <- seq(0, 1800, by = 120)      # in meters
# in order to get estimates for each donut seq it is necessary not to include the
# shift (dist_adj) to the  dataset (otherwise NA ).
run_donut <- function(gap_m, dat, yvar, xkm = "signed_distance_external_border_km",
                      one_sided_negative = FALSE) {
  to_km <- function(m) m / 1000
  d <- dat
  if (one_sided_negative) {
    d <- d %>% filter(!(get(xkm) >= -to_km(gap_m) & get(xkm) < 0)) %>%
      mutate(dist_adj = if_else(get(xkm) < 0,
                                get(xkm) + to_km(gap_m),
                                get(xkm)))
  } else {
    d <- d %>% filter(abs(get(xkm)) > to_km(gap_m))
  }
  # skip if too small or no variance
  if (nrow(d) < 200 || var(d[[xkm]], na.rm = TRUE) == 0 || var(d[[yvar]], na.rm = TRUE) == 0) {
    return(tibble(gap_m, est = NA, se = NA, p = NA, N = nrow(d)))
  }

  # safe rdrobust
  fit_safe <- purrr::safely(rdrobust, otherwise = NULL)
  res <- fit_safe(y = d[[yvar]], x = d[[xkm]],
                  #c = -(to_km(gap_m)),
                  p = 1,
                  h= 2, #to_km(gap_m),
                  #bwselect = "mserd",
                  #kernel = 'epanechnikov',
                  #cluster = cluster,
                  #bwcheck= 50000,
                  vce = "hc0")

  if (is.null(res$result)) {
    return(tibble(gap_m, est = NA, se = NA, p = NA, N = nrow(d)))
  }

  fit <- res$result
  tibble(
    gap_m = gap_m,
    est   = fit$Estimate[1, "tau.us"],
    se    = fit$Estimate[1, "se.us"],
    p     = 2 * pnorm(abs(est / se), lower.tail = FALSE),
    N     = sum(fit$Nh)
  )
}

donut_all <- purrr::map_dfr(vars_to_test, function(var) {
  res_one <- map_dfr(donuts_m, run_donut,
                     dat = analysis_data,
                     yvar = var,
                     one_sided_negative = TRUE) %>%
    mutate(one_sided = TRUE, variable = var)
  bind_rows( res_one)
})
donut_all <- donut_all %>%
  filter(!is.na(est)) %>%
  mutate(
    ci_lo = est - 1.96 * se,
    ci_hi = est + 1.96 * se,
    side_label = ifelse(one_sided, "Outside-only", "Symmetric"),
    variable_label = recode(variable,
                            "deforestation_area_perc" = "B. Deforestation (%)",
                            "pasture_area_perc" = "A. Pasture (%)",
                            "stocking_rate_ha" = "C. Cattle density (animals/pasture ha)",
                            "stocking_rate_ha_deforested_cell" = "D. Cattle density \n (animals/recently cleared pasture ha)"))
quad_fits <- donut_all %>%
  filter(gap_m < 1500) %>%
  group_by(variable_label) %>%
  group_modify(~{
    fit <- lm(est ~ poly(gap_m, 2, raw = TRUE), data = .x)
    coefs <- coef(fit)   # (Intercept), gap_m, gap_m^2

    # compute vertex of quadratic
    b1 <- coefs[2]
    b2 <- coefs[3]
    x_vertex <- -b1 / (2 * b2)

    tibble(xintercept = x_vertex)
  })
quad_fits
crit_points <- donut_all %>%filter(gap_m < 1500) %>%
  group_by(variable_label) %>%
  group_modify(~{
    fit <- lm(est ~ poly(gap_m, 3, raw = TRUE), data = .x)
    coefs <- coef(fit)   # b0, b1, b2, b3

    b1 <- coefs[2]
    b2 <- coefs[3]
    b3 <- coefs[4]

    # discriminant
    D <- 4*b2^2 - 12*b3*b1

    if (D < 0) {
      # no real critical points
      return(tibble(xintercept = NA_real_))
    }

    # two critical points
    x1 <- (-2*b2 + sqrt(D)) / (6*b3)
    x2 <- (-2*b2 - sqrt(D)) / (6*b3)

    tibble(xintercept = c(x1))
  })
# the cubic fit for the pasture area outputs an NA, so we take the quadratic fit.
crit_points$xintercept[crit_points$variable_label=='A. Pasture (%)']<-quad_fits$xintercept[quad_fits$variable_label=='A. Pasture (%)']
vline_df <- data.frame(
  variable_label = unique(quad_fits$variable_label),
  xintercept =  crit_points$xintercept  #c(1080, 1080, 240, 240)   # f
)
ggplot(donut_all %>% filter( gap_m<1600, side_label=='Outside-only'), aes(x = gap_m , y = est, color = side_label, linetype = side_label)) +
  geom_vline(data = vline_df,
             aes(xintercept = xintercept),
             color = "black", linetype = "dashed") +
  # quadratic fit
  geom_smooth(method = "lm", color='grey45',
              formula = y ~ poly(x, 3, raw = TRUE),
              se = FALSE,
              linewidth = 0.7) +
  geom_line() +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = side_label), alpha = 0.2, color = NA) +
  facet_wrap(~ variable_label, scales = "free_y") +
  scale_color_manual(values = c("Outside-only" = "#e31a1c", "Symmetric" = "#1f78b4")) +
  scale_fill_manual(values = c("Outside-only" = "#fb9a99", "Symmetric" = "#a6cee3")) +
  labs(
    x = "Donut exclusion (m) (minmum of cubic fit)",
    y = "RDD estimate (τ)",
    color = "Donut type",
    fill = "Donut type",
    linetype = "Donut type"
  ) +
  #xlim(c(0,1000))+
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(dir_fig_si, '4_donut_estimation_cubic_fit.png'), width = 7, height = 5.5, units = 'in')

# ---------------------------------------------------------------------------
# 6. Figure S5: pasture cohorts histogram
# ---------------------------------------------------------------------------
pasture_summary_pas <- analysis_data %>%
  filter(within_pa == 0) %>%
  group_by(pa_id) %>%
  summarise(mean_pasture = mean(pasture_area_perc, na.rm = TRUE))

pasture_summary <- pasture_summary_pas %>%
  mutate(pasture_group = cut(
    mean_pasture,
    breaks = quantile(mean_pasture,
                      probs = seq(0, 1, length.out = 7), # 6 groups = 7 breakpoints
                      na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("I: No to very low", "II: Very low", "III: Low", "IV: Medium", "V: High", "VI: Highest")
  ))

si_pasturegroups <- ggplot(pasture_summary, aes(x = mean_pasture, fill = pasture_group)) +
  geom_histogram(bins = 30, color = "white") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    #title = "Protected Areas grouped by Pasture Area in Buffer",
    x = "Mean pasture area (%)",
    y = "Count of PAs",
    fill = "Pasture %"
  ) +
  theme_minimal()

ggsave(file.path(dir_fig_si, "5_pasture_groups.jpg"), plot = si_pasturegroups,
       width = 170, height = 100, units = "mm", dpi = 300)


