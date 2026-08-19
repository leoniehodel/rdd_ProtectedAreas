
# Which rdrobust quantities to report. Set once, applies everywhere, so the
# estimate
#   "conventional" -> tau.us / se.us / pv[1]
#   "robust"       -> tau.bc / se.rb / pv[3]
# paper uses 'robust'
RDD_INFERENCE <- "robust"

rdd_pick <- function(fit, what) {
  if (RDD_INFERENCE == "robust") {
    switch(what, est = fit$Estimate[1, "tau.bc"],
           se  = fit$Estimate[1, "se.rb"],
           pv  = fit$pv[3])
  } else {
    switch(what, est = fit$Estimate[1, "tau.us"],
           se  = fit$Estimate[1, "se.us"],
           pv  = fit$pv[1])
  }
}

##One rdrobust object -> one tidy row --------------------------------
tidy_rd <- function(fit) {
  if (is.null(fit)) return(NULL)
  bl <- fit$beta_Y_p_l          # length p+1, so slopes exist only for p >= 1
  br <- fit$beta_Y_p_r
  tibble(
    est   = rdd_pick(fit, "est"),
    se    = rdd_pick(fit, "se"),
    pval  = rdd_pick(fit, "pv"),
    intercept_left  = bl[1],
    intercept_right = br[1],
    slope_left      = if (length(bl) > 1) bl[2] else NA_real_,
    slope_right     = if (length(br) > 1) br[2] else NA_real_,
    bw_left  = fit$bws["h", "left"],
    bw_right = fit$bws["h", "right"],
    N_left   = fit$N_h[1],
    N_right  = fit$N_h[2]
  )
}

## Data + outcome -> all (polynomial x bandwidth) fits, in long form --------
# donut_km is the donut hwidth in KM. The gap is cut on the outside only
# and the running variable is shifted to close it.
rdd_run <- function(data, yvar, donut_km, cluster_var = "pa_id",
                    include_2h = TRUE, polys = c(1, 0),
                    min_n = 500, bwcheck = 50000) {

  if (!yvar %in% names(data))
    stop("rdd_run: column '", yvar, "' not found in data", call. = FALSE)

  d <- data %>%
    filter(!(signed_distance_external_border_km >= -donut_km &
               signed_distance_external_border_km <= 0)) %>%
    mutate(dist_adj = if_else(signed_distance_external_border_km < 0,
                              signed_distance_external_border_km + donut_km,
                              signed_distance_external_border_km))

  y <- d[[yvar]]; x <- d$dist_adj; cl <- d[[cluster_var]]
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]; cl <- cl[ok]

  if (length(y) < min_n || var(y) == 0 || var(x) == 0) {
    warning("rdd_run: insufficient data for ", yvar, call. = FALSE)
    return(NULL)
  }

  safe_rd <- purrr::safely(rdrobust, otherwise = NULL)

  map_dfr(polys, function(p) {
    f1 <- safe_rd(y = y, x = x, p = p, bwselect = "mserd",
                  cluster = cl, bwcheck = bwcheck, vce = "hc0")$result
    if (is.null(f1)) return(NULL)

    rows <- tidy_rd(f1) %>% mutate(p = p, bandwidth_type = "h", .before = 1)

    if (include_2h) {
      f2 <- safe_rd(y = y, x = x, p = p, h = 2 * f1$bws[1],
                    cluster = cl, bwcheck = bwcheck, vce = "hc0")$result
      if (!is.null(f2))
        rows <- bind_rows(rows,
                          tidy_rd(f2) %>% mutate(p = p, bandwidth_type = "2h", .before = 1))
    }
    rows
  }) %>%
    mutate(variable = yvar, donut_m = donut_km * 1000, .before = 1)
}

## A grid of subsets x outcomes ---------------------------------------
# subsets: named list of data frames. donuts: named numeric vector (KM).
rdd_grid <- function(subsets, vars, donuts, subset_col = "subset", ...) {
  imap_dfr(subsets, function(d, subset_label) {
    map_dfr(vars, function(v) {
      message("RDD: ", subset_label, " / ", v)
      dk <- if (v %in% names(donuts)) unname(donuts[[v]]) else 0
      if (is.na(dk)) dk <- 0
      message("Donut km: ", dk, " / ")
      res <- rdd_run(d, v, donut_km = dk, ...)
      if (is.null(res)) return(NULL)
      res[[subset_col]] <- subset_label
      res
    })
  })
}

##Long results -> display-ready wide table ---------------------------
rdd_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
}

# est_digits / comp_digits let each table set its own precision, instead of a
# blanket round() that leaves the pasted columns ragged.
rdd_format <- function(res, var_labels, est_digits = 2, comp_digits = 3) {
  f <- function(v, d) formatC(v, format = "f", digits = d)

  res %>%
    mutate(
      variable_label = factor(unname(var_labels[variable]),
                              levels = unname(var_labels)),
      tau = paste0(f(est, est_digits), rdd_stars(pval),
                   " (", f(se, est_digits), ")"),
      slopes = ifelse(is.na(slope_left), NA_character_,
                      paste0(f(slope_left, comp_digits), " / ",
                             f(slope_right, comp_digits))),
      intercepts = paste0(f(intercept_left, comp_digits), " / ",
                          f(intercept_right, comp_digits)),
      h  = f(bw_left, comp_digits),
      Ns = paste0(format(N_left, big.mark = ","), " / ",
                  format(N_right, big.mark = ","))
    ) %>%
    select(variable_label, bandwidth_type, p, tau, slopes, intercepts, h, Ns,
           any_of(c("subset", "pasture_type", "exclusion"))) %>%
    pivot_wider(names_from = p, values_from = c(tau, slopes, intercepts, h, Ns),
                names_glue = "{.value}_p{p}") %>%
    select(-any_of(c("slopes_p0", "intercepts_p1", "Ns_p0"))) %>%
    rename(slopes = slopes_p1, intercepts = intercepts_p0, Ns = Ns_p1) %>%
    mutate(bandwidth_type = factor(bandwidth_type, levels = c("h", "2h"))) %>%
    arrange(variable_label,
            across(any_of(c("subset", "pasture_type", "exclusion"))),
            bandwidth_type)
}

## One gt styler for every table --------------------------------------
rdd_gt <- function(df, rowname_col = "bandwidth_type",
                   groupname_col = "variable_label", title = NULL,
                   footnote = NULL) {

  g <- df %>%
    gt(rowname_col = rowname_col, groupname_col = groupname_col) %>%
    tab_options(
      table.width = pct(100),
      table.font.names = "Times New Roman",
      table.font.size = 12,
      heading.align = "center",
      data_row.padding = px(3),
      row_group.as_column = TRUE,
      row_group.font.weight = "bold",
      row_group.padding = px(8),
      column_labels.font.weight = "bold",
      column_labels.background.color = "#f2f2f2",
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      table.border.top.width = px(1),
      table.border.bottom.width = px(1)
    ) %>%
    cols_align(align = "center",
               columns = any_of(c("tau_p1", "tau_p0", "intercepts", "slopes",
                                  "h_p1", "h_p0", "Ns"))) %>%
    tab_spanner(label = "\u03c4", columns = any_of(c("tau_p1", "tau_p0"))) %>%
    tab_spanner(label = "Model Components",
                columns = any_of(c("intercepts", "slopes"))) %>%
    tab_spanner(label = "Bandwidths & Sample Size",
                columns = any_of(c("h_p0", "h_p1", "Ns"))) %>%
    cols_label(
      tau_p1 = md("**(a) Linear**"),
      tau_p0 = md("**(b) Constant**"),
      intercepts = md("**Intercepts (out/in)**"),
      slopes = md("**Slopes (out/in)**"),
      h_p0 = md("**h Constant**"),
      h_p1 = md("**h Linear**"),
      Ns = md("**N (out/in)**")
    ) %>%
    # striping keyed on row position: bandwidth_type is consumed by rowname_col
    tab_style(style = cell_fill(color = "#fafafa"),
              locations = cells_body(rows = seq(1, nrow(df), by = 2)))

  if (!is.null(title))    g <- g %>% tab_header(title = md(title))
  if (!is.null(footnote)) g <- g %>% tab_source_note(source_note = md(footnote))
  g
}


##Write raw CSV and formatted HTML together --------------------------
save_rdd <- function(res, stem, raw_file = paste0(stem, "_raw.csv"),
                     rowname_col = "bandwidth_type",
                     groupname_col = "variable_label",
                     keep_h_only = FALSE, drop = NULL) {

  write_csv(res, file.path(dir_tables, raw_file))

  tab <- rdd_format(res, var_labels_rdd)
  if (keep_h_only) tab <- filter(tab, bandwidth_type == "h")
  if (!is.null(drop)) tab <- select(tab, -any_of(drop))

  g <- rdd_gt(tab, rowname_col = rowname_col, groupname_col = groupname_col)
  gtsave(g, file.path(dir_tables, paste0(stem, ".html")))
  g
}

