# Standalone regulatory activity / MLS regression plots.
# Dependencies: ggplot2; jsonlite only for loading data from the genome browser.
# Source this file; it does not fetch data or change the website on its own.
#
# d <- load_regulatory_activity("SIRT6", window_kb = 500, measurement = "count")
# result <- plot_regulatory_mls(d, "SIRT6_regulatory_MLS.pdf")
# plot_regulatory_mls(d, "SIRT6_mass_corrected.png", body_mass_correction = TRUE)
# plot_regulatory_mls(d, "SIRT6_raw.png", log_activity = FALSE, log_mls = FALSE)
# tests <- read.delim("peaks500kb_without_body_mass_all_tested.tsv.gz")
# plot_regulatory_mls(d, "SIRT6_FDR.png", fdr_results = tests)
# mass_tests <- read.delim("peaks500kb_with_body_mass_all_tested.tsv.gz")
# plot_regulatory_mls(d, "SIRT6_both_FDR.png", fdr_results = tests,
#                     body_mass_fdr_results = mass_tests)
#
# Input: one row per species; columns species, activity (raw species median),
# mls (years), and body_mass_g if correcting for mass. The browser loader also
# supplies available_fraction, region_eligible, and mls_lower_bound. If present,
# these fields are used to exclude unavailable regions / lower-bound lifespans.
# Default exclusions: human, macaque and bowhead whale.
# Pass exclude_species = character() to retain all three.
# No positive-ATAC-coverage filter is used. This is OLS, not PGLS.
# Log transforms have no pseudocount. Model: activity_axis ~ MLS_axis, optionally
# + log10(body_mass_g). For adjusted plots the line holds log mass at its mean;
# points are observed values on the selected axes, not residuals.
# All points are black, the fitted line is dashed red, and species labels are omitted.
# Default export dimensions are 6 x 6 inches (1:1 canvas).
# Exports any ggplot2-supported format, including PDF and PNG; returns the plot,
# fitted lm, included data, excluded data with reasons, and statistics.
# fdr_results must contain ALL tested genes (both signs, including nonsignificant
# genes) from the same window, measurement, scales, species policy and model.
# Columns gene and p_value are required; BH is computed with stats::p.adjust.
# The current gene's nominal P must match its regression here. Optional model
# columns are also checked. Never provide just significant genes for BH.
# No FDR is invented for a single-gene plot without a testing-family table.
# Gene title is inferred from data$gene, or supplied with gene = "SIRT6".
# On an uncorrected plot, body_mass_fdr_results adds the corrected model's BH P.
# It fits and validates that model separately; the red line remains uncorrected.

regulatory_bh_fdr <- function(results) {
  if (!is.data.frame(results) || !all(c("gene", "p_value") %in% names(results))) {
    stop("fdr_results must contain gene and p_value for all tested genes in one model/window.")
  }
  results$gene <- toupper(trimws(as.character(results$gene)))
  if (anyNA(results$gene) || any(!nzchar(results$gene)) || anyDuplicated(results$gene)) {
    stop("FDR input must contain one unique gene per test; do not combine windows or models.")
  }
  if (!is.numeric(results$p_value) || any(!is.na(results$p_value) &
      (!is.finite(results$p_value) | results$p_value < 0 | results$p_value > 1))) {
    stop("FDR P values must be numeric in [0, 1], or NA for untestable genes.")
  }
  valid <- !is.na(results$p_value)
  if (sum(valid) < 2L) stop("BH across genes needs the full testing family, not a single gene.")
  results$BH_FDR <- NA_real_
  results$BH_FDR[valid] <- stats::p.adjust(results$p_value[valid], method = "BH")
  attr(results, "n_tests") <- sum(valid)
  results
}

load_regulatory_activity <- function(gene, window_kb = 500,
                                    measurement = c("count", "signal"),
                                    api_url = "http://127.0.0.1:8765") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite first.")
  measurement <- match.arg(measurement)
  if (length(window_kb) != 1L || !window_kb %in% c(20, 50, 100, 200, 500)) {
    stop("window_kb must be 20, 50, 100, 200 or 500 (on each side of the TSS).")
  }
  if (length(gene) != 1L || is.na(gene) || !nzchar(trimws(gene))) stop("Supply one gene symbol.")
  base <- sub("/+$", "", api_url)
  strategy <- paste0(if (measurement == "count") "peaks" else "peak_signal", window_kb, "kb")
  payload <- jsonlite::fromJSON(paste0(base, "/api/gene-activity?gene=",
    utils::URLencode(toupper(trimws(gene)), reserved = TRUE), "&strategy=", strategy),
    simplifyVector = FALSE)
  if (!is.null(payload$error)) stop(payload$error)
  if (!length(payload$records)) stop("No gene measurements returned.")
  catalog <- jsonlite::fromJSON(paste0(base, "/api/catalog"), simplifyVector = FALSE)
  catalog <- stats::setNames(catalog, vapply(catalog, `[[`, character(1), "id"))
  number <- function(x) if (is.null(x) || length(x) != 1L) NA_real_ else as.numeric(x)
  rows <- lapply(payload$records, function(r) {
    sp <- catalog[[r$species]]
    if (is.null(sp)) stop("Species missing from browser catalog: ", r$species)
    life <- payload$lifespans[[sp$latin]]
    data.frame(gene = payload$gene, strategy = strategy, species = r$species,
      common_name = sp$name, scientific_name = sp$latin, locus = r$locus,
      activity = number(r$median), mls = number(life$years),
      body_mass_g = number(life$body_mass_g),
      mls_lower_bound = number(life$lower_bound_years),
      available_fraction = number(r$available_fraction),
      region_eligible = isTRUE(r$region_eligible),
      activity_label = if (measurement == "count") "Median peak count" else "Median signal in peaks (CPM x bp)",
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  # Select the highest-activity locus BEFORE applying region / model filters,
  # matching the browser; do not fall back to a lower-activity eligible locus.
  score <- ifelse(is.finite(d$activity), d$activity, -Inf)
  d <- d[order(d$species, -score, d$locus, method = "radix"), , drop = FALSE]
  d <- d[!duplicated(d$species), , drop = FALSE]
  rownames(d) <- NULL
  d
}

plot_regulatory_mls <- function(data, filename = NULL,
                                log_activity = TRUE, log_mls = TRUE,
                                body_mass_correction = FALSE,
                                exclude_species = c("Human", "Macaque", "Balaena_mysticetus",
                                  "Homo sapiens", "Homo_sapiens", "Macaca mulatta",
                                  "Macaca_mulatta", "Balaena mysticetus"),
                                activity_label = NULL, show_statistics = TRUE,
                                width = 6, height = 6, dpi = 300,
                                gene = NULL, fdr_results = NULL,
                                body_mass_fdr_results = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  required <- c("species", "activity", "mls", if (body_mass_correction) "body_mass_g")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("Input must be a data frame containing: ", paste(required, collapse = ", "))
  }
  if (!nrow(data) || anyNA(data$species) || any(!nzchar(as.character(data$species)))) {
    stop("Supply nonempty species identifiers.")
  }
  if (anyDuplicated(data$species)) stop("Use one species-median observation per species.")
  input_genes <- if ("gene" %in% names(data)) unique(toupper(trimws(as.character(data$gene)))) else character()
  if (length(input_genes) > 1L) stop("Plot one gene at a time.")
  if (is.null(gene) && length(input_genes)) gene <- input_genes
  if (!is.null(gene)) {
    if (length(gene) != 1L || is.na(gene) || !nzchar(trimws(gene))) stop("Supply one gene name.")
    gene <- toupper(trimws(gene))
    if (length(input_genes) && !identical(gene, input_genes)) stop("Gene title does not match data$gene.")
  }
  numeric_columns <- setdiff(required, "species")
  for (nm in intersect(c("available_fraction", "mls_lower_bound"), names(data))) {
    numeric_columns <- c(numeric_columns, nm)
  }
  if (!all(vapply(data[numeric_columns], is.numeric, logical(1)))) stop("Measurement columns must be numeric.")
  if (is.null(activity_label)) {
    labels <- if ("activity_label" %in% names(data)) unique(data$activity_label) else character()
    if (length(labels) > 1L) stop("Do not mix regulatory measurements in one plot.")
    activity_label <- if (length(labels) == 1L) labels else "Median regulatory activity"
  }
  reason <- rep("", nrow(data))
  flag <- function(bad, label) { reason[which(bad & reason == "")] <<- label }
  flag(data$species %in% exclude_species, "Excluded species")
  flag(!is.finite(data$activity), "Activity unavailable")
  if ("available_fraction" %in% names(data)) {
    flag(!is.finite(data$available_fraction) | data$available_fraction < 0.5,
         "Less than 50% region available, or availability unknown")
  }
  if ("region_eligible" %in% names(data)) {
    flag(is.na(data$region_eligible) | !data$region_eligible, "Region ineligible")
  }
  if ("mls_lower_bound" %in% names(data)) flag(is.finite(data$mls_lower_bound), "MLS is a lower bound")
  flag(!is.finite(data$mls) | data$mls <= 0, "MLS unavailable or nonpositive")
  if (log_activity) flag(data$activity <= 0, "Nonpositive activity: log10 undefined")
  if (body_mass_correction) {
    flag(!is.finite(data$body_mass_g) | data$body_mass_g <= 0, "Body mass unavailable or nonpositive")
  }
  excluded <- data[reason != "", , drop = FALSE]
  excluded$exclusion_reason <- reason[reason != ""]
  d <- data[reason == "", , drop = FALSE]
  d$.x <- if (log_mls) log10(d$mls) else d$mls
  d$.y <- if (log_activity) log10(d$activity) else d$activity
  if (nrow(d) < 3L + as.integer(body_mass_correction)) stop("Too few eligible species to fit and test this model.")
  if (length(unique(d$.x)) < 2L || length(unique(d$.y)) < 2L) stop("MLS and activity must vary.")
  if (body_mass_correction) d$.log_mass <- log10(d$body_mass_g)
  model <- stats::lm(if (body_mass_correction) .y ~ .x + .log_mass else .y ~ .x, data = d)
  if (model$rank < length(stats::coef(model))) stop("Predictors are collinear; cannot estimate the requested fit.")
  sm <- summary(model)
  cf <- sm$coefficients[".x", ]
  t_value <- unname(cf["t value"])
  r <- if (is.infinite(t_value)) sign(t_value) else sign(t_value) * sqrt(t_value^2 / (t_value^2 + model$df.residual))
  statistics <- data.frame(n = nrow(d), MLS_coefficient = unname(cf["Estimate"]),
    standard_error = unname(cf["Std. Error"]), p_value = unname(cf["Pr(>|t|)"]),
    r = unname(r), r_type = if (body_mass_correction) "Partial r, log10 mass adjusted" else "Pearson r",
    model_r_squared = sm$r.squared, df = model$df.residual,
    log10_activity = log_activity, log10_MLS = log_mls,
    body_mass_correction = body_mass_correction,
    BH_FDR = NA_real_, FDR_genes_tested = NA_integer_)
  fdr_table <- NULL
  if (!is.null(fdr_results)) {
    if (is.null(gene)) stop("Provide gene or data$gene to locate the gene in fdr_results.")
    fdr_table <- regulatory_bh_fdr(fdr_results)
    check_setting <- function(column, expected) {
      if (!column %in% names(fdr_table)) return(invisible(NULL))
      values <- unique(tolower(as.character(fdr_table[[column]])))
      if (length(values) != 1L || is.na(values) || values != tolower(as.character(expected))) {
        stop("FDR testing family has a different or mixed ", column, " setting.")
      }
    }
    check_setting("log10_activity", log_activity)
    check_setting("log10_MLS", log_mls)
    check_setting("body_mass_correction", body_mass_correction)
    if ("strategy" %in% names(data)) {
      strategy <- unique(data$strategy)
      if (length(strategy) != 1L) stop("Do not mix regulatory windows or measurements.")
      check_setting("strategy", strategy)
    }
    hit <- match(gene, fdr_table$gene)
    if (is.na(hit) || is.na(fdr_table$p_value[hit])) stop("Gene has no valid test in the FDR family.")
    reference_p <- fdr_table$p_value[hit]
    if (abs(reference_p - statistics$p_value) > 1e-6 * max(reference_p, statistics$p_value, 1e-300)) {
      stop("Gene's P value does not match the FDR family. Recompute all genes using the same data, filters and model.")
    }
    if ("n_species" %in% names(fdr_table) && fdr_table$n_species[hit] != nrow(d)) {
      stop("Gene's sample size does not match the FDR family.")
    }
    statistics$BH_FDR <- fdr_table$BH_FDR[hit]
    statistics$FDR_genes_tested <- attr(fdr_table, "n_tests")
  }
  body_mass_corrected <- NULL
  if (!is.null(body_mass_fdr_results)) {
    if (body_mass_correction) {
      stop("For a corrected fitted line, supply its testing family as fdr_results. body_mass_fdr_results adds a companion result to an uncorrected plot.")
    }
    # Reuse the same filtering, fitting and FDR validation for the companion
    # model. Missing mass can change its included species; report its n below.
    body_mass_corrected <- plot_regulatory_mls(data,
      log_activity = log_activity, log_mls = log_mls,
      body_mass_correction = TRUE, exclude_species = exclude_species,
      activity_label = activity_label, show_statistics = FALSE,
      gene = gene, fdr_results = body_mass_fdr_results)
  }
  prediction <- data.frame(.x = seq(min(d$.x), max(d$.x), length.out = 200))
  if (body_mass_correction) prediction$.log_mass <- mean(d$.log_mass)
  prediction$.y <- as.numeric(stats::predict(model, newdata = prediction))
  format_p <- function(value) format.pval(value, digits = 3, eps = 1e-300)
  p_text <- paste0("P = ", format_p(statistics$p_value))
  fdr_text <- if (is.finite(statistics$BH_FDR)) paste0("BH-adjusted P = ", format_p(statistics$BH_FDR),
    " (", format(statistics$FDR_genes_tested, big.mark = ",", trim = TRUE), " genes)") else "BH-adjusted P: testing family not supplied"
  subtitle <- if (show_statistics) paste0("n = ", nrow(d), "   ",
    if (body_mass_correction) "Partial r" else "Pearson r", " = ", sprintf("%.3f", r),
    "   ", p_text, "\n", fdr_text) else NULL
  if (show_statistics && !is.null(body_mass_corrected)) {
    mass_stats <- body_mass_corrected$statistics
    subtitle <- paste0("Uncorrected fit: n = ", nrow(d), "   r = ", sprintf("%.3f", r),
      "   ", p_text, "\n",
      "BH-adjusted P (uncorrected) = ",
      if (is.finite(statistics$BH_FDR)) format_p(statistics$BH_FDR) else "unavailable", "\n",
      "BH-adjusted P (body-mass corrected) = ", format_p(mass_stats$BH_FDR), "\n",
      if (isTRUE(statistics$FDR_genes_tested == mass_stats$FDR_genes_tested)) {
        paste0("BH across ", format(statistics$FDR_genes_tested, big.mark = ",", trim = TRUE), " genes per model")
      } else {
        paste0("BH tests: ", statistics$FDR_genes_tested, " uncorrected; ",
          mass_stats$FDR_genes_tested, " corrected")
      },
      if (mass_stats$n != nrow(d)) paste0("; corrected n = ", mass_stats$n) else "")
  } else if (show_statistics && body_mass_correction) {
    subtitle <- paste0("Body-mass-corrected fit\n", subtitle)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .x, y = .y)) +
    ggplot2::geom_line(data = prediction, colour = "red", linetype = "dashed", linewidth = 0.9) +
    ggplot2::geom_point(colour = "black", size = 2.5) +
    ggplot2::labs(x = if (log_mls) "log10(Maximum lifespan in years)" else "Maximum lifespan (years)",
      y = if (log_activity) paste0("log10(", activity_label, ")") else activity_label,
      title = gene, subtitle = subtitle) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E8E8E8", linewidth = 0.3),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(colour = "black", size = 11),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "none")
  if (!is.null(filename)) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename, p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  }
  invisible(list(plot = p, model = model, data = d, excluded = excluded,
                 statistics = statistics, fdr_results = fdr_table,
                 body_mass_corrected = body_mass_corrected))
}
