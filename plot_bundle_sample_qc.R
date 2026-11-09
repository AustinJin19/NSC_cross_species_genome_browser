#!/usr/bin/env Rscript

# Horizontal sample-QC box plots and assembly-count bars with a shared y axis.
# Dependency: install.packages("ggplot2")
#
# Run from a terminal:
#   Rscript plot_bundle_sample_qc.R
#   Rscript plot_bundle_sample_qc.R path/to/qc.tsv path/to/output_prefix
#   Rscript plot_bundle_sample_qc.R path/to/qc.tsv path/to/output_prefix path/to/assembly_counts.tsv
#
# Outputs: <output_prefix>.pdf (vector) and <output_prefix>.png (300 dpi).
# Default input: bundle_sample_qc_summary.tsv beside this script.
# Default output prefix: outputs/sample_qc/bundle_sample_qc beside the input.
# Default assembly table: assembly_sequence_counts.tsv beside the QC input.
# Assembly counts are extracted from the supplied reference workbook's
# "Local file evidence" sheet, column "FASTA sequences". The table retains the
# workbook, sheet, row and FASTA path for each count. To refresh the extract:
#   python3 output/sample_qc/extract_assembly_counts.py
# This counts FASTA records (chromosomes/scaffolds/contigs), not contigs obtained
# by splitting scaffolds at gaps. Each assembly appears once, not per sample.
# The fourth panel uses log10(count) bar lengths and labels the original counts.
# In RStudio, set the working directory to the input folder, then source this
# script. The ggplot object remains available as qc_plot.
#
# Each row is one sample. Batches/runs are pooled within each species.
# Common-name labels are ordered by median alignment rate, highest first.
# Missing values ("-", "NA", or blank) are omitted only for that metric.
# No samples are filtered based on QC thresholds or lifespan-analysis settings.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 first: install.packages(\"ggplot2\")", call. = FALSE)
}
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L || any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript plot_bundle_sample_qc.R [input.tsv] [output_prefix] [assembly_counts.tsv]\n")
  quit(status = if (length(args) > 3L) 1L else 0L)
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
input_file <- if (length(args) >= 1L) args[1] else {
  file.path(script_dir, "bundle_sample_qc_summary.tsv")
}
input_file <- normalizePath(input_file, mustWork = TRUE)
output_prefix <- if (length(args) >= 2L) args[2] else {
  file.path(dirname(input_file), "outputs", "sample_qc", "bundle_sample_qc")
}
assembly_file <- if (length(args) >= 3L) args[3] else {
  file.path(dirname(input_file), "assembly_sequence_counts.tsv")
}

qc <- read.delim(
  input_file, header = TRUE, sep = "\t", quote = "", comment.char = "",
  na.strings = c("-", "NA", ""), stringsAsFactors = FALSE,
  check.names = FALSE, strip.white = TRUE
)
metric_columns <- c("overall_align_pct", "final_read_pairs", "tss_enrichment")
required <- c("species", "sample", metric_columns)
missing_columns <- setdiff(required, names(qc))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (!nrow(qc) || anyNA(qc$species) || anyNA(qc$sample) ||
    any(!nzchar(trimws(qc$species))) || any(!nzchar(trimws(qc$sample)))) {
  stop("The input must contain samples with nonempty species/sample names.",
       call. = FALSE)
}
for (column in metric_columns) {
  original <- qc[[column]]
  numeric_values <- suppressWarnings(as.numeric(original))
  invalid <- !is.na(original) &
    (is.na(numeric_values) | !is.finite(numeric_values) | numeric_values < 0)
  if (any(invalid)) {
    stop("Invalid numeric values in ", column, " for samples: ",
         paste(qc$sample[invalid], collapse = ", "), call. = FALSE)
  }
  qc[[column]] <- numeric_values
}
if (any(qc$overall_align_pct > 100, na.rm = TRUE)) {
  stop("overall_align_pct must use the 0-100 percentage scale.", call. = FALSE)
}
if (anyDuplicated(qc$sample)) {
  warning("Repeated sample IDs found. Each input row will be plotted separately.")
}

# Common names follow the genome browser's species labels.
# Bramus lutescens: https://www.mammaldiversity.org/taxon/1002134/
# Cricetomys gambianus: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10085
common_names <- c(
  Acomys_cahirinus = "African spiny mouse",
  Artibeus_jamaicensis = "Jamaican fruit bat",
  Balaena_mysticetus = "Bowhead whale",
  Bramus_lutescens = "Transcaucasian mole vole",
  Castor_canadensis = "North American beaver",
  Cavia_porcellus = "Guinea pig",
  Chinchilla_lanigera = "Chinchilla",
  Cricetomys_gambianus = "Gambian pouched rat",
  Cricetulus_griseus = "Chinese hamster",
  Eonycteris_spelaea = "Cave nectar bat",
  Eptesicus_fuscus = "Big brown bat",
  Erethizon_dorsatum = "North American porcupine",
  Fukomys_damarensis = "Damaraland mole-rat",
  Heterocephalus_glaber = "Naked mole-rat",
  Hydrochoerus_hydrochaeris = "Capybara",
  Marmota_monax = "Woodchuck",
  Meriones_unguiculatus = "Mongolian gerbil",
  Mesocricetus_auratus = "Golden hamster",
  Microtus_pennsylvanicus = "Meadow vole",
  Molossus_molossus = "Pallas's mastiff bat",
  Mus_musculus = "Mouse",
  Myocastor_coypus = "Nutria",
  Nannospalax_galili = "Blind mole-rat",
  Octodon_degus = "Degu",
  Ondatra_zibethicus = "Muskrat",
  Peromyscus_maniculatus = "Deer mouse",
  Rattus_norvegicus = "Rat",
  Scalopus_aquaticus = "Eastern mole",
  Sciurus_carolinensis = "Eastern gray squirrel",
  Sciurus_niger = "Fox squirrel",
  Sylvilagus_floridanus = "Eastern cottontail",
  Tamias_striatus = "Eastern chipmunk",
  Tamiasciurus_hudsonicus = "American red squirrel"
)
unmapped_species <- setdiff(unique(qc$species), names(common_names))
if (length(unmapped_species)) {
  stop("Add common names to common_names for: ",
       paste(unmapped_species, collapse = ", "), call. = FALSE)
}

# One shared top-to-bottom order: descending species median alignment rate.
# Break ties by common name; species without alignment measurements go last.
median_alignment <- tapply(qc$overall_align_pct, qc$species, median, na.rm = TRUE)
species_order <- names(median_alignment)[order(
  -median_alignment, common_names[names(median_alignment)], na.last = TRUE
)]
assembly <- read.delim(
  assembly_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("species", "sequence_count") %in% names(assembly)) ||
    anyDuplicated(assembly$species)) {
  stop("Assembly table must have one species/sequence_count row per species.",
       call. = FALSE)
}
assembly <- assembly[match(species_order, assembly$species), ]
counts <- suppressWarnings(as.numeric(assembly$sequence_count))
if (anyNA(assembly$species) || anyNA(counts) || any(!is.finite(counts)) ||
    any(counts < 1 | counts != floor(counts))) {
  stop("A positive integer assembly count is required for every QC species.",
       call. = FALSE)
}
assembly$sequence_count <- counts
metric_labels <- c(
  "Alignment (%)", "Final read pairs (millions)", "TSS enrichment",
  "Assembly sequences\n(log10 count)"
)
metric_colors <- setNames(c("#3C5488", "#00A087", "#E64B35", "#66788A"), metric_labels)
assembly$species <- factor(assembly$species, levels = rev(species_order))
assembly$metric <- factor(metric_labels[4], levels = metric_labels)
assembly$value <- log10(assembly$sequence_count)
assembly$count_label <- format(assembly$sequence_count, big.mark = ",",
                               scientific = FALSE, trim = TRUE)
# Reserve space for exact count labels only in the assembly facet.
assembly_padding <- assembly[1, ]
assembly_padding$value <- max(assembly$value) + 1.45

# Base-R reshape keeps the only package dependency as ggplot2.
long_qc <- do.call(rbind, lapply(seq_along(metric_columns), function(i) {
  values <- qc[[metric_columns[i]]]
  if (metric_columns[i] == "final_read_pairs") values <- values / 1e6
  data.frame(species = qc$species, sample = qc$sample,
             metric = metric_labels[i], value = values)
}))
long_qc$species <- factor(long_qc$species, levels = rev(species_order))
long_qc$metric <- factor(long_qc$metric, levels = metric_labels)
sample_counts <- table(qc$species)
species_labels <- setNames(
  paste0(common_names[species_order], "  (n=", sample_counts[species_order], ")"),
  species_order
)
observed <- long_qc[!is.na(long_qc$value), ]
if (!nrow(observed)) stop("No QC measurements are available.", call. = FALSE)

# Retain all species in every facet, including those with no TSS measurements.
availability <- aggregate(
  list(n = as.integer(!is.na(long_qc$value))),
  long_qc[c("species", "metric")], sum
)
unavailable <- availability[availability$n == 0L, ]

qc_plot <- ggplot(observed, aes(x = value, y = species)) +
  geom_col(
    data = assembly, aes(x = value, y = species),
    inherit.aes = FALSE, orientation = "y", width = 0.56, fill = "#66788A"
  ) +
  geom_text(
    data = assembly, aes(x = value, y = species, label = count_label),
    inherit.aes = FALSE, hjust = -0.15, size = 3, colour = "#303846"
  ) +
  geom_blank(
    data = assembly_padding, aes(x = value, y = species),
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = data.frame(
      metric = factor("TSS enrichment", levels = metric_labels),
      threshold = 4
    ),
    aes(xintercept = threshold),
    linetype = "dashed", colour = "#555555", linewidth = 0.55
  ) +
  geom_boxplot(
    aes(fill = metric, colour = metric), orientation = "y",
    width = 0.56, alpha = 0.20, linewidth = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    aes(colour = metric),
    position = position_jitter(width = 0, height = 0.10, seed = 42),
    size = 1.75, alpha = 0.90
  ) +
  geom_text(
    data = unavailable, aes(x = -Inf, y = species, label = "NA"),
    inherit.aes = FALSE, hjust = -0.3, colour = "#777777", size = 3
  ) +
  # free_x allows metric-specific numeric scales; the species axis stays fixed.
  facet_grid(cols = vars(metric), scales = "free_x", drop = FALSE) +
  scale_y_discrete(labels = species_labels, drop = FALSE) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  scale_fill_manual(values = metric_colors, guide = "none") +
  scale_colour_manual(values = metric_colors, guide = "none") +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL, caption = NULL) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 11, margin = margin(9, 5, 9, 5)),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    axis.text.y = element_text(size = 10, colour = "#303846"),
    axis.text.x = element_text(size = 10, colour = "#303846"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E8ECF0", linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "#E8ECF0", linewidth = 0.3),
    panel.spacing.x = grid::unit(1.1, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
plot_height <- max(5, 0.30 * length(species_order) + 0.9)
ggsave(paste0(output_prefix, ".pdf"), qc_plot, device = grDevices::pdf,
       width = 17.5, height = plot_height, units = "in", bg = "white",
       useDingbats = FALSE)
ggsave(paste0(output_prefix, ".png"), qc_plot, device = "png",
       width = 17.5, height = plot_height, units = "in", dpi = 300, bg = "white")
if (interactive()) print(qc_plot)
message("Plotted ", nrow(qc), " samples across ", length(species_order), " species.")
message("Assembly counts: ", nrow(assembly), " species, extracted from reference workbook.")
message("Missing values by metric: ",
        paste(metric_columns, colSums(is.na(qc[metric_columns])),
              sep = "=", collapse = ", "))
message("Saved: ", output_prefix, ".pdf and ", output_prefix, ".png")
