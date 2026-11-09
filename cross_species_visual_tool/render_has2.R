# =====================================================================
# render_has2.R — static, non-interactive render of the HAS2 locus
# across NMR/BMR/DMR, using the exact same fetch_signal/fetch_models/
# fetch_sequence/build_plot logic as app.R's reactive pipeline, just
# driven by a fixed gene query instead of a live Shiny session.
#
#   Rscript render_has2.R
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
  library(data.table)
})

INDEX_DIR <- "index"
GENE      <- "HAS2"
WINDOW_KB <- 40      # half-window in kb, generous enough to see flanking peaks
BINS      <- 500

`%||%` <- function(a, b) if (is.null(a)) b else a

BAND <- list(
  mod_lo = -0.40, mod_hi = -0.10, mod_mid = -0.25,
  seq_lo = -0.86, seq_hi = -0.52,
  y_min  = -0.94, y_max  = 1.12
)
SEQ_GC_MAX <- 500000

manifest <- readRDS(file.path(INDEX_DIR, "manifest.rds"))
SP <- lapply(names(manifest), function(name) {
  prefix <- file.path(INDEX_DIR, manifest[[name]])
  meta <- readRDS(paste0(prefix, "_meta.rds"))
  meta$genes <- readRDS(paste0(prefix, "_genes.rds"))
  data.table::setkey(meta$genes, key)
  meta
})
names(SP) <- names(manifest)
SP_NAMES <- names(SP)

# ---- same fetch_* functions as app.R (kept identical on purpose) --------
fetch_signal <- function(meta, chrom, from, to, strand, n_bins) {
  seqname <- meta$map_bw[[chrom]]
  if (is.null(seqname) || is.na(seqname)) return(NULL)
  size <- meta$chrom_size[[chrom]]
  lo <- max(1L, as.integer(from)); hi <- min(as.integer(size), as.integer(to))
  if (lo > hi) return(NULL)
  gr <- GRanges(seqname, IRanges(lo, hi))
  vals <- try(as.numeric(
    rtracklayer::import(BigWigFile(meta$bigwig), which = gr, as = "NumericList")[[1]]
  ), silent = TRUE)
  if (inherits(vals, "try-error")) return(NULL)
  vals[is.na(vals)] <- 0
  full <- rep(NA_real_, as.integer(to) - as.integer(from) + 1L)
  full[(lo - as.integer(from) + 1L):(hi - as.integer(from) + 1L)] <- vals
  if (strand == "-") full <- rev(full)
  n <- length(full)
  bin <- cut(seq_len(n), breaks = n_bins, labels = FALSE)
  rel <- seq_len(n) - (n %/% 2L) - 1L
  data.table(rel = as.numeric(tapply(rel, bin, mean)),
             value = as.numeric(tapply(full, bin, mean, na.rm = TRUE)))
}

fetch_models <- function(meta, chrom, from, to) {
  gr <- GRanges(chrom, IRanges(max(1L, as.integer(from)), as.integer(to)))
  lines <- try(Rsamtools::scanTabix(TabixFile(meta$models), param = gr)[[1]], silent = TRUE)
  if (inherits(lines, "try-error") || !length(lines)) return(NULL)
  f <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  data.table(start = as.integer(f[[2]]) + 1L, end = as.integer(f[[3]]), tx = f[[4]],
             gname = f[[5]], strand = f[[6]], ex_s = f[[7]], ex_e = f[[8]])
}

# ---- anchor every species on its own HAS2 TSS ---------------------------
w <- WINDOW_KB * 1000
sig_list <- list(); mod_list <- list(); labels <- character(0); maxima <- numeric(0)

for (n in SP_NAMES) {
  hits <- SP[[n]]$genes[.("HAS2"), nomatch = 0L]
  if (!nrow(hits)) { message(n, ": HAS2 not found"); next }
  a <- list(chrom = hits$chrom[1], pos = as.numeric(hits$tss[1]), strand = hits$strand[1])
  centre <- a$pos
  from <- round(centre - w); to <- round(centre + w)

  sig <- fetch_signal(SP[[n]], a$chrom, from, to, a$strand, BINS)
  if (is.null(sig)) { message(n, ": no signal fetched"); next }
  sig[, species := n]
  sig_list[[n]] <- sig
  mx <- suppressWarnings(max(sig$value, na.rm = TRUE))
  maxima[n] <- if (is.finite(mx) && mx > 0) mx else 1
  labels[n] <- sprintf("%s    %s:%s–%s (%s)    max %.2f",
                       n, a$chrom, format(max(1, from), big.mark = ","),
                       format(to, big.mark = ","), a$strand, maxima[n])

  sgn <- if (a$strand == "-") -1 else 1
  mods <- fetch_models(SP[[n]], a$chrom, from, to)
  if (!is.null(mods) && nrow(mods)) {
    mods <- mods[order(gname, -(end - start))][, .SD[1], by = gname]
    rows <- lapply(seq_len(nrow(mods)), function(i) {
      es <- as.integer(strsplit(mods$ex_s[i], ",", fixed = TRUE)[[1]]) + 1L
      ee <- as.integer(strsplit(mods$ex_e[i], ",", fixed = TRUE)[[1]])
      x1 <- (es - centre) * sgn; x2 <- (ee - centre) * sgn
      data.table(species = n, tx = mods$tx[i], gname = mods$gname[i],
                 xmin = pmin(x1, x2), xmax = pmax(x1, x2),
                 body_lo = min((mods$start[i] - centre) * sgn, (mods$end[i] - centre) * sgn),
                 body_hi = max((mods$start[i] - centre) * sgn, (mods$end[i] - centre) * sgn))
    })
    mod_list[[n]] <- data.table::rbindlist(rows)
  }
  message(sprintf("%s: HAS2 at %s:%d-%d (%s), max signal %.3f",
                  n, a$chrom, from, to, a$strand, maxima[n]))
}

sig <- data.table::rbindlist(sig_list)
sig[, norm := value / maxima[species]]
sig[, species := factor(species, levels = SP_NAMES)]

mod <- if (length(mod_list)) data.table::rbindlist(mod_list) else NULL
if (!is.null(mod)) mod[, species := factor(species, levels = SP_NAMES)]

w_kb <- w / 1000
cols <- c(NMR = "#2E5E4E", BMR = "#7A4B72", DMR = "#B5622B")

p <- ggplot() +
  geom_area(data = sig, aes(x = rel / 1000, y = norm, fill = species),
            alpha = 0.9, colour = NA) +
  scale_fill_manual(values = cols, guide = "none") +
  facet_wrap(~species, ncol = 1)

if (!is.null(mod) && nrow(mod)) {
  lanes <- unique(mod[, .(species, tx, body_lo, body_hi)])
  p <- p +
    geom_segment(data = lanes, aes(x = body_lo / 1000, xend = body_hi / 1000,
                                    y = BAND$mod_mid, yend = BAND$mod_mid),
                 colour = "#5A6469", linewidth = 0.35) +
    geom_rect(data = mod, aes(xmin = xmin / 1000, xmax = xmax / 1000,
                               ymin = BAND$mod_lo + 0.04, ymax = BAND$mod_hi - 0.04),
              fill = "#3A4448", colour = NA)
}

lab <- data.table(species = factor(SP_NAMES, levels = SP_NAMES), text = labels[SP_NAMES])
p <- p +
  geom_text(data = lab, aes(x = -w_kb * 0.985, y = BAND$y_max - 0.06, label = text),
            hjust = 0, vjust = 1, size = 3.4, colour = "#14181A") +
  geom_vline(xintercept = 0, colour = "#000000", linewidth = 0.45) +
  geom_hline(yintercept = 0, colour = "#DCE1E3", linewidth = 0.3) +
  coord_cartesian(xlim = c(-w_kb, w_kb), ylim = c(BAND$y_min, BAND$y_max), expand = FALSE) +
  labs(title = "HAS2 — ATAC-seq accessibility across mole-rat species (liver)",
       subtitle = sprintf("Each panel anchored on that species' own HAS2 TSS, ±%g kb", WINDOW_KB),
       x = "Distance from HAS2 TSS (kb)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#EDF0F1", linewidth = 0.3),
        panel.spacing = grid::unit(0.6, "lines"),
        panel.border = element_rect(colour = "#E7EBED", fill = NA, linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "#5A6469", size = 10.5))

# cairo_pdf needs X11/cairo shared libs that aren't present on this machine;
# the default pdf() device has no such dependency.
ggsave("HAS2_cross_species.pdf", p, width = 10, height = 7)
ggsave("HAS2_cross_species.png", p, width = 10, height = 7, dpi = 200)
message("Saved HAS2_cross_species.pdf / .png")
