# =====================================================================
# Cross-species ATAC-seq browser  —  whole genome
#
#   Rscript build_index.R config.tsv index/     # once, per data set
#   shiny::runApp("app.R")
#
# Coordinate model
# ----------------
#   genomic position = anchor[sp] + nudge[sp] + pan + relative_x
#
#   anchor  per species; a gene's TSS, or the midpoint of a region you typed
#   nudge   per species; manual fine correction, in bp
#   pan     shared; moves every panel together
#   window  shared; half-width of the view
#
# Each panel is normalised onto the 0-1 plotting band using a per-species
# data range [range_lo, range_hi] (default 0-1, editable per species in the
# "Species" panel). Values outside that range saturate at the top/bottom of
# the band rather than distorting other panels. The true observed maximum is
# still printed on each panel as a reference for picking a sensible range.
# "Shared" mode uses one range (the min/max across active species' own
# ranges) for every panel instead of each species' own range, so relative
# magnitudes stay comparable across panels.
# =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
  library(data.table)
})

INDEX_DIR <- Sys.getenv("ATAC_INDEX_DIR", "index")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Vertical layout within each panel, in normalised units.
BAND <- list(
  mod_lo = -0.40, mod_hi = -0.10, mod_mid = -0.25,
  seq_lo = -0.86, seq_hi = -0.52,
  y_min  = -0.94, y_max  = 1.12
)

SEQ_LETTER_MAX <- 400      # view width below which individual bases are drawn
SEQ_GC_MAX     <- 500000   # view width below which GC content is computed
BASE_COLS <- c(A = "#2E7D5B", C = "#2A5C8A", G = "#B5622B",
               T = "#8C3B3B", N = "#9AA3A8")


# ---------------------------------------------------------------------
# Load the index once at process start, shared by every session
# ---------------------------------------------------------------------

if (!dir.exists(INDEX_DIR) || !file.exists(file.path(INDEX_DIR, "manifest.rds"))) {
  stop("No index in '", INDEX_DIR, "'. Run: Rscript build_index.R config.tsv ", INDEX_DIR)
}

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

# Search index: symbols ordered by how many species carry them.
#
# NOTE: data.table(key = ..., ...) collides with the data.table() constructor's
# own reserved `key` parameter (used to auto-set the new table's key to a
# vector of COLUMN NAMES). Naming a column "key" here makes it try to key the
# table by the gene symbols themselves, which fails with "some columns are
# not in the data.table: [<gene symbols>]" the moment two species are loaded.
# Build the column under a different name and rename it after construction.
GENE_KEYS <- local({
  tab <- data.table::rbindlist(lapply(SP_NAMES, function(n) {
    dt <- data.table(gene_key = unique(SP[[n]]$genes$key), species = n)
    data.table::setnames(dt, "gene_key", "key")
    dt
  }))
  counts <- tab[, .(n = .N), by = key][order(-n, key)]
  stats::setNames(counts$key,
                  sprintf("%s  (%d/%d)", counts$key, counts$n, length(SP_NAMES)))
})

message(sprintf("Loaded %d species, %s searchable symbols",
                length(SP), format(length(GENE_KEYS), big.mark = ",")))


# ---------------------------------------------------------------------
# Data access
# ---------------------------------------------------------------------

#' Binned ATAC signal over an interval, reversed for minus-strand anchors.
fetch_signal <- function(meta, chrom, from, to, strand, n_bins) {
  seqname <- meta$map_bw[[chrom]]
  if (is.null(seqname) || is.na(seqname)) return(NULL)

  size <- meta$chrom_size[[chrom]]
  lo <- max(1L, as.integer(from))
  hi <- min(as.integer(size), as.integer(to))
  if (lo > hi) return(NULL)

  gr <- GRanges(seqname, IRanges(lo, hi))
  vals <- try(as.numeric(
    rtracklayer::import(BigWigFile(meta$bigwig), which = gr, as = "NumericList")[[1]]
  ), silent = TRUE)
  if (inherits(vals, "try-error")) return(NULL)
  vals[is.na(vals)] <- 0

  # Pad off-chromosome regions with NA, not 0, so a truncated view does not
  # masquerade as a closed region.
  full <- rep(NA_real_, as.integer(to) - as.integer(from) + 1L)
  full[(lo - as.integer(from) + 1L):(hi - as.integer(from) + 1L)] <- vals
  if (strand == "-") full <- rev(full)

  n <- length(full)
  bin <- cut(seq_len(n), breaks = n_bins, labels = FALSE)
  rel <- seq_len(n) - (n %/% 2L) - 1L
  data.table(rel   = as.numeric(tapply(rel, bin, mean)),
             value = as.numeric(tapply(full, bin, mean, na.rm = TRUE)))
}


#' Transcript models overlapping an interval, from the tabix index.
fetch_models <- function(meta, chrom, from, to) {
  gr <- GRanges(chrom, IRanges(max(1L, as.integer(from)), as.integer(to)))
  lines <- try(Rsamtools::scanTabix(TabixFile(meta$models), param = gr)[[1]],
               silent = TRUE)
  if (inherits(lines, "try-error") || !length(lines)) return(NULL)

  f <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  data.table(start  = as.integer(f[[2]]) + 1L,
             end    = as.integer(f[[3]]),
             tx     = f[[4]],
             gname  = f[[5]],
             strand = f[[6]],
             ex_s   = f[[7]],
             ex_e   = f[[8]])
}


#' Reference sequence over an interval, reverse-complemented if needed.
fetch_sequence <- function(meta, chrom, from, to, strand) {
  seqname <- meta$map_fa[[chrom]]
  if (is.null(seqname) || is.na(seqname)) return(NULL)
  size <- meta$chrom_size[[chrom]]
  lo <- max(1L, as.integer(from)); hi <- min(as.integer(size), as.integer(to))
  if (lo > hi) return(NULL)

  s <- try(Rsamtools::scanFa(FaFile(meta$fasta),
                             GRanges(seqname, IRanges(lo, hi)))[[1]],
           silent = TRUE)
  if (inherits(s, "try-error")) return(NULL)
  if (strand == "-") s <- Biostrings::reverseComplement(s)
  s
}


parse_locus <- function(txt) {
  txt <- gsub("[, ]", "", trimws(txt))
  if (!nzchar(txt)) return(NULL)
  m <- regmatches(txt, regexec("^([^:]+):([0-9]+)(?:-([0-9]+))?$", txt))[[1]]
  if (!length(m)) return(NULL)
  start <- as.numeric(m[3])
  end   <- if (nzchar(m[4])) as.numeric(m[4]) else start
  if (end < start) { tmp <- start; start <- end; end <- tmp }
  list(chrom = m[2], start = start, end = end,
       mid = floor((start + end) / 2), width = end - start + 1)
}


# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------

css <- "
body { font-family:'Inter',-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
       background:#E9ECEE; color:#14181A; }
.container-fluid { max-width:1520px; }
h1 { font-size:19px; font-weight:620; letter-spacing:-.015em; margin:18px 0 2px; }
.sub { color:#5A6469; font-size:13px; margin:0 0 16px; }
.card { background:#fff; border:1px solid #D3D9DC; border-radius:3px;
        padding:14px 16px; margin-bottom:14px; }
.card h2 { font-size:12px; font-weight:600; color:#5A6469; margin:0 0 10px; }
.mono input { font-family:ui-monospace,'SF Mono',Menlo,Consolas,monospace; font-size:12px; }
.navrow { display:flex; gap:10px; align-items:flex-end; flex-wrap:wrap; }
.navrow .grow { flex:1 1 260px; min-width:200px; }
.sp-block { border-left:3px solid #D3D9DC; padding:2px 0 2px 10px; margin-bottom:14px; }
.sp-block .nm { font-size:13.5px; font-weight:600; display:flex; align-items:center; gap:6px; }
.sp-block .loc { font-family:ui-monospace,monospace; font-size:11px; color:#5A6469;
                 margin:1px 0 6px; }
.swatch { width:10px; height:10px; border-radius:2px; display:inline-block; }
.msg { font-size:12.5px; color:#5A6469; background:#fff; border:1px solid #D3D9DC;
       border-left:2px solid #B5622B; padding:8px 11px; border-radius:0 2px 2px 0;
       margin-bottom:12px; }
.form-group { margin-bottom:8px; }
label { font-weight:500; font-size:13px; color:#5A6469; }
.btn-xs { padding:3px 8px; font-size:11.5px; }
"

ui <- fluidPage(
  tags$head(tags$style(HTML(css))),

  h1("Cross-species ATAC-seq browser"),
  p(class = "sub",
    paste("Panels are locked to each species' own TSS. Pan and zoom move every panel",
          "together; type a region to move one panel on its own.")),

  div(class = "card",
    div(class = "navrow",
      div(class = "grow",
        selectizeInput("gene", "Gene", choices = NULL, width = "100%",
                       options = list(placeholder = "Search a symbol",
                                      maxOptions = 200))),
      div(numericInput("window_kb", "Half-window (kb)", value = 5,
                       min = 0.05, max = 2000, step = 1, width = "140px")),
      div(numericInput("bins", "Bins", value = 500, min = 100, max = 2000,
                       step = 100, width = "100px")),
      div(style = "padding-bottom:12px",
        actionButton("zoom_out", "\u2212",     class = "btn btn-default btn-xs"),
        actionButton("zoom_in",  "+",          class = "btn btn-default btn-xs"),
        actionButton("pan_l",    "\u2190",     class = "btn btn-default btn-xs"),
        actionButton("pan_r",    "\u2192",     class = "btn btn-default btn-xs"),
        actionButton("recentre", "Recentre",   class = "btn btn-default btn-xs"))
    )
  ),

  fluidRow(
    column(3,
      div(class = "card",
        h2("Species"),
        uiOutput("species_controls")),
      div(class = "card",
        h2("Display"),
        radioButtons("yaxis", NULL, inline = TRUE,
                     choices = c("Independent" = "free", "Shared" = "fixed"),
                     selected = "free"),
        checkboxInput("show_models", "Gene models", TRUE),
        checkboxInput("show_seq", "Sequence / GC", TRUE),
        checkboxInput("all_tx", "All transcripts", FALSE),
        checkboxInput("strand_fix", "Correct for strand", TRUE),
        tags$hr(style = "margin:10px 0"),
        downloadButton("dl_pdf", "PDF", class = "btn-xs"),
        downloadButton("dl_csv", "CSV", class = "btn-xs"))
    ),
    column(9,
      uiOutput("messages"),
      div(class = "card", style = "padding:8px 10px",
        plotOutput("browser", height = "660px",
                   brush = brushOpts(id = "zoom_brush", direction = "x",
                                     resetOnNew = TRUE),
                   dblclick = "zoom_reset")),
      p(class = "sub", style = "margin-top:-6px",
        "Drag across the plot to zoom to that span. Double-click to reset the window.")
    )
  )
)


# ---------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------

server <- function(input, output, session) {

  updateSelectizeInput(session, "gene", choices = GENE_KEYS,
                       selected = character(0), server = TRUE)

  view <- reactiveValues(
    anchor   = stats::setNames(vector("list", length(SP_NAMES)), SP_NAMES),
    nudge    = stats::setNames(as.list(rep(0, length(SP_NAMES))), SP_NAMES),
    scale    = stats::setNames(as.list(rep(1, length(SP_NAMES))), SP_NAMES),
    show     = stats::setNames(as.list(rep(TRUE, length(SP_NAMES))), SP_NAMES),
    range_lo = stats::setNames(as.list(rep(0, length(SP_NAMES))), SP_NAMES),
    range_hi = stats::setNames(as.list(rep(1, length(SP_NAMES))), SP_NAMES),
    pan    = 0,
    window = 5000,
    notes  = character(0)
  )

  # ---- picking a gene re-anchors every panel --------------------------
  observeEvent(input$gene, {
    # NOTE: this used to be `key <- input$gene`. Every species' genes table
    # has a column literally named "key" (that's what it's keyed on), and
    # data.table's `.()` inside `[i]` resolves a bare symbol against the
    # table's OWN columns before the calling environment — so `.(key)` was
    # silently joining the whole 30k+-row key column against itself instead
    # of looking up the single typed gene, producing a >1M-row cartesian
    # result and blowing up every search. Renamed to sidestep the collision
    # entirely rather than rely on data.table's `..` env-escape prefix.
    gene_query <- input$gene
    req(nzchar(gene_query))
    notes <- character(0)

    for (n in SP_NAMES) {
      hits <- SP[[n]]$genes[.(gene_query), nomatch = 0L]
      if (!nrow(hits)) {
        view$anchor[[n]] <- NULL
        notes <- c(notes, sprintf("%s: %s is not in this annotation.", n, gene_query))
        next
      }
      if (nrow(hits) > 1) {
        notes <- c(notes, sprintf("%s: %s maps to %d loci; using the longest, on %s.",
                                  n, gene_query, nrow(hits), hits$chrom[which.max(hits$end - hits$start)]))
        hits <- hits[which.max(end - start)]
      }
      view$anchor[[n]] <- list(chrom = hits$chrom[1], pos = as.numeric(hits$tss[1]),
                               strand = hits$strand[1], label = hits$gname[1],
                               manual = FALSE)
      updateTextInput(session, paste0("loc_", n), value = "")
    }
    for (n in SP_NAMES) view$nudge[[n]] <- 0
    view$pan <- 0
    view$notes <- notes
  })

  # ---- per-species controls -------------------------------------------
  output$species_controls <- renderUI({
    lapply(SP_NAMES, function(n) {
      a <- view$anchor[[n]]
      loc <- if (is.null(a)) "no anchor"
             else sprintf("%s:%s (%s)%s", a$chrom,
                          format(round(a$pos), big.mark = ","), a$strand,
                          if (isTRUE(a$manual)) "  manual" else "")
      div(class = "sp-block", style = paste0("border-left-color:", SP[[n]]$colour),
        div(class = "nm",
            tags$span(class = "swatch",
                      style = paste0("background:", SP[[n]]$colour)), n),
        div(class = "loc", loc),
        fluidRow(
          column(7, div(class = "mono",
            textInput(paste0("loc_", n), NULL, value = "",
                      placeholder = "chr:start-end"))),
          column(5, numericInput(paste0("scale_", n), NULL,
                                 value = isolate(view$scale[[n]]),
                                 min = 0, step = 0.1))),
        fluidRow(
          column(7, div(class = "mono",
            numericInput(paste0("nudge_", n), NULL,
                         value = isolate(view$nudge[[n]]), step = 100))),
          column(5, checkboxInput(paste0("show_", n), "show",
                                  value = isolate(view$show[[n]])))),
        fluidRow(
          column(6, numericInput(paste0("range_lo_", n), "Range min",
                                 value = isolate(view$range_lo[[n]]), step = 0.1)),
          column(6, numericInput(paste0("range_hi_", n), "Range max",
                                 value = isolate(view$range_hi[[n]]), step = 0.1)))
      )
    })
  })

  observe({
    for (n in SP_NAMES) {
      s <- input[[paste0("scale_", n)]]
      if (!is.null(s) && is.finite(s) && s >= 0) view$scale[[n]] <- s
      v <- input[[paste0("show_", n)]]
      if (!is.null(v)) view$show[[n]] <- isTRUE(v)
      g <- input[[paste0("nudge_", n)]]
      if (!is.null(g) && is.finite(g)) view$nudge[[n]] <- g
      lo <- input[[paste0("range_lo_", n)]]
      hi <- input[[paste0("range_hi_", n)]]
      if (!is.null(lo) && !is.null(hi) && is.finite(lo) && is.finite(hi) && hi > lo) {
        view$range_lo[[n]] <- lo
        view$range_hi[[n]] <- hi
      }
    }
  })

  # Typing a region re-anchors that one species; a span also sets the
  # shared window, since that is almost always what is meant.
  lapply(SP_NAMES, function(n) {
    observeEvent(input[[paste0("loc_", n)]], {
      txt <- input[[paste0("loc_", n)]]
      if (is.null(txt) || !nzchar(trimws(txt))) return()
      p <- parse_locus(txt)
      if (is.null(p)) {
        view$notes <- unique(c(view$notes,
          sprintf("%s: could not read '%s'. Expected chr:start-end.", n, txt)))
        return()
      }
      if (is.null(SP[[n]]$chrom_size[[p$chrom]])) {
        view$notes <- unique(c(view$notes,
          sprintf("%s: '%s' is not a sequence in this assembly.", n, p$chrom)))
        return()
      }
      view$anchor[[n]] <- list(chrom = p$chrom, pos = p$mid, strand = "+",
                               label = "region", manual = TRUE)
      view$nudge[[n]] <- 0
      if (p$width > 1) {
        view$window <- max(50, p$width / 2)
        updateNumericInput(session, "window_kb", value = round(view$window / 1000, 4))
        view$pan <- 0
      }
    }, ignoreInit = TRUE)
  })

  # ---- navigation -------------------------------------------------------
  set_window <- function(w) {
    view$window <- max(50, min(2e8, w))
    updateNumericInput(session, "window_kb", value = round(view$window / 1000, 4))
  }

  observeEvent(input$window_kb, {
    w <- input$window_kb
    if (is.finite(w) && w > 0) view$window <- w * 1000
  })
  observeEvent(input$zoom_in,  { set_window(view$window / 2) })
  observeEvent(input$zoom_out, { set_window(view$window * 2) })
  observeEvent(input$pan_l,    { view$pan <- view$pan - view$window * 0.4 })
  observeEvent(input$pan_r,    { view$pan <- view$pan + view$window * 0.4 })
  observeEvent(input$recentre, {
    view$pan <- 0
    for (n in SP_NAMES) {
      view$nudge[[n]] <- 0
      updateNumericInput(session, paste0("nudge_", n), value = 0)
    }
  })

  observeEvent(input$zoom_brush, {
    b <- input$zoom_brush
    req(b)
    lo <- b$xmin * 1000; hi <- b$xmax * 1000     # the plot is drawn in kb
    if (hi - lo < 20) return()
    view$pan <- view$pan + (lo + hi) / 2
    set_window((hi - lo) / 2)
    session$resetBrush("zoom_brush")
  })

  observeEvent(input$zoom_reset, { view$pan <- 0; set_window(5000) })

  # ---- gather everything the plot needs ---------------------------------
  frame <- reactive({
    active <- SP_NAMES[vapply(SP_NAMES, function(n)
      isTRUE(view$show[[n]]) && !is.null(view$anchor[[n]]), logical(1))]
    if (!length(active)) return(NULL)

    w <- view$window
    n_bins <- max(50, min(2000, as.integer(input$bins %||% 500)))
    span <- 2 * w
    notes <- character(0)

    sig_list <- list(); mod_list <- list(); seq_list <- list()
    labels <- character(0); maxima <- numeric(0)

    # "Shared" mode uses one data range for every panel: the widest span
    # covered by the active species' own range settings. "Independent" uses
    # each species' own [range_lo, range_hi] as configured in the Species
    # panel (default 0-1).
    shared_range <- identical(input$yaxis, "fixed")
    if (shared_range) {
      g_lo <- min(vapply(active, function(n) view$range_lo[[n]], numeric(1)))
      g_hi <- max(vapply(active, function(n) view$range_hi[[n]], numeric(1)))
    }

    for (n in active) {
      a <- view$anchor[[n]]
      strand <- if (isTRUE(input$strand_fix)) a$strand else "+"
      centre <- a$pos + view$nudge[[n]] +
                (if (strand == "-") -view$pan else view$pan)
      from <- round(centre - w); to <- round(centre + w)

      sig <- fetch_signal(SP[[n]], a$chrom, from, to, strand, n_bins)
      if (is.null(sig)) {
        notes <- c(notes, sprintf("%s: no signal for %s in the bigWig.", n, a$chrom))
        next
      }
      sig[, value := value * view$scale[[n]]]
      sig[, species := n]

      lo_n <- if (shared_range) g_lo else view$range_lo[[n]]
      hi_n <- if (shared_range) g_hi else view$range_hi[[n]]
      sig[, norm := (value - lo_n) / (hi_n - lo_n)]
      sig_list[[n]] <- sig

      mx <- suppressWarnings(max(sig$value, na.rm = TRUE))
      maxima[n] <- if (is.finite(mx) && mx > 0) mx else 1
      labels[n] <- sprintf("%s    %s:%s\u2013%s (%s)    max %.2f  range %.2g\u2013%.2g",
                           n, a$chrom, format(max(1, from), big.mark = ","),
                           format(to, big.mark = ","), strand, maxima[n], lo_n, hi_n)

      sgn <- if (strand == "-") -1 else 1

      if (isTRUE(input$show_models)) {
        mods <- fetch_models(SP[[n]], a$chrom, from, to)
        if (!is.null(mods) && nrow(mods)) {
          if (!isTRUE(input$all_tx)) {
            mods <- mods[order(gname, -(end - start))][, .SD[1], by = gname]
          } else if (nrow(mods) > 8) {
            mods <- mods[order(-(end - start))][1:8]
            notes <- c(notes, sprintf("%s: showing 8 of %d transcripts.", n, nrow(mods)))
          }
          rows <- lapply(seq_len(nrow(mods)), function(i) {
            es <- as.integer(strsplit(mods$ex_s[i], ",", fixed = TRUE)[[1]]) + 1L
            ee <- as.integer(strsplit(mods$ex_e[i], ",", fixed = TRUE)[[1]])
            x1 <- (es - centre) * sgn; x2 <- (ee - centre) * sgn
            data.table(species = n, tx = mods$tx[i], gname = mods$gname[i],
                       xmin = pmin(x1, x2), xmax = pmax(x1, x2),
                       body_lo = min((mods$start[i] - centre) * sgn,
                                     (mods$end[i] - centre) * sgn),
                       body_hi = max((mods$start[i] - centre) * sgn,
                                     (mods$end[i] - centre) * sgn))
          })
          mod_list[[n]] <- data.table::rbindlist(rows)
        }
      }

      if (isTRUE(input$show_seq) && span <= SEQ_GC_MAX) {
        s <- fetch_sequence(SP[[n]], a$chrom, from, to, strand)
        if (!is.null(s)) {
          chars <- toupper(strsplit(as.character(s), "", fixed = TRUE)[[1]])
          mid <- length(chars) %/% 2L
          if (span <= SEQ_LETTER_MAX) {
            seq_list[[n]] <- data.table(species = n, mode = "letters",
                                        rel = seq_along(chars) - mid - 1L,
                                        base = chars, value = NA_real_)
          } else {
            nb <- min(400L, max(20L, length(chars) %/% 50L))
            b <- cut(seq_along(chars), breaks = nb, labels = FALSE)
            seq_list[[n]] <- data.table(
              species = n, mode = "gc",
              rel = as.numeric(tapply(seq_along(chars) - mid - 1L, b, mean)),
              base = NA_character_,
              value = as.numeric(tapply(chars %in% c("G", "C"), b, mean)))
          }
        }
      }
    }

    if (!length(sig_list)) return(list(sig = NULL, notes = c(view$notes, notes)))

    sig <- data.table::rbindlist(sig_list)
    keep <- intersect(active, unique(sig$species))
    sig[, species := factor(species, levels = keep)]

    list(sig = sig,
         mod = if (length(mod_list)) data.table::rbindlist(mod_list) else NULL,
         seq = if (length(seq_list)) data.table::rbindlist(seq_list) else NULL,
         labels = labels, active = keep, window = w,
         notes = unique(c(view$notes, notes)))
  })

  # ---- plot --------------------------------------------------------------
  build_plot <- reactive({
    fr <- frame()
    validate(need(!is.null(fr) && !is.null(fr$sig),
                  "Search a gene, or type a region for one species."))

    w_kb <- fr$window / 1000
    cols <- stats::setNames(vapply(fr$active, function(n) SP[[n]]$colour,
                                   character(1)), fr$active)

    p <- ggplot() +
      geom_area(data = fr$sig, aes(x = rel / 1000, y = norm, fill = species),
                alpha = 0.9, colour = NA) +
      scale_fill_manual(values = cols, guide = "none") +
      facet_wrap(~species, ncol = 1)

    if (!is.null(fr$mod) && nrow(fr$mod)) {
      m <- data.table::copy(fr$mod)[species %in% fr$active]
      m[, species := factor(species, levels = fr$active)]
      lanes <- unique(m[, .(species, tx, body_lo, body_hi)])
      p <- p +
        geom_segment(data = lanes,
                     aes(x = body_lo / 1000, xend = body_hi / 1000,
                         y = BAND$mod_mid, yend = BAND$mod_mid),
                     colour = "#5A6469", linewidth = 0.35) +
        geom_rect(data = m,
                  aes(xmin = xmin / 1000, xmax = xmax / 1000,
                      ymin = BAND$mod_lo + 0.04, ymax = BAND$mod_hi - 0.04),
                  fill = "#3A4448", colour = NA)
    }

    if (!is.null(fr$seq) && nrow(fr$seq)) {
      s <- data.table::copy(fr$seq)[species %in% fr$active]
      s[, species := factor(species, levels = fr$active)]
      lt <- s[mode == "letters"]; gc <- s[mode == "gc"]
      if (nrow(lt)) {
        p <- p +
          geom_text(data = lt,
                    aes(x = rel / 1000, y = (BAND$seq_lo + BAND$seq_hi) / 2,
                        label = base, colour = base),
                    size = 2.5, show.legend = FALSE) +
          scale_colour_manual(values = BASE_COLS, guide = "none",
                              na.value = "#9AA3A8")
      }
      if (nrow(gc)) {
        p <- p +
          geom_hline(yintercept = BAND$seq_lo + 0.5 * (BAND$seq_hi - BAND$seq_lo),
                     colour = "#EDF0F1", linewidth = 0.3) +
          geom_line(data = gc,
                    aes(x = rel / 1000,
                        y = BAND$seq_lo + value * (BAND$seq_hi - BAND$seq_lo)),
                    colour = "#5A6469", linewidth = 0.35)
      }
    }

    lab <- data.table(species = factor(fr$active, levels = fr$active),
                      text = fr$labels[fr$active])
    p +
      geom_text(data = lab,
                aes(x = -w_kb * 0.985, y = BAND$y_max - 0.06, label = text),
                hjust = 0, vjust = 1, size = 3.1, colour = "#14181A") +
      geom_vline(xintercept = 0, colour = "#000000", linewidth = 0.45) +
      geom_hline(yintercept = 0, colour = "#DCE1E3", linewidth = 0.3) +
      coord_cartesian(xlim = c(-w_kb, w_kb),
                      ylim = c(BAND$y_min, BAND$y_max), expand = FALSE) +
      labs(x = "Distance from anchor (kb)", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(strip.text = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(colour = "#EDF0F1", linewidth = 0.3),
            panel.spacing = grid::unit(0.5, "lines"),
            panel.border = element_rect(colour = "#E7EBED", fill = NA,
                                        linewidth = 0.4))
  })

  output$browser <- renderPlot({ build_plot() }, res = 96)

  output$messages <- renderUI({
    fr <- try(frame(), silent = TRUE)
    msgs <- if (inherits(fr, "try-error") || is.null(fr)) view$notes else fr$notes
    msgs <- unique(msgs)
    if (!length(msgs)) return(NULL)
    div(class = "msg", lapply(msgs, function(x) div(x)))
  })

  # ---- downloads ---------------------------------------------------------
  output$dl_pdf <- downloadHandler(
    filename = function() sprintf("%s_%gkb.pdf", input$gene %||% "region",
                                  view$window / 1000),
    content = function(file) {
      fr <- frame()
      # cairo_pdf requires X11/cairo shared libs; not every deployment target
      # has them (confirmed missing in at least one environment this app was
      # run in), and it fails silently under ggsave rather than erroring
      # loudly. The default pdf() device has no such dependency and covers
      # the same glyphs used here (no special unicode beyond em dash/en dash).
      ggsave(file, build_plot(), width = 10,
             height = max(3, 1.7 * length(fr$active) + 1),
             limitsize = FALSE)
    })

  output$dl_csv <- downloadHandler(
    filename = function() sprintf("%s_%gkb.csv", input$gene %||% "region",
                                  view$window / 1000),
    content = function(file) {
      fr <- frame()
      data.table::fwrite(
        fr$sig[, .(species, rel_bp = rel, value, normalised = norm)], file)
    })
}

shinyApp(ui, server)
