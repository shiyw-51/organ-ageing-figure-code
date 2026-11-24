# Supplementary Figure S1. Study sample flow chart.
# Rscript --vanilla SupplementaryFigure/FigureS1.R [output_directory]
# Use Unicode paths even when launched from a non-UTF-8 Windows shell.
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script_dir <- dirname(normalizePath(script_path, mustWork=TRUE))
# FIGURE_DATA_DIR points to a local directory containing the plotting inputs.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1L) stop("Supply at most one output directory.", call.=FALSE)
input_root <- Sys.getenv("FIGURE_DATA_DIR")
if (!nzchar(input_root) || !dir.exists(path.expand(input_root))) {
  stop("Set FIGURE_DATA_DIR to your local plotting input directory.", call.=FALSE)
}
input_dir <- file.path(path.expand(input_root), "SupplementaryFigure")
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir, call.=FALSE)
output_dir <- if (length(args)) path.expand(args[[1L]]) else file.path(dirname(script_dir), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
library(grid)
out_dir <- output_dir
data_dir <- file.path(input_dir, "FigureS1")

read_table <- function(path) {
  read.delim(path, check.names=FALSE, fileEncoding="UTF-8", stringsAsFactors=FALSE)
}
counts <- jsonlite::fromJSON(file.path(data_dir, "ukb_flow_counts.json"))
fmt <- function(n) formatC(as.integer(n), format = "d", big.mark = ",")
counts_at <- function(table, row, columns) {
  i <- which(table$Characteristic == row)
  stopifnot(length(i) == 1L, all(columns %in% names(table)))
  values <- as.character(unlist(table[i, columns, drop = FALSE], use.names = FALSE))
  n <- as.integer(gsub(",", "", sub(" .*", "", values), fixed = TRUE))
  stopifnot(!anyNA(n))
  setNames(n, columns)
}

ukb_n <- counts$ukb_n
split <- counts$split
complete_cov <- counts$complete_cov
qc <- counts$protein_qc
samples <- counts_at(read_table(file.path(data_dir, "ukb_sample_sizes.tsv")), "Participants, n",
                     c("Clinical layer", "Proteomic layer", "Both layers"))
clinical <- samples[["Clinical layer"]]
proteomic <- samples[["Proteomic layer"]]
both <- samples[["Both layers"]]

# development / locked-test split within each analytic sample (written from the analysis frame)
dt <- read_table(file.path(data_dir, "flow_split_counts_v1.tsv")); rownames(dt) <- dt$sample
stopifnot(dt["clinical", "total"] == clinical, dt["proteomic", "total"] == proteomic, dt["both", "total"] == both,
          all(dt$D + dt$T == dt$total))
dt_line <- function(k) paste0("Development ", fmt(dt[k, "D"]), "  ·  locked test ", fmt(dt[k, "T"]))

cohorts <- c("ELSA", "UKHLS", "CHARLS", "NHANES")
external <- read_table(file.path(dirname(input_dir), "Figure7", "TableS_ext1_baseline_v1.tsv"))
baseline <- counts_at(external, "Participants, n", cohorts)
valid <- counts_at(external, "Valid lung function test", cohorts)
results <- read_table(file.path(dirname(input_dir), "Figure7", "forest_cohorts_v1.tsv"))
lung <- subset(results, exposure == "pulmonary" & group == "copd_type" & model == "B")
stopifnot(nrow(lung) == length(cohorts), setequal(lung$cohort, cohorts))
lung <- lung[match(cohorts, lung$cohort), ]
stopifnot(ukb_n == split$D + split$T,
          qc$n_protein_file - qc$drop_gt1000_missing - qc$drop_sex_mismatch ==
            qc$n_after_sample_qc,
          complete_cov <= ukb_n, clinical <= complete_cov, proteomic <= qc$n_after_sample_qc,
          both <= min(clinical, proteomic), all(lung$N <= valid), all(valid <= baseline),
          all(lung$events <= lung$N),
          all(lung$outcome == c(rep("chronic_lung_disease", 3), "CLRD_death")))

W <- 183
H <- 189
B0 <- 11   # downward shift of everything below the two lane cards
font <- "Times New Roman"
ink <- "#202A2E"
slate <- "#45565E"
rule <- "#C9D3D6"
connector <- "#74858B"
blue <- "#375889"
teal <- "#147B6C"
blue_fill <- "#EAF0F8"
teal_fill <- "#E6F3F0"
mm <- function(x) unit(x, "mm")
top <- function(y) mm(H - y)

label <- function(x, y, text, size = 8, bold = FALSE, colour = ink, align = "left") {
  grid.text(text, x = mm(x), y = top(y), just = c(align, "centre"),
            gp = gpar(fontfamily = font, fontsize = size, fontface = if (bold) 2 else 1,
                      col = colour))
}
rect <- function(x, y, w, h, fill = "white", border = NA, lwd = .55) {
  grid.rect(x = mm(x), y = top(y), width = mm(w), height = mm(h),
            just = c("left", "top"), gp = gpar(fill = fill, col = border, lwd = lwd))
}
line <- function(x, y, colour = rule, lwd = .65, head = FALSE) {
  grid.lines(x = mm(x), y = top(y),
             arrow = if (head) arrow(length = mm(1.25), type = "closed") else NULL,
             gp = gpar(col = colour, fill = colour, lwd = lwd, lineend = "butt", linejoin = "round"))
}
heading <- function(y, letter, title) {
  label(3.1, y, letter, 11.6, TRUE)
  label(10, y, title, 9.4, TRUE)
}
card <- function(x, y, w, h, fill = "white", border = rule) {
  grid.roundrect(x = mm(x), y = top(y), width = mm(w), height = mm(h),
                 just = c("left", "top"), r = mm(1.0),
                 gp = gpar(fill = fill, col = border, lwd = .65))
}
value_row <- function(left, right, y, text, value, colour = ink, emphasis = FALSE) {
  label(left, y, text, if (emphasis) 8.05 else 7.85, emphasis, colour)
  label(right, y, fmt(value), if (emphasis) 9.3 else 8.25, TRUE, colour, "right")
}
lane_header <- function(x, title, accent, tint) {
  card(x, 31, 84.2, 49)
  rect(x + .35, 31.35, 83.5, 8.2, tint)
  label(x + 3.3, 35.6, title, 8.9, TRUE, accent)
  line(c(x + 3.3, x + 80.9), c(61.7, 61.7))
}

draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  heading(5.2, "A", "UK Biobank")

  # The protein measurements are a subset of the baseline cohort before covariate exclusions.
  card(44, 9.6, 95, 14.4, "#F5F7F8")
  value_row(47.5, 135.5, 14.5, "Baseline participants, 2006–2010", ukb_n, emphasis = TRUE)
  label(47.5, 20.0, paste0("Development set ", fmt(split$D), "     ·     Locked test set ", fmt(split$T)),
        7.65, colour = slate)
  centres <- c(45.6, 137.4)
  line(c(91.5, 91.5), c(24, 27.5), connector)
  line(centres, c(27.5, 27.5), connector)
  for (x in centres) line(c(x, x), c(27.5, 30.8), connector, head = TRUE)

  lane_header(3.5, "Clinical clocks", blue, blue_fill)
  value_row(6.8, 84.4, 46.6, "Complete model B covariates", complete_cov)
  label(6.8, 52.5, paste0(fmt(ukb_n - complete_cov), " excluded for missing covariates"),
        7.25, colour = slate)
  value_row(6.8, 84.4, 66.3, "At least one clinical organ age", clinical, blue, TRUE)
  label(6.8, 71.1, paste0(fmt(complete_cov - clinical), " excluded for lack of a clinical organ age"),
        7.15, colour = slate)
  label(6.8, 76.3, dt_line("clinical"), 7.4, colour = blue)

  lane_header(95.3, "Proteomic clocks", teal, teal_fill)
  value_row(98.6, 176.2, 43.8, "Olink measurements available", qc$n_protein_file)
  value_row(98.6, 176.2, 49.4, "Passed quality control", qc$n_after_sample_qc)
  label(98.6, 54.1, paste0(fmt(qc$drop_gt1000_missing), " excluded for >1,000 missing proteins"),
        7.15, colour = slate)
  label(98.6, 58.2, paste0(fmt(qc$drop_sex_mismatch), " excluded for discordant reported and genetic sex"),
        7.15, colour = slate)
  value_row(98.6, 176.2, 66.3, "Proteomic age and complete covariates", proteomic, teal, TRUE)
  label(98.6, 71.1, paste0(fmt(qc$n_after_sample_qc - proteomic), " excluded for missing model B covariates"),
        7.15, colour = slate)
  label(98.6, 76.3, dt_line("proteomic"), 7.4, colour = teal)

  for (x in centres) line(c(x, x), c(80, 83.5), connector)
  line(centres, c(83.5, 83.5), connector)
  line(c(91.5, 91.5), c(83.5, 87.3), connector, head = TRUE)
  card(44, 87.5, 95, 18.4, "#F0F3F5")
  value_row(47.5, 135.5, 92.0, "Participants with both layers", both, emphasis = TRUE)
  label(47.5, 97.1, "All 14 organ ages and complete model B covariates", 7.6, colour = slate)
  label(47.5, 102.0, dt_line("both"), 7.4, colour = slate)

  line(c(3.5, 179.5), c(99.3, 99.3) + B0)
  heading(105.2 + B0, "B", "Lung clock validation in external cohorts")
  label(4.8, 113.0 + B0, "Cohort", 7.65, TRUE)
  label(68.1, 113.0 + B0, "Baseline, n", 7.65, TRUE, align = "right")
  label(77.3, 113.0 + B0, "Lung measure", 7.65, TRUE)
  label(130.0, 113.0 + B0, "Valid lung test, n", 7.65, TRUE, align = "right")
  label(155.1, 113.0 + B0, "Analysed, n", 7.65, TRUE, align = "right")
  label(178.0, 113.0 + B0, "Events, n", 7.65, TRUE, align = "right")
  # The full endpoint names distinguish new diagnoses from mortality without a glossary.
  group_bar <- function(y, text) {
    rect(3.5, y - 3.0, 176, 6.0, "#EDF1F3")
    label(5.0, y, text, 7.8, TRUE)
  }
  group_bar(120.0 + B0, "Incident chronic lung disease")

  selection <- c("Core members ≥50 y · nurse visit, wave 2 or 4", "≥40 y · nurse visit, wave 2",
                 "≥45 y · 2011 or 2015", "≥40 y · 1988–94 or 2007–12")
  measure <- c("Spirometry", "Spirometry", "Peak expiratory flow", "Spirometry")
  ys <- c(127.4, 138.2, 149.0, 166.4) + B0
  for (i in seq_along(cohorts)) {
    if (cohorts[i] == "NHANES") {
      group_bar(158.6 + B0, "Death from chronic lower respiratory disease through 2019")
    }
    y <- ys[i]
    label(5, y, cohorts[i], 8.05, TRUE)
    label(5, y + 4.2, selection[i], 7.0, colour = slate)
    label(68.1, y + 1.6, fmt(baseline[i]), 8.3, align = "right")
    label(77.3, y + 1.6, measure[i], 7.6)
    label(130.0, y + 1.6, fmt(valid[i]), 8.3, align = "right")
    label(155.1, y + 1.6, fmt(lung$N[i]), 8.35, TRUE, align = "right")
    label(178.0, y + 1.6, fmt(lung$events[i]), 8.3, align = "right")
  }
  for (y in c(134.6, 145.4, 174.4) + B0) line(c(3.5, 179.5), c(y, y), lwd = .5)
}

stem <- file.path(out_dir, "FigureS1_flowchart")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = W / 25.4, height = H / 25.4,
                     family = font, bg = "white")
draw()
invisible(dev.off())
# Align raster rows to four pixels for consistent Windows Cairo exports.
pixel_dimensions <- function(dpi) {
  c(width = 4L * round(W * dpi / (25.4 * 4)), height = round(H * dpi / 25.4))
}
size <- pixel_dimensions(300)
grDevices::png(paste0(stem, ".png"), width = size["width"], height = size["height"],
               res = 300, type = "cairo", bg = "white")
draw()
invisible(dev.off())
size <- pixel_dimensions(600)
grDevices::tiff(paste0(stem, ".tiff"), width = size["width"], height = size["height"],
                res = 600, compression = "lzw", type = "cairo", bg = "white")
draw()
invisible(dev.off())
