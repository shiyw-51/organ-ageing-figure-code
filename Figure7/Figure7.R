# Figure7
# Rscript --vanilla Figure7/Figure7.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure7")
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir, call.=FALSE)
output_dir <- if (length(args)) path.expand(args[[1L]]) else file.path(dirname(script_dir), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

write_figure_rasters <- function(stem, width_mm, height_mm, draw) {
  # Match the integer-point canvas of the Cairo PDF at each resolution.
  page_points <- floor(c(width_mm, height_mm) * 72 / 25.4)
  size <- function(dpi) ceiling(page_points * dpi / 72)
  px <- size(300)
  preview <- paste0(stem, ".png")
  grDevices::png(preview, width=px[1L], height=px[2L], res=300, type="cairo",
                 bg="white", antialias="gray")
  tryCatch(draw(), finally=invisible(grDevices::dev.off()))
  png::writePNG(png::readPNG(preview, native=TRUE), preview, dpi=c(300,300))

  px <- size(600)
  grDevices::tiff(paste0(stem, ".tiff"), width=px[1L], height=px[2L],
                  res=600, compression="lzw", type="cairo", bg="white", antialias="gray")
  tryCatch(draw(), finally=invisible(grDevices::dev.off()))
}

write_figure_svg <- function(stem, width_mm, height_mm, family, draw) {
  path <- paste0(stem, ".svg")
  grDevices::svg(path, width=width_mm/25.4, height=height_mm/25.4,
                 family=family, bg="white")
  tryCatch(draw(), finally=invisible(grDevices::dev.off()))
  lines <- readLines(path, encoding="UTF-8")
  tag <- grep("^<svg ", lines)
  stopifnot(length(tag)==1L)
  lines[tag] <- gsub('(width|height)="([0-9.]+)"', '\\1="\\2pt"', lines[tag])
  writeLines(lines, path, useBytes=TRUE)
}

library(grid)
out <- output_dir
data_dir <- input_dir

tab <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                                 fileEncoding = "UTF-8", na.strings = c("", "NA"))
cohort_est <- tab(file.path(data_dir, "forest_cohorts_v1.tsv"))
pooled_est <- tab(file.path(data_dir, "forest_pooled_v1.tsv"))
frozen <- tab(file.path(data_dir, "matched_ukb_reference.tsv"))
baseline <- tab(file.path(data_dir, "TableS_ext1_baseline_v1.tsv"))
accuracy <- tab(file.path(data_dir, "TableS_ext2_clock_accuracy_v1.tsv"))
era_inputs <- jsonlite::fromJSON(file.path(data_dir, "nhanes_clock_inputs_v1.json"), simplifyVector=FALSE)
one <- function(x) {
  stopifnot(nrow(x) == 1L)
  x
}
key <- function(x) paste(x$cohort, x$group, x$exposure, sep = "/")
equal_num <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10))
two <- function(x) sprintf("%.2f", x)

cohorts <- c("ELSA", "UKHLS", "CHARLS", "NHANES")
systems <- c("pulmonary", "cardiovascular", "renal", "hepatic", "immune",
             "musculoskeletal", "metabolic")
system_names <- setNames(c("Lung", "Cardiovascular", "Renal", "Liver", "Immune",
                          "Musculoskeletal", "Metabolic"), systems)
short_names <- c("Lung", "CV", "Renal", "Liver", "Immune", "MSK", "Metab.")
country <- c("England", "United Kingdom", "China", "United States")
years <- c("2004–05, 2008–09", "2010–12", "2011, 2015", "1988–94, 2007–12")
outcome1 <- c("CLD, asthma (2021–23)", "CLD, asthma, all-cause death",
              "CLD, asthma, all-cause death", "CLRD, influenza/pneumonia")
outcome2 <- c("All-cause death (2012)", "through 2023–24",
              "through 2021–23", "and all-cause deaths (2019)")
nrow0 <- one(baseline[baseline$Characteristic == "Participants, n", ])
n_baseline <- as.integer(gsub(",", "", unlist(nrow0[cohorts], use.names = FALSE)))
accuracy <- accuracy[match(system_names, accuracy$System), ]
ref <- "UK Biobank (reference)"
r_lung <- function(x) as.numeric(sub(" .*", "", x))
lung_r <- vapply(accuracy[1, cohorts, drop = FALSE], r_lung, numeric(1))
ref_r <- r_lung(accuracy[1, ref])
get_counts <- function(x) {
  if (x == "not measured") return(integer())
  hits <- regmatches(x, regexec("([0-9]+(?:/[0-9]+)*) inputs?\\b", x, perl = TRUE))[[1]]
  stopifnot(length(hits) == 2L)
  as.integer(strsplit(hits[2], "/", fixed = TRUE)[[1]])
}
ukb_n <- vapply(accuracy[[ref]], get_counts, integer(1))
coverage <- expand.grid(system = systems, cohort = cohorts, stringsAsFactors = FALSE)
coverage$ukb_total <- ukb_n[match(coverage$system, systems)]
coverage$counts <- lapply(seq_len(nrow(coverage)), function(i) {
  syst <- coverage$system[i]
  who <- coverage$cohort[i]
  s <- get_counts(accuracy[match(syst, systems), who])
  if (who == "NHANES") {
    by_era <- vapply(era_inputs[c("III", "2007-12")], function(a) length(a[[syst]]), integer(1))
    stopifnot(identical(sort(unique(unname(by_era))), sort(s)))
    return(unname(by_era))
  }
  s
})
stopifnot(!anyNA(accuracy$System), all(n_baseline > 0), all(is.finite(c(lung_r, ref_r))))
for (i in seq_len(nrow(coverage))) {
  stopifnot(all(coverage$counts[[i]] > 0 & coverage$counts[[i]] <= coverage$ukb_total[i]))
}

model_b <- cohort_est[cohort_est$model == "B", ]
available <- model_b[is.finite(model_b$exp_HR), ]
stopifnot(!anyDuplicated(key(model_b)), !anyDuplicated(key(frozen)), setequal(key(available), key(frozen)))
checked <- available[match(key(frozen), key(available)), ]
for (field in c("HR", "lo", "hi", "exp_HR")) stopifnot(equal_num(checked[[field]], frozen[[field]]))
stopifnot(equal_num(checked$HR / checked$exp_HR, frozen$ratio),
          all(checked$lo <= checked$HR & checked$HR <= checked$hi))
limits <- c(0.8, 1.25)
n_within <- sum(frozen$ratio >= limits[1] & frozen$ratio <= limits[2])
n_ratio <- nrow(frozen)
log_r <- cor(log(frozen$HR), log(frozen$exp_HR))

panels <- data.frame(system = c("pulmonary", "pulmonary", "cardiovascular", "renal", "hepatic"),
                     outcome = c("copd_type", "asthma", rep("copd_type", 3)),
                     pool = c(rep("full_spirometry", 2), rep("all", 3)))
widths <- c(31.8, 31.8, 26.4, 26.4, 26.4)
panels$left <- 24.2 + c(0, cumsum(head(widths + 3.2, -1)))
panels$right <- panels$left + widths
panels$mid <- (panels$left + panels$right) / 2
panels$min <- c(.9, .9, .85, .85, .85)
panels$max <- c(4.6, 4.6, 1.35, 1.35, 1.35)
pool_b <- pooled_est[pooled_est$model == "B", ]
forest <- lapply(seq_len(nrow(panels)), function(i) {
  p <- panels[i, ]
  rows <- model_b[model_b$exposure == p$system & model_b$group == p$outcome, ]
  summary <- one(pool_b[pool_b$exposure == p$system & pool_b$group == p$outcome & pool_b$pool == p$pool, ])
  ukb <- unique(as.numeric(rows$ukb_full_HR))
  stopifnot(length(ukb) == 1L, is.finite(ukb), nrow(rows) > 0, !anyDuplicated(rows$cohort),
            summary$k == sum(if (p$pool == "full_spirometry") rows$cohort != "CHARLS" else rep(TRUE, nrow(rows))))
  shown <- c(rows$lo, rows$hi, rows$exp_HR, summary$lo, summary$hi, ukb)
  stopifnot(all(is.finite(shown)), all(shown >= p$min & shown <= p$max))
  list(rows = rows, pool = summary, ukb = ukb)
})

columns <- do.call(rbind, lapply(cohorts, function(who) {
  groups <- if (who == "NHANES") c("copd_type", "death") else c("copd_type", "asthma", "death")
  data.frame(cohort = who, group = groups)
}))
hm_left <- 30.3
hm_right <- 179.8
hm_gap <- 2
cw <- (hm_right - hm_left - 3 * hm_gap) / nrow(columns)
columns$x <- hm_left + (seq_len(nrow(columns)) - 1L) * cw + (match(columns$cohort, cohorts) - 1L) * hm_gap
heat <- do.call(rbind, lapply(seq_len(nrow(columns)), function(j) {
  data.frame(cohort = columns$cohort[j], group = columns$group[j], exposure = systems,
             x = columns$x[j], y = 132 + (seq_along(systems) - 1) * 5.7)
}))
heat$ratio <- frozen$ratio[match(key(heat), key(frozen))]
heat$outside <- is.finite(heat$ratio) & (heat$ratio < limits[1] | heat$ratio > limits[2])
stopifnot(all(key(frozen) %in% key(heat)), sum(is.finite(heat$ratio)) == n_ratio, n_ratio == 68L,
          sum(heat$outside) == n_ratio - n_within)

W <- 183
H <- 187
font <- "Times New Roman"
ink <- "#20272B"
slate <- "#4B585E"
rule <- "#D4DCDD"
axis_col <- "#57636B"
clinical <- "#395688"
matched <- "#586B76"
pool_col <- "#34517F"
pooled_band <- "#EEF2F7"
na_fill <- "#E5E8E8"
xmm <- function(x) unit(x, "mm")
ymm <- function(y) unit(1, "npc") - unit(y, "mm")
txt <- function(x, y, text, size = 7, col = ink, face = 1, align = "left") {
  grid.text(text, xmm(x), ymm(y), just = c(align, "centre"),
            gp = gpar(fontfamily = font, fontface = face, fontsize = size, col = col, lineheight = .94))
}
line <- function(x1, y1, x2, y2, col = rule, width = .15) {
  grid.segments(xmm(x1), ymm(y1), xmm(x2), ymm(y2),
                gp = gpar(col = col, lwd = width/.2646, lineend = "butt"))
}
box <- function(x, y, w, h, fill, stroke = NA, width = .14) {
  grid.rect(xmm(x), ymm(y), xmm(w), xmm(h), just = c("left","top"),
            gp = gpar(fill = fill, col = stroke, lwd = width/.2646))
}
diamond <- function(left, mid, right, y, h, fill, stroke = NA, width = .19) {
  grid.polygon(xmm(c(left, mid, right, mid)), ymm(c(y, y-h, y, y+h)),
               gp = gpar(fill = fill, col = stroke, lwd = width/.2646, linejoin = "mitre"))
}
point_ci <- function(left, mid, right, y) {
  line(left, y, right, y, clinical, .29)
  line(left, y-.40, left, y+.40, clinical, .19)
  line(right, y-.40, right, y+.40, clinical, .19)
  grid.circle(xmm(mid), ymm(y), r = xmm(.72), gp = gpar(fill = clinical, col = NA))
}
tag <- function(y, letter) txt(3.12, y, letter, 10.35, face = 2)
group_header <- function(left, right, y, label) {
  txt((left+right)/2, y, label, 7.55, ink, 2, "centre")
  line(left+.3, y+2.25, right-.3, y+2.25, rule, .16)
}

draw_a <- function() {
  tag(4.8, "A")
  gx <- 69; step <- 10.8
  group_header(gx, gx + 7*step, 4.95, "Available clock inputs (UKB totals)")
  txt(3.15, 12, "Cohort", 7.05, face = 2)
  txt(17.95, 10.45, "Country", 7.05, face = 2)
  txt(17.95, 13.55, "Baseline", 6.5, slate)
  txt(49.0, 12.0, "N", 7.05, face = 2, align = "right")
  txt(59.2, 10.45, "Lung clock r", 7.0, face = 2, align = "centre")
  txt(59.2, 13.55, paste("UKB", two(ref_r)), 6.25, slate, align = "centre")
  txt(148.2, 10.45, "Outcomes", 7.05, face = 2)
  txt(148.2, 13.55, "End of follow-up", 6.5, slate)
  for (k in seq_along(systems)) {
    x <- gx + (k-.5)*step
    txt(x, 10.45, short_names[k], 6.6, face = 2, align = "centre")
    txt(x, 13.45, sprintf("(%d)", ukb_n[k]), 6.3, slate, align = "centre")
  }
  line(3.15, 15.65, 179.8, 15.65)
  for (i in seq_along(cohorts)) {
    y <- 19.9 + (i-1)*7.25
    txt(3.15, y, cohorts[i], 7.22, face = 2)
    txt(17.95, y-1.46, country[i], 6.85)
    txt(17.95, y+1.58, years[i], 6.40, slate)
    txt(49.0, y, formatC(n_baseline[i], format = "d", big.mark = ","), 7.1, align = "right")
    txt(59.2, y, two(lung_r[i]), 7.35, clinical, 2, "centre")
    txt(148.2, y-1.46, outcome1[i], 6.42)
    txt(148.2, y+1.58, outcome2[i], 6.42, slate)
    for (k in seq_along(systems)) {
      idx <- which(coverage$cohort == cohorts[i] & coverage$system == systems[k])
      n <- unique(coverage$counts[[idx]])
      full <- length(n) > 0 && all(n == ukb_n[k])
      fill <- if (!length(n)) na_fill else if (full) "#D4DFEE" else "#F0F3F6"
      label <- if (!length(n)) "–" else paste(n, collapse = ", ")
      if (cohorts[i] == "CHARLS" && systems[k] == "pulmonary") label <- paste(label, "(PEF)")
      x <- gx + (k-1)*step
      box(x+.15, y-2.53, step-.30, 5.06, fill)
      txt(x+step/2, y, label, 6.72, if (full) clinical else slate, if (full) 2 else 1, "centre")
    }
  }
  line(3.15, 46, 179.8, 46)
}

draw_b <- function() {
  tag(51.3, "B")
  group_header(panels$left[1], panels$right[2], 51.3, "Lung")
  for (i in 3:5) group_header(panels$left[i], panels$right[i], 51.3, system_names[panels$system[i]])
  txt(3.15, 56.8, "Cohort", 7.05, face = 2)
  line(3.15, 59.0, 179.8, 59.0)
  ys <- 62.6 + (0:3)*5.8
  for (j in seq_along(cohorts)) txt(3.15, ys[j], cohorts[j], 7.10)
  line(3.15, 83.15, 179.8, 83.15)
  txt(3.15, 86.35, "UKB full clock", 6.90, slate)
  box(3.15, 89.6, 176.65, 9.9, pooled_band)
  txt(4.0, 93.25, "Pooled HR", 7.10, clinical, 2)
  txt(4.0, 96.72, "(95% CI)", 6.40, slate)
  for (i in seq_len(nrow(panels))) {
    p <- panels[i, ]; f <- forest[[i]]
    left <- p$left + .95; right <- p$right - .95
    sx <- function(hr) left + (log(hr)-log(p$min))/(log(p$max)-log(p$min)) * (right-left)
    ticks <- if (p$system == "pulmonary") c(1,1.5,2,3,4) else c(.9,1,1.1,1.2,1.3)
    txt(p$mid, 56.8, if (p$outcome == "asthma") "Asthma" else "COPD-type", 6.95, align = "centre")
    line(sx(1), 59.65, sx(1), 82.70, "#BAC4C8", .16)
    line(sx(1), 90.15, sx(1), 94.70, "#BAC4C8", .16)
    for (j in seq_along(cohorts)) {
      row <- f$rows[f$rows$cohort == cohorts[j], ]
      y <- ys[j]
      if (!nrow(row)) {
        txt(p$mid, y, "–", 7.1, "#859197", align = "centre")
      } else {
        point_ci(sx(row$lo), sx(row$HR), sx(row$hi), y-.73)
        ex <- sx(row$exp_HR)
        diamond(ex-.74, ex, ex+.74, y+.94, .69, "white", matched, .21)
      }
    }
    ux <- sx(f$ukb)
    line(ux, 85.26, ux, 87.43, matched, .28)
    txt(ux+1.5, 86.35, two(f$ukb), 6.62, slate)
    z <- f$pool
    diamond(sx(z$lo), sx(z$HR), sx(z$hi), 92.40, 1.06, pool_col)
    txt(p$mid, 97.05, sprintf("%.2f (%.2f–%.2f)", z$HR, z$lo, z$hi), 7.04, clinical, 2, "centre")
    line(left, 100.8, right, 100.8, axis_col, .16)
    for (t in ticks) {
      line(sx(t), 100.8, sx(t), 101.55, axis_col, .14)
      txt(sx(t), 103.55, format(t, trim=TRUE), 6.1, slate, align = "centre")
    }
  }
  txt(3.15, 109.5, "Hazard ratio per SD (log scale)", 6.95, face = 2)
  point_ci(65.7, 68.2, 70.7, 109.5)
  txt(72.3, 109.5, "Observed HR (95% CI)", 6.85)
  diamond(108.9,109.7,110.5,109.5,.73,"white",matched,.22)
  txt(112.0, 109.5, "Matched UKB HR", 6.85)
  line(3.15, 113.1, 179.8, 113.1)
}

stops <- log(c(.5,.8,1,1.25,2))
colours <- c("#3066A0", "#C6DDE9", "#F8F8F6", "#F6D7C9", "#B7453C")
ramps <- lapply(seq_len(length(stops)-1L), function(j) grDevices::colorRamp(colours[j+c(0L,1L)], space="Lab"))
ratio_colour <- function(r) {
  lr <- max(stops[1], min(tail(stops,1), log(r)))
  j <- findInterval(lr, stops, all.inside=TRUE)
  rgb <- ramps[[j]]((lr-stops[j])/(stops[j+1]-stops[j]))
  grDevices::rgb(rgb[1],rgb[2],rgb[3],maxColorValue=255)
}
label_colour <- function(fill) {
  rgb <- as.numeric(col2rgb(fill))/255
  linear <- ifelse(rgb <= .04045, rgb/12.92, ((rgb+.055)/1.055)^2.4)
  if (sum(linear*c(.2126,.7152,.0722)) < .178) "white" else ink
}
draw_c <- function() {
  tag(118, "C")
  for (who in cohorts) {
    sub <- columns[columns$cohort == who, ]
    group_header(min(sub$x), max(sub$x)+cw, 118, who)
  }
  txt(3.15, 126.05, "Organ system", 7.10, face=2)
  for (j in seq_len(nrow(columns))) {
    cc <- columns[j, ]
    label <- if (cc$group == "asthma") c("Incident","asthma") else if (cc$group == "death") c("All-cause","death") else if (cc$cohort=="NHANES") c("CLRD","death") else c("Incident","CLD")
    txt(cc$x+cw/2, 124.15, label[1], 6.62, align="centre")
    txt(cc$x+cw/2, 127.05, label[2], 6.62, align="centre")
  }
  for (k in seq_along(systems)) txt(3.15, 132+(k-1)*5.7, system_names[k], 7.15)
  for (i in seq_len(nrow(heat))) {
    h <- heat[i, ]; missing <- !is.finite(h$ratio)
    fill <- if (missing) na_fill else ratio_colour(h$ratio)
    border <- if (h$outside) "#314754" else NA
    box(h$x+.12, h$y-2.58, cw-.24, 5.16, fill, border, .28)
    txt(h$x+cw/2, h$y, if (missing) "–" else two(h$ratio),
        7.0, if (missing) "#819097" else label_colour(fill), if (h$outside) 2 else 1, "centre")
  }
  txt(3.15, 173.45, "Observed / expected HR", 7.10, face=2)
  # Sample the same colour function as the cells, at equally spaced log ratios.
  gx <- 3.15; gw <- 47
  gradient_stops <- seq(min(stops), max(stops), length.out = 65L)
  gradient <- grid::linearGradient(vapply(exp(gradient_stops), ratio_colour, character(1)),
                                    stops = (gradient_stops-min(stops))/diff(range(stops)),
                                    x1=0, y1=.5, x2=1, y2=.5)
  box(gx, 177.1, gw, 2.05, gradient)
  for (r in c(.5,.8,1,1.25,2)) {
    x <- gx + gw*(log(r)-min(stops))/diff(range(stops))
    line(x,179.2,x,179.9,axis_col,.13)
    txt(x,182.05,format(r,trim=TRUE),6.32,slate,align=if (r==.5) "left" else if (r==2) "right" else "centre")
  }
  box(63.2,173.25,4.2,2.65,"white","#314754",.28)
  txt(69.2,174.58,"Outside 0.80–1.25",6.85)
  box(63.2,178.50,4.2,2.65,na_fill)
  txt(65.3,179.83,"–",6.6,slate,align="centre")
  txt(69.2,179.83,"Unavailable",6.85)
  txt(122.0,174.58,sprintf("%d of %d within range",n_within,n_ratio),7.25,ink,2)
  txt(122.0,179.83,sprintf("Pearson r for log HRs = %.2f",log_r),6.85,slate)
}

draw <- function() {
  grid.newpage()
  draw_a(); draw_b(); draw_c()
}
stem <- file.path(out,"Figure7")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=font,bg="white")
draw(); invisible(dev.off())
write_figure_rasters(stem, W, H, draw)
write_figure_svg(stem, W, H, font, draw)

cat(sprintf("Plotted %d cohort estimates in B and %d observed/expected ratios in C. %d within range, log-HR r = %.6f.\n",
            sum(vapply(forest,function(f)nrow(f$rows),integer(1))),n_ratio,n_within,log_r))
