# Figure 1B–E. Clock performance, cross-layer correlations and baseline disease contrasts.
# All drawing coordinates are millimetres from the upper-left corner.
suppressPackageStartupMessages(library(grid))

# Rscript --vanilla Figure1/Figure1.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure1")
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

out_dir <- output_dir
data_dir <- input_dir
validity_file <- file.path(data_dir, "clock_validity_v1.tsv")
stem <- file.path(out_dir, "Figure1")
read_tab <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
systems <- c("pulmonary", "cardiovascular", "renal", "hepatic", "immune", "musculoskeletal", "metabolic")
system_name <- c("Lung", "Cardiovascular", "Renal", "Liver", "Immune", "Musculoskeletal", "Metabolic")
short_name <- c("Lung", "CV", "Renal", "Liver", "Imm.", "MSK", "Met.")
disease_labels <- c(
  copd="COPD", asthma="Asthma", ild="ILD", hypertension="Hypertension", ihd="Ischaemic heart disease",
  stroke="Stroke", heart_failure="Heart failure", atrial_fibrillation="Atrial fibrillation",
  pad="Peripheral artery disease", aortic_aneurysm="Aortic aneurysm", ckd="Chronic kidney disease",
  cirrhosis="Cirrhosis", chronic_liver_disease="Chronic liver disease",
  rheumatoid_arthritis="Rheumatoid arthritis", osteoporosis="Osteoporosis",
  osteoarthritis="Osteoarthritis", diabetes="Diabetes"
)

accuracy <- read_tab(file.path(data_dir, "fig1_accuracy.tsv"))
cross <- read_tab(file.path(data_dir, "fig1_crosslayer_r.tsv"))
rownames(cross) <- cross[[1]]
correlation <- as.matrix(cross[systems, systems])
all_validity <- read_tab(validity_file)
# Each row of the main figure represents a disease with an assigned system.
# Within each system, follow the source table's clinical ordering.
target <- unique(all_validity[all_validity$matched %in% systems, c("disease","matched")])
target <- target[order(match(target$matched, systems), seq_len(nrow(target))), ]
diseases <- target$disease
undefined <- setdiff(diseases, names(disease_labels))
if (length(undefined)) stop("Missing display labels for ", paste(undefined, collapse=", "))
disease_name <- unname(disease_labels[diseases])
validity <- all_validity[all_validity$disease %in% diseases & all_validity$system %in% systems, ]
validity$outlined <- tolower(as.character(validity$is_matched)) == "true"
grid_keys <- c("layer","system","disease")
expected <- 2L * length(systems) * length(diseases)
stopifnot(!anyDuplicated(diseases), length(diseases) > 0L, nrow(accuracy) == 21L,
          all(is.finite(accuracy$r)), !anyDuplicated(accuracy[c("system","layer","sex")]),
          all(dim(correlation) == c(7,7)), all(is.finite(correlation)), nrow(validity) == expected,
          !anyDuplicated(validity[grid_keys]), all(is.finite(validity$diff)),
          all(validity$outlined == (validity$matched == validity$system)),
          setequal(validity$layer, c("clinical","protein")))

W <- 183
validity_top <- 137.4
validity_rh <- 3.80
validity_bottom <- validity_top + length(diseases)*validity_rh
canvas_top <- 56
H <- validity_bottom + 12 - canvas_top
family <- "Times New Roman"
ink <- "#20272B"; text2 <- "#46535A"; muted <- "#57666A"; rule <- "#D7DFDF"
clinical <- "#395688"; protein <- "#168372"
xx <- function(x) unit(x, "mm")
yy <- function(y) unit(H + canvas_top - y, "mm")
label <- function(x, y, value, size = 7, col = ink, face = 1, align = "left", angle = 0) {
  grid.text(value, x = xx(x), y = yy(y), just = c(align, "centre"), rot = angle,
            gp = gpar(fontfamily = family, fontsize = size, fontface = face, col = col, lineheight = 1.04))
}
line_mm <- function(x, y, col = rule, width = .16, end = "round", dash = 1) {
  grid.lines(xx(x), yy(y), gp = gpar(col = col, lwd = width / .2646, lineend = end, linejoin = "round", lty = dash))
}
box_mm <- function(x, y, w, h, fill = NA, stroke = NA, width = .15, radius = 0) {
  args <- list(x = xx(x), y = yy(y), width = xx(w), height = xx(h), just = c("left", "top"),
               gp = gpar(fill = fill, col = stroke, lwd = width / .2646))
  if (radius > 0) do.call(grid.roundrect, c(args, list(r = xx(radius)))) else do.call(grid.rect, args)
}
polygon_mm <- function(x, y, fill, stroke = NA, width = .13) {
  grid.polygon(xx(x), yy(y), gp = gpar(fill = fill, col = stroke, lwd = width / .2646, linejoin = "round"))
}
dot <- function(x, y, radius, fill, stroke = NA, width = .15) {
  grid.circle(xx(x), yy(y), r = xx(radius), gp = gpar(fill = fill, col = stroke, lwd = width / .2646))
}
diamond <- function(x, y, r, fill) polygon_mm(x + c(0, r, 0, -r), y + c(-r, 0, r, 0), fill)
tag <- function(x, y, letter) label(x, y, letter, size = 10.3, col = "#10181B", face = 2)
fmt2 <- function(x) sprintf("%.2f", ifelse(abs(x) < .005, 0, x))
# Interpolate on the published scale; clipped numeric values are printed in full.
colour_fun <- function(stops, colours) {
  ramp <- grDevices::colorRamp(colours, space = "Lab")
  positions <- seq(0, 1, length.out = length(stops))
  function(x) {
    position <- approx(stops, positions, xout = x, rule = 2)$y
    rgb <- ramp(position)
    grDevices::rgb(rgb[,1], rgb[,2], rgb[,3], maxColorValue = 255)
  }
}
cor_fill <- colour_fun(c(-.4, -.2, 0, .2, .4), c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"))
val_fill <- colour_fun(c(-1, -.3, 0, .3, 1, 2), c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B"))
cell_text <- function(fill) {
  s <- as.numeric(col2rgb(fill))/255
  lum <- sum(c(.2126, .7152, .0722) * ifelse(s <= .04045, s/12.92, ((s+.055)/1.055)^2.4))
  if (lum < .18) "white" else "#202427"
}
colour_key <- function(x, y, w, h, limits, pal, ticks, labels) {
  z <- seq(limits[1], limits[2], length.out = 65)
  for (i in seq_len(64)) box_mm(x + (i-1)*w/64, y, w/64+.015, h, pal(mean(z[i:(i+1)])))
  for (i in seq_along(ticks)) {
    at <- x + w*(ticks[i]-limits[1])/diff(limits)
    line_mm(c(at, at), c(y+h+.1, y+h+.65), muted, .12)
    label(at, y+h+2.35, labels[i], 5.55, text2, align = "centre")
  }
}

# B and C use identical frame dimensions and exactly the same y coordinates.
x_left <- 26.9; x_right <- 117.4; grid_width <- 60.9
cw <- grid_width/7
top_mid <- 71.6; mid_rh <- 5.6
mid_y <- top_mid + (seq_along(systems)-.5)*mid_rh
bottom_mid <- top_mid + length(systems)*mid_rh
accuracy_panel <- function() {
  tag(3.8,60.5,"B")
  x <- function(r) x_left + (r-.19)/(.81-.19)*grid_width
  for (i in seq_along(systems)) {
    y <- mid_y[i]
    line_mm(c(x_left,x_left+grid_width),c(y,y),"#E4E8E8",.13)
    label(x_left-1.3,y,system_name[i],7.05,text2,align="right")
    rows <- accuracy[accuracy$system == systems[i], ]
    a <- rows$r[rows$layer=="Clinical" & rows$sex=="female"]
    b <- rows$r[rows$layer=="Clinical" & rows$sex=="male"]
    p <- rows$r[rows$layer=="Protein" & rows$sex=="pooled"]
    stopifnot(length(a)==1,length(b)==1,length(p)==1)
    line_mm(x(c(a,b)),c(y,y),"#A7B3C8",.38)
    dot(x(a),y,.89,clinical)
    dot(x(b),y,.89,"white",clinical,.29)
    diamond(x(p),y,1.15,protein)
  }
  line_mm(c(x_left,x_left+grid_width),rep(bottom_mid,2),"#697579",.19)
  for (k in seq(.2,.8,by=.2)) {
    line_mm(rep(x(k),2),c(bottom_mid,bottom_mid+.85),"#697579",.16)
    label(x(k),bottom_mid+2.85,sprintf("%.1f",k),6.25,text2,align="centre")
  }
  dot(x_left-9.0,66.3,.81,clinical)
  label(x_left-7.3,66.3,"Clinical, women",6.15)
  dot(x_left+12.3,66.3,.82,"white",clinical,.27)
  label(x_left+14.1,66.3,"Clinical, men",6.15)
  diamond(x_left+31.4,66.3,1.00,protein)
  label(x_left+33.2,66.3,"Protein, pooled",6.15)
  label(x_left+grid_width/2,117.6,"Predicted vs chronological age (Pearson r)",6.9,align="centre")
  label(x_left+grid_width/2,121.2,"Healthy participants in the test set",6.15,muted,align="centre")
}
heat_cell <- function(x, y, w, h, value, fill, outlined = FALSE, size = 6.0) {
  gap <- .18
  box_mm(x+gap,y+gap,w-2*gap,h-2*gap,fill)
  if (outlined) box_mm(x+gap+.13,y+gap+.13,w-2*gap-.26,h-2*gap-.26,NA,"#303940",.25)
  label(x+w/2,y+h/2,fmt2(value),size,cell_text(fill),align="centre")
}
corr_panel <- function() {
  tag(94.3,60.5,"C")
  label(x_right+grid_width/2,62.2,"Protein age acceleration",7.05,align="centre")
  label(94.2,(top_mid+bottom_mid)/2,"Clinical age acceleration",6.45,angle=90,align="centre")
  for (i in seq_along(systems)) {
    label(x_right-1.1,mid_y[i],system_name[i],6.8,text2,align="right")
    label(x_right+(i-.5)*cw,68.9,short_name[i],6.75,text2,align="centre")
    for (j in seq_along(systems)) {
      v <- correlation[i,j]
      heat_cell(x_right+(j-1)*cw,top_mid+(i-1)*mid_rh,cw,mid_rh,v,cor_fill(v),i==j,6.8)
    }
  }
  label(x_right,116.2,"n = 19,448",6.65,text2)
  box_mm(x_right,118.9,2.05,1.85,"#F7F7F7","#303940",.21)
  label(x_right+3.05,119.9,"Same system",6.1,muted)
  label(x_right+25.7,116.55,"Pearson r",6.25,text2)
  colour_key(x_right+38,115.65,22.5,1.65,c(-.4,.4),cor_fill,c(-.4,0,.4),c("\u22120.4","0","0.4"))
}

validity_panel <- function(layer, x0, tag_x, letter) {
  tag(tag_x,128.6,letter)
  # The layer is an axis label, not a panel title.
  offset <- if (layer=="clinical") "Clinical systems" else "Proteomic systems"
  label(x0+grid_width/2,129.8,offset,7.55,align="centre")
  top <- validity_top; rh <- validity_rh
  for (j in seq_along(systems)) label(x0+(j-.5)*cw,134.7,short_name[j],6.80,text2,align="centre")
  d <- validity[validity$layer == layer, ]
  for (i in seq_along(diseases)) {
    y <- top+(i-1)*rh
    label(x0-1.05,y+rh/2,disease_name[i],6.72,text2,align="right")
    data_i <- d[d$disease==diseases[i], ]
    data_i <- data_i[match(systems, data_i$system), ]
    for (j in seq_along(systems)) {
      v <- data_i$diff[j]
      heat_cell(x0+(j-1)*cw,y,cw,rh,v,val_fill(v),data_i$outlined[j],6.78)
    }
  }
}
draw_figure <- function() {
  grid.newpage()
  box_mm(0,canvas_top,W,H,"white")
  accuracy_panel()
  corr_panel()
  line_mm(c(3.8,178.9),c(124.5,124.5),rule,.18)
  validity_panel("clinical",x_left,3.8,"D")
  validity_panel("protein",x_right,94.3,"E")
  label(4.4,validity_bottom+4.5,"Age- and sex-adjusted difference in age acceleration (SD)",6.35,text2)
  label(4.4,validity_bottom+8.1,"Participants with vs without each baseline disease",6.05,muted)
  box_mm(81.3,validity_bottom+4.0,2.6,2.3,"#F7F7F7","#303940",.23)
  label(85.4,validity_bottom+5.3,"Disease-matched system",6.2,text2)
  label(133.3,validity_bottom+4.9,"SD",6.3,text2)
  colour_key(140,validity_bottom+4.1,37.6,1.8,c(-1,2),val_fill,c(-1,0,1,2),c("\u2264\u22121","0","1","\u22652"))
}
stopifnot(dir.exists(out_dir))
grDevices::cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=family,bg="white")
draw_figure()
invisible(dev.off())
write_figure_rasters(stem, W, H, draw_figure)
write_figure_svg(stem, W, H, family, draw_figure)
cat(sprintf("B: %d estimates. C: %d correlations. D/E: %d diseases, %d contrasts.\n",
            nrow(accuracy), length(correlation), length(diseases), nrow(validity)))
