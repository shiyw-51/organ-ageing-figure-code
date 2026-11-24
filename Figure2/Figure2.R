# Figure2
# Rscript --vanilla Figure2/Figure2.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure2")
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
out_dir <- output_dir

tab <- function(file) read.delim(file.path(input_dir,file), check.names=FALSE,
                                 stringsAsFactors=FALSE, na.strings=c("", "NA"), fileEncoding="UTF-8")
results <- tab("Figure2_data.tsv")
accuracy <- tab("clock_accuracy.tsv")

systems <- c("pulmonary","cardiovascular","renal","hepatic","immune","musculoskeletal","metabolic")
system_labels <- c("Lung","Cardiovascular","Renal","Liver","Immune","Musculoskeletal","Metabolic")
endpoints <- c("copd","asthma","bronchiectasis","ild","ipf","pneumoconiosis",
               "pvd","pe","ph","osa","lung_cancer","pneumonia")
endpoint_labels <- c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD (J84.1)","Pneumoconiosis",
                     "Pulmonary vascular disease","Pulmonary embolism","Pulmonary hypertension",
                     "Sleep apnoea","Lung cancer","Pneumonia")
groups <- c("Airway","Parenchymal","Vascular","Other")
group_fills <- c("#E3E7F0","#DEEEE9","#F4E4DF","#EEE7DE")
nested <- c("ipf","pe","ph")
stopifnot(!anyDuplicated(results[c("layer","endpoint","system")]),
          !anyDuplicated(accuracy[c("layer","system")]))
panel_data <- function(layer, comparator) {
  cols <- c(systems,comparator)
  d <- results[results$layer==layer & results$system %in% cols & results$endpoint %in% endpoints,]
  idx <- match(as.vector(outer(endpoints,cols,paste,sep=":")), paste(d$endpoint,d$system,sep=":"))
  stopifnot(length(idx)==96L, nrow(d)==96L, !anyNA(idx), all(is.finite(d$HR) & d$HR>0))
  mat <- function(key) matrix(d[[key]][idx],nrow=length(endpoints),dimnames=list(endpoints,cols))
  qc <- mat("q_fdr")
  stopifnot(all(is.finite(qc[,systems])), all(qc[,systems]>=0 & qc[,systems]<=1),
            all(is.na(qc[,comparator])))
  clock <- accuracy[accuracy$layer==layer, ]
  r <- clock$r_T[match(systems,clock$system)]
  events <- mat("events")[,"pulmonary"]
  stopifnot(all(is.finite(r) & r>=0 & r<=.8), all(is.finite(events) & events>=0 & events==floor(events)))
  list(hr=mat("HR"), q=qc, r=r, events=events)
}
cl <- panel_data("clinical","phenoage")
pr <- panel_data("protein","organismal")

W <- 183
H <- 115
font <- "Times New Roman"
ink <- "#20272B"
muted <- "#4B585E"
axis_col <- "#616C70"
rule <- "#D8DEDF"
blue <- "#395688"
teal <- "#168372"
xx <- function(v) unit(v,"mm")
yy <- function(v) unit(H-v,"mm")
label <- function(x,y,s,size=7,col=ink,face=1,just="left",rot=0) {
  grid.text(s, x=xx(x), y=yy(y), just=c(just,"centre"), rot=rot,
            gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col))
}
line_mm <- function(x,y,col=rule,width=.15) {
  grid.lines(x=xx(x),y=yy(y),gp=gpar(col=col,lwd=width/.2646,lineend="butt"))
}
rect_mm <- function(x,y,w,h,fill) {
  grid.rect(x=xx(x),y=yy(y),width=xx(w),height=xx(h),just=c("left","top"),
            gp=gpar(col=NA,fill=fill))
}

# Interpolate the original diverging palette in Lab on a symmetric log-HR scale.
limit <- log(1.6)
palette <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B")
ramp <- grDevices::colorRamp(palette,space="Lab")
heat_colour <- function(hr) {
  pos <- pmin(1,pmax(0,(log(hr)+limit)/(2*limit)))
  rgb <- ramp(pos)
  grDevices::rgb(rgb[,1],rgb[,2],rgb[,3],maxColorValue=255)
}
cell_ink <- function(fill) {
  rgb <- grDevices::col2rgb(fill)/255
  lin <- ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
  lum <- drop(c(.2126,.7152,.0722) %*% lin)
  ifelse((lum+.05)/.05 >= 1.05/(lum+.05),"#000000","#FFFFFF")
}

x_cl <- 36.1
x_pr <- 110.15
cw <- 6.76
comp_gap <- 1.75
comp_w <- 7.65
event_offset <- 68.25
data_top <- 28.6
rh <- 4.65
group_gap <- 1.25
row_top <- data_top + (seq_along(endpoints)-1)*rh + (seq_along(endpoints)-1)%/%3*group_gap
row_mid <- row_top+rh/2
data_bottom <- tail(row_top,1)+rh

outcome_axis <- function() {
  label(8.35,24.9,"Respiratory outcome",7.25,face=2)
  for (i in seq_along(groups)) {
    y <- row_top[3*i-2]
    rect_mm(3.1,y+.09,3.75,3*rh-.18,group_fills[i])
    label(4.975,y+1.5*rh,groups[i],6.13,muted,just="centre",rot=90)
  }
  # Indent the three nested outcomes beneath their parent diagnosis.
  for (branch in list(c(4L,5L),c(7L,8L,9L))) {
    line_mm(c(8.65,8.65),c(row_mid[branch[1]]+1.45,row_mid[tail(branch,1)]),"#A2ADAF",.13)
    for (k in branch[-1]) line_mm(c(8.65,9.9),rep(row_mid[k],2),"#A2ADAF",.13)
  }
  for (i in seq_along(endpoints)) {
    sub <- endpoints[i] %in% nested
    label(if (sub) 10.55 else 8.35, row_mid[i],endpoint_labels[i],
          if (sub) 6.65 else 6.85,if (sub) muted else ink,if (sub) 3 else 1)
  }
}

accuracy_bars <- function(x0,r,colour) {
  top <- 11.85
  bottom <- 20.55
  height <- bottom-top
  for (v in c(0,.4,.8)) {
    y <- bottom-height*v/.8
    line_mm(c(x0-.85,x0+length(systems)*cw),c(y,y),if (v==0) "#A4AFB1" else "#E8ECEC",.12)
    label(x0-1.35,y,if(v==0) "0" else sprintf("%.1f",v),5.9,muted,just="right")
  }
  for (j in seq_along(r)) {
    h <- height*r[j]/.8
    rect_mm(x0+(j-.5)*cw-.27*cw,bottom-h,.54*cw,h,colour)
  }
}

association_panel <- function(x0,tag_x,tag,title,accuracy_label,comparator,colour,d) {
  label(tag_x,4.8,tag,10.35,ink,2)
  label(tag_x+5.25,4.9,title,9.05,ink,2)
  label(x0,9.0,accuracy_label,6.25,muted)
  accuracy_bars(x0,d$r,colour)
  xc <- x0+length(systems)*cw+comp_gap
  event_right <- x0+event_offset
  label(x0+length(systems)*cw/2,24.9,"Organ systems",7.05,face=2,just="centre")
  label(xc+comp_w/2,24.9,"Whole body",6.5,face=2,just="centre")
  label(event_right-3.4,23.45,"Events",7.0,face=2,just="centre")
  label(event_right-3.4,26.0,"(lung model)",5.95,muted,just="centre")

  xs <- c(x0+(seq_along(systems)-1)*cw,xc)
  ws <- c(rep(cw,length(systems)),comp_w)
  for (j in seq_along(xs)) {
    fills <- heat_colour(d$hr[,j])
    inks <- cell_ink(fills)
    for (i in seq_along(endpoints)) {
      rect_mm(xs[j]+.09,row_top[i]+.09,ws[j]-.18,rh-.18,fills[i])
      q <- d$q[i,j]
      number <- sprintf("%.2f%s",d$hr[i,j],if (is.finite(q) && q<.05) "*" else "")
      label(xs[j]+ws[j]/2,row_mid[i],number,6.7,inks[i],just="centre")
    }
  }
  all_names <- c(system_labels,comparator)
  for (j in seq_along(xs))
    label(xs[j]+ws[j]/2+.16,data_bottom+1.6,all_names[j],6.8,just="right",rot=48)
  for (i in seq_along(endpoints))
    label(event_right,row_mid[i],formatC(d$events[i],format="d",big.mark=","),6.95,muted,just="right")
}

shared_legend <- function() {
  label(3.1,105.25,"Hazard ratio per 1 SD (model B)",7.3,face=2)
  x <- 3.1; y <- 108.15; w <- 55.8; h <- 2.0
  edges <- seq(-limit,limit,length.out=145L)
  mid <- (head(edges,-1L)+tail(edges,-1L))/2
  fills <- heat_colour(exp(mid))
  for (i in seq_along(fills))
    rect_mm(x+(i-1)*w/length(fills),y,w/length(fills)+.009,h,fills[i])
  vals <- c(1/1.6,.8,1,1.25,1.6)
  at <- x+w*(log(vals)+limit)/(2*limit)
  labels <- c("\u22640.625","0.80","1.00","1.25","\u22651.60")
  for (j in seq_along(at)) {
    line_mm(rep(at[j],2),c(y+h,y+h+.65),axis_col,.12)
    align <- if (j==1L) "left" else if (j==length(at)) "right" else "centre"
    label(at[j],y+h+2.45,labels[j],6.15,muted,just=align)
  }
  label(68.6,106.75,"* FDR q < 0.05 for the seven organ systems.",7.0,ink)
  label(68.6,111.05,"Separate correction for nine main outcomes and three nested subtypes.",6.75,muted)
}

draw_figure <- function() {
  grid.newpage()
  rect_mm(0,0,W,H,"white")
  outcome_axis()
  label(8.35,14.35,"Clock accuracy",6.9,ink,2)
  label(8.35,18.15,"Healthy test-set Pearson r",6.40,muted)
  association_panel(x_cl,3.1,"A","Clinical age acceleration","Mean of sex-specific correlations",
                    "PhenoAge",blue,cl)
  association_panel(x_pr,108.1,"B","Proteomic age acceleration","Pooled correlation",
                    "Organismal",teal,pr)
  shared_legend()
}
stem <- file.path(out_dir,"Figure2")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=font,bg="white")
draw_figure()
invisible(dev.off())
write_figure_rasters(stem, W, H, draw_figure)
write_figure_svg(stem, W, H, font, draw_figure)
cat("Figure 2 exported with 192 hazard ratios and 24 pulmonary-model event counts.\n")
