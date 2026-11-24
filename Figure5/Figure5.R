# Figure5
# Rscript --vanilla Figure5/Figure5.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure5")
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
data_dir <- input_dir

read_tsv <- function(path) read.delim(path,stringsAsFactors=FALSE,check.names=FALSE,
                                      fileEncoding="UTF-8",na.strings=c("","NA"))
fits <- read_tsv(file.path(data_dir,"lung_interaction_results_v1.tsv"))
boot <- read_tsv(file.path(data_dir,"ic_bootstrap_summary_v1.tsv"))

layers <- c("clinical","protein")
systems <- c("cardiovascular","renal","hepatic","immune","musculoskeletal","metabolic")
system_names <- c("Cardiovascular","Renal","Liver","Immune","Musculoskeletal","Metabolic")
endpoints <- c("copd","asthma","bronchiectasis","ild","ipf","pneumoconiosis",
               "pvd","pe","ph","osa","lung_cancer","pneumonia")
endpoint_names <- c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD (J84.1)","Pneumoconiosis",
                    "Pulmonary vascular disease","Pulmonary embolism","Pulmonary hypertension",
                    "Sleep apnoea","Lung cancer","Pneumonia")
examples <- c("copd","pvd","ph","pneumonia")
display_names <- c(copd="COPD",pvd="PVD",ph="PH",pneumonia="Pneumonia",asthma="Asthma")

stopifnot(!anyDuplicated(fits[c("layer","endpoint","other")]),
          !anyDuplicated(boot$endpoint),setequal(boot$endpoint,examples),
          all(boot$n_boot_ok==1000L),all(fits$other %in% systems),
          all(fits$layer %in% layers),all(fits$endpoint %in% endpoints))
shown <- fits[fits$endpoint!="pneumoconiosis",]
stopifnot(nrow(shown)==2L*6L*11L,
          all(is.finite(as.matrix(shown[c("HR_int","HR_int_lo","HR_int_hi","F00","F10","F01",
                                           "F11","IC","IC_lo","IC_hi")]))),
          all(shown$HR_int_lo>0 & shown$HR_int_lo<=shown$HR_int & shown$HR_int<=shown$HR_int_hi),
          all(shown$IC_lo<=shown$IC & shown$IC<=shown$IC_hi),
          max(abs(shown$IC-(shown$F11-shown$F10-shown$F01+shown$F00)))<1e-8)
shown$mult_mark <- with(shown,HR_int_lo>1 | HR_int_hi<1)
shown$add_mark <- with(shown,IC_lo>0 | IC_hi<0)
excess <- shown$F11-shown$F00
shown$share <- ifelse(excess>.05,100*shown$IC/excess,NA_real_)

layer_col <- c(clinical="#395688",protein="#168372")
font <- "Times New Roman"
ink <- "#20272B"; muted <- "#4B585E"; rule <- "#D8DEDF"; axis_col <- "#5B666B"
grey_na <- "#E3E6E6"
group_names <- c("Airway","Parenchymal","Vascular","Other")
group_fill <- c("#E3E7F0","#DEEEE9","#F4E4DF","#EEE7DE")
W <- 183; H <- 244
right_tag <- 105.55
xx <- function(x) unit(x,"mm")
yy <- function(y) unit(H-y,"mm")
label <- function(x,y,s,size=7,col=ink,face=1,just="left",rot=0) {
  grid.text(s,x=xx(x),y=yy(y),just=c(just,"centre"),rot=rot,
            gp=gpar(fontfamily=font,fontface=face,fontsize=size,col=col,lineheight=.94))
}
line_mm <- function(x,y,col=rule,width=.15) {
  grid.lines(xx(x),yy(y),gp=gpar(col=col,lwd=width/.2646,lineend="butt",linejoin="round"))
}
box_mm <- function(x,y,w,h,fill) {
  grid.rect(xx(x),yy(y),width=xx(w),height=xx(h),just=c("left","top"),gp=gpar(fill=fill,col=NA))
}
point_mm <- function(x,y,r,col,filled=TRUE) {
  grid.circle(xx(x),yy(y),r=xx(r),gp=gpar(col=col,fill=if(filled) col else "white",lwd=.19/.2646))
}
header <- function(x,y,tag,title) {
  label(x,y,tag,10.35,face=2)
  label(x+5.25,y+.05,title,8.75,face=2)
}
scale_to <- function(x,domain,range) range[1]+(x-domain[1])*diff(range)/diff(domain)
number <- function(x,digits) sub("^-","−",sprintf(paste0("%.",digits,"f"),x))
ci_text <- function(v,lo,hi) sprintf("%s (%s–%s)",number(v,2),number(lo,2),number(hi,2))
one <- function(d) {stopifnot(nrow(d)==1L); d[1L,,drop=FALSE]}

# Both layers use the same neutral colour at the null and the same range for each measure.
palette <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B")
ramp <- grDevices::colorRamp(palette,space="Lab")
fill_value <- function(v,metric) {
  scaled <- if (metric=="HR_int") log(v)/log(1.12) else v/.45
  rgb <- ramp((pmax(-1,pmin(1,scaled))+1)/2)
  grDevices::rgb(rgb[,1],rgb[,2],rgb[,3],maxColorValue=255)
}
cell_ink <- function(fill) {
  rgb <- grDevices::col2rgb(fill)/255
  lin <- ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
  lum <- drop(c(.2126,.7152,.0722) %*% lin)
  ifelse((lum+.05)/.05>=1.05/(lum+.05),"#000000","#FFFFFF")
}
origins <- c(clinical=36.8,protein=110.95)
cw <- 11.44; rh <- 3.65; group_gap <- .91
grid_top <- c(HR_int=19.05,IC=79.50)
row_top <- function(metric) {
  k <- seq_along(endpoints)-1L
  grid_top[[metric]]+rh*k+group_gap*(k%/%3L)
}
outcomes <- function(metric) {
  top <- row_top(metric); mid <- top+rh/2
  mult <- metric=="HR_int"; head_y <- top[1]-9.0
  label(3.1,head_y,if(mult) "Multiplicative interaction" else "Additive interaction",7.15,face=2)
  label(3.1,head_y+3.85,if(mult) "Product-term HR" else "10-year IC (percentage points)",6.5,muted)
  for (g in seq_along(group_names)) {
    t <- top[3*g-2]
    box_mm(3.1,t+.08,3.75,3*rh-.16,group_fill[g])
    label(4.975,t+1.5*rh,group_names[g],5.65,muted,just="centre",rot=90)
  }
  for (branch in list(c(4L,5L),c(7L,8L,9L))) {
    line_mm(rep(8.65,2),c(mid[branch[1]]+1.05,mid[tail(branch,1)]),"#A2ADAF",.13)
    for (i in branch[-1]) line_mm(c(8.65,9.9),rep(mid[i],2),"#A2ADAF",.13)
  }
  for (i in seq_along(top)) {
    nested <- endpoints[i] %in% c("ipf","pe","ph")
    label(if(nested) 10.55 else 8.35,mid[i],endpoint_names[i],
          if(nested) 6.45 else 6.65,if(nested) muted else ink,if(nested) 3 else 1)
  }
}
atlas <- function(ly,metric) {
  x0 <- origins[[ly]]; top <- row_top(metric)
  z <- shown[shown$layer==ly,]
  for (j in seq_along(systems)) {
    x <- x0+(j-1)*cw
    label(x+cw/2,top[1]-5.85,system_names[j],6.75,just="centre",rot=38)
    series <- z[z$other==systems[j],]
    idx <- match(endpoints,series$endpoint)
    stopifnot(identical(which(is.na(idx)),match("pneumoconiosis",endpoints)))
    for (i in seq_along(endpoints)) {
      y <- top[i]
      if (is.na(idx[i])) {
        box_mm(x+.10,y+.10,cw-.20,rh-.20,grey_na)
        label(x+cw/2,y+rh/2,"–",6.3,"#7B878B",just="centre")
      } else {
        row <- series[idx[i],]
        v <- row[[metric]]
        sig <- if(metric=="HR_int") row$mult_mark else row$add_mark
        digits <- if(metric=="IC" || (sig && abs(v-1)<.005)) 3L else 2L
        s <- paste0(number(v,digits),if(sig) "*" else "")
        fill <- fill_value(v,metric)
        box_mm(x+.10,y+.10,cw-.20,rh-.20,fill)
        label(x+cw/2,y+rh/2,s,6.65,cell_ink(fill),just="centre")
      }
    }
  }
}
colour_key <- function(x,w,metric) {
  y <- 134.1
  is_hr <- metric=="HR_int"
  label(x,y-3.15,if(is_hr) "Product-term HR" else "10-year IC (percentage points)",7.0,face=2)
  stops <- seq(0,1,length.out=65)
  signed <- 2*stops-1
  vals <- if(is_hr) exp(signed*log(1.12)) else signed*.45
  gradient <- linearGradient(fill_value(vals,metric),stops=stops,x1=0,y1=.5,x2=1,y2=.5)
  box_mm(x,y,w,1.70,gradient)
  ticks <- if(is_hr) c("0.89","1.00","1.12") else c("−0.45","0","+0.45")
  for (j in seq_along(ticks)) {
    at <- x+(j-1)*w/2
    line_mm(rep(at,2),c(y+1.7,y+2.35),axis_col,.13)
    label(at,y+4.05,ticks[j],6.25,muted,just=c("left","centre","right")[j])
  }
}

# The label positions are fixed in physical units; the points are always at their fitted values.
scatter_labels <- list(
  clinical=data.frame(endpoint=c("pvd","ph","copd","pneumonia"),
                      x=c(63.5,45.1,44.7,64.8),y=c(164.7,168.2,179.9,173.9),
                      align=c("left","right","right","left")),
  protein=data.frame(endpoint=c("pvd","ph","copd","pneumonia","asthma"),
                     x=c(175.4,177.0,157.5,156.1,161.8),y=c(160.2,171.1,160.8,173.75,152.2),
                     align=c("centre","centre","right","right","left"))
)
scatter <- function(ly,l,r,tag_x,tag) {
  header(tag_x,146.0,tag,if(ly=="clinical") "Clinical" else "Proteomic")
  z <- shown[shown$layer==ly & is.finite(shown$share),]
  xd <- c(.905,1.13); yd <- c(-32,77)
  top <- 150.55; bottom <- 185.80
  sx <- function(v) scale_to(v,xd,c(l,r))
  sy <- function(v) scale_to(v,yd,c(bottom,top))
  stopifnot(all(z$HR_int>=xd[1] & z$HR_int<=xd[2] & z$share>=yd[1] & z$share<=yd[2]))
  line_mm(rep(sx(1),2),c(top,bottom),"#CBD2D3",.17)
  line_mm(c(l,r),rep(sy(0),2),"#CBD2D3",.17)
  line_mm(rep(l,2),c(top,bottom),axis_col,.22)
  line_mm(c(l,r),rep(bottom,2),axis_col,.22)
  for (v in c(.92,.96,1,1.04,1.08,1.12)) {
    x <- sx(v)
    line_mm(rep(x,2),c(bottom,bottom+.70),axis_col,.17)
    label(x,bottom+2.55,sprintf("%.2f",v),6.40,muted,just="centre")
  }
  for (v in seq(-20,60,20)) {
    y <- sy(v)
    line_mm(c(l-.82,l),rep(y,2),axis_col,.17)
    label(l-1.55,y,number(v,0),6.40,muted,just="right")
  }
  label((l+r)/2,bottom+6.4,"Product-term HR",7.0,just="centre")
  if (ly=="clinical") label(4.50,(top+bottom)/2,"Share of joint excess risk (%)",7.0,just="centre",rot=90)
  col <- layer_col[[ly]]
  annotations <- scatter_labels[[ly]]
  for (j in seq_len(nrow(annotations))) {
    a <- annotations[j,]
    row <- one(z[z$other=="cardiovascular" & z$endpoint==a$endpoint,])
    x_end <- a$x+switch(a$align,left=-.65,right=.65,centre=0)
    y_end <- if(a$align=="centre") a$y+sign(sy(row$share)-a$y)*1.28 else a$y
    line_mm(c(sx(row$HR_int),x_end),c(sy(row$share),y_end),
            if(ly=="clinical") "#9AAAC4" else "#88B8AE",.15)
  }
  z <- z[order(z$other=="cardiovascular"),]
  for (i in seq_len(nrow(z))) {
    feature <- z$other[i]=="cardiovascular"
    twin <- feature && ly=="protein" && z$endpoint[i] %in% c("pvd","ph")
    radius <- if(twin) .51 else if(feature) .79 else .47
    point_mm(sx(z$HR_int[i]),sy(z$share[i]),radius,
             if(feature) col else "#A9B3B6",z$mult_mark[i] && z$add_mark[i])
  }
  for (j in seq_len(nrow(annotations))) {
    a <- annotations[j,]
    label(a$x,a$y,display_names[[a$endpoint]],6.55,col,just=a$align)
  }
}
scatter_key <- function() {
  y <- 197.65
  point_mm(17.0,y,.69,layer_col[1]); point_mm(19.0,y,.69,layer_col[2])
  label(21.1,y,"Lung × cardiovascular",6.6)
  point_mm(57.9,y,.53,"#A9B3B6"); label(60.0,y,"Other lung pairs",6.6)
  point_mm(98.25,y,.65,ink); label(100.4,y,"Both 95% CIs exclude the null",6.6)
  point_mm(149.4,y,.65,ink,FALSE); label(151.45,y,"Other estimates",6.6)
}

# Absolute risks are percentages; the forest plots use percentage-point differences.
risk_table <- function() {
  header(3.1,205.95,"E","Clinical lung × cardiovascular")
  header(right_tag,205.95,"F","Proteomic lung × cardiovascular")
  ys <- c(221.0,226.0,231.0,236.0)
  label(3.5,216.3,"Outcome",6.95,face=2)
  for (i in seq_along(examples)) {
    if (i%%2==0L) box_mm(3.1,ys[i]-2.26,176.7,4.52,"#F2F5F5")
    nested <- examples[i]=="ph"
    label(if(nested) 5.85 else 3.52,ys[i],display_names[[examples[i]]],6.9,
          if(nested) muted else ink,if(nested) 3 else 1)
  }
  line_mm(c(3.76,3.76),c(ys[2]+1.40,ys[3]),"#A2ADAF",.13)
  line_mm(c(3.76,5.10),rep(ys[3],2),"#A2ADAF",.13)
  line_mm(c(3.1,179.8),rep(218.5,2),rule,.15)
  for (ly in layers) {
    gx <- if(ly=="clinical") 19.2 else 101.1
    cx <- gx+c(3.7,10.7,17.0,25.1,33.2)
    l <- gx+38.3; r <- gx+55.8
    endx <- gx+78.5
    col <- layer_col[[ly]]
    label(mean(range(cx)),211.65,"Standardised 10-year risk (%)",7.0,face=2,just="centre")
    label((l+endx)/2,211.65,"IC, percentage points (95% CI)",6.9,face=2,just="centre")
    line_mm(c(gx+.15,gx+37),rep(213.45,2),rule,.14)
    line_mm(c(l-.35,endx),rep(213.45,2),rule,.14)
    labs <- c("Neither","Lung","CV","Expected","Both")
    for (j in seq_along(labs)) label(cx[j],216.3,labs[j],6.55,if(j==5) col else ink,2,"centre")
    domain <- c(-.10,.56)
    sx <- function(v) scale_to(v,domain,c(l,r))
    line_mm(rep(sx(0),2),c(219.4,238.0),"#B1BDBE",.15)
    for (i in seq_along(examples)) {
      ep <- examples[i]
      v <- one(shown[shown$layer==ly & shown$other=="cardiovascular" & shown$endpoint==ep,])
      exp_risk <- with(v,F10+F01-F00)
      values <- c(v$F00,v$F10,v$F01,exp_risk,v$F11)
      for (j in seq_along(values)) label(cx[j],ys[i],number(values[j],2),6.75,
                                         if(j==5) col else ink,if(j==5) 2 else 1,"centre")
      if (ly=="protein") {
        b <- one(boot[boot$endpoint==ep,])
        stopifnot(abs(b$IC_delta-v$IC)<1e-7)
        lo <- b$IC_lo; hi <- b$IC_hi
      } else {
        lo <- v$IC_lo; hi <- v$IC_hi
      }
      stopifnot(domain[1]<=lo,lo<=v$IC,v$IC<=hi,hi<=domain[2])
      line_mm(sx(c(lo,hi)),rep(ys[i],2),col,.33)
      point_mm(sx(v$IC),ys[i],.57,col)
      label(endx,ys[i],ci_text(v$IC,lo,hi),6.7,just="right")
    }
    line_mm(c(l,r),rep(239.1,2),axis_col,.16)
    forest_ticks <- c(-.10,0,.25,.50)
    for (j in seq_along(forest_ticks)) {
      x <- sx(forest_ticks[j])
      line_mm(rep(x,2),c(239.1,239.75),axis_col,.15)
      label(x,241.40,c("−0.10","0","0.25","0.50")[j],6.1,muted,just=if(j==1L) "right" else "centre")
    }
  }
}
draw <- function() {
  grid.newpage()
  box_mm(0,0,W,H,"white")
  header(3.1,4.8,"A","Clinical age acceleration")
  header(right_tag,4.8,"B","Proteomic age acceleration")
  for (metric in names(grid_top)) {
    outcomes(metric)
    for (ly in layers) atlas(ly,metric)
  }
  colour_key(3.1,46.1,"HR_int"); colour_key(61.0,46.1,"IC")
  label(119.5,132.5,"* Nominal 95% CI excludes 1 (HR) or 0 (IC).",6.65)
  box_mm(119.5,136.7,3.0,1.95,grey_na)
  label(124.0,137.65,"Pneumoconiosis omitted (sparse events)",6.5,muted)
  line_mm(c(3.1,179.8),rep(141.9,2),rule,.17)
  scatter("clinical",16.8,87.5,3.1,"C")
  scatter("protein",108.0,178.7,right_tag,"D")
  scatter_key()
  line_mm(c(3.1,179.8),rep(201.9,2),rule,.17)
  risk_table()
}
stem <- file.path(out_dir,"Figure5")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=font,bg="white")
tryCatch(draw(),finally=invisible(dev.off()))
write_figure_rasters(stem, W, H, draw)
write_figure_svg(stem, W, H, font, draw)
