# Figure3
# Rscript --vanilla Figure3/Figure3.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure3")
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir, call.=FALSE)
output_dir <- if (length(args)) path.expand(args[[1L]]) else file.path(dirname(script_dir), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

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

read_tsv <- function(file) read.delim(file,stringsAsFactors=FALSE,check.names=FALSE,
                                      fileEncoding="UTF-8",na.strings=c("","NA"))
data_table <- function(name) read_tsv(file.path(data_dir,paste0(name,".tsv")))
weights <- data_table("clinical_lung_coef")
proteins <- data_table("single_protein_v1")
extra <- data_table("extra_models_v1")
joint <- extra[extra$model=="B + other-layer lung age",]
spline <- data_table("rcs_v1")
histogram <- data_table("hist_v1")
cif <- data_table("cif_v1")
risk_summary <- data_table("cif_summary_v1")
checks <- read_tsv(file.path(data_dir,"sensitivity_checks.tsv"))
stopifnot(!anyDuplicated(histogram[c("endpoint","layer","x0")]),!anyDuplicated(spline[c("endpoint","layer","x")]),
          !anyDuplicated(joint[c("endpoint","layer")]),!anyDuplicated(weights[c("input","sex")]),
          !anyDuplicated(cif[c("endpoint","layer","g","time")]),
          !anyDuplicated(risk_summary[c("endpoint","layer","g")]),
          !anyDuplicated(checks[c("panel","analysis")]),
          is.numeric(checks$p),all(is.finite(checks$p) & checks$p>=0 & checks$p<=1))
one <- function(d) {
  stopifnot(nrow(d)==1L)
  d[1L,,drop=FALSE]
}

font <- "Times New Roman"
ink <- "#20272B"; muted <- "#48565D"; rule <- "#CED6D8"; light <- "#EDF1F2"; stripe <- "#F3F6F6"
layer_col <- c(clinical="#395688",protein="#168372")
layers <- names(layer_col)
risk_col <- c("#4387A1","#68767B","#BB792F","#A43D35")
band_col <- function(col,alpha=.17) grDevices::adjustcolor(col,alpha.f=alpha)
int <- function(x) format(as.integer(x),big.mark=",",scientific=FALSE,trim=TRUE)
estimate <- function(h,lo,hi,d=2L) {
  pattern <- paste0("%.",d,"f (%.",d,"f\u2013%.",d,"f)")
  sprintf(pattern,h,lo,hi)
}
unit_x <- function(x) grid::unit(x,"mm")
unit_y <- function(y) grid::unit(1,"npc")-grid::unit(y,"mm")
text_at <- function(x,y,s,size=7,col=ink,face=1,align="left",angle=0) {
  grid::grid.text(s,x=unit_x(x),y=unit_y(y),just=c(align,"centre"),rot=angle,
                  gp=grid::gpar(fontfamily=font,fontface=face,fontsize=size,col=col))
}
path <- function(x,y,col=rule,mm=.17,lty=1) {
  grid::grid.lines(unit_x(x),unit_y(y),gp=grid::gpar(col=col,lwd=mm/.2646,lty=lty,
                                                     lineend="butt",linejoin="round"))
}
box_at <- function(x,y,w,h,fill) {
  grid::grid.rect(unit_x(x),unit_y(y),width=unit_x(w),height=unit_x(h),just=c("left","top"),
                  gp=grid::gpar(fill=fill,col=NA))
}
polygon_at <- function(x,y,fill) {
  grid::grid.polygon(unit_x(x),unit_y(y),gp=grid::gpar(fill=fill,col=NA))
}
dot <- function(x,y,col,shape="clinical",radius=.68,fill=col) {
  gp <- grid::gpar(col=col,fill=fill,lwd=.16/.2646)
  if (shape=="protein") {
    grid::grid.polygon(unit_x(x+c(0,radius,0,-radius)),unit_y(y+c(-radius,0,radius,0)),gp=gp)
  } else if (shape=="square") {
    grid::grid.rect(unit_x(x),unit_y(y),unit_x(1.65*radius),unit_x(1.65*radius),gp=gp)
  } else {
    grid::grid.circle(unit_x(x),unit_y(y),r=unit_x(.8*radius),gp=gp)
  }
}
linear <- function(v,domain,rng) rng[1]+diff(rng)*(v-domain[1])/diff(domain)
log_map <- function(v,domain,rng) linear(log(v),log(domain),rng)
haxis <- function(xs,y,labels,size=6.25,col=muted) {
  path(range(xs),rep(y,2),rule,.17)
  for (i in seq_along(xs)) {
    path(rep(xs[i],2),c(y,y+.72),col,.13)
    text_at(xs[i],y+2.45,labels[i],size,col,align="centre")
  }
}
hr_axis <- function(left,right,y,domain,ticks=c(1,1.5,2,3)) {
  path(c(left,right),rep(y,2),rule,.17)
  for (v in ticks) {
    x <- log_map(v,domain,c(left,right))
    path(c(x,x),c(y,y+.65),muted,.13)
    text_at(x,y+2.35,format(v,trim=TRUE),6.0,muted,align="centre")
  }
}
hr_mark <- function(d,yy,left,right,domain,layer,radius=.72) {
  xs <- log_map(c(d$lo,d$HR,d$hi),domain,c(left,right))
  stopifnot(all(xs>=left-.001 & xs<=right+.001))
  path(xs[c(1,3)],rep(yy,2),layer_col[layer],.31)
  dot(xs[2],yy,layer_col[layer],layer,radius)
}
divider <- function(y) path(c(3.1,179.8),rep(y,2),rule,.17)
left_origin <- 3.1; right_origin <- 95.0

panel_letter <- function(x,y,s) text_at(x,y,s,10.35,face=2)
swatch <- function(x,y,layer,label,draw_line=FALSE,size=6.65) {
  if (draw_line) path(x+c(-1.95,1.95),rep(y,2),layer_col[layer],.35)
  dot(x,y,layer_col[layer],layer,.83)
  text_at(x+if(draw_line) 3.0 else 1.8,y,label,size)
}
first_panel <- function() {
  panel_letter(3.1,4.8,"A")
  dot(15.7,9.55,layer_col[1],radius=.87)
  text_at(17.55,9.55,"Women",6.65)
  dot(31.15,9.55,layer_col[1],radius=.87,fill="white")
  text_at(33.0,9.55,"Men",6.65)
  names <- c("fev1","fvc","fev1_fvc","pef")
  labs <- c("FEV\u2081","FVC","FEV\u2081/FVC","PEF")
  ys <- 16.2+(0:3)*6.25
  l <- 15.0; r <- 45.4; dom <- c(-8.6,4.45)
  xx <- function(v) linear(v,dom,c(l,r))
  path(rep(xx(0),2),c(13.2,37.6),"#9DA7AD",.17)
  for (i in seq_along(names)) {
    z <- weights[weights$input==names[i],]
    stopifnot(nrow(z)==2L,setequal(z$sex,c("female","male")),all(z$weight>=dom[1]&z$weight<=dom[2]))
    text_at(3.24,ys[i],labs[i],6.85)
    x <- xx(z$weight[match(c("female","male"),z$sex)])
    path(x,rep(ys[i],2),"#98A9C5",.52)
    dot(x[1],ys[i],layer_col[1],radius=.91)
    dot(x[2],ys[i],layer_col[1],radius=.91,fill="white")
  }
  haxis(xx(c(-8,-4,0,4)),38.3,c("\u22128","\u22124","0","4"))
  text_at((l+r)/2,44.0,"Clock weight (years per SD)",6.8,align="centre")
}
protein_panel <- function() {
  panel_letter(49.2,4.8,"B")
  labs <- c(p_pulmonary="Lung clock",lamp3="LAMP3",scgb1a1="SCGB1A1 (CC16)",
            ccl18="CCL18",sftpd="SFTPD (SP-D)")
  ild <- proteins[proteins$endpoint=="ild",]
  d <- ild[match(names(labs),ild$protein),]
  stopifnot(identical(as.character(d$protein),names(labs)),!anyNA(d$HR),all(d$lo<=d$HR & d$HR<=d$hi),
            all(d$weight[-1L]>0 & d$weight[-1L]<=1.6))
  text_at(49.55,9.7,"Protein",6.75,face=2)
  text_at(79.15,9.7,"Weight (y/SD)",6.37,face=2,align="centre")
  text_at(107.85,9.7,"ILD HR (95% CI)",6.75,face=2,align="centre")
  y <- 15.55+(0:4)*4.86; l <- 90.1; r <- 105.7; dom <- c(.9,3.4)
  box_at(49.15,13.27,76.67,4.55,"#E4F0ED")
  for (i in c(3L,5L)) box_at(49.15,y[i]-2.26,76.67,4.52,"#F4F6F6")
  path(rep(log_map(1,dom,c(l,r)),2),c(13.35,37.65),"#99A6A9",.14,"32")
  for (i in seq_len(nrow(d))) {
    z <- d[i,]; composite <- i==1L
    text_at(49.55,y[i],unname(labs[i]),6.52,ink,if(composite) 2 else 1)
    if(composite) {
      text_at(87.55,y[i],"\u2013",6.6,muted,align="right")
    } else {
      box_at(70.2,y[i]-.65,linear(z$weight,c(0,1.6),c(0,9.9)),1.30,layer_col["protein"])
      text_at(87.55,y[i],sprintf("%.2f",z$weight),6.42,align="right")
    }
    xs <- log_map(c(z$lo,z$HR,z$hi),dom,c(l,r))
    path(xs[c(1,3)],rep(y[i],2),layer_col["protein"],.29)
    dot(xs[2],y[i],layer_col["protein"],if(composite) "protein" else "clinical",if(composite) .88 else .69)
    text_at(125.65,y[i],estimate(z$HR,z$lo,z$hi),6.58,ink,if(composite) 2 else 1,"right")
  }
  hr_axis(l,r,38.3,dom,c(1,2,3))
}
joint_panel <- function() {
  panel_letter(129.0,4.8,"C")
  swatch(139.0,9.55,"clinical","Clinical")
  swatch(155.0,9.55,"protein","Proteomic")
  e <- c("copd","ild","ipf"); x <- c(143.6,159.3,174.5)
  l <- 137.9; r <- 180.2; dom <- c(.87,3.23)
  yy <- function(v) log_map(v,dom,c(38.3,14.25))
  for (v in c(1,1.5,2,3)) {
    path(c(l,if(v==3) l+.65 else r),rep(yy(v),2),if(v==1) "#8F9C9F" else light,.14,if(v==1) "32" else 1)
    text_at(l-1.1,yy(v),format(v,trim=TRUE),6.1,muted,align="right")
  }
  text_at(129.9,26.4,"Mutually adjusted HR per SD",6.55,angle=90,align="centre")
  for (k in seq_along(layers)) {
    ly <- layers[k]; z <- joint[joint$layer==ly,]; z <- z[match(e,z$endpoint),]
    stopifnot(identical(as.character(z$endpoint),e),all(z$lo>dom[1]&z$hi<dom[2]))
    xp <- x+c(-.78,.78)[k]
    path(xp,yy(z$HR),band_col(layer_col[ly],.62),.36)
    for (i in seq_along(e)) {
      path(rep(xp[i],2),yy(c(z$lo[i],z$hi[i])),layer_col[ly],.30)
      for (v in c(z$lo[i],z$hi[i])) path(xp[i]+c(-.42,.42),rep(yy(v),2),layer_col[ly],.17)
      dot(xp[i],yy(z$HR[i]),layer_col[ly],ly,.87)
      above <- (i==1 && ly=="clinical") || (i>1 && ly=="protein")
      side <- i==1L && ly=="protein"
      yt <- if(side) yy(z$HR[i])+.65 else if(above) yy(z$hi[i])-1.60 else min(yy(z$lo[i])+1.45,yy(1)-1.4)
      text_at(xp[i]+if(side) 1.8 else 0,yt,sprintf("%.2f",z$HR[i]),6.3,layer_col[ly],
              align=if(side) "left" else "centre")
    }
  }
  path(c(l,r),c(38.3,38.3),rule,.17)
  for (j in seq_along(e)) text_at(x[j],41.5,c("COPD","ILD","Fibrotic ILD")[j],6.9,face=2,align="centre")
}
dose_cols <- c(clinical="#3C5488",protein="#00A087")
dose_ribbons <- c(clinical="#C9D1E6",protein="#C4E6DE")
axis_ink <- "#474A4B"

# Isolate the physical clip rectangle while drawing in the shared page coordinates.
in_plot_window <- function(l,top,r,b,draw) {
  grid::pushViewport(grid::viewport(clip="on"))
  on.exit(grid::popViewport())
  grid::grid.clip(x=unit_x(l),y=unit_y(top),width=unit_x(r-l),height=unit_x(b-top),
                  just=c("left","top"))
  force(draw)
}
solid_y_axis <- function(l,top,b,ys,labels) {
  path(c(l,l),c(top,b),axis_ink,.26)
  for (i in seq_along(ys)) {
    path(c(l-.85,l),rep(ys[i],2),axis_ink,.21)
    text_at(l-1.63,ys[i],labels[i],6.4,ink,align="right")
  }
}
solid_x_axis <- function(l,r,b,xs,labels) {
  path(c(l,r),rep(b,2),axis_ink,.26)
  for (i in seq_along(xs)) {
    path(rep(xs[i],2),c(b,b+.82),axis_ink,.21)
    text_at(xs[i],b+2.52,labels[i],6.45,ink,align="centre")
  }
}
fmt_p <- function(p,equality="= ") {
  if(p<.001) "< 0.001" else paste0(equality,sprintf("%.3f",p))
}
dose_panel <- function(x0,letter,ep) {
  panel_letter(x0,51.3,letter)
  d <- spline[spline$endpoint==ep,]
  l <- x0+12.35; r <- x0+84.5; top <- 52.05; b <- 83.85
  # Match the v2 coordinate window, including ggplot2's 1% x and 5% log-y expansion.
  xdomain <- c(-3,4); xd <- xdomain+c(-1,1)*.01*diff(xdomain)
  log_limits <- log(c(.25,16))
  yd <- exp(log_limits+c(-1,1)*.05*diff(log_limits))
  xx <- function(v) linear(v,xd,c(l,r))
  yy <- function(v) log_map(v,yd,c(b,top))
  stopifnot(all(d$x>=xdomain[1] & d$x<=xdomain[2]),all(d$lo>0 & d$lo<=d$HR & d$HR<=d$hi))
  in_plot_window(l,top,r,b,{
    for (ly in layers) {
      z <- d[d$layer==ly,]; z <- z[order(z$x),]
      polygon_at(c(xx(z$x),rev(xx(z$x))),c(yy(z$lo),rev(yy(z$hi))),band_col(dose_ribbons[ly],.55))
    }
    for (ly in layers) {
      z <- d[d$layer==ly,]; z <- z[order(z$x),]
      path(xx(z$x),yy(z$HR),dose_cols[ly],.53)
    }
  })
  yticks <- c(.25,.5,1,2,4,8,16)
  solid_y_axis(l,top,b,yy(yticks),c("0.25","0.5","1","2","4","8","16"))
  text_at(x0+1.7,(top+b)/2,paste0("HR for ",toupper(ep)," (95% CI)"),6.9,angle=90,align="centre")
  for (k in seq_along(layers)) {
    z <- d[d$layer==layers[k],]
    po <- unique(z$p_overall); pn <- unique(z$p_nonlinear)
    stopifnot(length(po)==1L,length(pn)==1L)
    row <- paste0("P overall ",fmt_p(po),"     P nonlinearity ",fmt_p(pn))
    text_at(l+1.6,yy(16/c(1,1.55)[k])+1.02,row,6.27,dose_cols[k])
  }

  # As in v2, each layer sums to one across the displayed 0.25-SD bins.
  h <- histogram[histogram$endpoint==ep & histogram$x1>xdomain[1] & histogram$x0<xdomain[2],]
  stopifnot(all(h$n>=0),all(abs(h$x1-h$x0-.25)<1e-8))
  h$prop <- h$n/ave(h$n,h$layer,FUN=sum)
  hist_top <- 85.0; floor <- 92.15
  hy <- function(p) linear(p,c(0,max(h$prop)),c(floor,hist_top))
  hc <- h[h$layer=="clinical",]
  for (i in seq_len(nrow(hc))) {
    box_at(xx(hc$x0[i])+.065,hy(hc$prop[i]),xx(hc$x1[i])-xx(hc$x0[i])-.13,
           floor-hy(hc$prop[i]),"#DCE1EE")
  }
  hp <- h[h$layer=="protein",]; hp <- hp[order(hp$x0),]
  path(xx(as.vector(rbind(hp$x0,hp$x1))),rep(hy(hp$prop),each=2),dose_cols["protein"],.32)
  solid_x_axis(l,r,floor,xx(c(-2,0,2,4)),c("\u22122","0","2","4"))
  text_at((l+r)/2,97.55,"Lung age acceleration (SD)",7.0,align="centre")
}
dose_legend <- function() {
  for (k in seq_along(layers)) {
    x <- c(64.0,98.15)[k]; y <- 101.20
    path(x+c(-2.25,2.25),rep(y,2),dose_cols[k],.53)
    text_at(x+3.25,y,c("Clinical lung age","Proteomic lung age")[k],7.05)
  }
}

step_path <- function(x,y) list(x=rep(x,each=2)[-1L],y=head(rep(y,each=2),-1L))
group_names <- c("<\u22121","\u22121 to <1","1 to <2","\u22652")
risk_panel <- function(x0,letter,ly,ep) {
  panel_letter(x0,107.0,letter)
  d <- cif[cif$layer==ly & cif$endpoint==ep & cif$time<=13,]
  s <- risk_summary[risk_summary$layer==ly & risk_summary$endpoint==ep,]
  s <- s[order(s$g),]
  ymax <- if(ep=="ild") 10 else 20
  stopifnot(identical(as.integer(s$g),0:3),min(d$lo)>=0,
            all(d$cif>=0 & 100*d$cif<=ymax),max(d$hi)<=1)
  l <- x0+12.35; r <- x0+84.5; top <- 112.0; b <- 141.65
  xx <- function(v) linear(v,c(0,13),c(l,r))
  yy <- function(v) linear(100*v,c(0,ymax),c(b,top))
  in_plot_window(l,top,r,b,{
    for (g in 0:3) {
      z <- d[d$g==g,]; z <- z[order(z$time),]
      low <- step_path(xx(z$time),yy(z$lo)); high <- step_path(xx(z$time),yy(z$hi))
      polygon_at(c(low$x,rev(high$x)),c(low$y,rev(high$y)),band_col(risk_col[g+1],.15))
    }
    for (g in 0:3) {
      z <- d[d$g==g,]; z <- z[order(z$time),]; col <- risk_col[g+1]
      curve <- step_path(xx(z$time),yy(z$cif))
      path(curve$x,curve$y,col,if(g==3L) .39 else .31)
      for (yr in c(5,10)) {
        p <- one(z[z$time==yr,])
        path(rep(xx(yr),2),yy(c(p$lo,p$hi)),col,.24)
        dot(xx(yr),yy(p$cif),col,radius=.46)
      }
      ten <- one(z[z$time==10,]); row <- s[g+1,]
      stopifnot(isTRUE(all.equal(c(ten$cif,ten$lo,ten$hi),c(row$risk10,row$risk10_lo,row$risk10_hi),tolerance=1e-10)))
    }
  })
  yticks <- if(ep=="ild") seq(0,10,2) else seq(0,20,5)
  solid_y_axis(l,top,b,yy(yticks/100),as.character(yticks))
  text_at(x0+1.7,(top+b)/2,paste0(toupper(ep)," cumulative incidence (%)"),6.85,angle=90,align="centre")
  solid_x_axis(l,r,b,xx(c(0,5,10,13)),c("0","5","10","13"))
  text_at((l+r)/2,147.6,"Follow-up (years)",7.0,align="centre")
  text_at(x0+.4,152.05,paste(if(ly=="clinical") "Clinical" else "Proteomic","lung age"),6.5,face=2)
  text_at(x0+.4,155.35,"acceleration (SD)",6.5,face=2)
  cx <- x0+c(26.4,39.2,52.05)
  text_at(cx[2],152.05,"Number at risk",6.5,muted,align="centre")
  for (j in 1:3) text_at(cx[j],155.35,c("0 y","5 y","10 y")[j],6.5,face=2,align="centre")
  text_at(x0+84.9,152.05,"10-year risk, %",6.6,face=2,align="right")
  text_at(x0+84.9,155.35,"(95% CI)",6.6,face=2,align="right")
  path(c(x0+.2,x0+85),rep(157.2,2),rule,.14)
  ys <- 160.0+(0:3)*4.22
  for (i in 1:4) {
    row <- s[i,]
    if (i%%2L==0L) box_at(x0+.2,ys[i]-1.98,84.8,3.96,stripe)
    path(x0+c(.55,3.65),rep(ys[i],2),risk_col[i],.5)
    text_at(x0+4.65,ys[i],group_names[i],6.75)
    for (j in seq_along(cx)) {
      nm <- c("at_risk_0","at_risk_5","at_risk_10")[j]
      text_at(cx[j],ys[i],int(row[[nm]]),6.7,align="centre")
    }
    text_at(x0+84.9,ys[i],estimate(100*row$risk10,100*row$risk10_lo,100*row$risk10_hi),6.85,align="right")
  }
}
h_labels <- c("Main model (B)","+ Whole-body age","+ Other lung age","Baseline FEV\u2081/FVC \u2265 0.70",
              "Clock without LAMP3 (post hoc)","+ Cystatin C-based eGFR","Never smokers","Follow-up \u22652 years","Test set")
h_analyses <- list(
  clinical=c("Model B (main)","+ PhenoAge","+ Protein lung age\u2020","FEV1/FVC \u2265 0.70 only",
             NA,NA,"Never smokers","Excluding first 2 years","Test set (T)"),
  protein=c("Model B (main)","+ Organismal age","+ Clinical lung age\u2020",NA,
            "Clock without LAMP3\u2021","+ eGFR (cystatin C)","Never smokers","Excluding first 2 years","Test set (T)")
)
robustness_panel <- function() {
  panel_letter(3.1,181.65,"H")
  text_at(3.5,187.2,"Analysis",6.95,face=2)
  nr <- seq_along(h_labels)
  ys <- 193.05+(nr-1)*4.65+.85*(nr>3)+.85*(nr>6)+.55*(nr>8)
  for (i in nr) {
    if (i==1L || i%%2L==0L)
      box_at(3.1,ys[i]-2.1,176.7,4.2,if(i==1L) "#EAF0F1" else "#F5F7F7")
    if (i %in% c(4,7,9)) path(c(3.1,179.8),rep(ys[i]-2.65,2),rule,.12)
    text_at(3.52,ys[i],h_labels[i],6.88,ink,if(i==1L) 2 else 1)
  }
  path(c(3.1,179.8),rep(189.9,2),rule,.15)
  domain <- c(.94,3.2)
  for (k in seq_along(layers)) {
    ly <- layers[k]; gx <- c(39.8,110.75)[k]
    e_col <- gx+7.4; l <- gx+10.6; r <- gx+35.2
    number <- gx+56.95; p_col <- gx+65.2
    text_at(e_col,187.2,"Events",6.9,face=2,align="right")
    dot(l+1.5,187.2,layer_col[ly],ly,.85)
    text_at(l+3.5,187.2,c("Clinical HR for COPD (95% CI)","Proteomic HR for ILD (95% CI)")[k],6.9,ink,2)
    text_at(p_col,187.2,"P",6.95,face=4,align="centre")
    path(rep(log_map(1,domain,c(l,r)),2),c(190.7,tail(ys,1)+2.3),"#99A6A8",.14,"32")
    block <- c("h_clinical_copd","h_protein_ild")[k]
    z <- checks[checks$panel==block,]; wanted <- h_analyses[[ly]]
    stopifnot(setequal(z$analysis,stats::na.omit(wanted)),nrow(z)==sum(!is.na(wanted)))
    for (i in nr) {
      if(is.na(wanted[i])) {
        text_at(number-8.5,ys[i],"\u2013",6.75,"#839093",align="centre")
        text_at(p_col,ys[i],"\u2013",6.75,"#839093",align="centre")
      } else {
        row <- one(z[z$analysis==wanted[i],])
        text_at(e_col,ys[i],int(row$events),6.85,align="right")
        hr_mark(row,ys[i],l,r,domain,ly,if(i==1L) .92 else .74)
        text_at(number,ys[i],estimate(row$HR,row$lo,row$hi),6.95,ink,if(i==1L) 2 else 1,"right")
        text_at(p_col,ys[i],fmt_p(row$p,equality=""),6.85,ink,if(i==1L) 2 else 1,"centre")
      }
    }
    hr_axis(l,r,tail(ys,1)+3.05,domain,c(1,1.5,2,3))
    text_at((l+r)/2,tail(ys,1)+8.30,"HR per SD",6.45,muted,align="centre")
  }
}
plot_figure <- function() {
  grid::grid.newpage()
  first_panel()
  protein_panel()
  joint_panel()
  divider(46.7)
  dose_panel(left_origin,"D","copd")
  dose_panel(right_origin,"E","ild")
  dose_legend()
  divider(103.85)
  risk_panel(left_origin,"F","clinical","copd")
  risk_panel(right_origin,"G","protein","ild")
  divider(177.25)
  robustness_panel()
}
height <- 245; width <- 183
stem <- file.path(out_dir,"Figure3")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=width/25.4,height=height/25.4,family=font,bg="white")
tryCatch(plot_figure(),finally=invisible(grDevices::dev.off()))
grDevices::png(paste0(stem,".png"),width=width,height=height,units="mm",res=450,type="cairo",bg="white",antialias="gray")
tryCatch(plot_figure(),finally=invisible(grDevices::dev.off()))
grDevices::tiff(paste0(stem,".tiff"),width=width,height=height,units="mm",res=600,type="cairo",bg="white",antialias="gray",compression="lzw")
tryCatch(plot_figure(),finally=invisible(grDevices::dev.off()))
write_figure_svg(stem, width, height, font, plot_figure)
