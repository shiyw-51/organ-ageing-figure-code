# Figure 4. Extrapulmonary associations and robustness.
# Rscript --vanilla Figure4/Figure4.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure4")
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


render_figure4 <- function(data, out_dir, view=c("full","main")) {
  stopifnot(dir.exists(out_dir))
  main_view <- match.arg(view)=="main"
  W <- 183; H <- 210
  family <- "Times New Roman"
  ink <- "#20272B"; muted <- "#4B585E"; rule <- "#D8DEDF"; axis_col <- "#768286"
  clinical <- "#395688"; protein <- "#168372"; missing_fill <- "#CBD1D1"
  systems <- c("cardiovascular","renal","hepatic","immune","musculoskeletal","metabolic")
  sys_names <- c("Cardiovascular","Renal","Liver","Immune","Musculoskeletal","Metabolic")
  eps <- c("copd","asthma","bronchiectasis","ild","ipf","pneumoconiosis",
           "pvd","pe","ph","osa","lung_cancer","pneumonia")
  ep_names <- c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD (J84.1)","Pneumoconiosis",
                "Pulmonary vascular disease","Pulmonary embolism","Pulmonary hypertension",
                "Sleep apnoea","Lung cancer","Pneumonia")
  rim_names <- c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD","Pneumoconiosis",
                 "PVD","PE","PH","Sleep apnoea","Lung cancer","Pneumonia")
  groups <- c("Airway","Parenchymal","Vascular","Other")
  group_fill <- c("#E3E7F0","#DEEEE9","#F4E4DF","#EEE7DE")

  required <- c("panel","endpoint","layer","system","analysis","HR","lo","hi")
  stopifnot(all(required %in% names(data)), setequal(tolower(data$panel),c("a","b","c")),
            all(data$endpoint %in% eps),
            all(is.finite(data$HR) & is.finite(data$lo) & is.finite(data$hi)),
            all(data$lo>0 & data$lo<=data$HR & data$HR<=data$hi))
  data$excludes <- data$lo>1 | data$hi<1
  forest <- data[tolower(data$panel)=="a",]
  stopifnot(setequal(forest$system,systems), setequal(forest$layer,c("clinical","protein")),
            all(forest$model=="C"), !anyDuplicated(forest[c("endpoint","system","layer")]))
  analyses <- function(third,ninth) c("Model B","+ Lung age (C)",third,"All seven (D)","+ Pack-years",
                                      "Never smokers","Excl. first 2 y","Fine\u2013Gray",ninth)
  make_ring <- function(code, layer, system, lev) {
    d <- data[tolower(data$panel)==code,]
    stopifnot(all(d$layer==layer),all(d$system==system),setequal(as.character(d$analysis),lev),
              !anyDuplicated(d[c("endpoint","analysis")]))
    keys <- cbind(match(d$endpoint,eps),match(d$analysis,lev))
    hr <- matrix(NA_real_,length(eps),length(lev),dimnames=list(eps,lev))
    sig <- matrix(FALSE,length(eps),length(lev),dimnames=list(eps,lev))
    hr[keys] <- d$HR; sig[keys] <- d$excludes
    shown <- if (main_view) c("Model B","+ Lung age (C)","All seven (D)",tail(lev,1L)) else lev
    list(hr=hr[,shown,drop=FALSE],sig=sig[,shown,drop=FALSE])
  }
  B <- make_ring("b","clinical","renal",analyses("+ PhenoAge","No cystatin C (C)"))
  C <- make_ring("c","protein","cardiovascular",analyses("+ Organismal","+ eGFR (C)"))
  omitted_c <- "pneumoconiosis"

  xx <- function(v) grid::unit(v,"mm")
  yy <- function(v) grid::unit(H-v,"mm")
  txt <- function(x,y,s,size=7,col=ink,face=1,just="left",rot=0) {
    grid::grid.text(s,x=xx(x),y=yy(y),just=c(just,"centre"),rot=rot,
                    gp=grid::gpar(fontfamily=family,fontsize=size,fontface=face,col=col,lineheight=.96))
  }
  line <- function(x,y,col=rule,width=.16,lty=1) {
    grid::grid.lines(xx(x),yy(y),gp=grid::gpar(col=col,lwd=width/.2646,lty=lty,lineend="butt",linejoin="round"))
  }
  rect <- function(x,y,w,h,fill) {
    grid::grid.rect(x=xx(x),y=yy(y),width=xx(w),height=xx(h),just=c("left","top"),
                    gp=grid::gpar(fill=fill,col=NA))
  }
  poly <- function(x,y,fill,stroke=NA,lw=.13,lty=1) {
    grid::grid.polygon(xx(x),yy(y),gp=grid::gpar(fill=fill,col=stroke,lwd=lw/.2646,lty=lty,linejoin="round"))
  }
  point <- function(x,y,layer,solid,col) {
    fill <- if (solid) col else "white"
    if (layer=="protein") {
      r <- .68
      poly(x+c(0,r,0,-r),y+c(-r,0,r,0),fill,col,.19)
    } else {
      grid::grid.circle(xx(x),yy(y),r=xx(.53),gp=grid::gpar(fill=fill,col=col,lwd=.19/.2646))
    }
  }
  header <- function(x,y,letter,title) {
    txt(x,y,letter,10.35,ink,2)
    txt(x+5.25,y+.1,title,9.05,ink,2)
  }

  # A. All six facets use the same log-HR coordinates and permanently ordered layer slots.
  first_x <- 36.3; facet_w <- 22.25; facet_gap <- 1.85
  row_h <- 4.68; group_gap <- 1.25
  top <- 16.1
  row_y <- top+(seq_along(eps)-1)*row_h+(seq_along(eps)-1)%/%3*group_gap
  mid <- row_y+row_h/2
  bottom <- tail(row_y,1)+row_h
  xlim <- range(c(.8,1.8,forest$lo,forest$hi))
  x_at <- function(v,x0) x0+facet_w*(log(v)-log(xlim[1]))/diff(log(xlim))

  draw_forest <- function() {
    header(3.1,4.8,"A","Extra-pulmonary age acceleration")
    txt(8.35,12.1,"Respiratory outcome",7.25,face=2)
    for (g in seq_along(groups)) {
      y <- row_y[3*g-2]
      rect(3.1,y,3.75,3*row_h,group_fill[g])
      txt(4.975,y+1.5*row_h,groups[g],6.13,muted,just="centre",rot=90)
    }
    for (branch in list(c(4L,5L),c(7L,8L,9L))) {
      line(c(8.68,8.68),c(mid[branch[1]]+1.45,mid[tail(branch,1)]),"#A2ADAF",.13)
      for (k in branch[-1]) line(c(8.68,9.95),rep(mid[k],2),"#A2ADAF",.13)
    }
    nest <- eps %in% c("ipf","pe","ph")
    for (i in seq_along(eps)) txt(if(nest[i]) 10.55 else 8.35,mid[i],ep_names[i],
                                  if(nest[i]) 6.65 else 6.85,if(nest[i]) muted else ink,if(nest[i]) 3 else 1)
    for (j in seq_along(systems)) {
      s <- systems[j]; x0 <- first_x+(j-1)*(facet_w+facet_gap)
      txt(x0+facet_w/2,12.1,sys_names[j],7.35,if(s=="renal") clinical else ink,2,just="centre")
      for (i in seq_along(eps)) {
        alternate <- i%%2L==0L
        fill <- if(s=="renal") {if(alternate) "#E5EBF3" else "#EFF3F8"} else {if(alternate) "#F5F7F7" else "white"}
        rect(x0,row_y[i],facet_w,row_h,fill)
      }
      for (g in 0:3) line(rep(x_at(1,x0),2),c(row_y[3*g+1],row_y[3*g+3]+row_h),"#8D989D",.15,"32")
      for (layer in c("clinical","protein")) {
        d <- forest[forest$system==s & forest$layer==layer,]
        ix <- match(eps,d$endpoint)
        ys <- mid+if(layer=="clinical") -.90 else .90
        col <- if(layer=="clinical") clinical else protein
        for (i in seq_along(eps)) {
          k <- ix[i]
          if (is.na(k)) {
            txt(x0+facet_w/2,ys[i],"\u2013",6.35,"#8B979A",just="centre")
          } else {
            line(x_at(c(d$lo[k],d$hi[k]),x0),rep(ys[i],2),col,.26)
            point(x_at(d$HR[k],x0),ys[i],layer,d$excludes[k],col)
          }
        }
      }
      axis_y <- bottom+.75
      line(c(x0,x0+facet_w),rep(axis_y,2),"#ABB5B8",.14)
      at <- c(.8,1,1.25,1.6)
      for (k in seq_along(at)) {
        xxk <- x_at(at[k],x0)
        line(rep(xxk,2),c(axis_y,axis_y+.65),axis_col,.13)
        txt(xxk,axis_y+2.25,c("0.8","1","1.25","1.6")[k],6.0,muted,just="centre")
      }
    }
    txt(107.6,82.5,"Hazard ratio per 1 SD after adjustment for lung age (model C)",7.1,just="centre")
    y <- 87.45
    line(c(3.3,7.4),rep(y,2),clinical,.26); point(5.35,y,"clinical",TRUE,clinical)
    txt(8.55,y,"Clinical",6.85)
    line(c(23.7,27.8),rep(y,2),protein,.26); point(25.75,y,"protein",TRUE,protein)
    txt(28.95,y,"Proteomic",6.85)
    point(51.1,y,"clinical",TRUE,ink); txt(53.2,y,"95% CI excludes 1",6.70)
    point(81.5,y,"clinical",FALSE,ink); txt(83.6,y,"95% CI includes 1",6.70)
    txt(111.7,y,"\u2013",7.0,muted,just="centre"); txt(114.1,y,"No estimate",6.70)
  }

  # Use the Figure 2 scale in its original range, with an extra dark-red stop for HRs above 1.60.
  anchors <- c(-log(1.8),seq(-log(1.6),log(1.6),length.out=7L),log(1.8))
  palette <- c("#0F477A","#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B","#721227")
  grad <- grDevices::colorRamp(palette,space="Lab")
  heat_colour <- function(hr) {
    pos <- stats::approx(anchors,seq(0,1,length.out=length(anchors)),xout=log(hr),rule=2)$y
    v <- grad(pos)
    grDevices::rgb(v[,1],v[,2],v[,3],maxColorValue=255)
  }
  symbol_colour <- function(fill) {
    v <- as.numeric(grDevices::col2rgb(fill))/255
    l <- sum(c(.2126,.7152,.0722)*ifelse(v<=.04045,v/12.92,((v+.055)/1.055)^2.4))
    if ((l+.05)/.05 >= 1.05/(l+.05)) "#000000" else "#FFFFFF"
  }
  # Angles are clockwise from twelve o'clock. Both rings share the same sector boundaries.
  ring_start <- 10
  slice <- (360-2*ring_start-3*2.9)/12
  a0 <- ring_start+(seq_along(eps)-1)*slice+(seq_along(eps)-1)%/%3*2.9
  a1 <- a0+slice; centre_angle <- (a0+a1)/2
  xy <- function(cx,cy,r,angle) {
    t <- angle*pi/180
    list(x=cx+r*sin(t),y=cy-r*cos(t))
  }
  tangent <- function(deg) ((-deg+90)%%180)-90
  wedge <- function(cx,cy,r0,r1,d0,d1,fill,stroke=NA,width=.13,lty=1) {
    a <- seq(d0,d1,length.out=max(3L,ceiling((d1-d0)/1.35)+1L))
    pts <- xy(cx,cy,c(rep(r1,length(a)),rep(r0,length(a))),c(a,rev(a)))
    poly(pts$x,pts$y,fill,stroke,width,lty)
  }
  ring_y <- 136.0; outer <- 31.15; inner <- 11.65
  draw_ring <- function(m,cx,omit=character()) {
    n_rings <- ncol(m$hr)
    pitch <- (outer-inner)/n_rings
    for (g in seq_along(groups)) {
      ia <- 3*g-2; ib <- 3*g
      wedge(cx,ring_y,8.8,inner-.4,a0[ia]+.16,a1[ib]-.16,group_fill[g])
      phi <- (a0[ia]+a1[ib])/2
      p <- xy(cx,ring_y,10.16,phi)
      txt(p$x,p$y,groups[g],5.75,muted,just="centre",rot=tangent(phi))
    }
    for (i in seq_along(eps)) {
      phi <- centre_angle[i]
      p <- xy(cx,ring_y,outer+2.0,phi)
      txt(p$x,p$y,rim_names[i],6.55,ink,if(eps[i] %in% c("ipf","pe","ph")) 3 else 1,
          just="centre",rot=tangent(phi))
      if (eps[i] %in% omit) {
        # The original panel C excluded this outcome for its low event count.
        # Preserve the display decision without conflating it with an unestimated model.
        wedge(cx,ring_y,inner+.08,outer-.08,a0[i]+.33,a1[i]-.33,"#F9FAFA","#8B989A",.16,"31")
        p <- xy(cx,ring_y,(inner+outer)/2,phi)
        txt(p$x,p$y,"Omitted\n73 events",5.95,muted,just="centre",rot=tangent(phi))
      } else {
        for (j in seq_len(ncol(m$hr))) {
          r1 <- outer-(j-1)*pitch; r0 <- r1-pitch
          v <- m$hr[i,j]
          fill <- if(is.na(v)) missing_fill else heat_colour(v)
          wedge(cx,ring_y,r0+.065,r1-.065,a0[i]+.21,a1[i]-.21,fill)
          if (isTRUE(m$sig[i,j])) {
            p <- xy(cx,ring_y,(r0+r1)/2,phi)
            txt(p$x,p$y+.22,"*",6.45,symbol_colour(fill),just="centre")
          }
        }
      }
    }
    for (j in seq_len(n_rings))
      txt(cx,ring_y-outer+(j-.5)*pitch,as.character(j),if(main_view) 6.6 else 5.35,muted,just="centre")
  }

  keys <- function() {
    txt(3.1,175.1,"Analyses from the outer to the inner ring",7.8,face=2)
    key <- function(x,y,n,s) {
      txt(x,y,as.character(n),6.8,ink,2)
      txt(x+3.25,y,s,6.85,ink)
    }
    y <- c(179.8,184.2,188.6)
    key(3.1,y[1],1,"Main model (B)")
    key(3.1,y[2],2,"Add lung age (model C)")
    if (main_view) {
      key(3.1,y[3],3,"All seven systems (model D)")
      key(80.0,y[1],4,"Panel-specific analysis under model C")
      key(83.25,y[2],"B","Renal clock retrained without cystatin C (post hoc)")
      key(83.25,y[3],"C","Additional adjustment for cystatin C-based eGFR")
    } else {
      key(3.1,y[3],3,"Model B + PhenoAge / organismal age")
      key(63.0,y[1],4,"All seven systems (model D)")
      key(63.0,y[2],5,"Model C + pack-years")
      key(63.0,y[3],6,"Never smokers (model C)")
      key(118.4,y[1],7,"2-year landmark (model C)")
      key(118.4,y[2],8,"Fine\u2013Gray subdistribution HR (model B)")
      key(118.4,y[3],9,"Renal clock without cystatin C (post hoc, model C)")
      txt(121.65,192.05,"Cardiovascular model C + cystatin C-based eGFR",6.5,muted)
    }
    line(c(3.1,179.2),rep(194.3,2),rule,.14)
    txt(3.1,197.0,"HR per 1 SD",7.1,face=2)
    bar_x <- 20.7; bar_w <- 43.0; bar_y <- 196.0; bar_h <- 1.95
    ticks <- c(.8,1,1.25,1.6,1.8)
    bounds <- log(range(c(ticks,B$hr,C$hr),na.rm=TRUE))
    edge <- seq(bounds[1],bounds[2],length.out=113)
    fills <- heat_colour(exp((head(edge,-1)+tail(edge,-1))/2))
    for (i in seq_along(fills)) rect(bar_x+(i-1)*bar_w/length(fills),bar_y,bar_w/length(fills)+.008,bar_h,fills[i])
    for (i in seq_along(ticks)) {
      at <- bar_x+bar_w*(log(ticks[i])-bounds[1])/diff(bounds)
      line(rep(at,2),c(bar_y+bar_h,bar_y+bar_h+.55),axis_col,.12)
      txt(at,bar_y+bar_h+2.2,c("0.80","1.00","1.25","1.60","1.80")[i],6.0,muted,
          just=if(i==1) "left" else "centre")
    }
    txt(71.4,197.8,"*",8.1,ink,2)
    txt(74.2,197.2,"95% CI excludes 1",6.85)
    if (anyNA(B$hr) || anyNA(C$hr[!rownames(C$hr) %in% omitted_c,,drop=FALSE])) {
      rect(100.95,196.3,2.05,1.9,missing_fill)
      txt(104.3,197.2,"No estimate",6.85)
    }
    txt(3.1,205.0,"PVD, pulmonary vascular disease. PE, pulmonary embolism. PH, pulmonary hypertension. eGFR, estimated glomerular filtration rate.",6.45,muted)
  }
  figure <- function() {
    grid::grid.newpage()
    rect(0,0,W,H,"white")
    draw_forest()
    line(c(3.1,179.2),rep(91.3,2),rule,.16)
    header(3.1,95.9,"B","Clinical renal ageing")
    header(94.0,95.9,"C","Proteomic cardiovascular (artery-dominated) ageing")
    draw_ring(B,46.2)
    draw_ring(C,137.0,omitted_c)
    keys()
  }
  stem <- file.path(out_dir,if(main_view) "Figure4" else "Figure4_full")
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=family,bg="white")
  figure(); invisible(grDevices::dev.off())
  grDevices::tiff(paste0(stem,".tiff"),width=W,height=H,units="mm",res=600,
                  compression="lzw",type="cairo",bg="white",antialias="gray")
  figure(); invisible(grDevices::dev.off())
  grDevices::png(paste0(stem,".png"),width=W,height=H,units="mm",res=450,
                 type="cairo",bg="white",antialias="gray")
  figure(); invisible(grDevices::dev.off())
  write_figure_svg(stem, W, H, family, figure)
  invisible(NULL)
}
read_figure4_data <- function(file) {
  read.delim(file,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8")
}
source_data <- read_figure4_data(file.path(input_dir, "associations_and_sensitivity.tsv"))
render_figure4(source_data, output_dir, view="main")
