# Shared drawing functions for Supplementary Figures S2 and S3.
library(grid)

family_name <- "Times New Roman"
page_width <- 183
ink <- "#20292D"
slate <- "#495963"
rule <- "#D0D9DC"
axis_ink <- "#65747B"
clinical_blue <- "#395688"
protein_teal <- "#168372"
joint_brown <- "#806448"
white <- "#FFFFFF"
endpoints <- c("copd","asthma","bronchiectasis","ild","ipf","pvd","pe","ph","osa","lung_cancer","pneumonia")
diseases <- setNames(c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD (J84.1)",
                       "Pulmonary vascular disease","Pulmonary embolism","Pulmonary hypertension",
                       "Sleep apnoea","Lung cancer","Pneumonia"), endpoints)
systems <- c("pulmonary","cardiovascular","renal","hepatic","immune","musculoskeletal","metabolic")
system_headers <- c("Lung","CV","Renal","Liver","Immune","MSK","Metab.")
models <- c("7 clinical","7 protein","All 14")
model_ids <- c("+clin7","+prot7","+14_ridge")
model_ink <- c(clinical_blue,protein_teal,joint_brown)
model_shapes <- c("circle","diamond","square")
cross_order <- c("copd","asthma","pneumonia","bronchiectasis","ild","ipf","pvd","pe","ph","osa","lung_cancer")

read_tsv <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE,
                                      fileEncoding="UTF-8", na.strings=c("","NA"))
lookup <- function(keys, available) {
  i <- match(keys,available)
  stopifnot(!anyDuplicated(available),!anyNA(i))
  i
}
same_numbers <- function(a,b) isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-9))
number <- function(x, digits=2L) {
  y <- sprintf(paste0("%.",digits,"f"),x)
  y <- sub("^-0([.]0+)?$","0\\1",y)
  sub("^-","−",y)
}
mm <- function(x) unit(x,"mm")
yy <- function(y) unit(1,"npc")-mm(y)
txt <- function(x,y,s,size=7.2,colour=ink,face=1,align="left",rot=0) {
  grid.text(s,mm(x),yy(y),just=c(align,"centre"),rot=rot,
            gp=gpar(fontfamily=family_name,fontsize=size,fontface=face,col=colour,lineheight=1.05))
}
stroke <- function(x,y,colour=rule,width=.15) {
  grid.lines(mm(x),yy(y),gp=gpar(col=colour,lwd=width/.2646,lineend="butt",linejoin="round"))
}
block <- function(x,y,w,h,fill) {
  grid.rect(mm(x),yy(y),mm(w),mm(h),just=c("left","top"),gp=gpar(fill=fill,col=NA))
}
symbol_at <- function(x,y,r,colour,shape="circle",fill=colour,width=.17) {
  gp <- gpar(col=colour,fill=fill,lwd=width/.2646,linejoin="round")
  if (shape=="diamond") {
    grid.polygon(mm(x+c(0,r,0,-r)),yy(y+c(-r,0,r,0)),gp=gp)
  } else if (shape=="square") {
    grid.rect(mm(x),yy(y),mm(1.65*r),mm(1.65*r),gp=gp)
  } else {
    grid.circle(mm(x),yy(y),r=mm(r),gp=gp)
  }
}
page <- function() {
  grid.newpage()
  grid.rect(gp=gpar(fill=white,col=NA))
}
tag <- function(x,y,s) txt(x,y,s,10.4,face=2)
scale_to <- function(z,lo,hi,a,b) a+(z-lo)*(b-a)/(hi-lo)

load_cross <- function(figure_dir,pred_dir) {
  archived <- read_tsv(file.path(figure_dir,"Figure6_data.tsv"))
  gain <- archived[archived$panel=="c",c("endpoint","model","src","est","lo","hi")]
  hr <- archived[archived$panel=="d",c("endpoint","system","HR_clin_both","p_clin_both",
                                        "HR_prot_both","p_prot_both")]
  boot <- read_tsv(file.path(pred_dir,"q4_prediction_bootstrap_v1.tsv"))
  boot <- boot[boot$from=="base" & boot$to %in% model_ids & boot$measure=="IPA points" &
                 as.logical(boot$meets_T100_allFU),]
  dev <- read_tsv(file.path(pred_dir,"q4_prediction_Dcv_v1.tsv"))
  dev <- dev[dev$from=="base" & dev$to %in% model_ids,]
  source_hr <- read_tsv(file.path(pred_dir,"q4_association_v1.tsv"))
  keys <- function(a,b) paste(a,b,sep="/")
  test <- gain$src=="Test set"
  stopifnot(nrow(gain)==33L,nrow(hr)==77L,nrow(source_hr)==77L,sum(test)==9L,
            setequal(keys(gain$endpoint,gain$model),as.vector(outer(endpoints,models,keys))),
            setequal(keys(hr$endpoint,hr$system),as.vector(outer(endpoints,systems,keys))),
            all(startsWith(gain$src[!test],"Cross-validated in D")),
            all(is.na(gain$lo[!test]) & is.na(gain$hi[!test])))
  ids <- model_ids[match(gain$model,models)]
  j <- lookup(keys(gain$endpoint[test],ids[test]),keys(boot$endpoint,boot$to))
  stopifnot(all(vapply(c("est","lo","hi"),function(n) same_numbers(gain[[n]][test],boot[[n]][j]),logical(1))))
  j <- lookup(keys(gain$endpoint[!test],ids[!test]),keys(dev$endpoint,dev$to))
  stopifnot(same_numbers(gain$est[!test],dev$gain_DCV[j]),all(dev$folds[j]==10L))
  j <- lookup(keys(hr$endpoint,hr$system),keys(source_hr$endpoint,source_hr$system))
  for (n in c("HR_clin_both","p_clin_both","HR_prot_both","p_prot_both")) {
    stopifnot(same_numbers(hr[[n]],source_hr[[n]][j]))
  }
  list(gain=gain, hr=hr)
}

hr_stops <- c(.8,1,1.25,1.6,2,2.6)
hr_colours <- c("#A8CBDF","#FAFAF7","#FAD7C3","#F39A7D","#D66052","#98263B")
hr_fill <- function(hr) {
  rgb <- col2rgb(hr_colours)/255
  z <- log(hr)
  channels <- lapply(1:3,function(i) approx(log(hr_stops),rgb[i,],xout=z,rule=2)$y)
  grDevices::rgb(channels[[1]],channels[[2]],channels[[3]])
}
cell_text_colour <- function(fill) {
  srgb <- col2rgb(fill)/255
  linear <- ifelse(srgb<=.04045,srgb/12.92,((srgb+.055)/1.055)^2.4)
  luminance <- as.vector(c(.2126,.7152,.0722) %*% linear)
  ifelse(luminance<.245,white,ink)
}

draw_cross <- function(data) {
  page()
  tag(3.1,5.35,"A")
  starts <- c(36.8,85.3,133.8)
  ends <- starts+44.1
  plot0 <- starts+.65
  plot1 <- starts+33.0
  for (j in seq_along(models)) {
    mid <- (starts[j]+ends[j])/2
    txt(mid,5.35,c("Clinical","Proteomic","Clinical + proteomic")[j],8.85,model_ink[j],2,"centre")
    txt(mid,10.45,c("7 organ ages","7 organ ages","14 organ ages · ridge")[j],7.05,slate,align="centre")
    stroke(c(starts[j],ends[j]),c(13.6,13.6))
  }
  txt(3.45,10.45,"Respiratory outcome",7.65,face=2)
  section <- function(y,name,detail) {
    block(3.15,y-2.5,176.6,5.0,"#EBF0F2")
    txt(4.35,y,name,7.75,face=2)
    txt(36.8,y,detail,7.0,slate)
  }
  section(18.0,"Locked test set","Point estimates and 95% bootstrap CIs")
  section(40.8,"Development set","Means from two repeats of fivefold cross-validation")
  cy <- c(23.6,29.2,34.8,46.5+5.25*0:7)
  subtype <- c("ipf","pe","ph")
  for (i in seq_along(cross_order)) {
    e <- cross_order[i]
    if (i %in% c(2,5,7,9,11)) block(3.15,cy[i]-2.35,176.6,4.7,"#F6F8F8")
    txt(if(e %in% subtype) 5.5 else 3.75,cy[i],diseases[e],
        if(e %in% subtype) 6.95 else 7.15,
        if(e %in% subtype) slate else ink,if(e %in% subtype) 3 else 1)
  }
  for (j in seq_along(models)) {
    rows <- data$gain[data$gain$model==models[j],]
    rows <- rows[lookup(cross_order,rows$endpoint),]
    stopifnot(all((rows$src=="Test set")==seq_along(cross_order)%in%1:3),
              all(rows$est>=-2 & rows$est<=6),all(rows$lo[1:3]>=-2 & rows$hi[1:3]<=6))
    xpos <- function(v) scale_to(v,-2,6,plot0[j],plot1[j])
    for (span in list(c(21.3,37.1),c(44.2,85.4))) stroke(rep(xpos(0),2),span,"#AAB8BD",.15)
    for (i in seq_along(cy)) {
      if (is.finite(rows$lo[i])) {
        ends_ci <- xpos(c(rows$lo[i],rows$hi[i]))
        stroke(ends_ci,rep(cy[i],2),model_ink[j],.28)
        for (x in ends_ci) stroke(rep(x,2),cy[i]+c(-.53,.53),model_ink[j],.18)
      }
      symbol_at(xpos(rows$est[i]),cy[i],.64,model_ink[j],model_shapes[j])
      txt(ends[j]-.2,cy[i],number(rows$est[i]),7.05,model_ink[j],align="right")
    }
    stroke(c(plot0[j],plot1[j]),rep(87.2,2),axis_ink,.17)
    for (v in seq(-2,6,2)) {
      x <- xpos(v)
      stroke(rep(x,2),c(87.2,88.05),axis_ink,.14)
      txt(x,90.15,number(v,0),6.4,slate,align="centre")
    }
  }
  txt(107.5,95.45,"10-year ΔIPA versus the base model (points) in the shared sample",7.7,align="centre")
  stroke(c(3.15,179.75),c(100.1,100.1),rule,.18)

  tag(3.1,105.5,"B")
  left <- c(37.0,109.4)
  cw <- 10.0
  for (j in 1:2) {
    txt(left[j]+3.5*cw,105.5,c("Clinical","Proteomic")[j],8.85,model_ink[j],2,"centre")
    stroke(c(left[j],left[j]+7*cw),rep(109.2,2))
    for(k in seq_along(systems)) txt(left[j]+(k-.5)*cw,113.0,system_headers[k],6.95,align="centre")
  }
  txt(8.25,113.0,"Respiratory outcome",7.4,face=2)
  heights <- c(3L,2L,3L,3L)
  groups <- c("Airway","Parenchymal","Vascular","Other")
  stripe <- c("#E3E7F0","#DEEEE9","#F4E4DF","#EEE7DE")
  rh <- 5.05
  top <- 116.25
  for (g in seq_along(heights)) {
    nr <- heights[g]
    start <- if(g==1L) 1L else 1L+sum(heights[seq_len(g-1L)])
    eps <- endpoints[start+seq_len(nr)-1L]
    block(3.15,top+.13,3.6,nr*rh-.26,stripe[g])
    glab <- if (g==2L) "Parench." else groups[g]
    txt(4.95,top+nr*rh/2,glab,6.05,slate,align="centre",rot=90)
    if(g %in% c(2,3)) {
      stroke(rep(8.65,2),c(top+rh/2+1.10,top+(nr-.5)*rh),"#A1AFB3",.12)
      for (i in 2:nr) stroke(c(8.65,9.75),rep(top+(i-.5)*rh,2),"#A1AFB3",.12)
    }
    for (i in seq_along(eps)) {
      e <- eps[i]; y <- top+(i-.5)*rh
      sub <- e %in% subtype
      txt(if(sub) 10.2 else 8.15,y,diseases[e],if(sub) 6.75 else 6.9,
          if(sub) slate else ink,if(sub) 3 else 1)
      d <- data$hr[data$hr$endpoint==e,]
      d <- d[lookup(systems,d$system),]
      for (j in 1:2) {
        hr <- if(j==1) d$HR_clin_both else d$HR_prot_both
        pv <- if(j==1) d$p_clin_both else d$p_prot_both
        stopifnot(all(hr>=min(hr_stops) & hr<=max(hr_stops)),all(is.finite(pv)))
        fill <- hr_fill(hr)
        fg <- cell_text_colour(fill)
        for(k in seq_along(systems)) {
          x <- left[j]+(k-1)*cw
          block(x+.08,y-rh/2+.10,cw-.16,rh-.20,fill[k])
          txt(x+cw/2,y,paste0(number(hr[k]),if(pv[k]<.05) "*" else ""),
              6.85,fg[k],align="centre")
        }
      }
    }
    top <- top+nr*rh+1.3
  }
  stroke(c(3.15,179.75),c(178.5,178.5))
  txt(3.3,183.0,"Mutually adjusted HR per SD",7.6,face=2)
  gx <- 45.15; gw <- 55.6
  fractions <- (log(hr_stops)-log(min(hr_stops)))/diff(log(range(hr_stops)))
  fill <- linearGradient(hr_colours,stops=fractions,x1=0,y1=.5,x2=1,y2=.5)
  block(gx,181.95,gw,1.95,fill)
  for (k in seq_along(hr_stops)) {
    x <- gx+fractions[k]*gw
    stroke(rep(x,2),c(184.0,184.65),axis_ink,.12)
    txt(x,186.50,format(hr_stops[k],trim=TRUE),6.45,slate,align="centre")
  }
  txt(113,183.0,"* Nominal two-sided P < 0.05",7.05)
}

load_calibration <- function(pred_dir, summary_path) {
  s <- read_tsv(summary_path)
  s <- s[as.logical(s$meets_T100_allFU) & s$endpoint %in% endpoints,]
  d <- read_tsv(file.path(pred_dir,"q3_calibration_v1.tsv"))
  key <- function(a,b) paste(a,b,sep="/")
  d <- d[d$model=="base+all7" & key(d$layer,d$endpoint) %in% key(s$layer,s$endpoint),]
  stopifnot(nrow(s)==20L,nrow(d)==200L,all(s$events_T>=100L),
            !anyDuplicated(key(s$layer,s$endpoint)),all(d$pred10>=0 & d$obs10>=0),
            all(d$pred10<=100 & d$obs10<=100),all(s$events_T_10y<=s$events_T))
  for (i in seq_len(nrow(s))) {
    z <- d[d$endpoint==s$endpoint[i] & d$layer==s$layer[i],]
    z <- z[order(z$decile),]
    stopifnot(identical(as.integer(z$decile),1:10),all(diff(z$pred10)>=0),
              sum(z$n)==s$N_T[i],all(z$n>0))
  }
  # Only the clinical layer meets the prespecified full-follow-up threshold for these two outcomes.
  stopifnot(setequal(s$endpoint[s$layer=="clinical"],endpoints),
            setequal(s$endpoint[s$layer=="protein"],setdiff(endpoints,c("ild","ipf"))))
  list(summary=s, deciles=d)
}

risk_ticks <- function(max_risk) {
  bounds <- c(1,2,4,6,8,10,16,20,30,40,60,80,100)
  upper <- bounds[which(bounds>=max_risk*1.018)[1L]]
  stopifnot(is.finite(upper))
  segments <- if (upper %in% c(6,30,60)) 3 else if (upper %in% c(10,20,100)) 5 else 4
  seq(0,upper,length.out=segments+1)
}

calibration_key <- function(x,y) {
  stroke(c(x,x+50.2),rep(y+3.8,2),rule,.16)
  txt(x,y+7.9,"Locked test set",8.1,face=2)
  stroke(c(x+.2,x+4.7),rep(y+13.0,2),clinical_blue,.24)
  symbol_at(x+2.45,y+13.0,.6,clinical_blue)
  txt(x+6.8,y+13.0,"Clinical clocks",7.25)
  stroke(c(x+.2,x+4.7),rep(y+17.55,2),protein_teal,.24)
  symbol_at(x+2.45,y+17.55,.78,protein_teal,"diamond",NA)
  txt(x+6.8,y+17.55,"Proteomic clocks",7.25)
  stroke(c(x+.2,x+4.7),c(y+22.95,y+21.15),"#9DABAF",.23)
  txt(x+6.8,y+22.05,"Equal observed and predicted risk",7.15)
  txt(x,y+26.85,"Points show deciles of predicted risk.",6.95,slate)
  txt(x,y+33.0,"O/E (events)",7.45,face=2)
  txt(x,y+37.35,"Observed-to-expected 10-year risk ratio",6.95,slate)
  txt(x,y+41.65,"Parentheses give the 10-year event count.",6.95,slate)
  txt(x,y+46.35,"—   Fewer than 100 events over full follow-up",6.85,slate)
}

draw_calibration <- function(data) {
  page()
  # The cells follow the airway, parenchymal, vascular and other outcome rows of Figure 2.
  tile_x <- c(8.0,65.6,123.2)
  tile_y <- c(1.8,53.5,105.2,156.9)
  slots <- c(1:5,7:12)
  for (y in tile_y[-1L]-.65) stroke(c(7.9,179.6),c(y,y),rule,.13)
  txt(3.5,102.5,"Observed 10-year risk (%)",8.2,align="centre",rot=90)
  txt(94.6,211.5,"Predicted 10-year risk (%)",8.2,align="centre")
  calibration_key(tile_x[3]+2.1,tile_y[2]+.2)
  for (i in seq_along(endpoints)) {
    e <- endpoints[i]
    slot <- slots[i]
    x <- tile_x[(slot-1)%%3+1L]
    y <- tile_y[(slot-1)%/%3+1L]
    tag(x,y+3.7,LETTERS[i])
    txt(x+5.0,y+3.7,diseases[e],8.3,face=2)
    z <- data$deciles[data$deciles$endpoint==e,]
    s <- data$summary[data$summary$endpoint==e,]
    ticks <- risk_ticks(max(z$obs10,z$pred10))
    high <- max(ticks)
    size <- 32.9; l <- x+10.25; t <- y+9.1; b <- t+size
    px <- function(v) l+size*v/high
    py <- function(v) b-size*v/high
    stroke(c(l,l+size),c(b,t),"#ACB8BD",.20)
    stroke(c(l,l+size),rep(b,2),axis_ink,.18)
    stroke(rep(l,2),c(t,b),axis_ink,.18)
    for (a in ticks) {
      lab <- format(a,trim=TRUE,scientific=FALSE)
      stroke(rep(px(a),2),c(b,b+.68),axis_ink,.14)
      txt(px(a),b+2.60,lab,6.25,slate,align="centre")
      stroke(c(l-.68,l),rep(py(a),2),axis_ink,.14)
      txt(l-1.35,py(a),lab,6.25,slate,align="right")
    }
    # Draw all polylines before symbols so one series cannot cover the other's symbols.
    layers <- c("clinical","protein")
    series <- lapply(layers,function(layer) {
      a <- z[z$layer==layer,]
      a[order(a$decile),]
    })
    for (j in seq_along(layers)) {
      a <- series[[j]]
      if (nrow(a)) stroke(px(a$pred10),py(a$obs10),c("#6D86AA","#55A393")[j],.32)
    }
    for (j in seq_along(layers)) {
      a <- series[[j]]
      colour <- c(clinical_blue,protein_teal)[j]
      for(k in seq_len(nrow(a))) {
        symbol_at(px(a$pred10[k]),py(a$obs10[k]),c(.48,.70)[j],colour,
                  c("circle","diamond")[j],if (j==1L) colour else NA,.20)
      }
      stat <- s[s$layer==layers[j],]
      sx <- x+c(3.2,29.5)[j]
      if(nrow(stat)) {
        symbol_at(sx,y+48.65,c(.51,.68)[j],colour,c("circle","diamond")[j],
                  if(j==1L) colour else NA,.18)
        label <- sprintf("%s (%s)",number(stat$O_E_all7),formatC(stat$events_T_10y,format="d",big.mark=","))
        txt(sx+1.8,y+48.65,label,7.05,colour)
      } else {
        symbol_at(sx,y+48.65,.68,"#9CA8AB","diamond",NA,.16)
        txt(sx+1.8,y+48.65,"—",7.05,slate)
      }
    }
  }
}

export_plot <- function(stem,height,draw) {
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=page_width/25.4,height=height/25.4,
                       family=family_name,bg=white)
  tryCatch(draw(),finally=invisible(grDevices::dev.off()))
  # Align RGB scan lines to four pixels for reproducible Windows Cairo output.
  px <- function(dpi) 4L*round(page_width*dpi/(25.4*4))
  grDevices::png(paste0(stem,".png"),width=px(300),height=round(height*300/25.4),
                 res=300,type="cairo-png",antialias="gray",bg=white)
  tryCatch(draw(),finally=invisible(grDevices::dev.off()))
  preview <- paste0(stem,".png")
  # Cairo's PNG writer omits physical resolution. Set it without resampling the pixels.
  pixels <- png::readPNG(preview,native=TRUE)
  png::writePNG(pixels,preview,dpi=c(300,300))
  grDevices::tiff(paste0(stem,".tiff"),width=px(600),height=round(height*600/25.4),
                  res=600,type="cairo",antialias="gray",compression="lzw",bg=white)
  tryCatch(draw(),finally=invisible(grDevices::dev.off()))
}

