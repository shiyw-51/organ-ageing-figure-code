# Figure6
# Rscript --vanilla Figure6/Figure6.R [output_directory]
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
input_dir <- file.path(path.expand(input_root), "Figure6")
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir, call.=FALSE)
output_dir <- if (length(args)) path.expand(args[[1L]]) else file.path(dirname(script_dir), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
library(grid)
out_dir <- output_dir
pred_dir <- input_dir

read_tab <- function(file) read.delim(file,stringsAsFactors=FALSE,check.names=FALSE,
                                     fileEncoding="UTF-8",na.strings=c("","NA"))
pred <- function(name) read_tab(file.path(pred_dir,paste0(name,"_v1.tsv")))
summary_t <- pred("q3_summary")
boot <- pred("q3_bootstrap")
archive <- read_tab(file.path(pred_dir,"Figure6_data.tsv"))

layers <- c("clinical","protein")
ep <- c("copd","asthma","bronchiectasis","ild","ipf","pneumoconiosis",
        "pvd","pe","ph","osa","lung_cancer","pneumonia")
diagnosis <- c("COPD","Asthma","Bronchiectasis","ILD","Fibrotic ILD (J84.1)","Pneumoconiosis",
               "Pulmonary vascular disease","Pulmonary embolism","Pulmonary hypertension",
               "Sleep apnoea","Lung cancer","Pneumonia")
pair_key <- function(layer,endpoint) paste(layer,endpoint,sep="/")
quantity_key <- function(layer,endpoint,qty) paste(layer,endpoint,qty,sep="/")
near <- function(a,b) isTRUE(all.equal(unname(a),unname(b),tolerance=1e-9,check.attributes=FALSE))
stopifnot(!anyDuplicated(summary_t[c("layer","endpoint")]),!anyDuplicated(boot[c("layer","endpoint","quantity")]))
ix <- match(as.vector(outer(ep,layers,function(e,l) pair_key(l,e))),pair_key(summary_t$layer,summary_t$endpoint))
stopifnot(length(ix)==24L,!anyNA(ix),nrow(summary_t)==length(ix))
dat <- summary_t[ix,]
dat$exploratory <- dat$events_T<100L
dat$available <- is.finite(dat$AUC_all7)
stopifnot(all(dat$events_T_10y<=dat$events_T),
          identical(as.logical(dat$meets_T100_allFU),!dat$exploratory),
          identical(dat$available,dat$events_T>=20L))
bk <- quantity_key(boot$layer,boot$endpoint,boot$quantity)
quant <- c("dAUC_lung","dAUC_others_after_lung","dAUC_all7","gain_lung","gain_others_after_lung","gain_all7")
for (q in quant) {
  j <- match(quantity_key(dat$layer,dat$endpoint,q),bk)
  stopifnot(identical(!is.na(j),dat$available),near(dat[[q]],boot$est[j]))
  dat[[paste0(q,"_lo")]] <- boot$lo[j]
  dat[[paste0(q,"_hi")]] <- boot$hi[j]
  stopifnot(all(boot$lo[j[dat$available]]<=boot$hi[j[dat$available]]))
}
stopifnot(near(dat$AUC_lung-dat$AUC_base,dat$dAUC_lung),
          near(dat$AUC_all7-dat$AUC_lung,dat$dAUC_others_after_lung),
          near(dat$dAUC_lung+dat$dAUC_others_after_lung,dat$dAUC_all7),
          near(dat$IPA_all7-dat$IPA_base,dat$gain_all7),
          near(dat$gain_lung+dat$gain_others_after_lung,dat$gain_all7))
old_a <- archive[archive$panel=="a",]
j <- match(pair_key(old_a$layer,old_a$endpoint),pair_key(dat$layer,dat$endpoint))
stopifnot(nrow(old_a)==sum(dat$available),!anyNA(j),
          near(old_a$AUC_base,dat$AUC_base[j]),near(old_a$AUC_lung,dat$AUC_lung[j]),
          near(old_a$AUC_all7,dat$AUC_all7[j]))

ink <- "#20272B"; muted <- "#4B585E"; axis_col <- "#546168"; rule <- "#D4DCDD"
dark <- c(clinical="#395688",protein="#168372")
light <- c(clinical="#AABDD6",protein="#9BC8BB")
groups <- c("Airway","Parenchymal","Vascular","Other")
group_fill <- c("#E3E7F0","#DEEEE9","#F4E4DF","#EEE7DE")
W <- 183; font <- "Times New Roman"
mm <- function(x) unit(x,"mm")
ym <- function(y) unit(1,"npc")-mm(y)
txt <- function(x,y,s,size=7,col=ink,face=1,just="left",rot=0) {
  grid.text(s,x=mm(x),y=ym(y),just=c(just,"centre"),rot=rot,
            gp=gpar(fontfamily=font,fontface=face,fontsize=size,col=col,lineheight=.94))
}
line <- function(x,y,col=rule,w=.16) {
  grid.lines(mm(x),ym(y),gp=gpar(col=col,lwd=w/.2646,lineend="butt",linejoin="round"))
}
box <- function(x,y,w,h,fill) {
  grid.rect(mm(x),ym(y),width=mm(w),height=mm(h),just=c("left","top"),gp=gpar(col=NA,fill=fill))
}
mark <- function(x,y,r,col,shape="circle",fill=col) {
  gp <- gpar(col=col,fill=fill,lwd=.15/.2646)
  if(shape=="diamond")
    grid.polygon(mm(x+c(0,r,0,-r)),ym(y+c(-r,0,r,0)),gp=gp)
  else if(shape=="square")
    grid.rect(mm(x),ym(y),width=mm(1.65*r),height=mm(1.65*r),gp=gp)
  else
    grid.circle(mm(x),ym(y),r=mm(r),gp=gp)
}
net_tip <- function(x,y) {
  grid.polygon(mm(x+c(0,-.60,.60)),ym(y+c(0,1.00,1.00)),gp=gpar(fill=ink,col=NA))
}
head_panel <- function(x,y,tag,s,col=ink) {
  txt(x,y,tag,10.35,face=2)
  txt(x+5.15,y+.06,s,8.8,col,2)
}
lin <- function(v,d,range) range[1]+(v-d[1])*diff(range)/diff(d)
num <- function(x,n) {
  s <- sprintf(paste0("%.",n,"f"),x)
  s <- sub("^-0([.]0+)?$","0\\1",s)
  sub("^-","−",s)
}
ci_str <- function(v,lo,hi,n) sprintf("%s (%s–%s)",num(v,n),num(lo,n),num(hi,n))
row_tops <- function(y,h,gap) {
  i <- seq_along(ep)-1L
  y+i*h+(i%/%3L)*gap
}
disease_axis <- function(tops,h) {
  centers <- tops+h/2
  for (g in seq_along(groups)) {
    t <- tops[3L*g-2L]
    box(3.15,t+.12,3.60,h*3-.24,group_fill[g])
    txt(4.95,t+h*1.5,groups[g],6.0,muted,just="centre",rot=90)
  }
  for (branch in list(c(4L,5L),c(7L,8L,9L))) {
    line(rep(8.5,2),c(centers[branch[1]]+1.25,centers[tail(branch,1)]),"#A1AEB1",.13)
    for (i in branch[-1]) line(c(8.5,9.78),rep(centers[i],2),"#A1AEB1",.13)
  }
  for (i in seq_along(ep)) {
    sub <- ep[i] %in% c("ipf","pe","ph")
    txt(if(sub) 10.55 else 8.25,centers[i],diagnosis[i],if(sub) 6.55 else 6.75,
        if(sub) muted else ink,if(sub) 3 else 1)
  }
}
shade_groups <- function(tops,h) {
  for (i in c(2L,4L)) {
    t <- tops[3*i-2L]
    box(7.25,t,172.55,3*h,"#F5F7F7")
  }
}

# The upper panels and the lower panels each have their own axes and outcome rows.
upper_height <- 99
auc_models <- c("AUC_base","AUC_lung","AUC_all7")
auc_ticks <- seq(.50,.90,by=.05)
auc_y_offsets <- c(-1.43,0,1.43)

auc_symbol <- function(x,y,model,layer) {
  shape <- c("circle","square","diamond")[model]
  radius <- c(.62,.69,.75)[model]
  stroke <- if(model==1L) muted else dark[[layer]]
  fill <- c("white",light[[layer]],dark[[layer]])[model]
  mark(x,y,radius,stroke,shape,fill)
}

absolute_panels <- function() {
  head_panel(3.1,5.05,"A","Clinical clocks",dark[[1L]])
  head_panel(106.4,5.05,"B","Proteomic clocks",dark[[2L]])
  row_height <- 5.55
  tops <- row_tops(15.4,row_height,1.10)
  yy <- tops+row_height/2
  bottom <- tail(tops,1)+row_height
  shade_groups(tops,row_height)
  disease_axis(tops,row_height)
  txt(8.25,11.55,"Respiratory outcome",7.05,face=2)

  for(ly in layers) {
    left <- if(ly=="clinical") 37.45 else 111.1
    right <- left+55.5
    events_x <- left+66.8
    scale_auc <- function(v) lin(v,range(auc_ticks),c(left,right))
    txt((left+right)/2,11.55,"Absolute 10-year AUC",7.05,face=2,just="centre")
    txt(events_x,10.05,"Events",6.85,face=2,just="right")
    txt(events_x,13.25,"(10 y)",6.30,muted,just="right")

    for(j in seq_along(auc_ticks)) {
      tick <- auc_ticks[j]
      x <- scale_auc(tick)
      grid_col <- if(j%%2L==1L) "#D6DEDF" else "#E5EAEA"
      line(rep(x,2),c(tops[1],bottom),grid_col,.13)
      line(rep(x,2),c(bottom+.85,bottom+1.54),axis_col,.15)
      txt(x,bottom+3.56,num(tick,2),6.45,muted,just="centre")
    }
    line(c(left,right),rep(bottom+.85,2),axis_col,.18)
    layer_data <- dat[dat$layer==ly,]
    stopifnot(identical(as.character(layer_data$endpoint),ep))
    for(i in seq_len(nrow(layer_data))) {
      a <- layer_data[i,]
      event_count <- formatC(a$events_T_10y,format="d",big.mark=",")
      txt(events_x,yy[i],event_count,6.60,if(a$available) ink else muted,just="right")
      if(!a$available) {
        txt(scale_auc(.70),yy[i],paste0("Not estimated (",a$events_T," events)†"),
            6.65,muted,just="centre")
        next
      }
      values <- as.numeric(unlist(a[auc_models],use.names=FALSE))
      stopifnot(all(is.finite(values)),all(values>=min(auc_ticks) & values<=max(auc_ticks)))
      xs <- scale_auc(values)
      ys <- yy[i]+auc_y_offsets
      line(xs,ys,"#8F9DA1",.22)
      for(model in seq_along(auc_models)) auc_symbol(xs[model],ys[model],model,ly)
      if(a$exploratory) txt(xs[3L]+1.35,ys[3L]-.43,"†",6.35,muted)
    }
  }
  legend_y <- bottom+8.54
  auc_symbol(38.0,legend_y,1L,"clinical")
  txt(40.0,legend_y,"Base model",6.85)
  auc_symbol(70.3,legend_y,2L,"clinical")
  auc_symbol(72.1,legend_y,2L,"protein")
  txt(74.1,legend_y,"Base + lung age",6.85)
  auc_symbol(115.3,legend_y,3L,"clinical")
  auc_symbol(117.2,legend_y,3L,"protein")
  txt(119.3,legend_y,"Base + all seven organ ages",6.85)
  line(c(3.15,179.8),rep(upper_height-1.8,2),rule,.15)
}

gain_panels <- function() {
  head_panel(3.1,5.05,"C","Clinical clocks",dark[1])
  head_panel(106.4,5.05,"D","Proteomic clocks",dark[2])
  tops <- row_tops(22.2,5.70,1.15)
  yy <- tops+5.70/2
  end <- tail(tops,1)+5.70
  shade_groups(tops,5.70)
  disease_axis(tops,5.70)
  txt(8.25,16.2,"Respiratory outcome",7.05,face=2)
  domain <- c(-.020,.085)
  for (ly in layers) {
    x0 <- if(ly=="clinical") 36.6 else 110.10
    l <- x0; r <- x0+35.50; cx <- x0+61.35; ipa <- x0+69.55
    sx <- function(v) lin(v,domain,c(l,r))
    txt((l+r)/2,13.55,"Change in 10-year AUC",6.95,face=2,just="centre")
    txt((l+r)/2,17.35,"Locked test set",6.45,muted,just="centre")
    txt(cx-10.85,13.55,"Total ΔAUC",6.95,face=2,just="centre")
    txt(cx-10.85,17.35,"(95% CI)",6.45,muted,just="centre")
    txt(ipa-2.35,13.55,"ΔIPA",6.95,face=2,just="centre")
    txt(ipa-2.35,17.35,"(points)",6.40,muted,just="centre")
    line(c(l,ipa),rep(20.15,2),rule,.14)
    line(rep(sx(0),2),c(tops[1],end),"#A7B3B7",.19)
    z <- dat[dat$layer==ly,]
    stopifnot(identical(as.character(z$endpoint),ep))
    for (i in seq_along(ep)) {
      a <- z[i,]; y <- yy[i]
      if(!a$available) {
        txt(sx(.004),y,paste0("Not estimated (",a$events_T," events)†"),6.35,muted)
        txt(cx,y,"–",6.8,muted,just="right")
        txt(ipa,y,"–",6.8,muted,just="right")
      } else {
        lung <- a$dAUC_lung; others <- a$dAUC_others_after_lung; total <- a$dAUC_all7
        origin_other <- if(lung*others<0) 0 else lung
        stopifnot(all(c(0,lung,origin_other+others,total)>=domain[1]),
                  all(c(0,lung,origin_other+others,total)<=domain[2]))
        # Keep very small increments at their actual width.
        box(min(sx(0),sx(lung)),y-1.12,abs(sx(lung)-sx(0)),2.12,dark[[ly]])
        box(min(sx(origin_other),sx(origin_other+others)),y-1.12,
            abs(sx(origin_other+others)-sx(origin_other)),2.12,light[[ly]])
        net_tip(sx(total),y+1.11)
        txt(cx,y,ci_str(total,a$dAUC_all7_lo,a$dAUC_all7_hi,3),6.55,just="right")
        txt(ipa,y,num(a$gain_all7,2),6.75,just="right")
      }
      if(a$available && a$exploratory) txt(cx+.27,y-.46,"†",6.45,muted)
    }
    line(c(l,r),rep(end+1.0,2),axis_col,.18)
    for (tick in c(-.02,0,.04,.08)) {
      x <- sx(tick)
      line(rep(x,2),c(end+1,end+1.75),axis_col,.16)
      txt(x,end+3.80,if(tick==0) "0" else num(tick,2),6.45,muted,just="centre")
    }
  }
  line(c(3.15,179.8),rep(101.08,2),rule,.15)
  y <- 105.8
  box(3.25,y-.87,2.18,1.74,dark[1]); box(5.43,y-.87,2.18,1.74,dark[2])
  txt(8.65,y,"Lung added first",6.75)
  box(37.4,y-.87,2.18,1.74,light[1]); box(39.58,y-.87,2.18,1.74,light[2])
  txt(42.80,y,"Other six added next",6.75)
  net_tip(81.3,y-.58); txt(83.4,y,"Net ΔAUC",6.75)
  txt(112.0,y,"† Fewer than 100 test-set events during full follow-up.",6.65,muted)
}

main_plot <- function() {
  grid.newpage()
  absolute_panels()
  pushViewport(viewport(y=ym(upper_height),height=mm(111),just=c("centre","top")))
  gain_panels()
  popViewport()
}

main_stem <- file.path(out_dir,"Figure6")
# Cairo PDF uses whole points for its page box. Use the same canvas on all devices.
main_size <- floor(c(W,upper_height+111)*72/25.4)/72
grDevices::cairo_pdf(paste0(main_stem,".pdf"),width=main_size[1],height=main_size[2],family=font,bg="white")
tryCatch(main_plot(),finally=invisible(dev.off()))
svg_path <- paste0(main_stem,".svg")
grDevices::svg(svg_path,width=main_size[1],height=main_size[2],family=font,bg="white")
tryCatch(main_plot(),finally=invisible(dev.off()))
# Cairo's SVG viewBox uses points; label the physical dimensions accordingly.
svg_lines <- readLines(svg_path,encoding="UTF-8")
svg_root <- grep("^<svg ",svg_lines)
stopifnot(length(svg_root)==1L)
svg_lines[svg_root] <- gsub('(width|height)="([0-9.]+)"','\\1="\\2pt"',svg_lines[svg_root])
writeLines(svg_lines,svg_path,useBytes=TRUE)
grDevices::png(paste0(main_stem,".png"),width=main_size[1],height=main_size[2],units="in",res=300,type="cairo",bg="white",antialias="gray")
tryCatch(main_plot(),finally=invisible(dev.off()))
grDevices::tiff(paste0(main_stem,".tiff"),width=main_size[1],height=main_size[2],units="in",res=600,type="cairo",compression="lzw",bg="white",antialias="gray")
tryCatch(main_plot(),finally=invisible(dev.off()))
