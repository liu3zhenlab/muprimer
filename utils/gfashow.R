args <- commandArgs(trailingOnly=T)
gfafile <- args[1] # gfa file, compatible to minigraph output
ref_x_gap_prop <- as.numeric(args[2]) # ratio of ref_seg gap to canvas x_length
nonref_x_gap_prop <- as.numeric(args[3]) # ratio of nonref_seg gap to canvas x_length
max_showseg_prop <- as.numeric(args[4]) # maximum length of seg to show in the actual length
main_label <- args[5] # main text
gtffile <- args[6] # gtf file of the gene
#if (is.null(gtffile)) { print("NULL - YES\n") }
outpdf <- args[7] # PDF output file

##############################################################################
# external parameters
##############################################################################
#ref_x_gap_prop <- 0.02
#max_showseg_prop <- 0.5
#nonref_x_gap_prop <- 0.05
#main_label <- "Graph of haplotypes of xxx"
#gfafile <- "/homes/liu3zhen/scripts2/homotools/homomg_dev/1_minigraph/all.gfa"
##############################################################################
# modules or packages
##############################################################################
if (! "RColorBrewer" %in% rownames(installed.packages())) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

##############################################################################
# parameter
##############################################################################
ref_y_bottom <- 0.1
ref_y_top <- 0.5
if (gtffile != "NULL") {
  ref_y_bottom <- 0.4
  ref_y_top <- 0.8
}
margin_prop <- 0.01

##############################################################################
# initiation
##############################################################################
seg_name <- NULL
seg_xstart <- NULL
seg_xend <- NULL
seg_ypos <- NULL

##############################################################################
# color generator
##############################################################################
col_generator <- function(ncol) {
  library(RColorBrewer)
  #set.seed(100)
  preset_cols <- c("darkseagreen4", "lightsalmon3", "mediumpurple3")
  if (ncol <= 3) {
    cols <- preset_cols[1:ncol]
  } else {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cols
  }
  cols
}

##############################################################################
# segments
##############################################################################
segplot <- function(gfa_data, xbase_pos=0, seg_name="segment",
                    pos_name="so", len_name="len", xmax=NULL,
                    ratio2max=NULL, seg_gap=0, y_bottom, y_top, ...) {
  # initiation
  out_seg_name <- NULL
  out_seg_xstart <- NULL
  out_seg_xend <- NULL
  out_seg_ypos <- NULL
  # plot segments
  for (eref in 1:nrow(gfa_data)) {
    gfa_data_seg_name <- gfa_data[eref, seg_name]
    gfa_data_seg_pos <- gfa_data[eref, pos_name]
    gfa_data_seg_len <- gfa_data[eref, len_name]
    gfa_data_seg_xleft <- xbase_pos + gfa_data_seg_pos + seg_gap * (eref - 1)
    gfa_data_seg_xright <-gfa_data_seg_xleft + gfa_data_seg_len
    
    # plot
    rect(gfa_data_seg_xleft, y_bottom, gfa_data_seg_xright, y_top, ...)
    
    # text label
    text_label <- gfa_data_seg_name
    if (!is.null(ratio2max)) {
      ratio <- gfa_data[eref, ratio2max]
      if (ratio > 1) {
        ratio <- round(ratio, 2)
        text_label <- paste0(gfa_data_seg_name, " (", ratio, "x)")
      }
    }
    
    # add labels
    if (! is.null(xmax) & gfa_data_seg_len > xmax*0.03) {
      #text(x=(gfa_data_seg_xleft + gfa_data_seg_xright)/2,
      #    y= y_bottom, pos=3, labels=text_label, col="red", cex=0.8)
      text(x=gfa_data_seg_xleft, y= (y_bottom + y_top) / 2, offset=0.2,
              pos=4, labels=text_label, col="orangered4", cex=0.8)
    }
    # seq start and end positions
    out_seg_name <- c(out_seg_name, gfa_data$segment[eref])
    out_seg_xstart <- c(out_seg_xstart, gfa_data_seg_xleft)
    out_seg_xend <- c(out_seg_xend, gfa_data_seg_xright)
    out_seg_ypos <- c(out_seg_ypos, (y_bottom + y_top) / 2)
  }
  # outupt
  refstart <- gfa_data[, pos_name]
  seglen <- gfa_data[, len_name]
  out <- list(name=out_seg_name, refstart=refstart, seglen=seglen,
              xstart=out_seg_xstart, xend=out_seg_xend, ypos=out_seg_ypos)
  invisible(out)
}

##############################################################################
# line transformation
##############################################################################
transform_curve <- function(startp, endp, midp=NULL, npoint=1000) {
  ### computer transform_curve coordinates
  if (is.null(midp)) {
    midp <- (startp + endp) / 2 
  }
  
  beizer_value <- sqrt(1:npoint) / sqrt(npoint)  # sqrt as default
  
  if (startp[1] == endp[1] & startp[1] == midp[1]) {
    curve.x <- rep(startp[1], 2*npoint)
  } else {
    curve.x1 <- seq(startp[1], midp[1], by = (midp[1] - startp[1])/(npoint-1))
    curve.x2 <- seq(midp[1], endp[1], by = (endp[1] - midp[1])/(npoint-1))
    curve.x <- c(curve.x1, curve.x2)
  }
  
  if (startp[2] == endp[2] & startp[2] == midp[2]) {
    curve.y <- rep(startp[2], 2*npoint)
  } else {
    curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
    curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
    curve.y <- c(curve.y1, curve.y2)
  }
  list(x=curve.x, y=curve.y)
}

##############################################################################
# links
##############################################################################
linker <- function(gfa_data, ...) {
  gfa_l <- gfa_data[gfa_data[,1] == "L", 2:5]
  colnames(gfa_l) <- c("seg1", "strand1", "seg2", "strand2")
  for (lrow in 1:nrow(gfa_l)) {
    seg1name <- gfa_l[lrow, "seg1"]
    seg1strand <- gfa_l[lrow, "strand1"]
    seg2name <- gfa_l[lrow, "seg2"]
    seg2strand <- gfa_l[lrow, "strand2"]
    
    seg1xstart <- seg_xstart[seg1name]
    seg1xend <- seg_xend[seg1name]
    seg1ypos <- seg_ypos[seg1name]
    
    seg2xstart <- seg_xstart[seg2name]
    seg2xend <- seg_xend[seg2name]
    seg2ypos <- seg_ypos[seg2name]
    
    # set y-axis adjustments for linker positions
    if (seg1ypos > seg2ypos) {
      y_adjust <- (ref_y_top - ref_y_bottom) / 2   
    } else if (seg1ypos < seg2ypos) {
      y_adjust <- (ref_y_bottom - ref_y_top) / 2
    } else {
      y_adjust <- 0
    }
    
    # plot lines
    startp <- c(seg1xend, seg1ypos - y_adjust)
    endp <- c(seg2xstart, seg2ypos + y_adjust)
    midp <- as.numeric((startp + endp)/2)
    if ((y_adjust==0) & (abs(seg2xstart - seg1xend - 1) > ref_x_gap)) {
      midp[2] <- midp[2] + (ref_y_top - ref_y_bottom) * 0.8
    }
    
    # plot linkers
    curve_out <- transform_curve(startp, endp, midp)
    lines(curve_out$x, curve_out$y, ...) 
  }
}

###########################################################
#' module to determine xaxis
###########################################################
smartaxis <- function(maxnum) {
# determine reasonable axis values
  numdigits <- nchar(maxnum)
  unit <- 10 ^ (numdigits - 1) / (2- round((maxnum / 10 ^ numdigits), 0)) # 1 or 5 e (numdigits - 1)
  subunit <- unit / 5 
  
  numsat <- unit * (0:10)
  numsat <- numsat[numsat < maxnum]
  
  if (numdigits >= 7) {
    numlabels <- numsat / 1000000
    label.scale <- "Mb"
  } else if (numdigits < 7 & numdigits >= 4) {
    numlabels <- numsat / 1000
    label.scale <- "kb"
  } else {
    numlabels <- numsat
    label.scale <- "bp"
  }
  
  subunits <- seq(0, maxnum, by = subunit)
  subunits <- subunits[!subunits %in% c(numsat, 0)] 
  # return
  list(numsat, numlabels, label.scale, subunits)
}

###########################################################
#' lifter
###########################################################
lifter <- function(coords, orimap, newmap, base) {
  # lift coords from orimap to the newmap
  if (sum(coords > orimap)>0) {
    last_start_order <- which.max(orimap[coords > orimap]) # last segment with the start less than coords
    last_start_coords <- max(orimap[coords > orimap]) # last start less than coords
    new_coords <- newmap[last_start_order] + coords - last_start_coords + 1
  } else {
    new_coords <- coords + base
  }
  new_coords
}

###########################################################
#' module: plotting gene structure
###########################################################
genegtf <- function(gtf, feature=c("exon", "CDS"), exoncol="lightskyblue3",
                    cdscol="lightskyblue4", xbase=0, ybottom=0.1, ytop=0.3) {
  # plot gene structure
  # GRMZM2G171650	ensembl	exon	2001	2423	.	+	.	gene_id "GRMZM2G171650"; gene_version "2";...
  stopifnot(nrow(gtf)>0)
  plot_gtf <- gtf[gtf[, 3] %in% c("exon", "CDS"), ]
  stopifnot(nrow(plot_gtf)>0)
  plot_gtf <- plot_gtf[order(plot_gtf[, 3], decreasing=T), ]
  plot_gtf$col <- exoncol
  plot_gtf$col[plot_gtf[,3] == "CDS"] <- cdscol

  # gene
  lines(c(min(plot_gtf[, c(4,5)]), max(plot_gtf[, c(4,5)])) + xbase, rep((ybottom+ytop)/2, 2),
        lwd=2, col="gray60")
  
  # exon and CDS
  for (i in 1:nrow(plot_gtf)) {
    rect(xleft=plot_gtf[i, 4] + xbase, xright=plot_gtf[i, 5] + xbase,
         ybottom=ybottom, ytop=ytop, col=plot_gtf[i, "col"], border=NA)
  }
}

##############################################################################
# input data
##############################################################################
gfa <- read.delim(gfafile, header=F, stringsAsFactors=F)
gfa_s <- gfa[gfa[,1] == "S", -c(1, 3)]
colnames(gfa_s) <- c("segment", "len", "sn", "so", "sr")
gfa_s$len <- as.numeric(gsub("LN\\:i\\:", "", gfa_s$len))
gfa_s$sn <- gsub("SN\\:Z\\:", "", gfa_s$sn)
gfa_s$so <- as.numeric(gsub("SO\\:i\\:", "", gfa_s$so))
gfa_s$sr <- as.numeric(gsub("SR:\\i\\:", "", gfa_s$sr))

slen <- gfa_s$len
sname <- gfa_s$sn
soffset <- gfa_s$so
srank <- gfa_s$sr

##############################################################################
### reference:
##############################################################################
ref_s <- gfa_s[gfa_s$sr==0, ]
total_reflen <- sum(slen[srank==0])  # total length of reference segments
ref_x_gap <- total_reflen * ref_x_gap_prop
x0_base <- total_reflen * margin_prop - 1
x_max <- x0_base + total_reflen * (1 + margin_prop) + ref_x_gap * (nrow(ref_s) - 1)

##############################################################################
# non-reference
##############################################################################
nlayer <- 1
max_showseg <- max_showseg_prop * total_reflen

if (sum(srank>0)>0) {
  nonref_s <- gfa_s[srank>0, ]
  nonref_s <- nonref_s[order(nonref_s$so), ]
  nonref_s$adj_len <- nonref_s$len
  nonref_s$adj_len[nonref_s$len > max_showseg] <- max_showseg
  nonref_s$ratio2max <- nonref_s$len / max_showseg
}

# layers
nonref_x_gap <- nonref_x_gap_prop * total_reflen
nonref_s$is_processed <- 0
nonref_s$layer <- 0
while (sum(nonref_s$is_processed == 0) > 0) {
  nlayer <- nlayer + 1
  cur_pos <- x0_base
  for (prow in 1:nrow(nonref_s)) {
    if (nonref_s[prow, "is_processed"]==0) {
      seq_pos <- nonref_s$so[prow]
      plot_len <- nonref_s$adj_len[prow]
      if (seq_pos > cur_pos - x0_base & (cur_pos + plot_len) <= x_max ) {
      # not overlap with previous one and not out of range
        nonref_s[prow, "is_processed"] <- 1
        nonref_s$layer <- nlayer
      }
      cur_pos <- cur_pos + plot_len + nonref_x_gap
    }
  }
}

##############################################################################
# plot
##############################################################################
pdf(outpdf, width=6, height=3.5 + nlayer*0.6  - 2)  # height adjusted by # layers

par(mgp=c(2,0.5,0), mar=c(3, 0, 2, 0))

# plot canvas
plot(NULL, NULL, xlim=c(0, x_max), ylim=c(0, nlayer-0.4), xlab="", ylab="",
     main=main_label, axes=F)

##############################################################################
# plot reference segments
##############################################################################
ref_seg_coordi <- segplot(gfa_data=ref_s, xbase_pos=x0_base, pos_name="so", len_name="len",
        seg_gap=ref_x_gap, y_bottom=ref_y_bottom, y_top=ref_y_top, col="azure4",
        border="azure4", xmax=x_max)

seg_name <- c(seg_name, ref_seg_coordi$name)
seg_xstart <- c(seg_xstart, ref_seg_coordi$xstart)
seg_xend <- c(seg_xend, ref_seg_coordi$xend)
seg_ypos <- c(seg_ypos, ref_seg_coordi$ypos)
names(seg_xstart) <- seg_name
names(seg_xend) <- seg_name
names(seg_ypos) <- seg_name

##############################################################################
# plot nonreference segments
##############################################################################
plot_nonref_sranks <- nonref_s$sr # non-redundant seg ranks
nranks <- length(unique(plot_nonref_sranks))
nonref_cols <- col_generator(nranks)
for (layer in 2:nlayer) {
  for (rank in 1:nranks) {
    plot_nonref_s <- nonref_s[nonref_s$layer==layer &
                              nonref_s$sr==rank, ]
    plot_nonref_s$adjusted_so <- sapply(plot_nonref_s$so, lifter, orimap=ref_seg_coordi$refstart,
                                        newmap=ref_seg_coordi$xstart, base=x0_base)
    nonref_seg_cooordi <- segplot(gfa_data=plot_nonref_s, xbase_pos=x0_base, pos_name="adjusted_so",
                                  len_name="adj_len", seg_gap=nonref_x_gap,
                                  y_bottom=layer+ref_y_bottom-1, y_top=layer+ref_y_top-1,
                                  col=nonref_cols[rank], border=nonref_cols[rank], xmax=x_max)
    seg_name <- c(seg_name, nonref_seg_cooordi$name)
    seg_xstart <- c(seg_xstart, nonref_seg_cooordi$xstart)
    seg_xend <- c(seg_xend, nonref_seg_cooordi$xend)
    seg_ypos <- c(seg_ypos, nonref_seg_cooordi$ypos)
    names(seg_xstart) <- seg_name
    names(seg_xend) <- seg_name
    names(seg_ypos) <- seg_name
  }
}

##############################################################################
# plot links
##############################################################################
linker(gfa, col="cornsilk3", lwd=3)

##############################################################################
# x-axis
##############################################################################
smart_axis_data <- smartaxis(total_reflen)
xtick_coords <- sapply(smart_axis_data[[1]], lifter, orimap=ref_seg_coordi$refstart,
       newmap=ref_seg_coordi$xstart, base=x0_base)
xtick_labels <- smart_axis_data[[2]]
xtick_labels_unit <- smart_axis_data[[3]]
xlabel <- paste("coordinates (", xtick_labels_unit, ")")
axis(side=1, at=xtick_coords, labels=xtick_labels)
mtext(text=xlabel, side=1, padj=3)

##############################################################################
# gtf plotting
##############################################################################
if (gtffile != "NULL") {
  # read GTF
  gtf <- read.delim(gtffile, stringsAsFactors=F, header=F)
  # gtf lifter
  gtf[,4] <- sapply(gtf[,4], lifter, orimap=ref_seg_coordi$refstart,
                    newmap=ref_seg_coordi$xstart, base=x0_base)
  gtf[,5] <- sapply(gtf[,5], lifter, orimap=ref_seg_coordi$refstart,
                    newmap=ref_seg_coordi$xstart, base=x0_base)
  genegtf(gtf=gtf, feature=c("exon", "CDS"), ybottom=0, ytop=0.25)
}
##############################################################################
# close pdf
##############################################################################
dev.off()

