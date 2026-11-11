library(grid)
library(svglite)
library(ragg)
library(jsonlite)

root <- normalizePath('.')
out <- file.path(root, '8-Visualization/refinements/fig1_global_revision_20261103_v3/compact_panels')
src <- file.path(root, '8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data')
pal <- read.csv(file.path(root, 'governance/data/subtype_label_color_dictionary.csv'), stringsAsFactors=FALSE)
study_pal <- read.csv(file.path(root, '0-Data/qc/study_color_dictionary.csv'), stringsAsFactors=FALSE)
colors <- setNames(pal$Color, pal$FinalSubtype)
states <- c('C1','C2','C3','C1_C2','C1_C3','C2_C3')
studies <- study_pal$Study[study_pal$AnalysisRole=='Discovery RNA-seq']
short <- setNames(c('Ahn','Dong','Deng','Jeon','Lee','Fan','TCGA'), studies)
u <- read.csv(file.path(src, 'Figure1C.csv'), stringsAsFactors=FALSE)
d <- read.csv(file.path(src, 'Figure1D.csv'), stringsAsFactors=FALSE)
six <- read.csv(file.path(src, 'SFig1K5.csv'), stringsAsFactors=FALSE)
k3 <- read.csv(file.path(src, 'SFig1K3.csv'), stringsAsFactors=FALSE)
k4 <- read.csv(file.path(src, 'SFig1K4.csv'), stringsAsFactors=FALSE)
stopifnot(nrow(u)==497, all(u$SubtypeClass=='Main'), !anyDuplicated(u$Sample),
          nrow(d)==21, sum(d$N)==497, nrow(six)==586, !anyDuplicated(six$Sample),
          identical(six[c('Sample','Study','FinalSubtype')],k3[c('Sample','Study','FinalSubtype')]),
          identical(six[c('Sample','Study','FinalSubtype')],k4[c('Sample','Study','FinalSubtype')]),
          all(six$Study %in% studies), all(six$FinalSubtype %in% states))
stopifnot(identical(as.integer(table(factor(u$FinalSubtype,levels=states[1:3]))), c(112L,223L,162L)))
counts <- as.data.frame(table(Study=factor(six$Study,levels=studies), FinalSubtype=factor(six$FinalSubtype,levels=states)), stringsAsFactors=FALSE)
names(counts)[3] <- 'N'
counts$Study <- as.character(counts$Study)
counts$FinalSubtype <- as.character(counts$FinalSubtype)
counts$Total <- as.integer(table(factor(six$Study,levels=studies))[counts$Study])
counts$Fraction <- counts$N / counts$Total
stopifnot(identical(as.integer(tapply(counts$N,factor(counts$FinalSubtype,levels=states),sum)),c(112L,223L,162L,41L,39L,9L)))
for (st in studies) {
  z <- d[d$Study==st,]
  stopifnot(length(unique(z$Total))==1L, sum(z$N)==z$Total[1], abs(sum(z$Fraction)-1)<1e-12)
  for (cl in states[1:3]) stopifnot(z$N[z$FinalSubtype==cl]==sum(u$Study==st & u$FinalSubtype==cl))
}
write.csv(u,file.path(out,'source_data/Figure1C_main_umap.csv'),row.names=FALSE)
write.csv(d,file.path(out,'source_data/Figure1D_study_composition.csv'),row.names=FALSE)
write.csv(counts,file.path(out,'source_data/SFig1H_sixstate_composition.csv'),row.names=FALSE)
write.csv(six[c('Sample','Study','FinalSubtype','SubtypeClass')],file.path(out,'source_data/SFig1H_sample_membership.csv'),row.names=FALSE)
write.csv(six,file.path(out,'source_data/SFig1G_sixstate_umap.csv'),row.names=FALSE)

# All visual marks and all preview/export devices use R. Text sizes are in points.
txt <- function(label,x,y,size=6.5,col='#222222',fontface='plain',just='centre',rot=0,lineheight=1.2) {
  grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,rot=rot,
            gp=gpar(fontfamily='Arial',fontsize=size,col=col,fontface=fontface,lineheight=lineheight))
}
line <- function(x0,y0,x1,y1,col='#555555',lwd=.45) {
  grid.segments(unit(x0,'mm'),unit(y0,'mm'),unit(x1,'mm'),unit(y1,'mm'),gp=gpar(col=col,lwd=lwd))
}
rect <- function(x,y,w,h,fill,col=NA,lwd=.4) {
  grid.rect(unit(x,'mm'),unit(y,'mm'),unit(w,'mm'),unit(h,'mm'),gp=gpar(fill=fill,col=col,lwd=lwd))
}
legend3 <- function(y,positions) {
  for (i in seq_along(positions)) {
    rect(positions[i]-3.8,y,1.8,1.8,colors[states[i]])
    txt(states[i],positions[i]-.8,y,6.5,just='left')
  }
}
draw_umap <- function() {
  grid.newpage()
  txt('Principal states (n = 497)',27,43.45,7.2,fontface='bold')
  lx <- c(15.2,29.3,43.4)
  for (i in 1:3) {
    rect(lx[i]-6.2,40.1,1.6,1.6,colors[states[i]])
    txt(paste0(states[i],' (',c(112,223,162)[i],')'),lx[i]-4.9,40.1,6.5,just='left')
  }
  xr <- c(-4.1,4.7); yr <- c(-3.2,3.4)
  scale <- 4.45
  x0 <- 11.8; y0 <- 8.3
  xp <- function(x) x0+(x-xr[1])*scale
  yp <- function(y) y0+(y-yr[1])*scale
  grid.points(unit(xp(u$UMAP1),'mm'),unit(yp(u$UMAP2),'mm'),pch=16,size=unit(.63,'mm'),
              gp=gpar(col=colors[u$FinalSubtype],alpha=.82))
  line(xp(xr[1]),yp(yr[1]),xp(xr[2]),yp(yr[1]))
  line(xp(xr[1]),yp(yr[1]),xp(xr[1]),yp(yr[2]))
  for (x in c(-4,-2,0,2,4)) {line(xp(x),y0,xp(x),y0-.9);txt(x,xp(x),y0-2.25)}
  for (y in c(-2,0,2)) {line(x0,yp(y),x0-.9,yp(y));txt(y,x0-1.8,yp(y),just='right')}
  txt('UMAP 1',(xp(xr[1])+xp(xr[2]))/2,2.2,6.5)
  txt('UMAP 2',3.4,(yp(yr[1])+yp(yr[2]))/2,6.5,rot=90)
}
draw_composition <- function() {
  grid.newpage()
  txt('Principal-state composition',27,38.3,7.2,fontface='bold')
  legend3(34.55,c(17,29,41))
  x0 <- 12.0; width <- 31.4
  for (j in seq_along(studies)) {
    z <- d[d$Study==studies[j],];z <- z[match(states[1:3],z$FinalSubtype),]
    y <- 29.9-(j-1)*3.4
    txt(short[studies[j]],.7,y,6.5,just='left')
    pos <- c(0,cumsum(z$Fraction))
    for (i in 1:3) {
      rect(x0+width*(pos[i]+pos[i+1])/2,y,width*z$Fraction[i],2.65,colors[states[i]],'#FFFFFF',.35)
      txt(z$N[i],x0+width*(pos[i]+pos[i+1])/2,y,6.5,col=if(i==1) '#FFFFFF' else '#111111')
    }
    txt(paste0('n=',z$Total[1]),53.5,y,6.5,just='right')
  }
  line(x0,7.1,x0+width,7.1,col='#999999')
  for (f in c(0,.5,1)) {line(x0+width*f,7.1,x0+width*f,6.4,col='#999999');txt(paste0(100*f,'%'),x0+width*f,4.9)}
  txt('Within-study proportion',x0+width/2,1.9,6.5)
}
draw_sixstate <- function() {
  grid.newpage()
  txt('All six states (n = 586)',30,38.35,7.2,fontface='bold')
  left <- 19.0; cw <- 6.6
  headers <- c('C1','C2','C3','C1/\nC2','C1/\nC3','C2/\nC3')
  txt('Study',1.4,34.3,6.5,fontface='bold',just='left')
  txt('n',17.0,34.3,6.5,fontface='bold',just='right')
  for (i in 1:6) {
    rect(left+(i-.5)*cw,34.3,cw-.2,4.7,colors[states[i]])
    txt(headers[i],left+(i-.5)*cw,34.3,6.5,col=if(i %in% c(1,4,5)) '#FFFFFF' else '#111111',lineheight=.92)
  }
  for (j in 1:7) {
    y <- 29.95-(j-1)*3.5
    z <- counts[counts$Study==studies[j],];z <- z[match(states,z$FinalSubtype),]
    txt(short[studies[j]],1.4,y,6.5,just='left')
    txt(z$Total[1],17.0,y,6.5,just='right')
    for (i in 1:6) {
      rect(left+(i-.5)*cw,y,cw-.2,3.32,adjustcolor(colors[states[i]],alpha.f=.12))
      txt(z$N[i],left+(i-.5)*cw,y,6.5)
    }
  }
  line(1.4,6.8,58.6,6.8,col='#AAAAAA')
  txt('Total',1.4,4.65,6.5,fontface='bold',just='left')
  txt('586',17,4.65,6.5,fontface='bold',just='right')
  totals <- c(112,223,162,41,39,9)
  for (i in 1:6) txt(totals[i],left+(i-.5)*cw,4.65,6.5,fontface='bold')
}
draw_sixstate_umap <- function() {
  grid.newpage()
  txt('Six molecular states (n = 586)',25,46.0,7.2,fontface='bold')
  for (i in 1:6) {
    j <- ((i-1)%%3)+1
    x <- c(10,24,38)[j]
    y <- if (i<=3) 42.3 else 38.9
    grid.points(unit(x-3,'mm'),unit(y,'mm'),pch=if(i<=3)16 else 1,size=unit(if(i<=3)1.25 else 1.35,'mm'),gp=gpar(col=colors[states[i]],lwd=.7))
    txt(gsub('_','/',states[i]),x-1.3,y,6.5,just='left')
  }
  xr <- c(-4.1,4.7); yr <- c(-3.2,3.4)
  scale <- 4.03
  x0 <- 11.0; y0 <- 8.0
  xp <- function(x) x0+(x-xr[1])*scale
  yp <- function(y) y0+(y-yr[1])*scale
  for (i in 1:6) {
    z <- six[six$FinalSubtype==states[i],]
    grid.points(unit(xp(z$UMAP1),'mm'),unit(yp(z$UMAP2),'mm'),pch=if(i<=3)16 else 1,
                size=unit(if(i<=3).6 else .85,'mm'),gp=gpar(col=colors[states[i]],alpha=if(i<=3).78 else 1,lwd=.65))
  }
  line(xp(xr[1]),yp(yr[1]),xp(xr[2]),yp(yr[1]))
  line(xp(xr[1]),yp(yr[1]),xp(xr[1]),yp(yr[2]))
  for (x in c(-4,-2,0,2,4)) {line(xp(x),y0,xp(x),y0-.9);txt(x,xp(x),y0-2.2)}
  for (y in c(-2,0,2)) {line(x0,yp(y),x0-.9,yp(y));txt(y,x0-1.8,yp(y),just='right')}
  txt('UMAP 1',(xp(xr[1])+xp(xr[2]))/2,2.0,6.5)
  txt('UMAP 2',3.1,(yp(yr[1])+yp(yr[2]))/2,6.5,rot=90)
}
export <- function(name,draw,w,h) {
  selected <- commandArgs(trailingOnly=TRUE)
  if(length(selected) && !name %in% selected) return(invisible(NULL))
  base <- file.path(out,'figures',name)
  svglite(paste0(base,'.svg'),width=w/25.4,height=h/25.4,bg='white');draw();dev.off()
  cairo_pdf(paste0(base,'.pdf'),width=w/25.4,height=h/25.4,family='Arial',bg='white');draw();dev.off()
  agg_png(paste0(base,'.png'),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
}
export('Figure1C_main_umap',draw_umap,54,45.83)
export('Figure1D_study_composition',draw_composition,54,40.12)
export('SFig1H_sixstate_composition',draw_sixstate,60,40)
export('SFig1G_sixstate_umap',draw_sixstate_umap,50,48)
manifest <- data.frame(panel=c('Figure1C_main_umap','Figure1D_study_composition','SFig1H_sixstate_composition','SFig1G_sixstate_umap'),
                       width_mm=c(54,54,60,50),height_mm=c(45.83,40.12,40,48),n=c(497,497,586,586),
                       minimum_font_pt=6.5,maximum_font_pt=7.2,png_dpi=600)
write.csv(manifest,file.path(out,'review/export_manifest.csv'),row.names=FALSE)
write_json(list(status='NUMERIC_CHECKS_PASS',source_rows=c(497,21,586),
                main_state_counts=c(C1=112,C2=223,C3=162),six_state_counts=setNames(c(112,223,162,41,39,9),states),
                umap_coordinate_changes=FALSE,umap_equal_axis_units=TRUE,
                sixstate_count_zero_cells=sum(counts$N==0),all_42_cells_displayed=TRUE,
                study_palette_loaded=TRUE,study_palette_used=FALSE,
                study_palette_reason='Studies are directly labelled; fills encode molecular states.'),
           file.path(out,'review/numerical_validation.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
