# Full-row synthesis: native editable R/grid marks; frozen scientific statistics.
F4SUMMARY <- '8-Visualization/refinements/fig4_summary_20261030_v2/'
f4summary_read <- function(name)read_source(paste0(F4SUMMARY,'source_data/',name,'.csv'))
prepare_fig4_summary_v2 <- function() {
  # Fig2 anchors define coverage; discordance/counterexamples are retained explicitly.
  plan<-data.table(Subtype=c(rep('C1',5),rep('C2',4),rep('C3',5)),
    Slot=c(1:5,1:4,1:5),
    Pathway=c('GOBP_AMINO_ACID_CATABOLIC_PROCESS','HALLMARK_FATTY_ACID_METABOLISM',
      'GOBP_XENOBIOTIC_METABOLIC_PROCESS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','REACTOME_TRANSLATION',
      'GOBP_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY','GOBP_CILIUM_ORGANIZATION',
      'REACTOME_TRANSLATION','REACTOME_RESPIRATORY_ELECTRON_TRANSPORT',
      'GOBP_GLUCOSE_CATABOLIC_PROCESS','HALLMARK_HYPOXIA',
      'GOBP_EPIDERMAL_CELL_DIFFERENTIATION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_E2F_TARGETS'),
    Label=c('Amino-acid catabolism','Fatty-acid metabolism','Xenobiotic metabolism',
      'Oxidative phosphorylation','Translation','WNT-related program','Cilium organization',
      'Translation','Respiratory chain\u2020','Glucose catabolism','Hypoxia','Epidermal differentiation',
      'EMT gene set','E2F targets'),
    SelectionRole=c(rep('Established metabolic anchor',4),'RNA/protein discordance',
      'Concordant C2 morphogenesis anchor','Concordant C2 morphogenesis anchor',
      'RNA/protein discordance','Full-RNA support; matched result not significant',
      'Established glycolytic anchor','Concordant hypoxia anchor','Concordant epithelial anchor',
      'Concordant remodeling anchor; not equivalent to epithelial differentiation',
      'Counterexample to general cross-omic cell-cycle activation'))
  inputs<-c(
    '4-Multiomic/refinements/fig4_matched_four_gsea_20261030/paired_gsea.csv',
    '4-Multiomic/publication/tables/step4_driver_snv_indel_study_stratified_tests.csv',
    '4-Multiomic/refinements/fig4_cnv_pairwise_20261030/pairwise_all.csv',
    '8-Visualization/refinements/fig2_sfig2_20261030_v14/tables/G_display_plan.csv',
    '8-Visualization/refinements/fig2_sfig2_20261030_v14/tables/G_evidence_audit.csv',
    '8-Visualization/scripts/render_fig4_summary_v2.R')
  evidence<-read_source(inputs[5]);inputs<-unique(c(inputs,evidence$SourceTable))
  paired<-read_source(inputs[1]);selected<-merge(unique(plan[,.(Pathway)]),paired,by.x='Pathway',by.y='pathway')
  stopifnot(nrow(selected)==39L,uniqueN(plan$Pathway)==13L)
  allcells<-rbindlist(lapply(c('RNA','Protein'),function(view) {
    z<-copy(selected);z[,`:=`(View=view,NES=get(paste0('NES_',view)),FDR=get(paste0('padj_',view)))]
    z[,.(Pathway,Subtype,View,NES,FDR,MeasuredSize,Collection)]
  }))
  allcells[,Significant:=FDR<.05]
  allcells[,Glyph:=ifelse(!Significant,'\u00b7',ifelse(NES>0,'\u2191','\u2193'))]
  display<-merge(plan,allcells,by=c('Subtype','Pathway'));setorder(display,Subtype,Slot,View)
  stopifnot(nrow(display)==28L)
  mutation<-read_source(inputs[2]);enriched<-mutation[CMHFDR<.05&CMHOddsRatio>1]
  stopifnot(nrow(enriched[Subtype=='C1'])==0L,
    setequal(enriched[Subtype=='C2',Gene],c('BAP1','IDH1','IDH2')),
    setequal(enriched[Subtype=='C3',Gene],c('KRAS','SMAD4','TP53')))
  cnv<-read_source(inputs[3]);cnv<-cnv[Direction=='Loss'&Contrast=='C3 vs C2'&FDR<.05]
  stopifnot(nrow(cnv)==8L,all(grepl('^1p36',cnv$Cytoband)),all(cnv$DifferencePP<0))
  prior<-read_source(inputs[4]);context<-prior[Layer!='Program']
  # Compact text is tied to the audited Fig2 evidence, including shared TF signals.
  cards<-data.table(Subtype=c('C1','C2','C3'),
    Title=c('Metabolic specialization','Morphogenesis / translation','Hypoxia / epithelial remodeling'),
    DNA=c('No prespecified driver enriched\nat CMH q < 0.05',
      'BAP1, IDH1/2 enriched\n1p36 loss > C3',
      'KRAS, SMAD4, TP53 enriched\n1p36 loss < C2'),
    TF=c('HNF4A / HNF1A / CEBPA / PPARA','HNF1B / ONECUT1','KLF5 / GRHL2 / HIF1A'),
    TME=c('Higher immune score\nMyeloid signals: method-dependent',
      'Lower immune / stromal scores',
      'Higher fibroblast / CAF signals\nMCP-counter / EPIC; xCell unclear'))
  respiration<-evidence[LabelID=='C2_P2'&Statistic=='NES'&grepl('RESPIRATORY_ELECTRON_TRANSPORT',Key)]
  stopifnot(nrow(respiration)==1L,respiration$Value>0,respiration$FDR<.05)
  cycle<-evidence[LabelID=='C3_P1'&Statistic=='NES'&grepl('MITOTIC_NUCLEAR_DIVISION',Key)]
  stopifnot(nrow(cycle)==1L,cycle$Value>0,cycle$FDR<.05)
  write<-function(x,name)fwrite(x,file.path(STEP8_ROOT,F4SUMMARY,'source_data',paste0(name,'.csv')),na='NA')
  write(plan,'display_plan');write(selected,'all_selected_pathways');write(allcells,'all_contrasts')
  write(display,'displayed_cells');write(cards,'subtype_cards');write(context,'fig2_context')
  write(evidence,'fig2_evidence_audit');write(respiration,'full_rna_respiration')
  write(cycle,'full_rna_cell_cycle')
  write(mutation,'mutation_tests');write(enriched,'mutation_enrichments');write(cnv,'cnv_context')
  write(data.table(Path=inputs,MD5=unname(tools::md5sum(file.path(STEP8_ROOT,inputs)))),'input_hashes')
}
f4summary_icon <- function(subtype,x,y,col) {
  # Conceptual icons carry no measurements, cell identity or mechanism claim.
  ln<-function(x,y)grid::grid.lines(grid::unit(x,'mm'),grid::unit(y,'mm'),gp=grid::gpar(col=col,lwd=.8))
  if(subtype=='C1') {
    a<-seq(0,2*pi,length.out=7);ln(x+2.4*cos(a),y+2.4*sin(a))
    for(k in 1:6)ln(x+c(2.4,3.5)*cos(a[k]),y+c(2.4,3.5)*sin(a[k]))
  } else if(subtype=='C2') {
    grid::grid.roundrect(x=grid::unit(x,'mm'),y=grid::unit(y,'mm'),width=grid::unit(8,'mm'),
      height=grid::unit(4.7,'mm'),r=grid::unit(2.2,'mm'),gp=grid::gpar(fill=NA,col=col,lwd=.8))
    ln(x+seq(-2.8,2.8,length.out=11),y+rep(c(-1.3,1.3),length.out=11))
  } else {
    for(k in -1:1) {
      fig4_v2_rect(x+k*2.7,y,2.5,4.3,'white',col=col,lwd=.8)
      grid::grid.circle(x=grid::unit(x+k*2.7,'mm'),y=grid::unit(y,'mm'),r=grid::unit(.65,'mm'),
        gp=grid::gpar(fill=NA,col=col,lwd=.65))
    }
    ln(x+c(-4.2,4.2),y+c(-2.5,-2.5))
  }
}
build_fig4_summary_v2 <- function() {
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  cards<-f4summary_read('subtype_cards');plan<-f4summary_read('display_plan');cells<-f4summary_read('displayed_cells')
  function() {
    grid::grid.newpage()
    for(i in 1:3) {
      st<-names(pal)[i];left<-4+(i-1)*56;center<-left+27;col<-pal[i];card<-cards[Subtype==st]
      # Three balanced columns, 174 mm overall: same full-row width as Figure4A.
      fig4_v2_rect(center,59,54,103,'white',col=grDevices::adjustcolor(col,alpha.f=.55),lwd=.65)
      fig4_v2_rect(center,100,54,21,grDevices::adjustcolor(col,alpha.f=.07))
      fig4_v2_text(st,left+3,105,11,TRUE,col=col,just='left')
      f4summary_icon(st,left+46,104,col)
      fig4_v2_text(card$Title,left+3,96.5,6.6,TRUE,just='left')
      fig4_v2_text('Selected program',left+3,87.3,6.1,TRUE,just='left')
      fig4_v2_text('RNA',left+42,87.3,5.9,TRUE)
      fig4_v2_text('Protein',left+49,87.3,5.9,TRUE)
      p<-plan[Subtype==st][order(Slot)]
      for(j in seq_len(nrow(p))) {
        yy<-80-(p$Slot[j]-1)*6
        if(j%%2==1)fig4_v2_rect(center,yy,50,5.6,'#F6F7F8')
        fig4_v2_text(p$Label[j],left+3,yy,6.1,just='left')
        for(v in c('RNA','Protein')) {
          z<-cells[Subtype==st&Pathway==p$Pathway[j]&View==v]
          stopifnot(nrow(z)==1L)
          cc<-if(!z$Significant)'#7D858B' else if(z$NES>0)'#AA573C' else '#33749D'
          fig4_v2_text(z$Glyph,left+if(v=='RNA')42 else 49,yy,if(z$Significant)9 else 10,TRUE,col=cc)
        }
      }
      if(st=='C2')fig4_v2_text('\u2020 Enriched in full RNA cohort (Fig. 2)',left+3,51.8,5.3,col='grey35',just='left')
      if(st=='C3')fig4_v2_text('RNA cohort: cell-cycle enrichment (Fig. 2)',left+3,51.8,5.3,col='grey35',just='left')
      grid::grid.lines(grid::unit(left+c(3,51),'mm'),grid::unit(c(48.5,48.5),'mm'),gp=grid::gpar(col='grey82',lwd=.45))
      fig4_v2_text('Genomic associations',left+3,44.8,6.4,TRUE,just='left')
      fig4_v2_text(card$DNA,left+3,38.3,6.1,just='left')
      fig4_v2_text('Inferred TF activity',left+3,30,6.4,TRUE,just='left')
      fig4_v2_text(card$TF,left+3,25.7,5.9,just='left')
      fig4_v2_text('Microenvironment context',left+3,18.2,6.4,TRUE,just='left')
      fig4_v2_text(card$TME,left+3,12.2,5.9,just='left')
    }
    fig4_v2_text('\u2191 / \u2193: matched RNA / protein GSEA q < 0.05     \u00b7: q \u2265 0.05',87,3.2,6,col='grey30')
  }
}
register_fig4_summary_v2 <- function() {
  m<-read_source('8-Visualization/tables/figure_manifest.csv');row<-copy(m[Figure=='Figure 4'&Panel=='F'][1])
  sources<-unique(c(f4summary_read('input_hashes')$Path,
    paste0(F4SUMMARY,'source_data/',c('display_plan','all_selected_pathways','all_contrasts','displayed_cells',
      'subtype_cards','fig2_context','fig2_evidence_audit','full_rna_respiration','full_rna_cell_cycle','mutation_tests',
      'mutation_enrichments','cnv_context','input_hashes'),'.csv'),
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R',
    '8-Visualization/scripts/render_fig4_genomics_v2.R','8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV))
  row[,`:=`(Role='Full-row subtype biological summary with matched RNA/protein evidence and separate genomic/TF/TME context',
    SourceTable=paste(sources,collapse=';'),WidthMM='174',HeightMM='114',Status='AnalysisReady',
    SourceMD5='',QAFile='',RenderedAt='',
    SubtypeScope='Matched RNA/protein Main248; primary RNA Main497 context; mutation Main308; CNV Main219; distinct source populations',
    Notes='Three columns at the same174mm width as4A.13 unique pathway anchors;14 subtype/program pairs;all78 matched contrasts in source,28 glyphs displayed as associated emphases,not exclusive identities. C2 respiration remains nonsignificant in matched data; dagger links full-RNA Fig2 support. Translation discordance and C3 E2F counterexample retained. Genomic events,TF inference and method-qualified TME are separate evidence types. No new tests/models; no causal or independent-validation claims.')]
  m[which(m$Figure=='Figure 4'&m$Panel=='F')]<-row
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
render_fig4_summary_v2 <- function() {
  row<-preflight_subfigure('Figure 4','F')
  save_panel(build_fig4_summary_v2(),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure=='Figure 4'&Panel=='F',`:=`(Status='Rendered',SourceMD5=hashes,QAFile=paste0(F4SUMMARY,'QA.md'),RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
