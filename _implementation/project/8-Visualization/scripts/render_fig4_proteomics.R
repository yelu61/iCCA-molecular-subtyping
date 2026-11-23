# Protein displays: frozen protein model plus separately versioned four-collection GSEA.
# R grid/ggplot2 and shared Cairo PDF export; Arial 7-9 pt; standalone PDFs.
F4P <- '8-Visualization/refinements/fig4_proteomics_20261030_v6/'
F4G <- '4-Multiomic/refinements/fig4_protein_gsea_20261030/'
f4p_program_plan <- function() {
  names <- c('FATTY_ACID_METABOLISM','BILE_ACID_METABOLISM','XENOBIOTIC_METABOLISM',
    'PEROXISOME','OXIDATIVE_PHOSPHORYLATION','CHOLESTEROL_HOMEOSTASIS',
    'E2F_TARGETS','G2M_CHECKPOINT','MYC_TARGETS_V1','MYC_TARGETS_V2','MITOTIC_SPINDLE',
    'EPITHELIAL_MESENCHYMAL_TRANSITION','TGF_BETA_SIGNALING','ANGIOGENESIS','APICAL_JUNCTION',
    'INFLAMMATORY_RESPONSE','INTERFERON_GAMMA_RESPONSE','COMPLEMENT','COAGULATION')
  labels <- c('Fatty acid metabolism','Bile acid metabolism','Xenobiotic metabolism','Peroxisome',
    'Oxidative phosphorylation','Cholesterol homeostasis','E2F targets','G2M checkpoint',
    'MYC targets v1','MYC targets v2','Mitotic spindle','Epithelial–mesenchymal transition',
    'TGF-beta signaling','Angiogenesis','Apical junction','Inflammatory response',
    'Interferon-gamma response','Complement','Coagulation')
  data.table(Pathway=paste0('HALLMARK_',names),Label=labels,
    Block=rep(c('Metabolism &\nhepatic function','Cell cycle\n& growth','Matrix &\nsignaling','Inflammation &\nhemostasis'),c(6,5,4,4)),RowOrder=seq_along(names))
}

prepare_fig4_proteomics <- function() {
  dir.create(file.path(STEP8_ROOT,F4P,'source_data'),recursive=TRUE,showWarnings=FALSE)
  gsea <- read_source(paste0(T4,'step4_mass_spec_ranked_gsea_main_subtypes.csv'))
  deg <- read_source(paste0(T4,'step4_mass_spec_study_adjusted_deg_main_subtypes.csv'))
  mat <- readRDS(file.path(STEP8_ROOT,'4-Multiomic/data/step4_mass_spec_within_study_standardized.rds'))
  master <- master_table(FALSE);roster <- master[match(colnames(mat),Sample)]
  stopifnot(!anyNA(roster$Sample),!anyDuplicated(roster$Sample),!anyNA(mat))
  roster <- roster[SubtypeClass=='Main',.(Sample,Study,FinalSubtype)]
  stopifnot(nrow(roster)==248L,all(roster$Study %in% c('2022Dong','2023Deng')))
  setorder(roster,FinalSubtype,Study,Sample);mat <- mat[,roster$Sample,drop=FALSE]
  counts <- roster[,.N,by=.(Study,FinalSubtype)]
  assay <- read_source(paste0(T4,'step4_assay_specific_inclusion.csv'))[Analysis=='Mass-spectrometry subtype features']
  totals <- counts[,.(Total=sum(N)),by=Study]
  stopifnot(all(totals$Total==assay[match(totals$Study,Study),MainStateN]),
    nrow(deg)==3L*nrow(mat),all(table(deg$Subtype)==nrow(mat)),setequal(deg$Protein,rownames(mat)))
  stopifnot(!anyDuplicated(gsea[,.(Subtype,pathway)]),all(gsea$padj>=0 & gsea$padj<=1),all(is.finite(gsea$NES)))
  plan <- f4p_program_plan();main <- merge(plan,gsea,by.x='Pathway',by.y='pathway',all.x=TRUE)
  stopifnot(nrow(main)==57L,!anyNA(main$Subtype),all(table(main$Pathway)==3L));setorder(main,RowOrder,Subtype)
  full <- gsea[Collection=='Hallmark'];stopifnot(all(table(full$Subtype)==45L))
  full[,`:=`(Pathway=pathway,Label=tools::toTitleCase(tolower(gsub('_',' ',sub('^HALLMARK_','',pathway)))),Block='')]
  full[,RowOrder:=match(Pathway,sort(unique(Pathway)))];setorder(full,RowOrder,Subtype)
  selected <- deg[FDR<.05 & StandardizedDifference>0][order(FDR,-StandardizedDifference,Protein),head(.SD,10),by=Subtype]
  selected[,SubtypeOrder:=match(Subtype,c('C1','C2','C3'))]
  setorder(selected,SubtypeOrder,FDR,-StandardizedDifference,Protein)
  genes <- unique(selected$Protein);stopifnot(length(genes)>0,all(genes%in%rownames(mat)))
  marker <- as.data.table(as.table(mat[genes,,drop=FALSE]));setnames(marker,c('Protein','Sample','StandardizedAbundance'))
  marker <- merge(marker,roster,by='Sample');marker[,ProteinOrder:=match(Protein,genes)]
  marker[,SampleOrder:=match(Sample,roster$Sample)]
  # Full statistics are retained; selection only controls descriptive protein labels.
  labels <- deg[FDR<.05][order(FDR,-abs(StandardizedDifference),Protein),head(.SD,4),by=.(Subtype,Direction)]
  inputs <- c(paste0(T4,c('step4_mass_spec_ranked_gsea_main_subtypes.csv','step4_mass_spec_study_adjusted_deg_main_subtypes.csv',
    'step4_protein_feature_coverage_by_study.csv','step4_assay_specific_inclusion.csv')),
    '4-Multiomic/data/step4_mass_spec_within_study_standardized.rds',SUBTYPE_TRUTH_CSV,
    '4-Multiomic/scripts/run_step4_core_analysis.R','8-Visualization/scripts/render_fig4_proteomics.R',
    paste0(F4G,c('all_collections_gsea.csv','input_hashes.csv','measured_membership.csv.gz')),
    '4-Multiomic/scripts/prepare_fig4_protein_gsea.R')
  write <- function(d,n)fwrite(d,file.path(STEP8_ROOT,F4P,'source_data',paste0(n,'.csv')),na='NA')
  write(plan,'display_plan');write(main,'main_programs');write(full,'all_hallmark');write(roster,'sample_roster')
  write(counts,'study_subtype_counts');write(selected,'selected_proteins');write(marker,'selected_abundance_long')
  write(deg,'all_protein_statistics');write(labels,'volcano_labels')
  summary <- deg[,.(Proteins=.N,HigherFDR05=sum(FDR<.05 & StandardizedDifference>0),
    LowerFDR05=sum(FDR<.05 & StandardizedDifference<0)),by=Subtype]
  write(summary,'protein_summary')
  prepare_fig4_protein_display_plan(gsea,deg,selected,labels,write)
  write(data.table(Path=inputs,MD5=unname(tools::md5sum(file.path(STEP8_ROOT,inputs)))),'input_hashes')
  print(counts);print(summary)
  invisible(list(Main=main,Full=full,Selected=selected,Marker=marker,Roster=roster))
}
f4p_read <- function(n)read_source(paste0(F4P,'source_data/',n,'.csv'))

f4p_pathway_label <- function(x) {
  x<-sub('^(HALLMARK_|GOBP_|REACTOME_|KEGG_MEDICUS_REFERENCE_)','',x)
  x<-tools::toTitleCase(tolower(gsub('_',' ',x)))
  aliases<-c('Mrna Processing'='mRNA processing','Metabolism of Rna'='RNA metabolism',
    'Cell Surface Receptor Protein Serine Threonine Kinase Signaling Pathway'='Cell-surface Ser/Thr kinase receptor signaling',
    'Ecm Proteoglycans'='ECM proteoglycans',
    'Processing of Capped Intron Containing Pre Mrna'='Processing of capped intron-containing pre-mRNA',
    'Itga b Talin Vinculin Signaling Pathway'='Integrin–talin–vinculin signaling',
    'Itga b Rhog Rac Signaling Pathway'='Integrin–RhoG–Rac signaling',
    'Arno Arf Actb g Signaling Pathway'='ARNO–ARF–actin signaling',
    'Rab7 Regulated Microtubule Minus End Directed Transport'='RAB7-regulated microtubule minus-end transport',
    'Mitochondrial Complex Ucp1 in Thermogenesis'='Mitochondrial UCP1 in thermogenesis',
    '26s Proteasome Mediated Protein Degradation'='26S proteasome-mediated protein degradation')
  aliases<-c(aliases,'Dna Strand Elongation'='DNA strand elongation',
    'Dna Replication'='DNA replication','Dna Metabolic Process'='DNA metabolic process',
    'Orc1 Removal from Chromatin'='ORC1 removal from chromatin',
    'Signaling by Rho Gtpases Miro Gtpases and Rhobtb3'='Rho/MIRO/RHOBTB3 GTPase signaling',
    'Signaling by Tgf Beta Receptor Complex'='TGF-beta receptor complex signaling',
    'Keap1 Nrf2 Signaling Pathway'='KEAP1–NRF2 signaling',
    'Electron Transfer in Complex i'='Electron transfer in complex I',
    'Nuclear Export of Mrna'='Nuclear export of mRNA',
    'Arl8 Regulated Microtubule Plus End Directed Transport'='ARL8-regulated microtubule plus-end transport',
    'Gf Rtk Ras Erk Signaling Pathway'='Growth factor–RTK–RAS–ERK signaling',
    'Itga b Fak Rac Signaling Pathway'='Integrin–FAK–RAC signaling',
    'Itga b Fak Cdc42 Signaling Pathway'='Integrin–FAK–CDC42 signaling',
    'Itga b Rhogap Rhoa Signaling Pathway'='Integrin–RhoGAP–RhoA signaling',
    'Itga b Rhogef Rhoa Signaling Pathway'='Integrin–RhoGEF–RhoA signaling',
    'Egf Egfr Actin Signaling Pathway'='EGF–EGFR–actin signaling',
    'Ccr2 Gnb g Pi3k Nfkb Signaling Pathway'='CCR2–G protein–PI3K–NF-kappaB signaling',
    'Cxcr4 Gnb g Plcb Pkc Signaling Pathway'='CXCR4–G protein–PLC-beta–PKC signaling')
  replace(x,x%in%names(aliases),unname(aliases[x[x%in%names(aliases)]]))
}
prepare_fig4_protein_display_plan <- function(old,deg,selected,labels,write) {
  g<-read_source(paste0(F4G,'all_collections_gsea.csv'))
  g[,`:=`(Pathway=pathway,padj=FDRAllCollections,Label=f4p_pathway_label(pathway),
    Side=ifelse(NES>0,'Positive','Negative'))]
  membership<-read_source(paste0(F4G,'measured_membership.csv.gz'))
  members<-split(membership$Gene,membership$Pathway)
  # Full set membership, not convenient label proteins, controls display redundancy.
  pick<-function(candidates,n,used=character()) {
    keep<-character()
    for(p in candidates) {
      if(p%in%c(used,keep))next
      overlap<-vapply(keep,function(q)length(intersect(members[[p]],members[[q]]))/
        length(union(members[[p]],members[[q]])),numeric(1))
      if(!length(overlap)||max(overlap)<.5)keep<-c(keep,p)
      if(length(keep)>=n)break
    };keep
  }
  # Visualization-stage representative themes, including the contrary C3 E2F result.
  # This curated list is not a preregistered top-p rule; all raw results are exported.
  mainplan<-data.table(Subtype=rep(c('C1','C2','C3'),each=6),Side=rep(rep(c('Negative','Positive'),each=3),3),
    Pathway=c('GOBP_ACTIN_FILAMENT_BASED_PROCESS','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','HALLMARK_INFLAMMATORY_RESPONSE',
      'HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_BILE_ACID_METABOLISM','HALLMARK_XENOBIOTIC_METABOLISM',
      'GOBP_ORGANIC_ACID_METABOLIC_PROCESS','REACTOME_BIOLOGICAL_OXIDATIONS','KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION',
      'HALLMARK_TGF_BETA_SIGNALING','GOBP_CELL_JUNCTION_ORGANIZATION','KEGG_MEDICUS_REFERENCE_RAB7_REGULATED_MICROTUBULE_MINUS_END_DIRECTED_TRANSPORT',
      'REACTOME_METABOLISM_OF_RNA','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_E2F_TARGETS',
      'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_COMPLEMENT','HALLMARK_COAGULATION'),
    MainLabel=c('Actin filament process','Extracellular matrix','Inflammatory response',
      'Fatty acid metabolism','Bile acid metabolism','Xenobiotic metabolism',
      'Organic acid metabolism','Biological oxidations','Translation initiation',
      'TGF-beta signaling','Cell junction organization','RAB7 minus-end transport',
      'RNA metabolism','Fatty acid metabolism','E2F targets',
      'Epithelial–mesenchymal transition','Complement','Coagulation'))
  extra<-data.table(Subtype=rep(c('C1','C2','C3'),each=6),Side=rep(rep(c('Negative','Positive'),each=3),3),
    Pathway=c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_ANGIOGENESIS',
      'HALLMARK_PEROXISOME','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
      'HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_CHOLESTEROL_HOMEOSTASIS','HALLMARK_INTERFERON_GAMMA_RESPONSE',
      'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','REACTOME_ECM_PROTEOGLYCANS','HALLMARK_MITOTIC_SPINDLE',
      'HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','HALLMARK_PEROXISOME',
      'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_APICAL_JUNCTION','HALLMARK_INTERFERON_GAMMA_RESPONSE'),
    MainLabel=c('EMT program','TGF-beta signaling','Angiogenesis',
      'Peroxisome','Oxidative phosphorylation','Cholesterol homeostasis',
      'Fatty acid metabolism','Cholesterol homeostasis','Interferon-gamma response',
      'EMT program','ECM proteoglycans','Mitotic spindle',
      'G2M checkpoint','MYC targets v1','Peroxisome',
      'Inflammatory response','Apical junction','Interferon-gamma response'))
  mainplan<-rbindlist(list(mainplan,extra))
  mainplan[Pathway=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',MainLabel:='EMT program']
  mainplan[,DisplayOrder:=seq_len(.N),by=.(Subtype,Side)]
  main<-merge(mainplan,g,by=c('Subtype','Side','Pathway'))
  main[,Label:=MainLabel];stopifnot(nrow(main)==36L,all(main$padj<.05))
  main<-rbindlist(lapply(c('C1','C2','C3'),function(st)rbindlist(lapply(c('Negative','Positive'),function(side) {
    z<-main[Subtype==st&Side==side][order(DisplayOrder)]
    z[,RepresentativeProteins:=vapply(seq_len(.N),function(i) {
      lead<-strsplit(z$leadingEdge[i],';',fixed=TRUE)[[1]]
      hits<-deg[Subtype==st&Protein%in%lead&FDR<.05&
        sign(StandardizedDifference)==sign(z$NES[i])][order(FDR,-abs(StandardizedDifference),Protein)]
      paste(head(hits$Protein,2),collapse=', ')
    },character(1))];z
  }))))
  # Theme coverage, not a fixed collection quota: map original Hallmark anchors
  # to measured-membership neighbors; also retain all main-panel terms.
  plan<-f4p_program_plan()
  quality<-g[,.(MinQ=min(padj)),by=.(Collection,Pathway)]
  mappings<-rbindlist(lapply(c('GO_BP','Reactome'),function(collection) {
    candidates<-quality[Collection==collection&MinQ<.05,Pathway]
    rbindlist(lapply(plan$Pathway,function(anchor) {
      scores<-rbindlist(lapply(candidates,function(p) {
        overlap<-length(intersect(members[[anchor]],members[[p]]))
        data.table(Anchor=anchor,Pathway=p,Collection=collection,SharedProteins=overlap,
          Cosine=overlap/sqrt(length(members[[anchor]])*length(members[[p]])))
      }))
      scores<-scores[SharedProteins>=5 & Cosine>=.2]
      if(!nrow(scores))return(NULL)
      scores[,MinQ:=quality$MinQ[match(Pathway,quality$Pathway)]]
      scores[order(-Cosine,-SharedProteins,MinQ,Pathway)][1]
    }))
  }))
  # The smaller KEGG reference library is retained in full when significant in
  # any subtype; variant/pathogen sets remain in full source data, not this display.
  kegg<-quality[Collection=='KEGG'&MinQ<.05&grepl('^KEGG_MEDICUS_REFERENCE_',Pathway),Pathway]
  keep<-unique(c(plan$Pathway,mappings$Pathway,main$Pathway,kegg))
  supp<-g[Pathway%in%keep]
  supp[Pathway%in%plan$Pathway,Label:=plan$Label[match(Pathway,plan$Pathway)]]
  supp[,Original4C:=Pathway%in%plan$Pathway]
  # A common circular ordering of the actual three-NES profile. This is a
  # descriptive display coordinate, never a sample trajectory or biological time.
  profiles<-dcast(supp,Pathway+Collection~Subtype,value.var='NES')
  profiles[,ProfileX:=C1-(C2+C3)/2]
  profiles[,ProfileY:=sqrt(3)/2*(C2-C3)]
  profiles[,ProfileAngle:=(atan2(ProfileY,ProfileX)+pi/6)%%(2*pi)]
  profiles[,CollectionOrder:=match(Collection,c('Hallmark','GO_BP','Reactome','KEGG'))]
  setorder(profiles,CollectionOrder,ProfileAngle,Pathway)
  profiles[,RowOrder:=seq_len(.N)]
  supp<-merge(supp,profiles[,.(Pathway,ProfileAngle,RowOrder)],by='Pathway')
  supp[,Block:=Collection];setorder(supp,RowOrder,Subtype)
  write(mappings,'pathway_theme_mapping');write(profiles,'pathway_profile_order')
  write(quality[,.(Collection,Pathway,MinQ,Selected=Pathway%in%keep)],'pathway_selection_audit')
  overlap<-rbindlist(lapply(c('C1','C2','C3'),function(st) {
    a<-selected[Subtype==st,Protein];b<-labels[Subtype==st&Direction=='Subtype higher',Protein]
    data.table(Subtype=st,HeatmapPositiveN=length(a),VolcanoPositiveLabels=length(b),
      SharedN=length(intersect(a,b)),SharedProteins=paste(intersect(a,b),collapse=';'))
  }))
  main[,AbsNES:=abs(NES)]
  setorderv(main,c('Subtype','Side','AbsNES','Pathway'),c(1,1,-1,1))
  main[,DisplayOrder:=seq_len(.N),by=.(Subtype,Side)]
  write(g,'extended_gsea_all_results');write(main,'main_gsea_selection')
  write(supp,'supplement_gsea_selection');write(overlap,'protein_label_overlap')
  stopifnot(nrow(supp)==3L*uniqueN(supp$Pathway),nrow(main)==36L,
    all(plan$Pathway%in%supp$Pathway),all(overlap$SharedN==4L))
}

f4p_wrap_term <- function(text,width=25,size=6.8) {
  words<-strsplit(text,' ',fixed=TRUE)[[1]];lines<-character();line<-''
  for(word in words) {
    candidate<-paste(c(if(nzchar(line))line,word),collapse=' ')
    grob<-grid::textGrob(candidate,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=size,fontface='plain'))
    if(nzchar(line)&&grid::convertWidth(grid::grobWidth(grob),'mm',valueOnly=TRUE)>width) {
      lines<-c(lines,line);line<-word
    } else line<-candidate
  }
  paste(c(lines,line),collapse='\n')
}
f4p_integrated_programs <- function() {
  programs<-f4p_read('main_gsea_selection');d<-f4p_read('all_protein_statistics')
  labels<-f4p_read('volcano_labels');counts<-f4p_read('sample_roster')[,.N,by=FinalSubtype]
  d[,LogQ:=-log10(pmax(FDR,1e-300))];labels[,LogQ:=-log10(pmax(FDR,1e-300))]
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')];lower<-'#526779'
  return(function() {
    grid::grid.newpage();grobs<-list()
    metrics<-copy(programs[,.(Subtype,Side,Pathway,Label)])
    metrics[,TextWidthMM:=vapply(Label,function(z)grid::convertWidth(grid::grobWidth(
      grid::textGrob(z,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.8))),
      'mm',valueOnly=TRUE),numeric(1))]
    metrics[,AvailableMM:=31.5];stopifnot(all(metrics$TextWidthMM<=metrics$AvailableMM))
    fwrite(metrics,file.path(STEP8_ROOT,F4P,'term_text_fit.csv'))
    for(k in seq_along(pal)) {
      st<-names(pal)[k];left<-3+(k-1)*70
      fig4_v2_rect(left+33,116,66,4.5,pal[[st]])
      fig4_v2_text(paste0(st,'   n = ',counts[FinalSubtype==st,N]),left+33,116,8.3,TRUE)
      for(side in c('Negative','Positive')) {
        x<-left+if(side=='Negative')0 else 34.5;accent<-if(side=='Negative')lower else pal[[st]]
        fig4_v2_text(paste(side,'NES'),x,109,7,TRUE,col=accent,just='left')
        z<-programs[Subtype==st&Side==side][order(DisplayOrder)]
        stopifnot(all(diff(abs(z$NES))<=0),all(abs(z$NES)<=3.5))
        for(j in seq_len(nrow(z))) {
          y<-104-(j-1)*5.8
          lab<-z$Label[j]
          fig4_v2_text(lab,x,y,6.8,just=c('left','top'))
          # Identical physical baseline and scale in all six direction columns.
          fig4_v2_rect(x+15.75,y-3.4,31.5,.55,'#ECEFF1')
          w<-31.5*abs(z$NES[j])/3.5
          fig4_v2_rect(x+w/2,y-3.4,w,.55,accent)
        }
      }
      dd<-d[Subtype==st];ll<-labels[Subtype==st]
      p<-ggplot(dd,aes(StandardizedDifference,LogQ))+
        geom_point(data=dd[Direction=='Not significant'],color='grey82',size=.28,alpha=.5)+
        geom_point(data=dd[Direction!='Not significant'],aes(color=Direction),size=.43,alpha=.65)+
        geom_vline(xintercept=0,color='grey72',linewidth=.2)+
        geom_hline(yintercept=-log10(.05),linetype=2,color='grey65',linewidth=.2)+
        ggrepel::geom_text_repel(data=ll,aes(label=Protein),size=2.25,family=STEP8_FONT$family,
          max.overlaps=Inf,box.padding=.3,point.padding=.2,min.segment.length=0,
          ylim=c(3,23.4),force=8,force_pull=.15,
          seed=20261030,max.time=5,max.iter=50000,segment.size=.18,segment.color='grey50')+
        scale_color_manual(values=c('Subtype higher'=pal[[st]],'Subtype lower'=lower),guide='none')+
        scale_x_continuous(limits=c(-1.4,1.8),breaks=c(-1,0,1),expand=expansion(mult=.04))+
        scale_y_continuous(limits=c(0,26),breaks=c(0,10,20),expand=expansion(mult=0))+
        annotate('text',x=-1.3,y=25.4,label=sum(dd$Direction=='Subtype lower'),hjust=0,size=2.25,color=lower)+
        annotate('text',x=1.7,y=25.4,label=sum(dd$Direction=='Subtype higher'),hjust=1,size=2.25,color=pal[[st]])+
        labs(x=NULL,y=if(k==1)expression(-log[10]*'FDR') else ' ')+theme_step8()+
        theme(legend.position='none',axis.title=element_text(size=7.3),axis.text=element_text(size=6.5),
          axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.2),
          plot.margin=margin(.8,.8,.3,.8,'mm'))
      grobs[[k]]<-ggplotGrob(p)
    }
    widths<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
    for(k in seq_along(grobs)) {
      grobs[[k]]$widths<-widths
      grobs[[k]]$vp<-grid::viewport(x=grid::unit(3+(k-1)*70,'mm'),y=grid::unit(12,'mm'),
        width=grid::unit(66,'mm'),height=grid::unit(56,'mm'),just=c('left','bottom'))
      grid::grid.draw(grobs[[k]])
    }
    fig4_v2_text('Study-adjusted standardized difference',106,7.2,7.3)
    fig4_v2_text('NES bars: shared |NES| scale 0–3.5',106,2.4,6.2,col='grey35')
  })
}

f4p_supplement_height <- function() {
  ceiling(uniqueN(f4p_read('supplement_gsea_selection')$Pathway)*3.2+59)
}
f4p_extended_heatmap <- function() {
  d<-f4p_read('supplement_gsea_selection')
  plan<-unique(d[,.(Pathway,Label,RowOrder,Collection)])[order(RowOrder)]
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  counts<-f4p_read('sample_roster')[,.N,by=FinalSubtype]
  # Figure-local categorical annotation colors, unrelated to governed Study colors.
  sourcecols<-c(Hallmark='#807694',GO_BP='#6E91AB',Reactome='#6F9B86',KEGG='#B19A6A')
  sourcelabel<-c(Hallmark='Hallmark',GO_BP='GO BP',Reactome='Reactome',KEGG='KEGG MEDICUS')
  return(function() {
    grid::grid.newpage();height<-f4p_supplement_height();top<-height-25;rh<-3.2
    fig4_v2_text('Proteomic pathway programs',3,height-4,9,TRUE,just='left')
    colfun<-circlize::colorRamp2(c(-3.5,0,3.5),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))
    for(k in seq_along(pal)) {
      x<-88+(k-1)*15;fig4_v2_rect(x,height-13,13,5,pal[[k]])
      fig4_v2_text(names(pal)[k],x,height-13,8,TRUE)
      fig4_v2_text(paste0('n = ',counts[FinalSubtype==names(pal)[k],N]),x,height-19,7)
    }
    for(i in seq_len(nrow(plan))) {
      gap<-(match(plan$Collection[i],names(sourcecols))-1)*1.3
      y<-top-(i-.5)*rh-gap
      fig4_v2_rect(78,y,2.0,rh-.2,sourcecols[[plan$Collection[i]]])
      label<-paste(strwrap(plan$Label[i],width=61),collapse='\n')
      fig4_v2_text(label,74,y,6.7,just='right')
      for(k in seq_along(pal)) {
        z<-d[Pathway==plan$Pathway[i]&Subtype==names(pal)[k]];x<-88+(k-1)*15
        fig4_v2_rect(x,y,13,rh-.25,colfun(pmax(-3.5,pmin(3.5,z$NES))))
        if(z$padj<.05)fig4_v2_text(fdr_label(z$padj),x,y,7,TRUE,col=if(abs(z$NES)>2)'white' else 'black')
      }
    }
    # All legends form one compact footer; none competes with subtype headers.
    fig4_v2_text('Gene-set source',5,26.5,6.8,TRUE,just='left')
    for(k in seq_along(sourcecols)) {
      x<-7+((k-1)%%2)*33;y<-21-floor((k-1)/2)*7
      fig4_v2_rect(x,y,2.5,2.5,sourcecols[k])
      fig4_v2_text(sourcelabel[k],x+2.5,y,6.6,just='left')
    }
    vals<-seq(-3.5,3.5,length.out=100)
    fig4_v2_text('Protein NES',101.5,26.5,6.8,TRUE)
    for(i in seq_along(vals))fig4_v2_rect(80+(i-.5)*.43,21,.44,2.5,colfun(vals[i]))
    for(i in 0:2)fig4_v2_text(c('\u2264 -3.5','0','\u2265 3.5')[i+1],80+i*21.5,16,6.7)
    fig4_v2_text('* q < 0.05    ** q < 0.01    *** q < 0.001',65,4.5,6.5)
  })
}

f4p_program_heatmap <- function(full=FALSE) {
  d <- f4p_read(if(full)'all_hallmark' else 'main_programs')
  plan <- unique(d[,.(Pathway,Label,Block,RowOrder)])[order(RowOrder)]
  counts <- f4p_read('sample_roster')[,.N,by=FinalSubtype]
  pal <- load_subtype_palette()$colors[c('C1','C2','C3')]
  height <- if(full)225 else 142;top<-height-25;rh<-if(full)3.55 else 4.4
  return(function() {
    grid::grid.newpage()
    fig4_v2_text(if(full)'Protein Hallmark enrichment' else 'Proteomic pathway programs',3,height-4,9,TRUE,just='left')
    colfun<-circlize::colorRamp2(c(-3.5,0,3.5),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))
    for(k in seq_along(pal)) {
      xx<-108+(k-1)*20
      fig4_v2_rect(xx,height-13,18,5,pal[[k]])
      fig4_v2_text(names(pal)[k],xx,height-13,8,TRUE)
      fig4_v2_text(paste0('n = ',counts[FinalSubtype==names(pal)[k],N]),xx,height-19,7)
    }
    if(!full)for(b in unique(plan$Block)) {
      idx<-which(plan$Block==b);ymid<-top-(mean(idx)-.5)*rh
      fig4_v2_rect(15,ymid,27,length(idx)*rh-.8,'#F1F1F1')
      fig4_v2_text(b,3,ymid,7,just='left')
    }
    for(i in seq_len(nrow(plan))) {
      yy<-top-(i-.5)*rh;fig4_v2_text(plan$Label[i],96,yy,if(full)7 else 7.5,just='right')
      for(k in seq_along(pal)) {
        z<-d[Pathway==plan$Pathway[i]&Subtype==names(pal)[k]]
        stopifnot(nrow(z)==1)
        xx<-108+(k-1)*20
        fig4_v2_rect(xx,yy,18,rh-.3,colfun(z$NES))
        if(z$padj<.05)fig4_v2_text(fdr_label(z$padj),xx,yy,7,bold=TRUE,col=if(abs(z$NES)>2)'white' else 'black')
      }
    }
    # Shared NES color scale; source FDR remains from the full owner family.
    vals<-seq(-3.5,3.5,length.out=100)
    for(i in seq_along(vals))fig4_v2_rect(105+(i-.5)*.48,16,.49,2.7,colfun(vals[i]))
    fig4_v2_text('Protein NES',101,16,7,just='right')
    for(i in 0:2)fig4_v2_text(c('-3.5','0','3.5')[i+1],105+i*24,11,7)
    fig4_v2_text('* q < 0.05    ** q < 0.01    *** q < 0.001',87,4,7)
  })
}

build_fig4_proteomics <- function(id) {
  if(id=='C')return(f4p_integrated_programs())
  if(id=='S10P')return(f4p_extended_heatmap())
  if(id=='S10N') {
    d<-f4p_read('selected_abundance_long');samples<-f4p_read('sample_roster')
    genes<-unique(d[order(ProteinOrder),Protein])
    return(function() {
      grid::grid.newpage();pal<-load_subtype_palette()$colors;study<-load_study_palette()$colors
      fig4_v2_text('Representative protein abundance',3,156,9,TRUE,just='left')
      left<-24;right<-168;top<-141;bottom<-34;gap<-1.5;rh<-(top-bottom)/length(genes)
      colfun<-circlize::colorRamp2(c(-2.5,0,2.5),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))
      for(st in c('C1','C2','C3')) {
        ss<-samples[FinalSubtype==st]
        # Use a fixed full-cohort drawing width across all three blocks.
        w<-(168-24-2*gap)*nrow(ss)/nrow(samples);cw<-w/nrow(ss)
        z<-d[FinalSubtype==st]
        grid::grid.rect(grid::unit(left+(match(z$Sample,ss$Sample)-.5)*cw,'mm'),
          grid::unit(top-(match(z$Protein,genes)-.5)*rh,'mm'),grid::unit(cw+.01,'mm'),grid::unit(rh,'mm'),
          gp=grid::gpar(fill=colfun(pmax(-2.5,pmin(2.5,z$StandardizedAbundance))),col=NA))
        fig4_v2_rect(left+w/2,146,w,5,pal[[st]]);fig4_v2_text(paste0(st,'  n = ',nrow(ss)),left+w/2,146,8,TRUE)
        grid::grid.rect(grid::unit(left+(seq_len(nrow(ss))-.5)*cw,'mm'),grid::unit(30,'mm'),
          grid::unit(cw+.01,'mm'),grid::unit(2.5,'mm'),gp=grid::gpar(fill=study[ss$Study],col=NA))
        left<-left+w+gap
      }
      for(i in seq_along(genes))fig4_v2_text(genes[i],22,top-(i-.5)*rh,7,just='right')
      fig4_v2_text('Study',22,30,7,just='right')
      for(k in 1:2) {
        st<-c('2022Dong','2023Deng')[k];xx<-30+(k-1)*37
        fig4_v2_rect(xx,21,2.5,2.5,study[[st]]);fig4_v2_text(st,xx+2,21,7,just='left')
      }
      vals<-seq(-2.5,2.5,length.out=100)
      for(i in seq_along(vals))fig4_v2_rect(105+(i-.5)*.5,21,.51,2.5,colfun(vals[i]))
      fig4_v2_text('Within-study standardized abundance',130,12,7)
      for(i in 0:2)fig4_v2_text(c('\u2264 -2.5','0','\u2265 2.5')[i+1],105+i*25,17,7)
    })
  }
  if(id=='S10O') {
    d<-f4p_read('all_protein_statistics');labels<-f4p_read('volcano_labels')
    d[,LogQ:=-log10(pmax(FDR,1e-300))];labels[,LogQ:=-log10(pmax(FDR,1e-300))]
    cols<-c('Subtype higher'=SCALE_EFFECT_HIGH,'Subtype lower'=SCALE_EFFECT_LOW,'Not significant'='#BDBDBD')
    return(ggplot(d,aes(StandardizedDifference,LogQ,color=Direction))+
      geom_point(size=.45,alpha=.5)+geom_vline(xintercept=0,color='grey65',linewidth=.3)+
      geom_hline(yintercept=-log10(.05),linetype=2,color='grey60',linewidth=.3)+
      ggrepel::geom_text_repel(data=labels,aes(label=Protein),size=2.3,max.overlaps=Inf,show.legend=FALSE,
        box.padding=.2,point.padding=.1,min.segment.length=0,seed=20261030)+
      facet_wrap(~Subtype,nrow=1)+scale_color_manual(values=cols)+theme_step8()+
      labs(x='Study-adjusted standardized difference',y='-log10(FDR)',color=NULL)+
      theme(legend.position='top'))
  }
  stop('Unknown protein display: ',id)
}

register_fig4_proteomics <- function() {
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  base<-copy(m[Figure=='Figure 4'&Panel=='A'][1])
  inputs<-f4p_read('input_hashes')
  common<-c(inputs$Path,paste0(F4P,'source_data/',c('input_hashes','sample_roster','display_plan','main_programs',
    'all_hallmark','selected_proteins','selected_abundance_long','all_protein_statistics','volcano_labels',
    'main_gsea_selection','supplement_gsea_selection','extended_gsea_all_results','protein_label_overlap','pathway_theme_mapping',
    'pathway_profile_order','pathway_selection_audit'),'.csv'),
    '8-Visualization/scripts/render_fig4_genomics_v2.R','8-Visualization/scripts/fig_style.R',
    '8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/figure_config.json',
    'governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
  ids<-c('C','S10N','S10P')
  files<-c(C='main/Figure4C_proteomic_programs.pdf',S10N='supplement/Supplement_S10N_protein_abundance.pdf',
    S10O='supplement/Supplement_S10O_protein_volcano.pdf',S10P='supplement/Supplement_S10P_protein_gsea_programs.pdf')
  roles<-c(C='Differential proteins and representative ranked GSEA programs',S10N='Selected protein abundance (descriptive)',S10O='All protein contrast statistics',
    S10P='Four-collection representative protein GSEA programs')
  for(id in ids) {
    row<-copy(base)
    src<-paste(unique(common),collapse=';')
    row[,`:=`(Figure=if(startsWith(id,'S'))'Supplement' else 'Figure 4',Panel=id,Role=roles[[id]],AnalysisOwner='Step4',
      SourceTable=src,SourceFigure='',TargetOutput=paste0('8-Visualization/figures/',files[[id]]),Status='AnalysisReady',
      WidthMM=if(id=='C')'212' else if(id=='S10P')'130' else '174',HeightMM=as.character(c(C=120,S10N=160,S10O=105,S10P=f4p_supplement_height())[[id]]),
      Notes='20261030 v6: frozen model/GSEA; C36 terms sorted by absolute NES with shared-scale bars; 212x120mm single-line terms and coordinated plot geometry; S10N top10 positive proteins per subtype; S10P73 source-annotated pathways sorted by actual three-NES profiles; 130mm width and grouped footer legends. Full expanded-family FDR/results retained. Profile order is descriptive, not progression.',
      SourceMD5='',QAFile='',RenderedAt='')]
    row[,SubtypeScope:='Mass-spec Main n=248 (C1=51,C2=126,C3=71); Dong184/Deng64']
    idx<-which(m$Figure==row$Figure&m$Panel==id)
    if(length(idx))m[idx]<-row else m<-rbindlist(list(m,row),use.names=TRUE)
  }
  m<-m[!(Figure=='Figure 4'&Panel=='D') & !(Figure=='Supplement'&Panel=='S10O')]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}

render_fig4_proteomics <- function(id) {
  figure<-if(startsWith(id,'S'))'Supplement' else 'Figure 4'
  row<-preflight_subfigure(figure,id)
  save_panel(build_fig4_proteomics(id),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure==figure&Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
  message('Rendered protein display: ',row$TargetOutput)
}
