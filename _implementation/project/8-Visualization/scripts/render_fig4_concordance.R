# Matched RNA/protein concordance; frozen owner statistics, R/Cairo PDF only.
F4D <- '8-Visualization/refinements/fig4_concordance_20261030_v6/'
F4DG <- '4-Multiomic/refinements/fig4_matched_four_gsea_20261030/'
f4d_mode_order<-c('Both higher','Both lower','RNA higher / protein lower',
  'RNA lower / protein higher','RNA only','Protein only')
f4d_select_by_rules <- function(pathway) {
  d<-copy(pathway)
  d[,`:=`(RNA_Sig=padj_RNA<.05,Protein_Sig=padj_Protein<.05)]
  d[,Mode:=ifelse(RNA_Sig&Protein_Sig,
    ifelse(NES_RNA>0&NES_Protein>0,'Both higher',ifelse(NES_RNA<0&NES_Protein<0,'Both lower',
      ifelse(NES_RNA>0,'RNA higher / protein lower','RNA lower / protein higher'))),
    ifelse(RNA_Sig,'RNA only',ifelse(Protein_Sig,'Protein only','Neither')))]
  d[,ModeOrder:=match(Mode,f4d_mode_order)]
  d[,EvidenceQ:=ifelse(RNA_Sig&Protein_Sig,pmax(padj_RNA,padj_Protein),pmin(padj_RNA,padj_Protein))]
  d[,JointNES:=ifelse(RNA_Sig&Protein_Sig,pmin(abs(NES_RNA),abs(NES_Protein)),
    ifelse(RNA_Sig,abs(NES_RNA),abs(NES_Protein))) ]
  candidates<-d[Mode!='Neither'];setorder(candidates,Subtype,ModeOrder,EvidenceQ,-JointNES,Pathway)
  members<-read_source(paste0(F4DG,'membership.csv.gz'))
  ids<-sort(unique(candidates$Pathway));genes<-sort(unique(members$Gene))
  z<-members[Pathway%in%ids]
  mat<-Matrix::sparseMatrix(i=match(z$Pathway,ids),j=match(z$Gene,genes),x=1,
    dims=c(length(ids),length(genes)))
  shared<-as.matrix(Matrix::tcrossprod(mat));sizes<-as.numeric(Matrix::rowSums(mat))
  names(sizes)<-ids;dimnames(shared)<-list(ids,ids)
  audits<-list();n<-0L
  for(st in c('C1','C2','C3'))for(mode in f4d_mode_order) {
    z<-copy(candidates[Subtype==st&Mode==mode]);if(!nrow(z))next
    kept<-character();z[,`:=`(ClusterRepresentative='',OverlapJaccard=0,OverlapContainment=0)]
    for(i in seq_len(nrow(z))) {
      p<-z$Pathway[i];hits<-character();jac<-contain<-numeric()
      if(length(kept)) {
        ov<-shared[p,kept];jac<-ov/(sizes[p]+sizes[kept]-ov);contain<-ov/pmin(sizes[p],sizes[kept])
        names(jac)<-names(contain)<-kept
        hits<-kept[jac>=.5|contain>=.8]
      }
      if(length(hits)) {
        q<-hits[1];z[i,`:=`(ClusterRepresentative=q,OverlapJaccard=jac[q],OverlapContainment=contain[q])]
      } else {kept<-c(kept,p);z[i,ClusterRepresentative:=p]}
    }
    z[,Nonredundant:=Pathway==ClusterRepresentative]
    z[,MainSelected:=Pathway==kept[1]]
    n<-n+1L;audits[[n]]<-z
  }
  audit<-rbindlist(audits)
  audit[,`:=`(MainSelected=FALSE,MainOverlapException=FALSE)]
  # Cover modes in a fixed order, while also preventing redundant main entries
  # across modes of the same subtype. Never deduplicate across subtypes.
  for(st in c('C1','C2','C3')) {
    chosen<-character()
    for(mode in f4d_mode_order) {
      rows<-which(audit$Subtype==st&audit$Mode==mode&audit$Nonredundant)
      if(!length(rows))next
      valid<-vapply(rows,function(i) {
        if(!length(chosen))return(TRUE)
        p<-audit$Pathway[i];ov<-shared[p,chosen]
        !any(ov/(sizes[p]+sizes[chosen]-ov)>=.5|ov/pmin(sizes[p],sizes[chosen])>=.8)
      },logical(1))
      i<-if(any(valid))rows[which(valid)[1]] else rows[1]
      audit[i,`:=`(MainSelected=TRUE,MainOverlapException=!any(valid))]
      chosen<-c(chosen,audit$Pathway[i])
    }
  }
  main<-copy(audit[MainSelected==TRUE])
  setorder(main,Subtype,ModeOrder);main[,ID:=seq_len(.N)]
  main[,SourceCode:=c(Hallmark='H',GO_BP='GO',Reactome='R',KEGG='K')[Collection]]
  main[,Label:=f4d_rule_label(Pathway)]
  main[,`:=`(BothSignificant=RNA_Sig&Protein_Sig,Direction=ifelse(NES_RNA>0&NES_Protein>0,'Both higher',
    ifelse(NES_RNA<0&NES_Protein<0,'Both lower',ifelse(NES_RNA>0,'RNA higher / protein lower','RNA lower / protein higher'))))]
  coverage<-CJ(Subtype=c('C1','C2','C3'),Mode=f4d_mode_order,sorted=FALSE)
  coverage<-merge(coverage,audit[,.(Candidates=.N,Representatives=sum(Nonredundant),
    MainEntries=sum(MainSelected)),by=.(Subtype,Mode)],by=c('Subtype','Mode'),all.x=TRUE)
  for(col in c('Candidates','Representatives','MainEntries'))set(coverage,which(is.na(coverage[[col]])),col,0L)
  full<-merge(d,audit[,.(Subtype,Pathway,ClusterRepresentative,Nonredundant,MainSelected,
    OverlapJaccard,OverlapContainment,MainOverlapException)],by=c('Subtype','Pathway'),all.x=TRUE)
  full[Mode=='Neither',`:=`(Nonredundant=FALSE,MainSelected=FALSE)]
  stopifnot(all(coverage[Candidates>0,MainEntries]==1L),all(coverage[Candidates==0,MainEntries]==0L))
  fwrite(audit,file.path(STEP8_ROOT,F4D,'source_data/selection_candidates.csv'))
  fwrite(coverage,file.path(STEP8_ROOT,F4D,'source_data/mode_coverage.csv'))
  fwrite(full,file.path(STEP8_ROOT,F4D,'source_data/selection_audit.csv'))
  main
}
f4d_rule_label <- function(x) {
  x<-gsub('_',' ',sub('^(HALLMARK_|GOBP_|REACTOME_|KEGG_MEDICUS_REFERENCE_)','',x))
  x<-tools::toTitleCase(tolower(x))
  # Readability aliases are applied only after rule-based selection; they never
  # enter candidate eligibility, clustering or ranking.
  aliases<-c('Organic Acid Metabolic Process'='Organic acid metabolism',
    'Cell Projection Organization'='Cell projection organization',
    'Actin Filament Based Process'='Actin filament processes',
    'Formation of Fibrin Clot Clotting Cascade'='Fibrin clot formation',
    'Synapse Organization'='Synapse organization','Cytoplasmic Translation'='Cytoplasmic translation',
    'Negative Regulation of Nucleobase Containing Compound Metabolic Process'='Nucleobase metabolism (negative regulation)',
    'Response to Elevated Platelet Cytosolic Ca2'='Platelet cytosolic Ca2+ response',
    'Metabolism of Amino Acids and Derivatives'='Amino acid and derivative metabolism',
    'Ribosome Associated Quality Control'='Ribosome-associated quality control',
    'Metabolism of Rna'='RNA metabolism','Rna Processing'='RNA processing',
    'Nervous System Development'='Nervous system development')
  replace(x,x%in%names(aliases),unname(aliases[x[x%in%names(aliases)]]))
}

prepare_fig4_concordance <- function() {
  dir.create(file.path(STEP8_ROOT,F4D,'source_data'),recursive=TRUE,showWarnings=FALSE)
  gene<-read_source(paste0(T4,'step4_matched_rna_protein_subtype_contrasts_long.csv'))
  summary<-read_source(paste0(T4,'step4_matched_rna_protein_subtype_concordance.csv'))
  pathway<-read_source(paste0(F4DG,'paired_gsea.csv'))
  ps<-pathway[,.(Pathways=.N,SpearmanRho=cor(NES_RNA,NES_Protein,method='spearman')),by=.(Subtype,Collection)]
  e<-new.env(parent=emptyenv());load(file.path(STEP8_ROOT,'1-Clustering/data/combined_expr_after_correction.rds'),envir=e)
  rna<-e$corrected_expr_limma
  protein<-readRDS(file.path(STEP8_ROOT,'4-Multiomic/data/step4_mass_spec_within_study_standardized.rds'))
  shared_samples<-intersect(colnames(protein),colnames(rna));shared_genes<-intersect(rownames(protein),rownames(rna))
  master<-master_table(FALSE)
  roster<-master[Sample%in%shared_samples&SubtypeClass=='Main',.(Sample,Study,FinalSubtype)]
  setorder(roster,FinalSubtype,Study,Sample)
  stopifnot(nrow(roster)==248L,length(shared_genes)==4706L,nrow(gene)==14118L,
    all(table(gene$Subtype)==4706L),!anyDuplicated(gene[,.(Subtype,Gene)]),
    setequal(gene$Gene,shared_genes),all(is.finite(gene$RNAModeratedT)),all(is.finite(gene$ProteinModeratedT)),
    nrow(pathway)==8940L,!anyDuplicated(pathway[,.(Subtype,pathway)]),
    all(table(pathway[Collection=='Hallmark',Subtype])==45L),
    all(table(pathway[Collection=='Reactome',Subtype])==671L),
    all(table(pathway[Collection=='GO_BP',Subtype])==2223L),
    all(table(pathway[Collection=='KEGG',Subtype])==41L),
    all(is.finite(pathway$NES_RNA)),all(is.finite(pathway$NES_Protein)),
    all(pathway$padj_RNA>=0&pathway$padj_RNA<=1),
    all(pathway$padj_Protein>=0&pathway$padj_Protein<=1))
  # Descriptive verification of the saved coefficient correlations; no new test.
  audit<-gene[,.(ReproducedRho=cor(RNAModeratedT,ProteinModeratedT,method='spearman'),
    SameSign=sum(sign(RNAModeratedT)==sign(ProteinModeratedT)),Genes=.N),by=Subtype]
  audit<-merge(audit,summary[,.(Subtype,SavedRho=SpearmanRho)],by='Subtype')
  stopifnot(all(abs(audit$ReproducedRho-audit$SavedRho)<1e-12))
  pathway[,`:=`(Pathway=pathway,BothSignificant=padj_RNA<.05&padj_Protein<.05,
    Direction=ifelse(NES_RNA>0&NES_Protein>0,'Both higher',
      ifelse(NES_RNA<0&NES_Protein<0,'Both lower',
        ifelse(NES_RNA>0,'RNA higher / protein lower','RNA lower / protein higher'))))]
  selected<-f4d_select_by_rules(pathway)
  plan<-selected[,.(Subtype,Pathway,Label,SourceCode,Mode,ModeOrder,ID,EvidenceQ,JointNES)]
  setorder(selected,Subtype,ID)
  hall<-copy(pathway[Collection=='Hallmark']);hall[,Selected:=paste(Subtype,Pathway)%in%paste(plan$Subtype,plan$Pathway)]
  hall[,ID:=plan$ID[match(paste(Subtype,Pathway),paste(plan$Subtype,plan$Pathway))]]
  write<-function(x,n)fwrite(x,file.path(STEP8_ROOT,F4D,'source_data',paste0(n,'.csv')),na='NA')
  write(gene,'gene_contrasts');write(summary,'gene_summary');write(pathway,'all_pathway_results')
  write(ps,'pathway_summary');write(roster,'matched_sample_roster');write(audit,'correlation_audit')
  write(plan,'pathway_display_plan');write(selected,'selected_pathways');write(hall,'all_hallmark_display')
  write(pathway[,.N,by=.(Subtype,Direction,BothSignificant)],'pathway_quadrant_counts')
  inputs<-c(paste0(T4,c('step4_matched_rna_protein_subtype_contrasts_long.csv',
    'step4_matched_rna_protein_subtype_concordance.csv','step4_matched_rna_protein_pathway_concordance.csv',
    'step4_matched_rna_protein_pathway_concordance_summary.csv')),
    '1-Clustering/data/combined_expr_after_correction.rds',
    '4-Multiomic/data/step4_mass_spec_within_study_standardized.rds',SUBTYPE_TRUTH_CSV,
    '4-Multiomic/scripts/run_step4_core_analysis.R','8-Visualization/scripts/render_fig4_concordance.R',
    paste0(F4DG,c('paired_gsea.csv','input_hashes.csv','coverage.csv','membership.csv.gz')),
    '4-Multiomic/scripts/prepare_fig4_matched_gsea.R')
  write(data.table(Path=inputs,MD5=unname(tools::md5sum(file.path(STEP8_ROOT,inputs)))),'input_hashes')
  print(roster[,.N,by=.(Study,FinalSubtype)]);print(audit)
}
f4d_read<-function(n)read_source(paste0(F4D,'source_data/',n,'.csv'))
f4d_place <- function(grobs,bottom,height) {
  widths<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
  for(k in seq_along(grobs)) {
    grobs[[k]]$widths<-widths
    grobs[[k]]$vp<-grid::viewport(x=grid::unit(2+(k-1)*58,'mm'),y=grid::unit(bottom,'mm'),
      width=grid::unit(54,'mm'),height=grid::unit(height,'mm'),just=c('left','bottom'))
    grid::grid.draw(grobs[[k]])
  }
}
f4d_wrap_rule_label <- function(text) {
  words<-strsplit(text,' ',fixed=TRUE)[[1]];lines<-character();line<-''
  for(word in words) {
    trial<-paste(c(if(nzchar(line))line,word),collapse=' ')
    w<-grid::convertWidth(grid::grobWidth(grid::textGrob(trial,gp=grid::gpar(
      fontfamily=STEP8_FONT$family,fontsize=7.8))),'mm',valueOnly=TRUE)
    if(nzchar(line)&&w>51){lines<-c(lines,line);line<-word}else line<-trial
  }
  paste(c(lines,line),collapse='\n')
}
build_fig4_concordance <- function(id) {
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  counts<-f4d_read('matched_sample_roster')[,.N,by=FinalSubtype]
  if(id=='D1') {
    d<-f4d_read('gene_contrasts');s<-f4d_read('gene_summary')
    xr<-ceiling(max(abs(d$RNAModeratedT)));yr<-ceiling(max(abs(d$ProteinModeratedT)))
    return(function() {
      grid::grid.newpage();grobs<-list()
      fig4_v2_text('RNA–protein subtype contrast concordance',3,86,9,TRUE,just='left')
      fig4_v2_text('248 matched tumors  |  4,706 shared genes',3,80,7,col='grey35',just='left')
      for(k in seq_along(pal)) {
        st<-names(pal)[k];x<-29+(k-1)*58
        fig4_v2_rect(x,72,54,5,pal[[st]])
        fig4_v2_text(paste0(st,'   n = ',counts[FinalSubtype==st,N]),x,72,8.3,TRUE)
        fig4_v2_text(sprintf('Spearman rho = %.3f',s[Subtype==st,SpearmanRho]),x,65.5,7)
        p<-ggplot(d[Subtype==st],aes(RNAModeratedT,ProteinModeratedT))+
          geom_hline(yintercept=0,color='grey70',linewidth=.25)+
          geom_vline(xintercept=0,color='grey70',linewidth=.25)+
          geom_point(color=pal[[st]],size=.42,alpha=.24,stroke=0)+
          scale_x_continuous(limits=c(-xr,xr),breaks=scales::breaks_pretty(5),expand=expansion(mult=.03))+
          scale_y_continuous(limits=c(-yr,yr),breaks=scales::breaks_pretty(5),expand=expansion(mult=.03))+
          labs(x=NULL,y=if(k==1)'Protein moderated t' else ' ')+theme_step8()+
          theme(axis.title=element_text(size=7.3),axis.text=element_text(size=6.5),
            axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.2),
            plot.margin=margin(.8,.8,.5,.8,'mm'))
        grobs[[k]]<-ggplotGrob(p)
      }
      f4d_place(grobs,14,47)
      fig4_v2_text('RNA moderated t',87,7,7.3)
    })
  }
  if(id=='D1_vertical') {
    d<-f4d_read('gene_contrasts');s<-f4d_read('gene_summary')
    xr<-ceiling(max(abs(d$RNAModeratedT)));yr<-ceiling(max(abs(d$ProteinModeratedT)))
    return(function() {
      grid::grid.newpage();grobs<-list()
      fig4_v2_text('RNA–protein concordance',4,212,9,TRUE,just='left')
      fig4_v2_text('248 matched tumors | 4,706 shared genes',4,206,6.8,col='grey35',just='left')
      for(k in seq_along(pal)) {
        st<-names(pal)[k];y<-197-(k-1)*65
        fig4_v2_rect(42.5,y,77,5,pal[[st]])
        fig4_v2_text(paste0(st,'   n = ',counts[FinalSubtype==st,N]),42.5,y,8.3,TRUE)
        fig4_v2_text(sprintf('Spearman rho = %.3f',s[Subtype==st,SpearmanRho]),42.5,y-6,7)
        p<-ggplot(d[Subtype==st],aes(RNAModeratedT,ProteinModeratedT))+
          geom_hline(yintercept=0,color='grey70',linewidth=.25)+
          geom_vline(xintercept=0,color='grey70',linewidth=.25)+
          geom_point(color=pal[[st]],size=.42,alpha=.24,stroke=0)+
          scale_x_continuous(limits=c(-xr,xr),breaks=scales::breaks_pretty(5),expand=expansion(mult=.03))+
          scale_y_continuous(limits=c(-yr,yr),breaks=scales::breaks_pretty(5),expand=expansion(mult=.03))+
          labs(x=NULL,y='Protein moderated t')+theme_step8()+
          theme(axis.title=element_text(size=7.3),axis.text=element_text(size=6.5),
            axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.2),
            plot.margin=margin(.8,.8,.5,.8,'mm'))
        grobs[[k]]<-ggplotGrob(p)
      }
      widths<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
      for(k in seq_along(grobs)) {
        grobs[[k]]$widths<-widths
        grobs[[k]]$vp<-grid::viewport(x=grid::unit(4,'mm'),y=grid::unit(140-(k-1)*65,'mm'),
          width=grid::unit(77,'mm'),height=grid::unit(47,'mm'),just=c('left','bottom'))
        grid::grid.draw(grobs[[k]])
      }
      fig4_v2_text('RNA moderated t',42.5,3.8,7.3)
    })
  }
  if(id=='D2') {
    d<-f4d_read('selected_pathways');plan<-f4d_read('pathway_display_plan')
    lim<-max(4,ceiling(max(abs(c(d$NES_RNA,d$NES_Protein)))*2)/2)
    d[,PointFill:=ifelse(BothSignificant,unname(pal[Subtype]),'white')]
    return(function() {
      grid::grid.newpage()
      fig4_v2_text('RNA–protein pathway concordance',3,116,9.5,TRUE,just='left')
      fig4_v2_text(sprintf('Matched cohort: 248 tumors  |  %d pathway comparisons',nrow(d)),3,109.5,7.5,col='grey35',just='left')
      p<-ggplot(d,aes(NES_RNA,NES_Protein))+
        annotate('rect',xmin=0,xmax=lim,ymin=0,ymax=lim,fill='#F0F5EF',color=NA)+
        annotate('rect',xmin=-lim,xmax=0,ymin=-lim,ymax=0,fill='#F0F2F4',color=NA)+
        annotate('rect',xmin=-lim,xmax=0,ymin=0,ymax=lim,fill='#FCF7F2',color=NA)+
        annotate('rect',xmin=0,xmax=lim,ymin=-lim,ymax=0,fill='#FCF7F2',color=NA)+
        geom_hline(yintercept=0,color='grey60',linewidth=.25)+
        geom_vline(xintercept=0,color='grey60',linewidth=.25)+
        annotate('text',x=lim*.08,y=lim*.96,label='Both higher',hjust=0,vjust=1,size=2.45,
          color='grey40',family=STEP8_FONT$family)+
        annotate('text',x=-lim*.95,y=-lim*.94,label='Both lower',hjust=0,vjust=0,size=2.45,
          color='grey40',family=STEP8_FONT$family)+
        annotate('text',x=-lim*.95,y=lim*.94,label='RNA− / protein+',hjust=0,vjust=1,size=2.45,
          color='grey40',family=STEP8_FONT$family)+
        annotate('text',x=lim*.95,y=-lim*.94,label='RNA+ / protein−',hjust=1,vjust=0,size=2.45,
          color='grey40',family=STEP8_FONT$family)+
        geom_point(aes(color=Subtype,fill=PointFill),shape=21,size=2.6,stroke=.5)+
        scale_color_manual(values=pal,guide='none')+scale_fill_identity()+
        ggrepel::geom_text_repel(aes(label=sprintf('%02d',ID)),family=STEP8_FONT$family,
          size=2.65,seed=20261030,max.overlaps=Inf,box.padding=.25,point.padding=.2,ylim=c(-3.35,3.35),
          min.segment.length=0,segment.size=.18,segment.color='grey50',force=5,
          max.iter=90000,max.time=7)+
        scale_x_continuous(limits=c(-lim,lim),breaks=seq(-lim,lim,2),expand=expansion(mult=.02))+
        scale_y_continuous(limits=c(-lim,lim),breaks=seq(-lim,lim,2),expand=expansion(mult=.02))+
        coord_fixed()+labs(x='RNA NES',y='Protein NES')+theme_step8()+
        theme(axis.title=element_text(size=8),axis.text=element_text(size=7.3),
          axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.2),
          plot.margin=margin(.8,.8,.5,.8,'mm'))
      g<-ggplotGrob(p);g$vp<-grid::viewport(x=grid::unit(65,'mm'),y=grid::unit(12,'mm'),
        width=grid::unit(102,'mm'),height=grid::unit(95,'mm'),just=c('left','bottom'))
      grid::grid.draw(g)
      # Align the complete significance key with the actual square panel center.
      grid::grid.force()
      vnames<-grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name
      panelvp<-vnames[grepl('^panel\\.',vnames)][1]
      stopifnot(!is.na(panelvp));grid::seekViewport(panelvp)
      left<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)$x*25.4
      right<-grid::deviceLoc(grid::unit(1,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)$x*25.4
      bottompanel<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)$y*25.4
      toppanel<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(1,'npc'),valueOnly=TRUE)$y*25.4
      grid::upViewport(0);center<-(left+right)/2
      # Equal top/bottom alignment gives both columns the same visual extent.
      gp<-grid::gpar(fontfamily=STEP8_FONT$family,fontsize=7.8,lineheight=1.05,col='grey15')
      labelsleft<-vapply(seq_len(nrow(d)),function(i)
        f4d_wrap_rule_label(paste0(d$Label[i],' [',d$SourceCode[i],']')),character(1))
      heights<-vapply(labelsleft,function(x)grid::convertHeight(grid::grobHeight(
        grid::textGrob(x,gp=gp)),'mm',valueOnly=TRUE),numeric(1))
      widthsleft<-vapply(labelsleft,function(x)grid::convertWidth(grid::grobWidth(
        grid::textGrob(x,gp=gp)),'mm',valueOnly=TRUE),numeric(1))
      groupgap<-4.2;rowgap<-(toppanel-bottompanel-sum(heights)-2*groupgap)/(nrow(d)-3)
      stopifnot(rowgap>1.5,all(widthsleft<=51))
      rowtop<-toppanel;metrics<-list()
      for(k in seq_along(pal)) {
        st<-names(pal)[k];indices<-which(d$Subtype==st);grouptop<-rowtop
        for(i in indices) {
          grid::grid.text(labelsleft[i],grid::unit(11,'mm'),grid::unit(rowtop,'mm'),
            just=c('left','top'),gp=gp)
          fig4_v2_text(sprintf('%02d',d$ID[i]),5,rowtop-1.4,7.8,TRUE,col=pal[[st]],just='left')
          bottom<-rowtop-heights[i]
          metrics[[i]]<-data.table(Subtype=st,ID=d$ID[i],PlotText=labelsleft[i],
            WidthMM=widthsleft[i],HeightMM=heights[i],TopMM=rowtop,BottomMM=bottom)
          rowtop<-bottom+if(i==tail(indices,1))-groupgap else -rowgap
        }
        fig4_v2_rect(2.5,(grouptop+bottom)/2,1.1,grouptop-bottom,pal[[st]])
      }
      metrics<-rbindlist(metrics);stopifnot(min(metrics$BottomMM)>6.5)
      fwrite(metrics,file.path(STEP8_ROOT,F4D,'D2_text_metrics.csv'))
      labels<-c('Both q < 0.05','One layer q < 0.05')
      tw<-vapply(labels,function(x)grid::convertWidth(grid::grobWidth(grid::textGrob(x,
        gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=7.2))),'mm',valueOnly=TRUE),numeric(1))
      widths<-3+tw;gap<-8;legendwidth<-sum(widths)+gap;start<-center-legendwidth/2
      for(k in 1:2) {
        x<-start+if(k==1)0 else widths[1]+gap
        grid::grid.points(grid::unit(x,'mm'),grid::unit(7.5,'mm'),pch=21,size=grid::unit(2.5,'mm'),
          gp=grid::gpar(col='grey30',fill=if(k==1)'grey30' else 'white'))
        fig4_v2_text(labels[k],x+3,7.5,7.2,just='left')
      }
      fwrite(data.table(PanelLeftMM=left,PanelRightMM=right,PanelTopMM=toppanel,PanelBottomMM=bottompanel,
        SidebarTopMM=max(metrics$TopMM),SidebarBottomMM=min(metrics$BottomMM),RowGapMM=rowgap,PanelCenterMM=center,
        LegendCenterMM=start+legendwidth/2,LegendWidthMM=legendwidth),
        file.path(STEP8_ROOT,F4D,'D2_layout_metrics.csv'))
    })
  }
  stop('Unknown concordance panel: ',id)
}
register_fig4_concordance <- function() {
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  inputs<-f4d_read('input_hashes')
  sources<-unique(c(inputs$Path,paste0(F4D,'source_data/',
    c('input_hashes','gene_contrasts','gene_summary','all_pathway_results','pathway_summary',
      'matched_sample_roster','correlation_audit','pathway_display_plan','selected_pathways',
      'all_hallmark_display','pathway_quadrant_counts','selection_candidates','mode_coverage','selection_audit'),'.csv'),
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R',
    '8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/render_fig4_genomics_v2.R','8-Visualization/scripts/figure_config.json',
    SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV))
  for(id in c('D1','D1_vertical','D2')) {
    idx<-which(m$Figure=='Figure 4'&m$Panel==id)
    if(!length(idx)) {
      row<-copy(m[Figure=='Figure 4'&Panel=='D1'][1]);row[,Panel:=id]
      m<-rbindlist(list(m,row),use.names=TRUE)
    }
    filename<-c(D1='Figure4D1_rna_protein_concordance.pdf',
      D1_vertical='Figure4D1_rna_protein_concordance_vertical.pdf',D2='Figure4D2_pathway_quadrant.pdf')[[id]]
    m[Figure=='Figure 4'&Panel==id,`:=`(
      Role=if(id=='D2')'Matched four-collection pathway concordance, single quadrant' else if(id=='D1_vertical')'Matched gene-level concordance (vertical alternative)' else 'Matched gene-level contrast concordance',
      SourceFigure='',SourceTable=paste(sources,collapse=';'),
      TargetOutput=paste0('8-Visualization/figures/main/',filename),
      WidthMM=if(id=='D1_vertical')'85' else if(id=='D2')'168' else '174',
      HeightMM=as.character(c(D1=90,D1_vertical=216,D2=120)[[id]]),Status='AnalysisReady',SourceMD5='',QAFile='',RenderedAt='',
      SubtypeScope='248 matched Main tumors (C1=51,C2=126,C3=71); 4706 common genes',
      Notes='D1 unchanged frozen gene-level contrasts; vertical alternative C1/C2/C3 top-to-bottom. D2 one shared quadrant, subtype colors,rule-selected entries from four-collection candidates; no manual functional theme or database quotas; D2 balanced168x120mm with7.8pt inline-bracket terms, matched sidebar/quadrant vertical extents, stripe-only groups and panel-centered significance legend; matched2980-pathway BH per subtype/view. Full8940 paired entries retained separately; not full-cohort single-omics GSEA.')]
  }
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
render_fig4_concordance <- function(id) {
  row<-preflight_subfigure('Figure 4',id)
  save_panel(build_fig4_concordance(id),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure=='Figure 4'&Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
