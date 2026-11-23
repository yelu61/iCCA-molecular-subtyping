# Native R/grid alluvial: every ribbon is a genuine same-sample three-label path.
F4ALL <- '8-Visualization/refinements/fig4_alluvial_20261030/'
F4DUAL <- '4-Multiomic/refinements/fig4_mofa_dual_20261030/'
f4all_read<-function(name)read_source(paste0(F4ALL,'source_data/',name,'.csv'))
prepare_fig4_alluvial <- function() {
  dir.create(file.path(STEP8_ROOT,F4ALL,'source_data'),recursive=TRUE,showWarnings=FALSE)
  dual<-jsonlite::fromJSON(file.path(STEP8_ROOT,F4DUAL,'run_manifest.json'))
  stopifnot(isTRUE(dual$converged),dual$fit_samples==204L,dual$main_samples==182L)
  inputs<-read_source(paste0(F4DUAL,'input_hashes.csv'))
  outputs<-read_source(paste0(F4DUAL,'output_hashes.csv'))
  stopifnot(identical(unname(tools::md5sum(file.path(STEP8_ROOT,inputs$Path))),inputs$MD5),
    identical(unname(tools::md5sum(file.path(STEP8_ROOT,F4DUAL,outputs$File))),outputs$MD5))
  d<-read_source(paste0(F4DUAL,'main_joined_assignments.csv'));d[,Sample:=as.character(Sample)]
  master<-master_table(FALSE);m<-master[match(d$Sample,Sample)]
  stopifnot(nrow(d)==182L,!anyDuplicated(d$Sample),all(m$SubtypeClass=='Main'),
    all(m$Study=='2022Dong'),identical(d$FinalSubtype,m$FinalSubtype),
    identical(as.integer(table(d$FinalSubtype)),c(42L,92L,48L)))
  paths<-read_source(paste0(F4DUAL,'alluvial_paths.csv'))
  stopifnot(nrow(paths)==27L,sum(paths$N)==182L)
  agreement<-read_source(paste0(F4DUAL,'agreement.csv'))
  primary_order<-c('C1','C2','C3')
  permutations<-list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  order_for<-function(field) {
    ids<-sort(unique(d[[field]]));stopifnot(length(ids)==3L)
    scores<-vapply(permutations,function(p)sum(vapply(1:3,function(k)
      sum(d$FinalSubtype==primary_order[k]&d[[field]]==ids[p[k]]),numeric(1))),numeric(1))
    ids[permutations[[which.max(scores)]]]
  }
  rp_order<-order_for('RPCluster');rpc_order<-order_for('RPCCluster')
  nodes<-list();tables<-list();unit_height<-(85-21-2*4)/182
  for(axis in c('RP','Primary','RPC')) {
    field<-c(RP='RPCluster',Primary='FinalSubtype',RPC='RPCCluster')[[axis]]
    ids<-switch(axis,RP=rp_order,Primary=primary_order,RPC=rpc_order)
    count<-d[,.N,by=c(field)];setnames(count,field,'ClusterID');count[,ClusterID:=as.character(ClusterID)]
    count[,Order:=match(ClusterID,as.character(ids))];setorder(count,Order)
    top<-85
    for(k in seq_len(nrow(count))) {
      height<-count$N[k]*unit_height;bottom<-top-height
      nodes[[length(nodes)+1L]]<-data.table(Axis=axis,ClusterID=count$ClusterID[k],Order=k,N=count$N[k],
        X=c(RP=23,Primary=87,RPC=151)[[axis]],Width=if(axis=='Primary')18 else 6,
        Top=top,Bottom=bottom)
      top<-bottom-4
    }
    temp<-copy(count);setnames(temp,'ClusterID','Cluster')
    temp[,Configuration:=axis];tables[[length(tables)+1L]]<-temp
  }
  nodes<-rbindlist(nodes)
  flow<-copy(paths[N>0]);flow[,`:=`(PrimaryOrder=match(FinalSubtype,primary_order),
    RPOrder=match(RPCluster,rp_order),RPCOrder=match(RPCCluster,rpc_order))]
  for(axis in c('RP','Primary','RPC')) {
    field<-c(RP='RPCluster',Primary='FinalSubtype',RPC='RPCCluster')[[axis]]
    ordercols<-switch(axis,RP=c('PrimaryOrder','RPCOrder','PathID'),
      Primary=c('RPOrder','RPCOrder','PathID'),RPC=c('PrimaryOrder','RPOrder','PathID'))
    for(k in seq_len(nrow(nodes[Axis==axis]))) {
      node<-nodes[Axis==axis][k];z<-copy(flow[as.character(get(field))==node$ClusterID])
      setorderv(z,ordercols)
      upper<-node$Top-c(0,head(cumsum(z$N)*unit_height,-1));lower<-upper-z$N*unit_height
      stopifnot(abs(tail(lower,1)-node$Bottom)<1e-10)
      idx<-match(z$PathID,flow$PathID)
      set(flow,i=idx,j=paste0(axis,'Top'),value=upper)
      set(flow,i=idx,j=paste0(axis,'Bottom'),value=lower)
    }
  }
  # Adjacent margins are descriptive checks; drawing uses the full triple table.
  countrows<-list()
  for(config in c('RP','RPC')) {
    field<-paste0(config,'Cluster');ids<-if(config=='RP')rp_order else rpc_order
    tab<-CJ(FinalSubtype=primary_order,Cluster=as.character(ids))
    actual<-d[,.N,by=c('FinalSubtype',field)];setnames(actual,field,'Cluster');actual[,Cluster:=as.character(Cluster)]
    tab<-merge(tab,actual,by=c('FinalSubtype','Cluster'),all.x=TRUE);tab[is.na(N),N:=0L]
    tab[,`:=`(Configuration=config,ClusterOrder=match(Cluster,as.character(ids)))]
    tab[,SubtypeN:=sum(N),by=FinalSubtype];tab[,RowPercent:=100*N/SubtypeN]
    countrows[[length(countrows)+1L]]<-tab
  }
  counts<-rbindlist(countrows)
  source_paths<-unique(c(inputs$Path,paste0(F4DUAL,c('run_manifest.json','input_hashes.csv','output_hashes.csv',
    'main_joined_assignments.csv','full_fit_assignments.csv','alluvial_paths.csv','agreement.csv',
    'factor_scores.csv','variance_explained_per_factor.csv','training_statistics.csv','mofa_dual_model.hdf5')),
    '8-Visualization/scripts/render_fig4_alluvial.R'))
  write<-function(x,name)fwrite(x,file.path(STEP8_ROOT,F4ALL,'source_data/',paste0(name,'.csv')),na='NA')
  write(d,'main_joined_assignments');write(paths,'all_possible_paths');write(flow,'path_layout');write(nodes,'node_layout')
  write(counts,'correspondence_tables');write(agreement,'agreement')
  write(data.table(WidthMM=174,HeightMM=116,GraphTopMM=85,GraphBottomMM=21,NodeGapMM=4,MMPerSample=unit_height),'layout_spec')
  write(data.table(Path=source_paths,MD5=unname(tools::md5sum(file.path(STEP8_ROOT,source_paths)))),'input_hashes')
}
f4all_ribbon <- function(x0,x1,top0,top1,bottom0,bottom1,fill) {
  t<-seq(0,1,length.out=70);dx<-x1-x0
  x<-(1-t)^3*x0+3*(1-t)^2*t*(x0+.45*dx)+3*(1-t)*t^2*(x1-.45*dx)+t^3*x1
  curve<-function(y0,y1)(1-t)^3*y0+3*(1-t)^2*t*y0+3*(1-t)*t^2*y1+t^3*y1
  grid::grid.polygon(grid::unit(c(x,rev(x)),'mm'),
    grid::unit(c(curve(top0,top1),rev(curve(bottom0,bottom1))),'mm'),
    gp=grid::gpar(fill=fill,col=NA))
}
build_fig4_alluvial <- function(id) {
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  stats<-f4all_read('agreement')
  if(id=='G') {
    paths<-f4all_read('path_layout');nodes<-f4all_read('node_layout')
    return(function() {
      grid::grid.newpage()
      fig4_v2_text('MOFA clusters across two and three omics',4,111,9.5,TRUE,just='left')
      fig4_v2_text('Same 182 Main tumors  |  Flow colors: primary RNA subtype',4,105,7.2,col='grey35',just='left')
      fig4_v2_text('RNA + protein\nMOFA factor k-means',23,95,7.5,TRUE)
      fig4_v2_text('Primary RNA\nC1 / C2 / C3',87,95,7.5,TRUE)
      fig4_v2_text('RNA + protein + CNV\nMOFA factor k-means',151,95,7.5,TRUE)
      z<-copy(paths);setorder(z,-N,PathID)
      for(i in seq_len(nrow(z))) {
        row<-z[i];fill<-grDevices::adjustcolor(pal[[row$FinalSubtype]],alpha.f=.46)
        f4all_ribbon(26,78,row$RPTop,row$PrimaryTop,row$RPBottom,row$PrimaryBottom,fill)
        f4all_ribbon(96,148,row$PrimaryTop,row$RPCTop,row$PrimaryBottom,row$RPCBottom,fill)
      }
      for(i in seq_len(nrow(nodes))) {
        row<-nodes[i];yy<-(row$Top+row$Bottom)/2
        fill<-if(row$Axis=='Primary')grDevices::adjustcolor(pal[[row$ClusterID]],alpha.f=.24) else '#E7E7E7'
        border<-if(row$Axis=='Primary')pal[[row$ClusterID]] else 'grey50'
        fig4_v2_rect(row$X,yy,row$Width,row$Top-row$Bottom,fill,col=border,lwd=.45)
        label<-paste0(if(row$Axis=='Primary')row$ClusterID else paste0(row$Axis,row$ClusterID),'\nn = ',row$N)
        fig4_v2_text(label,if(row$Axis=='RP')19 else if(row$Axis=='RPC')155 else row$X,yy,7.4,
          bold=row$Axis=='Primary',just=if(row$Axis=='RP')'right' else if(row$Axis=='RPC')'left' else 'centre')
      }
      fig4_v2_text(sprintf('ARI vs primary: %.3f',stats$ARI[1]),23,12.5,7.4,TRUE)
      fig4_v2_text(sprintf('ARI vs primary: %.3f',stats$ARI[2]),151,12.5,7.4,TRUE)
      fig4_v2_text('Both fits: same 204 cases; paths show actual three-label assignments of the 182 Main cases.',87,4,6.7,col='grey35')
    })
  }
  if(id=='S10S') {
    d<-f4all_read('correspondence_tables')
    d[,`:=`(ClusterLabel=paste0(Configuration,Cluster),FinalSubtype=factor(FinalSubtype,levels=c('C3','C2','C1')),
      Label=sprintf('%d\n%.1f%%',N,RowPercent))]
    levels<-d[order(Configuration,ClusterOrder),unique(ClusterLabel)]
    d[,ClusterLabel:=factor(ClusterLabel,levels=levels)]
    d[,Configuration:=factor(Configuration,levels=c('RP','RPC'),labels=c(
      sprintf('RNA + protein | ARI = %.3f',stats$ARI[1]),sprintf('RNA + protein + CNV | ARI = %.3f',stats$ARI[2])))]
    return(ggplot(d,aes(ClusterLabel,FinalSubtype,fill=RowPercent))+
      geom_tile(color='white',linewidth=.6)+geom_text(aes(label=Label),size=2.8,family=STEP8_FONT$family,lineheight=1.15)+
      facet_wrap(~Configuration,nrow=1,scales='free_x')+
      scale_fill_gradient(low='#F5F7FA',high='#AAC2D7',limits=c(0,100),guide='none')+
      scale_y_discrete(labels=c('C3 (n = 48)','C2 (n = 92)','C1 (n = 42)'),expand=expansion(add=.03))+
      scale_x_discrete(expand=expansion(add=.03))+
      labs(title='Exact correspondence by modality set',x='MOFA factor k-means cluster (original ID)',y='Primary RNA subtype',
        caption='182 shared Main cases. Cells show n and within-subtype percentage.')+
      theme_step8()+theme(plot.title=element_text(size=9,face='bold'),axis.title=element_text(size=7.5),
        axis.text=element_text(size=7.1),axis.line=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=7.4),
        plot.caption=element_text(size=6.7,hjust=.5,color='grey35')))
  }
  stop('Unknown alluvial panel: ',id)
}
register_fig4_alluvial <- function() {
  m<-read_source('8-Visualization/tables/figure_manifest.csv');base<-copy(m[Figure=='Figure 4'&Panel=='E'][1])
  sources<-unique(c(f4all_read('input_hashes')$Path,
    paste0(F4ALL,'source_data/',c('input_hashes','main_joined_assignments','all_possible_paths','path_layout',
      'node_layout','correspondence_tables','agreement','layout_spec'),'.csv'),
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R',
    '8-Visualization/scripts/render_fig4_genomics_v2.R','8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV))
  for(id in c('G','S10S')) {
    row<-copy(base);row[,`:=`(Figure=if(id=='G')'Figure 4' else 'Supplement',Panel=id,SourceFigure='',
      Role=if(id=='G')'Same-sample dual MOFA / primary RNA / triple MOFA alluvial' else 'Exact two-view and three-view MOFA correspondence counts',
      SourceTable=paste(sources,collapse=';'),
      TargetOutput=if(id=='G')'8-Visualization/figures/main/Figure4G_mofa_alluvial.pdf' else '8-Visualization/figures/supplement/Supplement_S10S_mofa_correspondence.pdf',
      Status='AnalysisReady',WidthMM='174',HeightMM=if(id=='G')'116' else '90',SourceMD5='',QAFile='',RenderedAt='',
      SubtypeScope='Same204 fitted Dong cases and same182 Main displayed cases (42/92/48); exact locked RNA/Protein inputs',
      Notes='New additive RNA+Protein MOFA removes only CNV from locked triple-view inputs; frozen triple model/labels unchanged. One fixed seed20260715,10 initial factors,Gaussian/ARD/spike-slab settings; derived kmeans3/n_init100. True sample-level triple paths,no marginal-table splicing; flows colored by central primary subtype and original cluster IDs retained. Single-fit correspondence,not independent validation or optimal-k/multi-seed evidence.')]
    idx<-which(m$Figure==row$Figure&m$Panel==id)
    if(length(idx))m[idx]<-row else m<-rbindlist(list(m,row),use.names=TRUE)
  }
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
render_fig4_alluvial <- function(id) {
  figure<-if(id=='G')'Figure 4' else 'Supplement';row<-preflight_subfigure(figure,id)
  save_panel(build_fig4_alluvial(id),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure==figure&Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile=paste0(F4ALL,'QA.md'),RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
