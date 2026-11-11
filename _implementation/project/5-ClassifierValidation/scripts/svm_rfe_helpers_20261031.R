# Pure linear-SVM coefficient recovery and training-only elimination.
svm_pair_weights <- function(model){
  stopifnot(model$kernel==0,model$nclasses==3)
  starts<-c(1L,1L+cumsum(model$nSV));res<-list()
  for(i in 1:2)for(j in (i+1):3){
    ii<-starts[i]:(starts[i+1]-1L);jj<-starts[j]:(starts[j+1]-1L)
    w<-colSums(model$SV[ii,,drop=FALSE]*model$coefs[ii,j-1])+colSums(model$SV[jj,,drop=FALSE]*model$coefs[jj,i])
    res[[paste(model$levels[model$labels[c(i,j)]],collapse='/')]]<-w
  }
  do.call(cbind,res)
}
rfe_path <- function(train_z,y,cost,counts){
  stopifnot(ncol(train_z)>=max(counts),min(counts)>=2,all(is.finite(train_z)))
  genes<-colnames(train_z);sets<-list();drops<-list();errors<-numeric()
  repeat{
    n<-length(genes)
    fit<-e1071::svm(x=train_z[,genes,drop=FALSE],y=y,kernel='linear',cost=cost,scale=FALSE,probability=FALSE)
    w<-svm_pair_weights(fit);score<-rowMeans(w^2)
    if(n %in% counts)sets[[as.character(n)]]<-genes
    if(n==ncol(train_z)||n==min(counts)){
      zz<-train_z[seq_len(min(7,nrow(train_z))),genes,drop=FALSE]
      ref<-attr(predict(fit,zz,decision.values=TRUE),'decision.values')
      calc<-sweep(zz%*%w,2,fit$rho,'-');stopifnot(identical(colnames(ref),colnames(calc)))
      errors<-c(errors,max(abs(ref-calc)))
    }
    if(n==min(counts))break
    gone<-names(score)[order(score,names(score))][1]
    drops[[length(drops)+1]]<-data.table(NBefore=n,Dropped=gone,Score=unname(score[gone]))
    genes<-genes[genes!=gone]
  }
  stopifnot(setequal(names(sets),as.character(counts)),max(errors)<1e-8)
  list(sets=sets,drops=rbindlist(drops),max_decision_error=max(errors),initial=colnames(train_z),cost=cost)
}
