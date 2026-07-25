#!/usr/bin/env Rscript

root <- normalizePath(".")
out <- file.path(root,"results","comparisons","validation_graphs")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
no_dir <- file.path(root,"results","baseline","external_validation")
px_dir <- file.path(root,"results","proxy","external_validation")

datasets <- c("pcp_aac","pcp_aac_ctd","pcp_aac_ctd_dpc")
dataset_labels <- c("PCP + AAC","PCP + AAC + CTD","PCP + AAC + CTD + DPC"); names(dataset_labels)<-datasets
models <- c("svm_rbf","svm_linear","logistic_regression","random_forest","knn","xgboost")
model_labels <- c("SVM-RBF","SVM-linear","LogReg","Random forest","kNN","XGBoost"); names(model_labels)<-models
cols <- c("Internal - no proxy"="#244D70","External - no proxy"="#70A1C3",
          "Internal - proxy"="#B75412","External - proxy"="#F2A65A")
pchs <- c("Internal - no proxy"=16,"External - no proxy"=17,"Internal - proxy"=16,"External - proxy"=17)

read_tab <- function(p) read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
clean <- function(x) sub("^s1_proxy_no_overlap_","",x)
shared <- intersect(read_tab(file.path(no_dir,"combined_external_validation_set.tsv"))$name,
                    read_tab(file.path(px_dir,"combined_external_validation_set.tsv"))$name)

load_condition <- function(folder,proxy=FALSE) {
  int <- read_tab(file.path(folder,"internal_metrics_by_seed.tsv")); int$dataset<-clean(int$dataset)
  pred <- read_tab(file.path(folder,"combined_external_predictions_by_seed.tsv")); pred<-pred[pred$name%in%shared,];pred$dataset<-clean(pred$dataset)
  ext <- aggregate(predicted_label~dataset+model+seed,pred,mean);names(ext)[4]<-"external"
  merge(int,ext,by=c("dataset","model","seed"))
}
no <- load_condition(no_dir); px <- load_condition(px_dir,TRUE)

means <- function(d) aggregate(cbind(test_recall,external,test_mcc)~dataset+model,d,mean)
n <- means(no); p <- means(px)

draw_direct <- function() {
  layout(matrix(1:6,2,3,byrow=TRUE),heights=c(1.35,0.8))
  par(family="sans",oma=c(4.5,4.5,5.2,1.5))
  offsets <- c(-0.24,-0.08,0.08,0.24)
  series <- names(cols)
  for (di in seq_along(datasets)) {
    ds<-datasets[di]; a<-n[n$dataset==ds,];b<-p[p$dataset==ds,]
    a<-a[match(models,a$model),];b<-b[match(models,b$model),]
    par(mar=c(7,if(di==1)4.2 else 1.5,3.8,0.8))
    plot(NA,xlim=c(0.5,6.5),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
    abline(h=seq(0,1,.2),col="#E1E6E9",lwd=.8)
    vals<-list(a$test_recall,a$external,b$test_recall,b$external)
    for(si in 1:4) points(seq_len(6)+offsets[si],vals[[si]],pch=pchs[series[si]],col=cols[series[si]],cex=1.15)
    axis(1,at=1:6,labels=model_labels[models],las=2,tick=FALSE,cex.axis=.73,line=-.3)
    if(di==1) axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),"%"),las=1,col=NA,cex.axis=.8)
    box(bty="l",col="#9EADB3")
    title(dataset_labels[ds],font.main=2,cex.main=1.08)
    if(di==1)mtext("Recall / positive detection",side=2,line=3,cex=.82)
  }
  for(di in seq_along(datasets)) {
    ds<-datasets[di]; a<-n[n$dataset==ds,];b<-p[p$dataset==ds,]
    a<-a[match(models,a$model),];b<-b[match(models,b$model),]
    par(mar=c(7,if(di==1)4.2 else 1.5,3.2,.8))
    plot(NA,xlim=c(.5,6.5),ylim=c(.35,.95),axes=FALSE,xlab="",ylab="")
    abline(h=seq(.4,.9,.1),col="#E1E6E9",lwd=.8)
    for(i in 1:6){segments(i-.08,a$test_mcc[i],i+.08,b$test_mcc[i],col="#89969B",lwd=1.5);points(i-.08,a$test_mcc[i],pch=16,col=cols[1],cex=1.05);points(i+.08,b$test_mcc[i],pch=16,col=cols[3],cex=1.05)}
    axis(1,at=1:6,labels=model_labels[models],las=2,tick=FALSE,cex.axis=.73,line=-.3)
    if(di==1)axis(2,at=seq(.4,.9,.1),labels=sprintf("%.1f",seq(.4,.9,.1)),las=1,col=NA,cex.axis=.8)
    box(bty="l",col="#9EADB3");title("Internal MCC",font.main=2,cex.main=.95)
    if(di==1)mtext("MCC",side=2,line=3,cex=.82)
  }
  mtext("Internal and external validation: no-proxy versus proxy-integrated training",outer=TRUE,side=3,line=3,cex=1.45,font=2)
  mtext("Means across 10 seeds; both external conditions use the same 11 peptides",outer=TRUE,side=3,line=1.4,cex=.83,col="#52636B")
  legend("bottom",inset=-.07,xpd=NA,horiz=TRUE,bty="n",legend=series,col=cols,pch=pchs,cex=.82,pt.cex=1.1)
}

png(file.path(out,"05_direct_internal_external_proxy_comparison.png"),width=15,height=9.5,units="in",res=300,bg="white")
draw_direct();dev.off()
pdf(file.path(out,"05_direct_internal_external_proxy_comparison.pdf"),width=15,height=9.5,bg="white",useDingbats=FALSE)
draw_direct();dev.off()
cat("Created direct comparison figure\n")
