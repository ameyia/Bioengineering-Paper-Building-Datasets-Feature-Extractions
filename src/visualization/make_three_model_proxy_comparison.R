#!/usr/bin/env Rscript
root<-normalizePath(".");out<-file.path(root,"results","comparisons","validation_graphs");dir.create(out,recursive=TRUE,showWarnings=FALSE)
nd<-file.path(root,"results","baseline","external_validation");pd<-file.path(root,"results","proxy","external_validation")
rd<-function(p)read.delim(p,stringsAsFactors=FALSE,check.names=FALSE);cl<-function(x)sub("^s1_proxy_no_overlap_","",x)
shared<-intersect(rd(file.path(nd,"combined_external_validation_set.tsv"))$name,rd(file.path(pd,"combined_external_validation_set.tsv"))$name)
mods<-c("svm_rbf","knn","logistic_regression");mlabs<-c("SVM-RBF","kNN","Logistic regression")
dsets<-c("pcp_aac","pcp_aac_ctd","pcp_aac_ctd_dpc");dlabs<-c("PCP + AAC","+ CTD","+ CTD + DPC")
loadc<-function(folder){
 i<-rd(file.path(folder,"internal_metrics_by_seed.tsv"));i$dataset<-cl(i$dataset);i<-i[i$model%in%mods,]
 p<-rd(file.path(folder,"combined_external_predictions_by_seed.tsv"));p<-p[p$name%in%shared & p$model%in%mods,];p$dataset<-cl(p$dataset)
 e<-aggregate(predicted_label~dataset+model+seed,p,mean);names(e)[4]<-"external"
 merge(i,e,by=c("dataset","model","seed"))
}
means<-function(x)aggregate(cbind(test_recall,external)~dataset+model,x,mean)
n<-means(loadc(nd));p<-means(loadc(pd))
cols<-c("Internal - no proxy"="#204D73","External - no proxy"="#72A6C8","Internal - with proxy"="#B95618","External - with proxy"="#F2A65A")
draw<-function(){
 par(mfrow=c(1,3),family="sans",oma=c(5,4.5,5.5,1.2),mar=c(5,ifelse(par("mfg")[2]==1,4.2,1.2),4,0.8))
 for(mi in 1:3){
  m<-mods[mi];a<-n[n$model==m,];b<-p[p$model==m,];a<-a[match(dsets,a$dataset),];b<-b[match(dsets,b$dataset),]
  vals<-rbind(a$test_recall,a$external,b$test_recall,b$external)*100
  par(mar=c(5,if(mi==1)4.2 else 1.2,4,.8))
  bp<-barplot(vals,beside=TRUE,col=cols,border=NA,ylim=c(0,105),names.arg=dlabs,axes=FALSE,cex.names=.9)
  abline(h=seq(0,100,20),col="#E0E6E9",lwd=.8)
  if(mi==1)axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),las=1,col=NA,cex.axis=.85)
  for(i in 1:4)for(j in 1:3)text(bp[i,j],vals[i,j]+3,paste0(round(vals[i,j]),"%"),cex=.65,font=2,col="#263238")
  box(bty="l",col="#9EADB3");title(mlabs[mi],font.main=2,cex.main=1.2)
  if(mi==1)mtext("Recall / positive detection",side=2,line=3,cex=.85)
 }
 mtext("Internal and external validation with and without proxy peptides",outer=TRUE,side=3,line=3.3,cex=1.55,font=2)
 mtext("Means across 10 seeds | identical 11-peptide external set",outer=TRUE,side=3,line=1.5,cex=.88,col="#52636B")
 legend("bottom",inset=-.12,xpd=NA,horiz=TRUE,bty="n",legend=names(cols),fill=cols,border=NA,cex=.82)
}
png(file.path(out,"07_three_model_internal_external_proxy_comparison.png"),width=13.333,height=7.5,units="in",res=300,bg="white");draw();dev.off()
pdf(file.path(out,"07_three_model_internal_external_proxy_comparison.pdf"),width=13.333,height=7.5,bg="white",useDingbats=FALSE);draw();dev.off()
cat("Created three-model comparison\n")
