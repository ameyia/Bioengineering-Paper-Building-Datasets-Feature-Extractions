#!/usr/bin/env Rscript
root<-normalizePath(".");out<-file.path(root,"results","comparisons","validation_graphs");dir.create(out,recursive=TRUE,showWarnings=FALSE)
rows<-data.frame(
 feature=c("PCP + AAC","PCP + AAC + CTD","PCP + AAC + CTD + DPC"),
 model=c("XGBoost","SVM-RBF","SVM-RBF"),
 frequency=c("5 of 10 seeds","7 of 10 seeds","6 of 10 seeds"),
 mcc=c("0.724","0.805","0.845"),stringsAsFactors=FALSE)

draw<-function(){
 par(mar=c(0,0,0,0),family="sans");plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 text(.05,.91,"Best Models Used for SHAP Analysis",adj=c(0,.5),cex=2,font=2,col="#173A57")
 text(.05,.84,"Most frequently selected model by validation MCC across 10 seeds",adj=c(0,.5),cex=1.05,col="#52636B")
 left<-.05;right<-.95;top<-.75;rh<-.115
 widths<-c(.34,.24,.22,.20);xs<-c(left,left+cumsum(widths)*(right-left));headers<-c("Feature dataset","Most selected model","Selection frequency","Mean test MCC*")
 rect(left,top-rh,right,top,col="#173A57",border=NA)
 for(j in 1:4)text(xs[j]+.012,top-rh/2,headers[j],adj=c(0,.5),cex=.95,font=2,col="white")
 for(i in 1:3){
  y1<-top-rh*i;y0<-y1-rh;fill<-if(i==3)"#DDF1EE" else if(i%%2)"#F4F7F8" else "white"
  rect(left,y0,right,y1,col=fill,border="#D5DEE2",lwd=1)
  vals<-unname(unlist(rows[i,]))
  for(j in 1:4)text(xs[j]+.012,(y0+y1)/2,vals[j],adj=c(0,.5),cex=1.03,font=if(i==3)2 else 1,col="#18364D")
 }
 ycall<-top-rh*4-.07
 rect(left,ycall-.05,right,ycall+.05,col="#EAF6F4",border="#2A9D8F",lwd=1.5)
 text(left+.018,ycall,"Overall SHAP choice: PCP + AAC + CTD + DPC with SVM-RBF",adj=c(0,.5),cex=1.05,font=2,col="#155E59")
 text(left,0.075,"*Mean held-out test MCC of the validation-selected model in each seed; SHAP was calculated after model selection.",adj=c(0,.5),cex=.70,col="#59686E")
 text(left,0.035,"Other selected models: PCP + AAC - SVM-RBF (3), kNN (1), RF (1); CTD - LogReg (2), XGBoost (1); DPC - LogReg (2), RF (2).",adj=c(0,.5),cex=.62,col="#59686E")
}
png(file.path(out,"08_best_models_for_shap_table.png"),width=13.333,height=6.2,units="in",res=300,bg="white");draw();dev.off()
pdf(file.path(out,"08_best_models_for_shap_table.pdf"),width=13.333,height=6.2,bg="white",useDingbats=FALSE);draw();dev.off()
cat("Created SHAP model table\n")
