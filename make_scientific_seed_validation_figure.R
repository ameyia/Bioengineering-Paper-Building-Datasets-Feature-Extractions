#!/usr/bin/env Rscript

root <- normalizePath(".")
out <- file.path(root, "outputs", "validation_graphs")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

no_dir <- file.path(root, "external_validation_no_overlap")
px_dir <- file.path(root, "external_validation_proxy_no_overlap")
datasets <- c("pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc")
dataset_labels <- c("PCP + AAC", "PCP + AAC + CTD", "PCP + AAC + CTD + DPC")
names(dataset_labels) <- datasets
models <- c("svm_rbf", "svm_linear", "logistic_regression", "random_forest", "knn", "xgboost")
model_labels <- c("SVM-RBF", "SVM-linear", "Logistic regression", "Random forest", "kNN", "XGBoost")
names(model_labels) <- models
condition_cols <- c("No proxy"="#35618D", "Proxy-integrated"="#D87522")
dataset_cols <- c("#3B6FB6", "#E07B25", "#3C9141"); names(dataset_cols) <- datasets
model_pch <- c(21,22,23,24,25,8); names(model_pch) <- models

read_tab <- function(p) read.delim(p, stringsAsFactors=FALSE, check.names=FALSE)
clean_dataset <- function(x) sub("^s1_proxy_no_overlap_", "", x)

shared <- intersect(
  read_tab(file.path(no_dir,"combined_external_validation_set.tsv"))$name,
  read_tab(file.path(px_dir,"combined_external_validation_set.tsv"))$name
)

load_condition <- function(folder, condition) {
  internal <- read_tab(file.path(folder,"internal_metrics_by_seed.tsv"))
  internal$dataset <- clean_dataset(internal$dataset)
  pred <- read_tab(file.path(folder,"combined_external_predictions_by_seed.tsv"))
  pred <- pred[pred$name %in% shared,]
  pred$dataset <- clean_dataset(pred$dataset)
  ext <- aggregate(predicted_label ~ dataset + model + seed, pred, mean)
  names(ext)[4] <- "external_detection"
  d <- merge(internal, ext, by=c("dataset","model","seed"))
  d$condition <- condition
  d$gap <- d$external_detection - d$test_recall
  d
}

all <- rbind(load_condition(no_dir,"No proxy"), load_condition(px_dir,"Proxy-integrated"))
all$combo <- paste(unname(dataset_labels[all$dataset]), unname(model_labels[all$model]), sep=" | ")
combo_order <- as.vector(sapply(datasets, function(d) paste(dataset_labels[d], model_labels[models], sep=" | ")))
all$combo_id <- match(all$combo, combo_order)

summarize_metric <- function(metric) {
  a <- aggregate(all[[metric]], list(combo=all$combo, combo_id=all$combo_id,
                                     dataset=all$dataset, model=all$model,
                                     condition=all$condition),
                 function(x)c(mean=mean(x), sd=sd(x), n=length(x)))
  vals <- if (is.matrix(a$x)) a$x else do.call(rbind,a$x)
  a$mean <- vals[,"mean"]; a$ci <- 1.96*vals[,"sd"]/sqrt(vals[,"n"])
  a$x <- NULL
  a
}

draw_interval_panel <- function(metric, xlim, title, subtitle, zero=FALSE, percent=FALSE) {
  s <- summarize_metric(metric)
  n <- length(combo_order); ybase <- rev(seq_len(n))
  plot(NA,xlim=xlim,ylim=c(0.4,n+0.6),axes=FALSE,xlab="",ylab="")
  if (zero) abline(v=0,col="#444444",lwd=1.2,lty=2)
  abline(v=pretty(xlim,6),col="#E2E7EA",lwd=0.8)
  offsets <- c("No proxy"=0.16,"Proxy-integrated"=-0.16)
  set.seed(42)
  for (cond in names(condition_cols)) {
    raw <- all[all$condition==cond,]
    for (i in seq_len(nrow(raw))) {
      val <- raw[[metric]][i]
      if (percent) val <- val*100
      points(val,ybase[raw$combo_id[i]]+offsets[cond]+runif(1,-0.07,0.07),
             pch=16,cex=0.43,col=adjustcolor(condition_cols[cond],alpha.f=0.30))
    }
    ss <- s[s$condition==cond,]
    for (i in seq_len(nrow(ss))) {
      mean <- ss$mean[i]; ci <- ss$ci[i]
      if (percent) {mean<-mean*100;ci<-ci*100}
      yy <- ybase[ss$combo_id[i]]+offsets[cond]
      segments(mean-ci,yy,mean+ci,yy,col=condition_cols[cond],lwd=2.2)
      segments(mean-ci,yy-0.07,mean-ci,yy+0.07,col=condition_cols[cond],lwd=1.2)
      segments(mean+ci,yy-0.07,mean+ci,yy+0.07,col=condition_cols[cond],lwd=1.2)
      points(mean,yy,pch=21,bg=condition_cols[cond],col="white",cex=1.15,lwd=0.8)
    }
  }
  ticks <- pretty(xlim,6)
  axis(1,at=ticks,labels=if(percent) paste0(ticks,"%") else format(ticks,nsmall=1),
       col=NA,col.axis="#43545B",cex.axis=0.8)
  axis(2,at=ybase,labels=combo_order,las=1,tick=FALSE,cex.axis=0.55,line=-0.4)
  box(bty="l",col="#AAB7BC")
  title(title,adj=0,font.main=2,cex.main=1.08,line=1.7)
  mtext(subtitle,side=3,adj=0,line=0.25,cex=0.68,col="#5B6A70")
}

draw_change_panel <- function() {
  keys <- c("dataset","model","seed")
  no <- all[all$condition=="No proxy",c(keys,"test_mcc","external_detection")]
  px <- all[all$condition=="Proxy-integrated",c(keys,"test_mcc","external_detection")]
  d <- merge(no,px,by=keys,suffixes=c("_no","_px"))
  d$dmcc <- d$test_mcc_px-d$test_mcc_no
  d$dext <- (d$external_detection_px-d$external_detection_no)*100
  s <- aggregate(cbind(dmcc,dext)~dataset+model,d,mean)
  plot(s$dmcc,s$dext,xlim=c(-0.18,0.12),ylim=c(-70,20),axes=FALSE,xlab="",ylab="",
       pch=model_pch[s$model],bg=dataset_cols[s$dataset],col="white",cex=1.45,lwd=1)
  abline(h=0,v=0,col="#555555",lty=2,lwd=1.1); grid(col="#E2E7EA")
  points(s$dmcc,s$dext,pch=model_pch[s$model],bg=dataset_cols[s$dataset],
         col="white",cex=1.45,lwd=1)
  axis(1,at=seq(-0.15,0.1,0.05),labels=sprintf("%+.2f",seq(-0.15,0.1,0.05)),col=NA,cex.axis=0.8)
  axis(2,at=seq(-60,20,20),labels=sprintf("%+d pp",seq(-60,20,20)),las=1,col=NA,cex.axis=0.8)
  box(bty="l",col="#AAB7BC")
  title("D  Proxy-induced change",adj=0,font.main=2,cex.main=1.08,line=1.7)
  mtext("Upper-right = improvement in internal MCC and external detection",side=3,adj=0,line=0.25,cex=0.68,col="#5B6A70")
  mtext("Change in internal test MCC",side=1,line=2.5,cex=0.8)
  mtext("Change in external detection",side=2,line=3.4,cex=0.8)
  legend("bottomleft",legend=dataset_labels[datasets],pch=21,pt.bg=dataset_cols[datasets],
         col="white",bty="n",cex=0.63,pt.cex=1.1)
  legend("topright",legend=model_labels[models],pch=model_pch[models],col="#555555",
         bty="n",cex=0.56,pt.cex=0.9)
}

draw_figure <- function() {
  layout(matrix(1:4,2,2,byrow=TRUE),widths=c(1.18,1),heights=c(1,1))
  par(family="sans",mar=c(3.5,15.7,3.7,1),oma=c(1,1,4.2,1))
  draw_interval_panel("test_mcc",c(0.35,0.95),"A  Internal discrimination (MCC)",
                      "Dots = individual seeds; large points and bars = mean and 95% CI")
  par(mar=c(3.5,15.7,3.7,1))
  draw_interval_panel("external_detection",c(0,100),"B  External positive detection",
                      "Calculated on the same 11 external positives in both conditions",percent=TRUE)
  par(mar=c(3.5,15.7,3.7,1))
  draw_interval_panel("gap",c(-100,40),"C  Generalization gap",
                      "External detection minus internal recall; zero indicates agreement",zero=TRUE,percent=TRUE)
  par(mar=c(4.2,4.7,3.7,1))
  draw_change_panel()
  mtext("Seed-level internal and external validation: effect of proxy-integrated training",
        outer=TRUE,side=3,line=2.5,cex=1.5,font=2)
  mtext("All summaries use 10 seeds; external comparisons use an identical 11-peptide set",
        outer=TRUE,side=3,line=0.9,cex=0.82,col="#52636B")
  legend("bottom",inset=-0.02,xpd=NA,horiz=TRUE,bty="n",legend=names(condition_cols),
         pch=21,pt.bg=condition_cols,col="white",pt.cex=1.25,cex=0.8)
}

png(file.path(out,"04_scientific_seed_level_validation_summary.png"),width=18,height=13,units="in",res=300,bg="white")
draw_figure(); dev.off()
pdf(file.path(out,"04_scientific_seed_level_validation_summary.pdf"),width=18,height=13,bg="white",useDingbats=FALSE)
draw_figure(); dev.off()

cat("Created scientific seed-level validation figure\n")
