#!/usr/bin/env Rscript

root <- normalizePath(".")
no_dir <- file.path(root, "external_validation_no_overlap")
px_dir <- file.path(root, "external_validation_proxy_no_overlap")
out_dir <- file.path(root, "outputs", "validation_graphs")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_names <- c(
  pcp_aac="PCP + AAC",
  pcp_aac_ctd="PCP + AAC + CTD",
  pcp_aac_ctd_dpc="PCP + AAC + CTD + DPC"
)
model_names <- c(
  svm_rbf="SVM-RBF", svm_linear="SVM-linear",
  logistic_regression="Logistic regression", random_forest="Random forest",
  knn="kNN", xgboost="XGBoost"
)
dataset_order <- unname(dataset_names)
model_order <- unname(model_names)
dataset_cols <- c("#4C78A8", "#F58518", "#54A24B")
names(dataset_cols) <- dataset_order

normalize_dataset <- function(x) {
  x <- sub("^s1_proxy_no_overlap_", "", x)
  unname(dataset_names[x])
}

read_tab <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)

shared_names <- intersect(
  read_tab(file.path(no_dir, "combined_external_validation_set.tsv"))$name,
  read_tab(file.path(px_dir, "combined_external_validation_set.tsv"))$name
)

build_summary <- function(folder) {
  internal <- read_tab(file.path(folder, "internal_metrics_summary_ci.tsv"))
  pred <- read_tab(file.path(folder, "combined_external_predictions_by_seed.tsv"))
  pred <- pred[pred$name %in% shared_names, ]
  seed_means <- aggregate(predicted_label ~ dataset + model + seed, pred, mean)
  ext_mean <- aggregate(predicted_label ~ dataset + model, seed_means, mean)
  names(ext_mean)[3] <- "external_detection"
  result <- merge(internal, ext_mean, by=c("dataset", "model"))
  result$dataset_label <- normalize_dataset(result$dataset)
  result$model_label <- unname(model_names[result$model])
  result$order <- match(result$dataset_label, dataset_order) * 10 + match(result$model_label, model_order)
  result[order(result$order), ]
}

no <- build_summary(no_dir)
px <- build_summary(px_dir)

open_png <- function(name, width, height) {
  png(file.path(out_dir, paste0(name, ".png")), width=width, height=height,
      units="in", res=300, bg="white")
}
open_pdf <- function(name, width, height) {
  pdf(file.path(out_dir, paste0(name, ".pdf")), width=width, height=height,
      bg="white", useDingbats=FALSE)
}

draw_paired <- function() {
  par(mfrow=c(1,2), mar=c(4.5, ifelse(par("mfg")[2] == 1, 17, 2), 4.5, 1),
      oma=c(2,0,4.5,0), family="sans")
  for (idx in seq_along(list(no, px))) {
    dat <- list(no, px)[[idx]]
    par(mar=c(4.5, if (idx == 1) 17 else 2, 4.5, 1))
    y <- rev(seq_len(nrow(dat)))
    plot(NA, xlim=c(0,100), ylim=c(0.5,nrow(dat)+0.5), axes=FALSE,
         xlab="Positive detection / recall (%)", ylab="")
    abline(v=seq(0,100,20), col="#D9E0E3", lwd=1)
    for (i in seq_len(nrow(dat))) {
      col <- dataset_cols[dat$dataset_label[i]]
      segments(dat$mean_test_recall[i]*100, y[i], dat$external_detection[i]*100, y[i],
               col=adjustcolor(col, alpha.f=0.6), lwd=2.2)
      points(dat$mean_test_recall[i]*100, y[i], pch=21, bg="#263238", col="white", cex=1.15)
      points(dat$external_detection[i]*100, y[i], pch=21, bg=col, col="white", cex=1.35)
    }
    axis(1, at=seq(0,100,20), labels=paste0(seq(0,100,20), "%"), col=NA, col.axis="#43545B")
    if (idx == 1) {
      labs <- paste(dat$dataset_label, dat$model_label, sep="  |  ")
      axis(2, at=y, labels=labs, las=1, tick=FALSE, cex.axis=0.72, line=-0.5)
    }
    box(bty="l", col="#AAB7BC")
    title(c("No-proxy training", "Proxy-integrated training")[idx], font.main=2, cex.main=1.25)
  }
  mtext("Internal recall versus external positive detection", outer=TRUE, side=3, line=2.5, cex=1.65, font=2)
  mtext("Both panels use the same 11 external peptides; means across 10 seeds",
        outer=TRUE, side=3, line=0.8, cex=0.9, col="#52636B")
  legend("bottom", inset=-0.18, xpd=NA, horiz=TRUE, bty="n",
         legend=c("Internal test recall", "External detection"),
         pt.bg=c("#263238", "#4C78A8"), pch=21, pt.cex=1.3, cex=0.9)
}

draw_scatter <- function() {
  par(mfrow=c(1,2), mar=c(4.5,4.8,4,1), oma=c(3.5,0,4.5,0), family="sans")
  pchs <- c(21,22,23,24,25,8); names(pchs) <- model_order
  for (idx in seq_along(list(no, px))) {
    dat <- list(no, px)[[idx]]
    plot(dat$mean_test_mcc, dat$external_detection*100, xlim=c(0.4,0.9), ylim=c(0,90),
         xlab="Internal test MCC", ylab=if(idx==1) "External positive detection (%)" else "",
         pch=pchs[dat$model_label], bg=dataset_cols[dat$dataset_label],
         col=dataset_cols[dat$dataset_label], cex=1.45, lwd=1.5)
    grid(col="#D9E0E3")
    points(dat$mean_test_mcc, dat$external_detection*100,
           pch=pchs[dat$model_label], bg=dataset_cols[dat$dataset_label],
           col="white", cex=1.45, lwd=1.2)
    title(c("No-proxy training", "Proxy-integrated training")[idx], font.main=2, cex.main=1.25)
  }
  mtext("Internal discrimination versus external detection", outer=TRUE, side=3, line=2.5, cex=1.65, font=2)
  mtext("Upper-right indicates stronger balance; both panels use the same 11 external peptides",
        outer=TRUE, side=3, line=0.8, cex=0.9, col="#52636B")
  legend("bottom", inset=-0.26, xpd=NA, horiz=TRUE, bty="n", ncol=3,
         legend=dataset_order, pch=21, pt.bg=dataset_cols, col=dataset_cols, cex=0.85)
  legend("bottom", inset=-0.39, xpd=NA, horiz=TRUE, bty="n", ncol=6,
         legend=model_order, pch=pchs, col="#59666B", pt.bg="white", cex=0.75)
}

draw_proxy_effect <- function() {
  dat <- merge(no[,c("dataset_label","model_label","external_detection")],
               px[,c("dataset_label","model_label","external_detection")],
               by=c("dataset_label","model_label"), suffixes=c("_no","_px"))
  dat$order <- match(dat$dataset_label, dataset_order)*10 + match(dat$model_label, model_order)
  dat <- dat[order(dat$order),]
  y <- rev(seq_len(nrow(dat)))
  par(mar=c(5,17,5,1), family="sans")
  plot(NA, xlim=c(0,100), ylim=c(0.5,nrow(dat)+0.5), axes=FALSE,
       xlab="External positive detection (%)", ylab="")
  abline(v=seq(0,100,20), col="#D9E0E3")
  for (i in seq_len(nrow(dat))) {
    a <- dat$external_detection_no[i]*100; b <- dat$external_detection_px[i]*100
    col <- if (b >= a) "#2A9D8F" else "#D95F59"
    segments(a,y[i],b,y[i],col=adjustcolor(col,alpha.f=0.75),lwd=2.8)
    points(a,y[i],pch=21,bg="#6C7A80",col="white",cex=1.15)
    points(b,y[i],pch=21,bg=col,col="white",cex=1.35)
  }
  axis(1,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),col=NA,col.axis="#43545B")
  axis(2,at=y,labels=paste(dat$dataset_label,dat$model_label,sep="  |  "),las=1,tick=FALSE,cex.axis=0.74,line=-0.5)
  box(bty="l",col="#AAB7BC")
  title("Effect of proxy-integrated training on external detection",font.main=2,cex.main=1.5,line=2)
  mtext("Same 11 external peptides in both conditions; means across 10 seeds",side=3,line=0.7,cex=0.9,col="#52636B")
  legend("topright",inset=c(0.01,0.01),xpd=FALSE,horiz=FALSE,bty="n",
         legend=c("No proxy","Proxy improved","Proxy decreased"),pch=21,
         pt.bg=c("#6C7A80","#2A9D8F","#D95F59"),col="white",pt.cex=1.3,cex=0.9)
}

for (device in c("png","pdf")) {
  if (device == "png") open_png("01_internal_recall_vs_external_detection",15,10.5) else open_pdf("01_internal_recall_vs_external_detection",15,10.5)
  draw_paired(); dev.off()
  if (device == "png") open_png("02_internal_mcc_vs_external_detection",14,7.2) else open_pdf("02_internal_mcc_vs_external_detection",14,7.2)
  draw_scatter(); dev.off()
  if (device == "png") open_png("03_proxy_effect_on_external_detection",11,10.5) else open_pdf("03_proxy_effect_on_external_detection",11,10.5)
  draw_proxy_effect(); dev.off()
}

writeLines(c(
  paste("Shared external peptides:", length(shared_names)),
  "External detection for both training conditions was calculated on the same shared peptide set.",
  "Internal metrics and external detection are means across 10 seeds."
), file.path(out_dir, "README.txt"))

cat("Created validation graphs in", out_dir, "\n")
