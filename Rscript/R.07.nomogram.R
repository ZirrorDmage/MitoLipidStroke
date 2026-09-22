

print(Sys.time())
rm(list = ls());gc()
####---- 项目根目录定位（脚本位于 <项目根>/Rscript/，脚本内全部使用相对路径）----####
.script_dir <- function() {
  a <- commandArgs(FALSE)
  i <- grep("^--file=", a)
  if (length(i) > 0) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  for (k in seq_len(sys.nframe())) {
    f <- tryCatch(sys.frame(k)$ofile, error = function(e) NULL)
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f, mustWork = FALSE)))
  }
  if (basename(getwd()) == "Rscript") getwd() else file.path(getwd(), "Rscript")
}
SCRIPT_DIR <- .script_dir()
ROOT_CAND <- unique(c(SCRIPT_DIR, dirname(SCRIPT_DIR), getwd(), dirname(getwd())))
ROOT_HIT <- ROOT_CAND[file.exists(file.path(ROOT_CAND, "00.rawdata"))]
root <- if (nzchar(Sys.getenv("PROJECT_ROOT"))) {
  normalizePath(Sys.getenv("PROJECT_ROOT"), mustWork = FALSE)
} else if (length(ROOT_HIT) > 0) {
  normalizePath(ROOT_HIT[1], mustWork = FALSE)
} else if (basename(SCRIPT_DIR) == "Rscript") {
  normalizePath(dirname(SCRIPT_DIR), mustWork = FALSE)
} else {
  normalizePath(SCRIPT_DIR, mustWork = FALSE)
}
options(project.root = root)   # 供脚本中 rm(list = ls()) 之后继续使用
cat("项目根目录：", root, "\n")

wd <- file.path(root,'07.nomogram')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)



rm(list = ls())

library(magrittr)
library(ggplot2)
library(car)
library(rms)
library(lance)
library(pROC)



hubgene <- read.csv('../06.expressionROC//keygenes.csv')[,1]
dat <- read.csv('../00.rawdata/GSE16561/01.data_GSE16561.csv',row.names = 1) 
group <- read.csv('../00.rawdata/GSE16561/02.group_GSE16561.csv',row.names = 1)
hub_exp <- t(dat[hubgene,])%>%as.data.frame()

colnames(group)

hub_exp$group <- group$group

ddist <- datadist(hub_exp)

options(datadist='ddist')

lrm_data_prog <- as.formula(paste0('group~',
                                   paste(hubgene,
                                         sep = "",
                                         collapse = '+')))

lrm <-lrm(lrm_data_prog, data=hub_exp, x=TRUE, y=TRUE,maxit=1000,penalty=0.1)#,penalty=0.1

print(lrm)
hl2 <- stats::resid(lrm,"gof")
pval <- signif(hl2[5], 3)


library(regplot)

nom3<-regplot(lrm,
              observation=hub_exp[28,], ## 
              center=TRUE,
              title="Nomogram",
              point=TRUE, ## point
              odds=FALSE,##
              showP=FALSE,## 
              rank="sd", ## 
              linesize = -4,
              clickable = F) ## 

save.p <- recordPlot()

pdf("01.nomogram.pdf", width = 12, height = 6, family = 'Times')
replayPlot(save.p)
dev.off()

png("01.nomogram.png", width = 12, height = 6, family = 'Times', units = 'in', res = 600)
replayPlot(save.p)
dev.off()





{cal1 <- calibrate(lrm, cmethod='KM', method="boot", B=1000)#B
  
  pdf("02.calibrate.pdf",width=7,height=6)
  par(mar = c(6,5,2,2))
  plot(cal1, lwd=2, lty=1,
       cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2,
       xlim=c(0, 1), ylim= c(-0.1, 1),
       xlab="Nomogram-Predicted Probability of Biomarker",
       ylab="Actual disease (proportion)",
       col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
       legend=FALSE)
  lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
  abline(0, 1, lty=3, lwd=2)
  legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"),
         lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
  if(pval>0.5){
  text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
  text(x = 0.21, y = 0.74, as.expression(bquote(italic('p')==.(pval))))}
  dev.off()
  png("02.calibrate.png",width=7,height=6, units = "in", res = 600)
  par(mar = c(6,5,2,2))
  plot(cal1, lwd=2, lty=1,
       cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2,
       xlim=c(0, 1), ylim= c(-0.1, 1),
       xlab="Nomogram-Predicted Probability of disease risk",
       ylab="Actual disease (proportion)",
       col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
       legend=FALSE)
  lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
  abline(0, 1, lty=3, lwd=2)
  legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"),
         lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
  if(pval>0.5){
    text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
    text(x = 0.21, y = 0.74, as.expression(bquote(italic('p')==.(pval))))}
  dev.off()
}

##ROC plot
library(rms)
library(pROC)
predicted<-predict(lrm,newdata = hub_exp)
lrm_predict<-ifelse(predicted >0.5, "disease", "control") %>% as.factor()

roc_curve = roc(hub_exp$group,predicted)
roc_curve
{png(file='04.roc.png', width = 4, height = 4, res = 600, units = "in")
  par(pin = c(4,4), mar = c(6,6,6,1))
  plot(roc_curve,col="#FF7F00",
       print.auc=TRUE,
       legacy.axes=T,
       asp = 1 ,
       cex.axis=1,
       cex.lab=1.0,
       cex.main=1,
       main='',
       font.lab = 2,
       font.main = 2,
       font.sub =2)
  dev.off()
  
  pdf(file='04.roc.pdf', width = 4, height = 4)
  par(pin = c(4,4), mar = c(6,6,6,1))
  plot(roc_curve,col="#FF7F00",
       print.auc=TRUE,
       legacy.axes=T,
       asp = 1 ,
       cex.axis=1,
       cex.lab=1.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.main=1,   ##标题的缩放倍数
       main='',
       font.lab = 2,
       font.main = 2,
       font.sub =2)
  dev.off()
}


##DCA plot
hub_exp2 <- hub_exp
table(hub_exp2$group)
hub_exp2$group <- ifelse(hub_exp2$group=='Disease',1,0)

library(rmda)
hubgene

no.variates <- hubgene

dca.formula <- sapply(no.variates, function(x) {as.formula(paste("group~", x))})


dca.model <- lapply(dca.formula, function(x) {decision_curve(x, data = hub_exp2,
                                                             family = binomial(link ='logit'),
                                                             thresholds = seq(0,1, by = 0.01),
                                                             confidence.intervals= 0.95,
                                                             study.design = 'case-control',
                                                             population.prevalence= 0.3)})

no.dc <- as.formula(paste("group~", paste0(no.variates, collapse = "+")))
dca.model[["nomogram"]] <- decision_curve(no.dc, data = hub_exp2,
                                          family = binomial(link ='logit'),
                                          thresholds = seq(0,1, by = 0.01),
                                          confidence.intervals= 0.95,
                                          study.design = 'case-control',
                                          population.prevalence= 0.3)


dca.list <- dca.model
saveRDS(dca.model,"dca.model.rds")
dca.list <- dca.model[["nomogram"]]
names(dca.list)


pdf("03.nomogram_dca.pdf", width = 9, height = 9)
plot_decision_curve(dca.list,
                    curve.names = names(dca.list),#
                    cost.benefit.axis = FALSE, col = c( '#039f89',"#0099DD","#FF9933",'#b883d3',"#a65628","#FF4858","#BEBAD8","gray20"),
                    confidence.intervals = FALSE,
                    standardize = FALSE)
dev.off()
png("03.nomogram_dca.png", width = 9, height = 9,res = 600,units = 'in')
plot_decision_curve(dca.list,
                    curve.names = names(dca.list),#
                    cost.benefit.axis = FALSE, col = c( '#039f89',"#0099DD","#FF9933",'#b883d3',"#a65628","#FF4858","#BEBAD8","gray20"),
                    confidence.intervals = FALSE,
                    standardize = FALSE)
dev.off()





