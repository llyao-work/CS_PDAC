library(randomForest)
library(tidyverse)
library(glmnet)
library(xgboost)
library(caret)
library(ggplot2)
load("gene_intersection.Rdata")
load("exp_cl.Rdata")
exp_lastgene <- t(exp_unicox[last_gene,]) 
exp_lastgene <- as.data.frame(exp_lastgene)
dt <- cbind(cl[,c(2,3)],exp_lastgene)
dt$status <- as.numeric(dt$status)
train_cox <- coxph(as.formula(paste0("Surv(OS,status) ~ ",paste(last_gene,collapse = " + "))),data =dt)
coef <- coef(train_cox)
myFun=function(x){crossprod(as.numeric(x),coef)}
trainScore=apply(exp_lastgene,1,myFun)
risk <- as.vector(ifelse(trainScore>median(trainScore),"high","low"))
exp <- cbind(dt,score=as.vector(trainScore),risk)
fit <- survfit(Surv(OS, status) ~ risk, data = exp)
ggsurvplot(
  fit,
  data = exp,
  censor = T, 
  censor.shape = "|", censor.size = 4,
  conf.int = T, 
  conf.int.style = "ribbon",
  conf.int.alpha = 0.3,
  pval = TRUE, 
  pval.size = 5,
  legend = "top", 
  legend.title = 'Risk Score',
  legend.labs = c("High Risk","Low Risk"),
  xlab = "Years",
  ylab = "Survival probablity",
  title = "Cohort",
  palette = c('#ed0000','#00468b'),
  ggtheme = theme_classic(), 
  risk.table = TRUE, 
  risk.table.col = "strata", 
  risk.table.title = 'Number at risk',
  fontsize = 4,
  risk.table.y.text = FALSE, 
  risk.table.height = 0.2,
)

dt_ROC <- exp[,c(1,2,length(colnames(exp))-1)]

timeROC <- with(dt_ROC, 
                     timeROC(T = OS,
                             delta = status,
                             marker = score, 
                             cause = 1, 
                             times = c(1,2,3), 
                             ROC = TRUE, 
                             iid = TRUE, 
                             weighting = "marginal")) 
print(timeROC)

###timeROC

pdf(file="AUC.pdf",width=5,height=5)
plot(timeROC,time=1,col="#92C5DE",title=FALSE,lwd=2)
plot(timeROC,time=2,col='#F4A582',add=TRUE,title=FALSE,lwd=2)
plot(timeROC,time=3,col='#66C2A5',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",GSE79668_timeROC$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",GSE79668_timeROC$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",GSE79668_timeROC$AUC[3]))),
       col=c("#92C5DE", "#F4A582", "#66C2A5"),lwd=2,bty = 'n')
dev.off()
