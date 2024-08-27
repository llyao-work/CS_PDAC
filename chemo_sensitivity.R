library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
load("cl_exp.Rdata")
exp2 <- apply(exp, 2, as.numeric)
rownames(exp2) <- rownames(exp)
exp <- as.data.frame(exp2)
exp <- log2(exp+1)
testExpr <- exp
testExpr[1:4,1:4]  
colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)