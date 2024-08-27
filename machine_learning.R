library(randomForest)
library(tidyverse)
library(glmnet)
library(xgboost)
library(caret)
setwd(workspace)
train_exp <- cbind(cl[,2:3],t(exp))
train_exp$OS <- as.numeric(train_exp$OS)
x1 <- as.matrix(train_exp[,3:ncol(train_exp)])
x2 <- as.matrix(Surv(train_exp$OS, train_exp$status))
lasso_cvfit = cv.glmnet(x1, x2,
                        nfold=10,
                        family = "cox"
)
lasso_fit <- glmnet(x1, x2,
                    family = "cox")
lasso_coef=coef(lasso_fit, s = lasso_cvfit$lambda.min)
lasso_index=which(lasso_coef != 0)
lassoGene=row.names(lasso_coef)[lasso_index]

train_exp$status <- as.factor(train_exp$status)
x = as.matrix(train_exp[,3:ncol(train_exp)]) 
y = train_exp$status
train_forest <- randomForest(x = x,
                             y = y,
                             importance = TRUE,
                             ntree = 500,
                             proximity = TRUE )
optionTrees=which.min(train_forest$err.rate[,1])
optionTrees
rf_forest <- randomForest(x = x,
                          y = y,
                          importance = TRUE, 
                          ntree = optionTrees)
importance=importance(rf_forest)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]

TrainControl <- trainControl( method = "repeatedcv")
model<- train(x=train_exp[,3:ncol(train_exp)],y=as.factor(train_exp$status),  method = "xgbTree", trControl = TrainControl)
plot(varImp(model))
importance <- varImp(model)
head(importance)
important <- as.data.frame(importance$importance)