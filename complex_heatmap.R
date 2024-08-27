load("comple.RData")
icb_merge <- cbind(im_ssgsea,im_quantiseq,im_xcell,im_estimate)
icb_num <- apply(icb_merge, 2, as.numeric)
rownames(icb_num) <- rownames(icb_merge)
icb_num <- t(icb_num)
icb_scaled <- t(scale(t(icb_num)))
icb_scaled <- t(icb_scaled)
icb_final <- cbind(icb_scaled,out_exp)
ssgsea <- icb_final[,1:14]
quantiseq <- icb_final[,15:25]
xcell <- icb_final[,26:53]
load("PAAD_cl.Rdata")
table(cl$Age)
dt <- cl
dt$Estimate <- icb_final$ESTIMATEScore_estimate
dt$ImmuneScore <- icb_final$ImmuneScore_estimate
dt$StromalScore <- icb_final$StromalScore_estimate
dt$Tumorpurity <- icb_final$TumorPurity_estimate
dt$Age <- as.numeric(dt$Age)
table(dt$Gender)
table(dt$Stage)
dt$Stage <- ifelse(dt$Stage %in% c("stage ia","stage ib"),"Stage I",
                   ifelse(dt$Stage %in%c("stage iia","stage iib"),"Stage II",
                          ifelse(dt$Stage=="stage iii","Stage III","Stage IIII")))
table(dt$T_stage)
dt$T_stage <- ifelse(dt$T_stage %in% c("T1","T2"),"T1/T2","T3/T4")
table(dt$M_stage)
table(dt$N_stage)
dt <- dt[dt$N_stage %in% c("N0","N1","N1b"),]
dt$N_stage <- ifelse(dt$N_stage=="N0","N0","N1")
table(dt$Histologic_grade)
dt <- dt[!(dt$Histologic_grade == "GX"),]
dt$Histologic_grade <- ifelse(dt$Histologic_grade %in% c("G1","G2"),"G1/2","G3/4")
table(dt$num_of_positive_node)

library(ComplexHeatmap)
library(dplyr)

out_exp <- out_exp[rownames(dt),]
dt <- cbind(dt,out_exp[,1:2])
cdata <- dt %>% dplyr::select(OS,everything())
gdata <- ssgsea
gdata <- gdata[rownames(cdata),]
gdata <- t(gdata)
Heatmap(gdata)

gdata1 <- quantiseq
gdata1 <- gdata1[rownames(cdata),]
gdata1 <- t(gdata1)
rownames(gdata1) <- gsub('_quantiseq','  ',rownames(gdata1))
Heatmap(gdata1)

gdata2 <- xcell
gdata2 <- gdata2[rownames(cdata),]
gdata2 <- t(gdata2)
rownames(gdata2) <- gsub('_xcell','   ',rownames(gdata2))
Heatmap(gdata2)

status = as.character(cdata$status) 
gender =  as.character(cdata$Gender)
stage =  as.character(cdata$Stage)
age <- as.numeric(cdata$Age)
group <- as.character(cdata$Group)
grade <- as.character(cdata$Histologic_grade)
T_stage <- as.character(cdata$T_stage)
N_stage <- as.character(cdata$N_stage)
Tumorpurity <- as.numeric(cdata$Tumorpurity)
ImmuneScore <- as.numeric(cdata$ImmuneScore)
Estimate <- as.numeric(cdata$Estimate)
StromalScore <- as.numeric(cdata$StromalScore)

ha = HeatmapAnnotation(group=group,stage = stage, gender = gender, age=age,T_stage=T_stage,N_stage=N_stage,Grade=grade,status=status,Tumorpurity=Tumorpurity,Estimate=Estimate,StromalScore=StromalScore,ImmuneScore=ImmuneScore,col = list(group=c("high"="#df536b","low"="#2297e6"),stage=c("Stage I"="#98d0a5","Stage II"="#969CCC","Stage III"="#efe593","Stage IIII"="#ee6060"),gender=c("female"="#4d8b5f","male"="#4b5c6f"),T_stage=c("T1/T2"="#e2514d","T3/T4"="#8aa0bf"),N_stage=c("N0"="#3E0F3A","N1"="#8C5F2F"),Grade=c("G1/2"="#9FACC4","G3/4"="#28414A"),status=c("0"="#C98EBD","1"="#B4B8B8")),border = T,simple_anno_size = unit(0.4, "cm"))

library(circlize)
f = colorRamp2(c(quantile(gdata,.05),median(gdata),quantile(gdata,.95)), c("#8fc9e2", "white", "#eb9184"))
f1 = colorRamp2(c(quantile(gdata,.05),median(gdata1),quantile(gdata1,.95)), c("#8fc9e2", "white", "#eb9184"))
f2 = colorRamp2(c(quantile(gdata,.05),median(gdata2),quantile(gdata2,.95)), c("#8fc9e2", "white", "#eb9184"))

hey=rowAnnotation(foo = anno_text(rownames(gdata), location = 0.5, just = "center",gp = gpar(fill = rep(2:4, each = 4), col = "white", border = "black")))

ht1 = Heatmap(gdata, name = "ssGSEA", row_title = "ssGSEA",top_annotation = ha,column_split = group,row_names_side = "left",show_column_names = F,show_row_names = T,col = f,cluster_columns = T,cluster_rows = T,clustering_method_columns = "complete",column_dend_height = unit(0.2, "cm"), row_dend_width = unit(0.2, "cm"),row_km = 2,border = T)
ht2 = Heatmap(gdata1, name = "Quantiseq", row_title = "Quantiseq",show_column_names = F,show_row_names = T,column_split = group,row_names_side = "left",col = f1,cluster_columns = T,cluster_rows = T,clustering_method_columns = "complete",column_dend_height = unit(0.2, "cm"), row_dend_width = unit(0.2, "cm"),row_km = 2,border = T)
ht3 = Heatmap(gdata2, name = "xcell", row_title = "xcell",show_column_names = F,show_row_names = T,column_split = group,row_names_side = "left",col = f2,cluster_columns = T,cluster_rows = T,clustering_method_columns = "complete",column_dend_height = unit(0.2, "cm"), row_dend_width = unit(0.2, "cm"),row_km = 2,border = T)

ht_list = ht1 %v% ht2 %v% ht3
pdf("heatmap.pdf",width = 15,height = 10)
draw(ht_list)
dev.off()