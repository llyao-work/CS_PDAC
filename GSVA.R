library(tidyverse)
library(clusterProfiler)
library(msigdbr) 
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(dplyr)
library(BiocParallel)
library(org.Hs.eg.db)
load("count_exp_nonlog.Rdata")
load("exp.Rdata")
load("cl.Rdata")
dat <- exp
meta <- data.frame(sample=colnames(dat),group=cl$risk)
GO_df_all <- msigdbr(species = "Homo sapiens", 
                     category = "C5")
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
go_list <- split(GO_df$gene_symbol, GO_df$gs_name)
data <- as.matrix(dat)
gsva_mat <- gsva(expr=data, 
                 gset.idx.list=go_list, 
                 kcdf="Gaussian", 
                 verbose=T, 
                 parallel.sz = parallel::detectCores())
gsva_mat <- t(gsva_mat)
gsva_mat <- as.data.frame(gsva_mat)
go <- t(go) %>% as.data.frame()
identical(rownames(exp),rownames(go))
GO_merge <- cbind(go,exp$score)
cor_list <- list()
for (i in colnames(GO_merge)) {
  tar <- as.numeric(GO_merge$`exp$score`)
  cor_res <- cor(x=tar,y=as.numeric(GO_merge[,i]),method = "spearman")
  cor_pal <- cor.test(x=tar,y=as.numeric(GO_merge[,i]))$p.value
  final_res <- data.frame(tar_name="score",
                          gene_name=i,
                          cor_res=cor_res,
                          cor_p=cor_pal)
  cor_list[[i]] <- final_res
}
gene_corres <- do.call("rbind",cor_list)
gene_corres <- gene_corres[gene_corres$cor_p<0.05,]
GO_relation <- gene_corres
GO_matrix <- go[c(rownames(GO_relation)),] %>% t()
merge_matrix <- cbind(exp$score,GO_matrix)

library(tidyverse)
library(Hmisc)
merge_matrix <- merge_matrix[,-21]
colnames(merge_matrix) <- substring(colnames(merge_matrix),6)
colnames(merge_matrix) <- c("Score","REGULATION_OF_CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN", "MITOTIC_G2_DNA_DAMAGE_CHECKPOINT_SIGNALING", "INTERLEUKIN_2_MEDIATED_SIGNALING_PATHWAY", "REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY", "MITOTIC_CELL_CYCLE_CHECKPOINT_SIGNALING", "NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION", "REGULATORY_T_CELL_DIFFERENTIATION", "INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY", "INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY", "T_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION", "NK_T_CELL_PROLIFERATION", "T_CELL_TOLERANCE_INDUCTION" , "POSITIVE_REGULATION_OF_MACROPHAGE_MIGRATION", "MACROPHAGE_PROLIFERATION", "REGULATION_OF_MACROPHAGE_CHEMOTAXIS", "PROXIMAL_TUBULE_BICARBONATE_RECLAMATION", "PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM", "TGF_BETA_SIGNALING_PATHWAY", "PRIMARY_IMMUNODEFICIENCY","DNA_REPLICATION")
cor_res <- rcorr(merge_matrix,type = "spearman")
cor_mat <- cor_res$r
p_mat <- cor_res$P
cor_mat_long <- cor_mat %>% 
  as.data.frame() %>% 
  mutate(x = factor(rownames(cor_mat), levels = rownames(cor_mat))) %>% 
  pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>% 
  mutate(y = factor(y, levels = rownames(cor_mat)))
head(cor_mat_long)
ggplot() +
  geom_point(data = cor_mat_long, aes(x, y, size = cor, fill = cor), 
             shape = 21, color = "#c4bcba")+
  scale_fill_gradient2(low = "#72bcd5", high = "#e76254")+
  scale_x_discrete(expand = c(0.025, 0.025))+
  scale_y_discrete(expand = c(0.025, 0.025))+
  geom_vline(aes(xintercept =seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggsave("p1.pdf", height = 8, width = 8)


cor_mat_up <- cor_mat
cor_mat_up[lower.tri(cor_mat_up, diag = T)] <- NA

cor_mat_up_long <- cor_mat_up %>% 
  as.data.frame() %>% 
  mutate(x = factor(rownames(cor_mat_up), levels = rownames(cor_mat_up))) %>% 
  pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>% 
  mutate(y = factor(y, levels = rownames(cor_mat_up)))

ggplot() +
  geom_point(data = cor_mat_up_long, aes(x, y, size = cor, fill = cor), 
             shape = 21, color = "#c4bcba")+
  scale_fill_gradient2(low = "#4178a7", mid = "#ffffff", high = "#e2201c")+
  scale_x_discrete(expand = c(0.025, 0.025))+
  scale_y_discrete(expand = c(0.025, 0.025))+
  geom_vline(aes(xintercept =seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggsave("p2.pdf", height = 8, width = 8)

p_mat_lower <- p_mat
p_mat_lower[!lower.tri(p_mat_lower, diag = F)] <- 1

p_mat_lower_long <- p_mat_lower %>% 
  as.data.frame() %>% 
  mutate(x = factor(rownames(p_mat_lower), levels = rownames(p_mat_lower))) %>% 
  pivot_longer(cols = !x, names_to = "y", values_to = "p") %>% 
  mutate(y = factor(y, levels = rownames(p_mat_lower)))


p_mat_lower_long$p_sig <- as.character(symnum(p_mat_lower_long$p, 
                                              cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                              symbols = c("***", "**", "*", "")))
ggplot() +
  geom_point(data = p_mat_lower_long, 
             aes(x, y, color = -log10(p)),
             size = 10,
             shape = 15)+
  geom_text(data = p_mat_lower_long, 
            aes(x, y, label = p_sig))+
  scale_color_gradientn(colours = c("#ffffff", "#72bcd5", "#e3eeef", "#e76254"))+
  scale_x_discrete(expand = c(0.025, 0.025))+
  scale_y_discrete(expand = c(0.025, 0.025))+
  geom_vline(aes(xintercept =seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggsave("p3.pdf", height = 8, width = 8.5)

ggplot() +
  geom_point(data = p_mat_lower_long, 
             aes(x, y, color = -log10(p)),
             size = 10,
             shape = 15)+
  geom_text(data = p_mat_lower_long, 
            aes(x, y, label = p_sig))+
  scale_color_gradientn(colours = c("#ffffff", "#72bcd5", "#e3eeef", "#e76254"))+
  geom_point(data = cor_mat_up_long, aes(x, y, size = cor, fill = cor), 
             shape = 21, color = "#c4bcba")+
  scale_fill_gradient2(low = "#4178a7", mid = "#ffffff", high = "#e2201c")+
  scale_x_discrete(expand = c(0.025, 0.025))+
  scale_y_discrete(expand = c(0.025, 0.025), position = "right")+
  geom_vline(aes(xintercept =seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90),
        legend.position = "top")+
  guides(fill = guide_colorbar("correlation", title.position = "top",
                               title.theme = element_text(face = "bold")),
         color = guide_colorbar("-log10(P value)", title.position = "top",
                                title.theme = element_text(face = "bold")),
         size = "none")
ggsave("GSVA.pdf", height = 13.5, width = 13.5)