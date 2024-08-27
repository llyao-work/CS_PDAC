library(IOBR)
load("cl_exp.Rdata")
exp <- read.table("PDA_exp_FPKm.txt",header = T,row.names = 1,sep = "\t")
expr_coad <- log2(exp+1)
expr_coad[1:4,1:4]
tme_deconvolution_methods

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr_coad,
                             method = "mcpcounter"
)

# EPIC
im_epic <- deconvo_tme(eset = expr_coad,
                       method = "epic",
                       arrays = F
)


# xCell
im_xcell <- deconvo_tme(eset = expr_coad,
                        method = "xcell",
                        arrays = F
)

# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr_coad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)


# IPS
im_ips <- deconvo_tme(eset = expr_coad,
                      method = "ips",
                      plot = F
)
## 
## >>> Running Immunophenoscore

# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr_coad,
                            method = "quantiseq",
                            scale_mrna = T
)

# ESTIMATE
im_estimate <- deconvo_tme(eset = expr_coad,
                           method = "estimate"
)

# TIMER
im_timer <- deconvo_tme(eset = expr_coad
                        ,method = "timer"
                        ,group_list = rep("coad",dim(expr_coad)[2])
)


#ssgsea
load(file = "ssGSEA28.Rdata")
im_ssgsea <- calculate_sig_score(eset = expr_coad
                                 , signature = cellMarker 
                                 , method = "ssgsea"
)
