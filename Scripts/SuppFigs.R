##============================-- 1. Set up Environment ----
## 1.1 load packages ----
rm(list=ls())

pkgs <- c("Seurat", "tidyverse", "pals", "immunarch",
          "ggpubr", "org.Hs.eg.db", "ggsci", "clusterProfiler", 
          "pheatmap", "devtools", "R.utils", "ggh4x", "ggalluvial", 
          "ggpubr", "ggprism")

lapply(pkgs, require, character.only = TRUE)

## 1.2 load auxiliary script ----
source("Scripts/Utility/RNA_Utility.R")

## 1.4 set system colors ----
colset8 <- readRDS("RDS/color_set_CD8.rds")
colset4 <- readRDS("RDS/color_set_CD4.rds")

## 1.5 load CAR4 and CAR8 seurat objects
CAR8 <- readRDS("RDS/CAR8_for_manuscript_Apr262023.rds")
CAR4 <- readRDS("RDS/CAR4_for_manuscript_Apr262023.rds")


##============================-- Sup2 Percentage of CAR-T cells: CD4 vs CD8 ----
require(gtools)
CD48Col <- c("#F2AE30", "#024959")
names(CD48Col) <- c("CAR4", "CAR8")

dfx <- read.csv("Data/WBC_Percent.csv") %>%
  mutate(CAR_Count = CAR_abundance * WBC_Count * 1000) %>%
  mutate(CAR_perc = CAR_abundance*100) %>%
  rename(time=date) %>%
  mutate(patientShort = plyr::mapvalues(patient,
                                        from = c("P031", "P033", "P034", "P039", "P048", "P049", "P116"),
                                        to = c("P04", "P01", "P03", "P05", "P06", "P02", "P07")))
  

selectMeta <- c("annoCoarse", "stageFine","time", "patientShort")
meta4 <- FetchData(CAR4, vars = selectMeta) %>% mutate(annoCoarse = paste0("CAR4_", annoCoarse))
meta8 <- FetchData(CAR8, vars = selectMeta) %>% mutate(annoCoarse = paste0("CAR8_", annoCoarse))
metaAll <- bind_rows(meta4, meta8) %>% 
  mutate(celltype = str_split(annoCoarse, "_", simplify = T)[,1]) %>%
  group_by(patientShort, time, stageFine) %>%
  count(celltype) %>%
  mutate(PropRaw = n/sum(n)) %>%
  rename(Freq = n) %>%
  ungroup() %>%
  left_join(dfx, by = c("time", "patientShort")) %>%
  drop_na() %>%
  select(-patient) %>%
  rename(patient = patientShort)

metaD0 <- data.frame(patient = rep(paste0("P0", 1:7), each = 2), 
                     time = "d00", celltype = rep(c("CAR4", "CAR8"), 7)) %>%
  mutate(stageFine="Exp", Freq=0, PropRaw=0, T_percent=0, CAR_abundance=0, 
         CAR_perc=0, WBC_Count=0, CAR_Count=0)

metaAllFinal <- bind_rows(metaD0, metaAll) %>%
  mutate(date = as.numeric(gsub("d", "", time))) %>%
  mutate(patient = factor(patient, levels = mixedsort(unique(patient)))) %>%
  mutate(celltype = factor(celltype, levels = c("CAR4", "CAR8"))) %>%
  mutate(Prop = PropRaw * CAR_perc)


### panel2b: area plot ----
rects <- data.frame(xstart = c(0,14), xend = c(14,30), stage=c("Exp", "Per"))
pts <- paste0("P0", 1:7)
rects_final <- bind_rows(lapply(pts, function(px){
  dfx <- rects
  dfx$patient <- px
  dfx
})) %>%
  mutate(patient=factor(patient, levels = mixedsort(unique(patient)))) 

ps2b <- ggplot() +
  geom_rect(data = rects_final, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                                    fill = stage), alpha = 0.4) +
  geom_area(data=metaAllFinal, aes(x=date, y=Prop, fill=celltype),
            alpha=0.6, size=1, color="black") +
  geom_vline(xintercept = 14, linetype="dashed") +
  facet_wrap(~patient, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = c(0,7, 14, 21, 28), 
                     labels = paste0("d", c(0,7, 14, 21, 28)),
                     limits = c(0,31)) +
  scale_fill_manual(values=c(CD48Col, colset8$stageCoarse[c(2:3)]))+
  labs(x="Timepoint", y="%CAR") +
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.spacing.x = unit(0.5,"cm")) +
  guides(fill=guide_legend(nrow=1))
saveFig(ps2b, prefix = "Plots/FigS2b_CAR_Percent_Dynamics_CD4_vs_CD8", 
        wd = 2400, ht = 2400)

### panel2a: bar plot ----
metaAllFinal <- metaAll %>%
  mutate(date = as.numeric(gsub("d", "", time))) %>%
  mutate(patient = factor(patient, levels = mixedsort(unique(patient)))) %>%
  mutate(celltype = factor(celltype, levels = c("CAR4", "CAR8"))) %>%
  mutate(Prop = PropRaw)

ps2a <- ggplot(metaAllFinal, aes(x=time, y=Prop, fill=celltype)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~patient, scales = "free", ncol = 2) +
  scale_fill_manual(values=c(CD48Col))+
  labs(x="Timepoint", y="Proportion of CAR-T") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0.5),
        strip.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.position = "top",
        legend.spacing.x = unit(0.5,"cm")) +
  guides(fill=guide_legend(nrow=1))
saveFig(ps2a, prefix = "Plots/FigS2a_CAR_Percent_Barplot_CD4_vs_CD8", 
        wd = 1500, ht = 2400)

##============================-- Sup3 Clonality Further Details ----
### panel3a: overlap between CAR and Endo ----
repAll <- readRDS("RDS/AggRepList_b123_PI_AllT_Only_CR_wMeta_Celltype_stageFine.rds")
overlapAll <- repOverlap(repAll$data, .method = "morisita")
ps3a <- vis(overlapAll)
saveFig(ps3a, 
        prefix = "Plots/FigS3a_Overlap_Coefficients_TileMap_Endo_and_CAR", 
        wd = 2400, ht = 2000)

### panel3b: overlap between CAR and Endo: barplot ----
repAll <- readRDS("RDS/AggRepList_b123_PI_AllT_Only_CR_wMeta_Patient_Celltype_stageFine.rds")
overlapAll <- repOverlap(repAll$data, .method = "morisita")

dfovAll <- overlapAll %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = !Sample1, names_to = "Sample2", values_to = "Overlap_Coeff") %>%
  filter(Overlap_Coeff >= 0) %>%
  mutate(s1Pt = str_split(Sample1, "-", simplify = T)[,1]) %>%
  mutate(s2Pt = str_split(Sample2, "-", simplify = T)[,1]) %>%
  mutate(s1Type = str_split(Sample1, "-", simplify = T)[,2]) %>%
  mutate(s2Type = str_split(Sample2, "-", simplify = T)[,2]) %>%
  mutate(s1Stage = str_split(Sample1, "-", simplify = T)[,3]) %>%
  mutate(s2Stage = str_split(Sample2, "-", simplify = T)[,3]) %>%
  filter(s1Pt == s2Pt & s1Type != s2Type) %>%
  mutate(s1Group = paste0(s1Type, "-", s1Stage)) %>%
  mutate(s2Group = paste0(s2Type, "-", s2Stage)) %>%
  mutate(group = unlist(Map(function(x,y){
    paste0(sort(c(x,y)), collapse = "_")
  }, s1Group, s2Group))) %>%
  select(s1Pt, group, Overlap_Coeff) %>%
  distinct()

ps3b <- ggplot(dfovAll, aes(x=group, y=Overlap_Coeff, fill=group)) +
  stat_summary(geom = "col", fun = mean, width = 0.7) +
  stat_summary(geom = "errorbar",
               fun = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.3) +
  theme_prism() +
  scale_fill_npg() +
  labs(x="", y="Clonotype Overlap Coefficient") +
  theme(
    strip.text = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=12, angle=60, vjust=1, hjust=1),
    axis.text = element_text(size=10),
    legend.text = element_text(size=10)
  )
saveFig(ps3b, 
        prefix = "Plots/FigS3b_Overlap_Coefficients_barPlot_Endo_and_CAR", 
        wd = 2400, ht = 2000)

##============================-- Sup6 Clonality of each CD8 and CD4 Clusters ----
### S6a: CD8 Clonality
repAnno8 <- readRDS("RDS/AggRepList_b123_PI_CAR8_Only_CR_wMeta_Patient_AnnoCoarse.rds")
imm_rare <- repClonality(repAnno8$data, .method = "rare")
ps6a <- vis(imm_rare) +
  labs(x="", title = "CAR8 Clonotype Count Distribution", subtitle = "") +
  theme_prism() +
  theme(plot.subtitle = element_blank(),
        axis.text.x = element_text(size=12, angle=60, vjust=1, hjust=1)) +
  scale_fill_manual(values = rev(c("#D94436", "#733F3F", "#CC8D1A", "#0897B4", "#B7BF99","grey90")))
saveFig(ps6a, "Plots/FigS6a_CAR8_TCR_clonotype_distribution_bySize", wd = 2400, ht = 2000)

### S6b: CD4 Clonality
repAnno4 <- readRDS("RDS/AggRepList_b123_PI_CAR4_Only_CR_wMeta_Patient_AnnoCoarse.rds")
imm_rare <- repClonality(repAnno4$data, .method = "rare")
ps6b <- vis(imm_rare) +
  labs(x="", title = "CAR4 Clonotype Count Distribution", subtitle = "") +
  theme_prism() +
  theme(plot.subtitle = element_blank(),
        axis.text.x = element_text(size=12, angle=60, vjust=1, hjust=1)) +
  scale_fill_manual(values = rev(c("#D94436", "#733F3F", "#CC8D1A", "#0897B4", "#B7BF99","grey90")))
saveFig(ps6b, "Plots/FigS6b_CAR4_TCR_clonotype_distribution_bySize", wd = 2400, ht = 2000)





##============================-- Sup7 Stage Comparison Further Details ----
### panel 7a: volcano plot of stageDEG
mStageSC <- readRDS("RDS/SingleCell_DEG_S1_vs_S2_Raw_May112023.rds")
require(ggrepel)
volcano_col <- c(unname(colset8$stageCoarse[2:3]), 'grey80')
Exp_genes <- c("CD69", "CD44", "HLA-DRA", "HLA-DRB1", "CD74", "HIF1A", "PDE4D",
               "PRDM1", "RGS1", "CTLA4", "NR4A2", "NFKBIA", "NFKB1", "NFKBIE", 
               "DUSP4", "BCL3", "DUSP5", "CD83")
Per_genes <- c("KLRB1", "KLRD1", "FCGR3A", "PRF1", "GZMB", "KIR3DL1", "ZEB2",
               "CX3CR1", "S1PR1", "S1PR5", "SMAD9", "CISH",
               "GZMM", "NKG7", "IRF7", "STAT1", "IFITM1", "IFITM2", "IRF2")
gene_highlight <- c(Exp_genes, Per_genes)

## preprocessing
res <- mStageSC %>%
  rownames_to_column(var = "gene") %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(logp = ifelse(logp>300, 300, logp)) %>%
  mutate(avg_log2FC = clipMtx(avg_log2FC, 1.2)) %>%
  mutate(category = case_when(
    logp > 3 & avg_log2FC > 0.25 ~ "Up in Exp",
    logp > 3 & avg_log2FC < -0.25 ~ "Up in Per",
    TRUE ~ "Non-DEGs"
  )) %>%
  mutate(gene_highlight = ifelse(gene %in% gene_highlight, gene, "")) %>%
  mutate(gene_highlight = ifelse(category=="Non-DEGs", "", gene_highlight)) %>%
  mutate(category = factor(category, levels = c("Up in Exp", "Up in Per", "Non-DEGs")))

res_Highlight <- subset(res, gene_highlight!="") 

ps7a<- ggplot(res, aes(x=avg_log2FC, y=logp, color=category, label=gene_highlight))+
  geom_point(size=1) +
  geom_point(data = res_Highlight, size=2, color="black") +
  geom_text_repel(color="black", min.segment.length = 0, max.overlaps = Inf, size=3, force = 1.5) +
  geom_vline(xintercept = 0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_vline(xintercept = -0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_hline(yintercept = 3, linetype="dashed", color="black", size=0.5) +
  theme_pubr()+
  labs(x=bquote(~Log[2]~ 'Fold Change'), y=bquote('-'~Log[10]~P['adj']))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_blank())+
  scale_color_manual(values = volcano_col)+
  scale_y_continuous(expand = expansion(add=30))+
  guides(color = guide_legend(override.aes = list(size=6))) 
saveFig(ps7a, 
        prefix = "Plots/FigS7a_volcanoPlot_DEG_Exp_vs_Per_CAR8", 
        wd = 3000, ht = 2000)

### panel 7c: Xgboost Accuracy
xgbx <- xgblist$model
bpx_pred <- as.data.frame(predict(xgbx, newdata = xgblist$scoreDf, type = 'prob'))
row.names(bpx_pred) <- row.names(xgblist$scoreDf)
metaDf <- FetchData(CAR8, vars = c("stageCoarse", "patientShort", "time"))
metaDf$time <- factor(metaDf$time, levels = gtools::mixedsort(unique(metaDf$time)))
dfx <- bind_cols(metaDf, bpx_pred) %>%
  mutate(Pred = ifelse(Exp>0.5, "Exp", "Per")) %>%
  group_by(patientShort, time) %>%
  count(Pred) %>%
  mutate(Freq = n/sum(n)) %>%
  ungroup() %>%
  rename(patient = patientShort) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient)))) %>%
  mutate(Pred= factor(Pred, levels = c("Exp", "Per")))

ps7c <- ggplot(dfx, aes(x=time, y=Freq, fill=Pred)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~patient, scales = "free", ncol = 2) +
  scale_fill_manual(values=colset8$stageCoarse[c(2:3)])+
  labs(x="Timepoint", y="Proportion of CAR-T") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0.5),
        strip.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.position = "top",
        legend.spacing.x = unit(0.5,"cm")) +
  guides(fill=guide_legend(nrow=1))
saveFig(ps7c, prefix = "Plots/FigS7c_CAR_Stage_Prediction_Accuracies", 
        wd = 1500, ht = 2400)



##============================-- Sup9 IP pre group Further Details ----
## panel a: UMAP showing the precursor distribution ----
IP8$preGroup <- ifelse(IP8$inExp=="Yes", "Exp_Pre", 
                       ifelse(IP8$inExp=="No" & IP8$inPer=="Yes", "Per_Pre", "Others"))
IP8$isNA <- ifelse(is.na(IP8$preGroup),"NA", "nonNA")
subIP8 <- subset(IP8, isNA=="nonNA")
subIP8$preGroup <- factor(subIP8$preGroup, levels = c("Exp_Pre", "Per_Pre", "Others"))
gpCols <- c(colset8$stageCoarse[c(2:3)],"grey90")
names(gpCols) <-  c("Exp_Pre", "Per_Pre", "Others")


ps9a_p1 <- DimPlot(subIP8, reduction = "umap", group.by = "preGroup",
                  pt.size = 1.5, raster = F, 
                  order = c("Exp_Pre", "Per_Pre", "Others"))+
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = gpCols) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(ps9a_p1, prefix = "Plots/FigS9a_p1_UMAP_IP8_Precursor_Group_combined", 
        wd = 3600, ht = 3200)

ps9a_p2 <- DimPlot(subIP8, reduction = "umap", group.by = "preGroup", split.by = "preGroup",
                  pt.size = 1.5, raster = F, 
                  order = c("Others","Per_Pre", "Exp_Pre"))+
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=10,color = 'black'),
        axis.title = element_text(size=15,face = 'bold'),
        plot.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=15),
        strip.background = element_rect(fill = "transparent", colour = "transparent"))+
  scale_color_manual(values = c(unname(colset8$stageCoarse[c(2:3)]),"grey90"))
saveFig(ps9a_p2, prefix = "Plots/FigS9a_p2_UMAP_IP8_Precursor_Group_split", 
        wd = 3600, ht = 1600)


## panel b: CAR expression in each cluster ----
ps9b <- VlnPlot(IP8, "YESCAR", group.by = "annoCoarse", pt.size = 0.5) +
  labs(x="", y="Normalized Expression", title="CAR Transgene") +
  theme_pubr() +
  theme(axis.text.x = element_text(size=15,color = 'black'),
        axis.title = element_text(size=15,face = 'bold'),
        axis.text.y = element_text(size=10),
        plot.title = element_text(size=15, hjust=0.5),
        legend.position = "none")+
  scale_fill_manual(values = colset8$annoIP8)
saveFig(ps9b, prefix = "Plots/FigS9b_VlnPlot_IP8_CAR_RNA_Expr", 
        wd = 2000, ht = 1200)




##============================-- Sup10 scATAC ----
## load data
CAR8 <- readRDS("RDS/CAR8ATAC_for_manuscript_May142023.rds")

## panel a: UMAPs
### annoCoarse; WNNUMAP
ps10a_p1 <- DimPlot(CAR8, reduction = "wnn.umap", group.by = "annoCoarse",
                    pt.size = 1.5, raster = F) +
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(ps10a_p1, prefix = "Plots/FigS10a_p1_UMAP_CAR8ATAC_by_annoCoarse_WNNUMAP", 
        wd = 3600, ht = 3200)

### annoCoarse; RNAUMAP
ps10a_p2 <- DimPlot(CAR8, reduction = "umap.rna", group.by = "annoCoarse",
                    pt.size = 1.5, raster = F) +
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(ps10a_p2, prefix = "Plots/FigS10a_p2_UMAP_CAR8ATAC_by_annoCoarse_RNAUMAP", 
        wd = 3600, ht = 3200)

### annoCoarse; ATACUMAP
ps10a_p3 <- DimPlot(CAR8, reduction = "umap.atac", group.by = "annoCoarse",
                    pt.size = 1.5, raster = F) +
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(ps10a_p3, prefix = "Plots/FigS10a_p3_UMAP_CAR8ATAC_by_annoCoarse_ATACUMAP", 
        wd = 3600, ht = 3200)

### stageCoarse; WNNUMAP
ps10a_p4 <- DimPlot(CAR8, reduction = "wnn.umap", group.by = "stageFine",
                    pt.size = 1.5, raster = F) +
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$stageFine[c(2:4)]) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(ps10a_p4, prefix = "Plots/FigS10a_p4_UMAP_CAR8ATAC_by_stageFine_WNNUMAP", 
        wd = 3600, ht = 3200)


## panel b: Density Plots
require(Nebulosa)

keyMarkers <- unique(c("KLRG1", "TBX21", 
                "MKI67","IL7R", "CCR7", 
                "TCF7", "GZMB", "PRF1",
                "TOX", "TIGIT", "NR4A2", "PDCD1","CTLA4",
                "HAVCR2", "LAG3", "ENTPD1", "IRF4", "CD27", 
                "B3GAT1", "IL7R", "LEF1", 
                "CX3CR1", "NKG7", "GNLY", "IFNG"))
for (gx in keyMarkers){
  ps10b <- plot_density(CAR8, gx, reduction = "wnn.umap")
  saveFig(ps10b, prefix = paste0("Plots/Density_Plots/FigS10b_CAR8ATAC_RNA_density_markers_", gx),
          wd = 1800, ht = 1500)
}


## panel c: peak plots
library(Signac)
library(colorspace)
dfPeaks <- read.csv("Data/Ex_associated_pk_list_in_Human.csv")
DefaultAssay(CAR8) <- "ATAC"
NewCovPlot <- function(sobj=CAR8, rx, gx, gpVar="stageCoarse",
                       gpCols=darken(colset8$stageCoarse[c(2:3)], amount = 0)){
  p1x <- CoveragePlot(sobj, annotation = F, peaks = F,
                      group.by = gpVar, region = rx) +
    scale_fill_manual(values = gpCols)
  
  p2x <- AnnotationPlot(sobj, region = rx)
  # p2g <- ggplot_build(p2)
  # p2g[["data"]][[5]]["size"] <- 3
  # p2n <- ggplot_gtable(p2g)
  
  p3x <- ExpressionPlot(sobj, features = gx, assay = "RNA", group.by = gpVar) +
    scale_fill_manual(values = gpCols)
  
  pA <- CombineTracks(
    plotlist = list(p1x, p2x),
    expression.plot = p3x,
    heights = c(10, 1),
    widths = c(10,1))
  
  return(pA)
}


CoveragePlot(object = CAR8, region = "ENTPD1", annotation = T, peaks = F, group.by = "stageCoarse", 
             region.highlight = "chr10-95800000-95810000")

NewCovPlot(rx = "chr10-95700000-95880000", gx = "ENTPD1")

rList <- dfPeaks$peak
gList <- dfPeaks$gene

### part 1: by stageCoarse
for (ii in 1:length(rList)){
  rr <- rList[ii]
  gg <- gList[ii]
  p1 <- NewCovPlot(rx=rr, gx = gg, gpVar = "stageCoarse", gpCols = colset8$stageCoarse[c(2:3)])
  saveFig(p1, prefix = paste0("Plots/CovPlot/FigS10c_p1_CAR8ATAC_CovPlot_byStage_", gg),
          wd = 2400, ht = 1800)
}


### part 2: by annoCoarse
for (ii in 1:length(rList)){
  rr <- rList[ii]
  gg <- gList[ii]
  p1 <- NewCovPlot(rx=rr, gx = gg, gpVar = "annoCoarse", gpCols = colset8$annoCoarse)
  saveFig(p1, prefix = paste0("Plots/CovPlot/FigS10c_p2_CAR8ATAC_CovPlot_annoCoarse_", gg),
          wd = 2400, ht = 1800)
}


### part 3: get the peaks from CD39 paper
require(GenomicRanges)
require(IRanges)
JWList <- readRDS("RDS/Ref_CAR8ATAC_JWHumanAtlasATAC_exPeaks.rds")
dpList1 <- readRDS("RDS/temp_CAR8ATAC_EMexh_DEP_May16.rds")
dpList2 <- readRDS("RDS/temp_CAR8ATAC_CM_DEP_May16.rds")
dpList <- c(dpList1, dpList2)

xx <- gsub("EMEX", "EM_Exh-like", names(dpList))
names(dpList) <- xx


dpDf <- do.call("rbind", lapply(names(dpList), function(labx){
  dfx <- dpList[[labx]]
  upLab <- str_split(labx, "_vs_")[[1]][1]
  dnLab <- str_split(labx, "_vs_")[[1]][2]
  dfxSig <- dfx %>% 
    mutate(upIn = ifelse(avg_log2FC>0, upLab, dnLab)) %>%
    mutate(upIn = factor(upIn, levels = unique(upIn))) %>%
    filter(p_val_adj<0.05) %>%
    count(upIn, .drop = F) %>%
    mutate(group = labx) %>%
    mutate(n = ifelse(upIn==upLab, n, -n))
  dfxSig
}))

ps10cp3 <- ggplot(dpDf, aes(x=group, y=n, fill = upIn)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = colset8$annoCoarse) +
  labs(y="Number of DAR", x="", fill="") +
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=15,color = 'black'),
        axis.title = element_text(size=15,face = 'bold'),
        legend.text = element_text(size=12))
saveFig(ps10cp3, prefix = "Plots/FigS10c_p3_CAR8ATAC_N_Diff_Peaks",
        wd = 2400, ht = 1500)
