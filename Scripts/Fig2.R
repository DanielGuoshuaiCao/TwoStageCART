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

## 1.5 load CAR4 and CAR8 seurat objects
CAR8 <- readRDS("RDS/CAR8_for_manuscript_Apr262023.rds")

##============================-- Panel 2a CAR8 UMAP ----
## main annotations
p2a_p1 <- DimPlot(CAR8, reduction = "umap", group.by = "annoCoarse",
               pt.size = 1.5, raster = F)+
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p2a_p1, prefix = "Plots/Fig2a_p1_UMAP_CAR8_annoCoarse", 
        wd = 3600, ht = 3200)

## stage annotations
p2a_p2 <- DimPlot(CAR8, reduction = "umap", group.by = "stageFine",
                  pt.size = 1.2, raster = F)+
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20))+
  scale_color_manual(values = colset8$stageFine[2:4]) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p2a_p2, prefix = "Plots/Fig2a_p2_UMAP_CAR8_stageCoarse", 
        wd = 3600, ht = 3200)


##============================-- Panel 2b CAR8 Markers VlnPlot ----
Idents(CAR8) <- factor(CAR8$annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                   "EM_Eff-like", "CM", "Prolif"))
## RNA markers
keyFeatures <- c("CD8A", "CD8B", "YESCAR",
                 "ZEB2", "TBX21", "CX3CR1","KLRG1",
                 "GZMB", "PRF1", "NKG7", 
                 "CXCR3", "CD38", "GZMK",
                 "TOX", "NR4A2", "PDCD1", "LAG3", "TIGIT",
                 "HAVCR2", "TCF7", "SELL", "MKI67", "TOP2A")

px <- VlnPlot(CAR8, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset8$annoCoarse)
saveFig(px, prefix = "Plots/Fig2b_p1_RNA_Markers_CAR8", 
        wd = 2400, ht = 3600)

## Mixed markers
row.names(CAR8@assays$ADT@data) <- gsub("-TotalSeqC", "", row.names(CAR8@assays$ADT@data))
row.names(CAR8@assays$ADT@scale.data) <- gsub("-TotalSeqC", "", row.names(CAR8@assays$ADT@scale.data))

keyFeatures <- c("CD8A", "ADT_CD8A","YESCAR",
                 "CX3CR1", "TBX21", "KLRB1", "KLRG1","GZMB","GNLY",
                 "CXCR3", "CD38", "ADT_CD45RO", "ADT_CD45RA", 
                 "IL7R", "GZMK", "CXCR6", "ADT_CD127",
                 "PDCD1", "TOX", "NR4A2", "ADT_PD-1",
                 "ADT_TIM-3", "ADT_LAG-3", "ADT_CD39",
                  "ADT_CCR7", "IRF7",
                 "MKI67", "TOP2A")

p2b <- VlnPlot(CAR8, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset8$annoCoarse)
saveFig(p2b, prefix = "Plots/Fig2b_p2_mixed_Markers_CAR8", 
        wd = 3600, ht = 5000)

##============================-- Panel 2c CAR8 phenotype Distribution ----
## overall distribution
metaDf <- FetchData(CAR8, vars = c("annoCoarse", "stageFine")) %>%
  rename(anno=annoCoarse, time=stageFine)

dfx <- metaDf %>%
  group_by(time) %>%
  dplyr::count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = rev(c("Exp", "Per1", "Per2")))) %>%
  mutate(anno = factor(anno, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", "EM_Eff-like", "CM", "Prolif")))

p2c_p1 <- ggplot(dfx, aes(x=time, y=Proportion, fill=anno)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values=colset8$annoCoarse) +
  theme_bw() + 
  labs(x="time") +
  coord_flip() +
  scale_y_continuous(position = "right") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth=1, color = "black")) +
  guides(fill=guide_legend(override.aes = list(size=6), nrow = 2))
saveFig(p2c_p1, 
        prefix = "Plots/Fig2c_p1_annoCoarse_distrib_overall", 
        wd = 2000, ht = 1200)

## per patient distribution
metaDf <- FetchData(CAR8, vars = c("annoCoarse", "stageFine", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageFine, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = rev(c("Exp", "Per1", "Per2")))) %>%
  mutate(anno = factor(anno, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                        "EM_Eff-like", "CM", "Prolif"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))

p2c_p2 <- ggplot(dfx, aes(x=time, y=Proportion, fill=anno)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values=colset8$annoCoarse) +
  theme_bw() + 
  labs(x="time") +
  coord_flip() +
  scale_y_continuous(position = "left") +
  facet_wrap2(~patient, ncol = 2) +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(size=1, color = "black"),
        strip.background = element_rect(fill = "transparent", colour = "transparent"),
        strip.text = element_text(size=12)) +
  guides(fill=guide_legend(override.aes = list(size=6), nrow = 2))
saveFig(p2c_p2, 
        prefix = "Plots/Fig2c_p2_annoCoarse_distrib_perPatient", 
        wd = 2800, ht = 2000)

##============================-- Panel 2d CAR8 phenotype Proportions ----
### 3.4.1 stageFine
metaDf <- FetchData(CAR8, vars = c("annoCoarse", "stageFine", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageFine, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("Exp", "Per1", "Per2"))) %>%
  mutate(anno = factor(anno, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                        "EM_Eff-like", "CM", "Prolif"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))


p2d_p1 <- ggboxplot(dfx, x="time", y="Proportion", add = "jitter", facet.by = "anno", 
          fill = "time", scales = "free_y", nrow = 2) +
  stat_compare_means(comparisons = list(c("Exp", "Per1"), c("Per1", "Per2"), c("Exp", "Per2")), 
                     label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colset8$stageFine[2:4]) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p2d_p1, 
        prefix = "Plots/Fig2d_p1_Pheno_Boxplot_annoCoarse_StageFine", 
        wd = 3000, ht = 2000)
  
### 3.4.2 StageCoarse
metaDf <- FetchData(CAR8, vars = c("annoCoarse", "stageCoarse", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageCoarse, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("Exp", "Per"))) %>%
  mutate(anno = factor(anno, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                        "EM_Eff-like", "CM", "Prolif"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))


p2d_p2 <- ggboxplot(dfx, x="time", y="Proportion", add = "jitter", facet.by = "anno", 
                fill = "time", scales = "free_y", nrow = 2) +
  stat_compare_means(label = "p.signif", comparisons = list(c("Exp", "Per"))) +
  scale_y_continuous(expand = expansion(mult = 0.05, add = 0.1)) +
  scale_fill_manual(values = colset8$stageCoarse[2:3]) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p2d_p2, 
        prefix = "Plots/Fig2d_p2_Pheno_Boxplot_annoCoarse_StageCoarse",  
        wd = 3000, ht = 2000)