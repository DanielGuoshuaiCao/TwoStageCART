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
colset4 <- readRDS("RDS/color_set_CD4.rds")

## 1.5 load CAR4 and CAR4 seurat objects
CAR4 <- readRDS("RDS/CAR4_for_manuscript_Apr262023.rds")

##============================-- Panel 2a CAR4 UMAP ----
## main annotations
p3a_p1 <- DimPlot(CAR4, reduction = "umap", group.by = "annoCoarse",
               pt.size = 1.5, raster = F)+
  xlab("UMAP1")+ylab("UMAp3")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset4$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p3a_p1, prefix = "Plots/Fig3a_p1_UMAP_CAR4_annoCoarse", 
        wd = 3600, ht = 3200)

## stage annotations
p3a_p2 <- DimPlot(CAR4, reduction = "umap", group.by = "stageFine",
                  pt.size = 1.5, raster = F)+
  xlab("UMAP1")+ylab("UMAp3")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20))+
  scale_color_manual(values = colset4$stageFine[2:4]) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p3a_p2, prefix = "Plots/Fig3a_p2_UMAP_CAR4_stageCoarse", 
        wd = 3600, ht = 3200)


##============================-- Panel 2b CAR4 Markers VlnPlot ----
Idents(CAR4) <- factor(CAR4$annoCoarse, levels = c("Th1","Treg","MEM","Prolif","ISG"))
## RNA markers
keyFeatures1 <- c("CD4", "YESCAR",
                 "TBX21", "CX3CR1","NKG7","CXCR6", 
                 "PRF1", "GNLY",
                 "TCF7", "GZMK", "KLRB1", "SELL", "CD27",
                 "FOXP3", "IL2RA", "IKZF2", 
                 "MKI67", "STMN1", 
                 "ISG15", "MX1", "OAS1", "IRF7")

px <- VlnPlot(CAR4, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset4$annoCoarse)
saveFig(px, prefix = "Plots/Fig3b_p1_RNA_Markers_CAR4", 
        wd = 2400, ht = 3600)

## Mixed markers
row.names(CAR4@assays$ADT@data) <- gsub("-TotalSeqC", "", row.names(CAR4@assays$ADT@data))
row.names(CAR4@assays$ADT@scale.data) <- gsub("-TotalSeqC", "", row.names(CAR4@assays$ADT@scale.data))

keyFeatures2 <- c("YESCAR", "TBX21", "FOXP3", "IL2RA", "IKZF2", "CXCR5", "CCR7",
                 "IL7R", "GZMK", "KLRB1", "IRF7", "CTLA4", "LAG3", "MKI67",
                 "ADT_CD4", "ADT_CD25", "ADT_CD45RO", "ADT_PD-1")
keyFeatures <- unique(c(keyFeatures1, keyFeatures2))

p3b <- VlnPlot(CAR4, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset4$annoCoarse)
saveFig(p3b, prefix = "Plots/Fig3b_p2_mixed_Markers_CAR4", 
        wd = 3600, ht = 5000)

##============================-- Panel 2c CAR4 phenotype Distribution ----
## overall distribution
metaDf <- FetchData(CAR4, vars = c("annoCoarse", "stageFine")) %>%
  rename(anno=annoCoarse, time=stageFine)

dfx <- metaDf %>%
  group_by(time) %>%
  dplyr::count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = rev(c("Exp", "Per1", "Per2")))) %>%
  mutate(anno = factor(anno, levels = c("Treg", "Th1","Prolif","MEM", "ISG")))

p3c_p1 <- ggplot(dfx, aes(x=time, y=Proportion, fill=anno)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values=colset4$annoCoarse) +
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
saveFig(p3c_p1, 
        prefix = "Plots/Fig3c_p1_annoCoarse_distrib_overall", 
        wd = 2000, ht = 1200)

## per patient distribution
metaDf <- FetchData(CAR4, vars = c("annoCoarse", "stageFine", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageFine, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = rev(c("Exp", "Per1", "Per2")))) %>%
  mutate(anno = factor(anno, levels = c("Treg", "Th1","Prolif","MEM", "ISG"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))

p3c_p2 <- ggplot(dfx, aes(x=time, y=Proportion, fill=anno)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values=colset4$annoCoarse) +
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
saveFig(p3c_p2, 
        prefix = "Plots/Fig3c_p2_annoCoarse_distrib_perPatient", 
        wd = 2800, ht = 2000)

##============================-- Panel 2d CAR4 phenotype Proportions ----
### 3.4.1 stageFine
metaDf <- FetchData(CAR4, vars = c("annoCoarse", "stageFine", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageFine, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("Exp", "Per1", "Per2"))) %>%
  mutate(anno = factor(anno, levels = c("Treg", "Th1","Prolif","MEM", "ISG"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))


p3d_p1 <- ggboxplot(dfx, x="time", y="Proportion", add = "jitter", facet.by = "anno", 
          fill = "time", scales = "free_y", nrow = 2) +
  stat_compare_means(comparisons = list(c("Exp", "Per1"), c("Per1", "Per2"), c("Exp", "Per2")), 
                     label = "p.signif", method = "t.test") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colset4$stageFine[2:4]) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p3d_p1, 
        prefix = "Plots/Fig3d_p1_Pheno_Boxplot_annoCoarse_StageFine", 
        wd = 3000, ht = 2000)
  
### 3.4.2 StageCoarse
metaDf <- FetchData(CAR4, vars = c("annoCoarse", "stageCoarse", "patientShort")) %>%
  rename(anno=annoCoarse, time=stageCoarse, patient = patientShort)

dfx <- metaDf %>%
  group_by(time, patient) %>%
  count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("Exp", "Per"))) %>%
  mutate(anno = factor(anno, levels = c("Treg", "Th1","Prolif","MEM", "ISG"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient))))


p3d_p2 <- ggboxplot(dfx, x="time", y="Proportion", add = "jitter", facet.by = "anno", 
                fill = "time", scales = "free_y", nrow = 2) +
  stat_compare_means(label = "p.signif", comparisons = list(c("Exp", "Per"))) +
  scale_y_continuous(expand = expansion(mult = 0.05, add = 0.1)) +
  scale_fill_manual(values = colset4$stageCoarse[2:3]) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p3d_p2, 
        prefix = "Plots/Fig3d_p2_Pheno_Boxplot_annoCoarse_StageCoarse",  
        wd = 3000, ht = 2000)