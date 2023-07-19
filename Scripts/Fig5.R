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
colset8 <- readRDS("RDS/color_set_overall.rds")

## 1.5 load CAR4 and CAR8 seurat objects
IP8 <- readRDS("RDS/IP8_for_manuscript_Apr262023.rds")


##============================-- Panel 5a IP CD8 Annotation ----
## main annotations
p5a_p1 <- DimPlot(IP8, reduction = "umap", group.by = "annoCoarse",
                  pt.size = 1.5, raster = F)+
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoIP8) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p5a_p1, prefix = "Plots/Fig5a_p1_UMAP_IP8_annoCoarse", 
        wd = 3600, ht = 3200)

## annotation distribution barplots: by patient
metaDf <- FetchData(IP8, vars = c("annoCoarse", "patientShort")) %>%
  rename(anno=annoCoarse, patient = patientShort)

dfx <- metaDf %>%
  group_by(patient) %>%
  dplyr::count(anno) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(patient = factor(patient, levels = rev(gtools::mixedsort(unique(patient))))) %>%
  mutate(anno = factor(anno, levels = c("EFF_TBX21", "EFF_GATA3", "Naive", 
                                        "Prolif")))

p5a_p2 <- ggplot(dfx, aes(x=patient, y=Proportion, fill=anno)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values=colset8$annoIP8) +
  theme_bw() + 
  labs(x="patient") +
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
  guides(fill=guide_legend(override.aes = list(size=6), nrow = 1))
saveFig(p5a_p2, 
        prefix = "Plots/FigS8c_annoCoarse_distrib_inIP", 
        wd = 2000, ht = 1200)



##============================-- Panel 5b IP8 Markers VlnPlot ----
Idents(IP8) <- factor(IP8$annoCoarse, levels = c("EFF_TBX21", "EFF_GATA3", "Naive", 
                                                   "Prolif"))
## RNA markers
keyFeatures <- c("CD8A", "CD8B", "YESCAR",
                 "GZMB", "PRF1", "NKG7", 
                 "TBX21","KLRG1","GNLY",
                 "GZMK", "GZMM", "GZMH",
                 "GATA3", "CCR4", "IL4R",
                 "FOXP3", "IL10RA", "TNFRSF4",
                 "LAG3", "TIGIT",  "HAVCR2",
                 "CCR7", "TCF7", "SELL",
                 "LEF1", "LRRN3", "ADT_CD45RA",
                 "MKI67", "TOP2A")

px <- VlnPlot(IP8, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset8$annoIP8)
saveFig(px, prefix = "Plots/Fig5b_p1_RNA_Markers_IP8", 
        wd = 2400, ht = 3600)

## Mixed markers
row.names(IP8@assays$ADT@data) <- gsub("-TotalSeqC", "", row.names(IP8@assays$ADT@data))
row.names(IP8@assays$ADT@scale.data) <- gsub("-TotalSeqC", "", row.names(IP8@assays$ADT@scale.data))

keyFeatures <- unique(c("CD8A", "ADT_CD8A","YESCAR",
                 "IFITM1", "IFITM2", "ISG20", "IRF7",
                 "CXCR3", "CD38", "ADT_CD45RO", "ADT_CD45RA", 
                 "IL7R", "CXCR6", "ADT_CD127", 
                 "CXCR4","ADT_CXCR3","ADT_CCR4", "CCR4",
                 "PDCD1", "TOX", "NR4A2", "ADT_PD-1",
                 "ADT_TIM-3", "ADT_LAG-3", "ADT_CD39",
                 "ADT_CCR7","CD8A", "CD8B", "YESCAR",
                 "GZMB", "PRF1", "NKG7", 
                 "TBX21","KLRG1","GNLY",
                 "GZMK", "GZMM", "GZMH",
                 "GATA3", "CCR4", "IL4R",
                 "FOXP3", "IL10RA", "TNFRSF4",
                 "LAG3", "TIGIT",  "HAVCR2",
                 "CCR7", "TCF7", "SELL",
                 "LEF1", "LRRN3", "ADT_CD45RA",
                 "MKI67", "TOP2A"))

p2b <- VlnPlot(IP8, keyFeatures, stack = T, flip = T, fill.by = "ident") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 330, hjust = 0, size=15),
        axis.text.y.right = element_text(hjust = 0, size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y.left = element_text(size=20),
        strip.text = element_text(size=15, hjust = 0)) +
  scale_fill_manual(values = colset8$annoIP8)
saveFig(p2b, prefix = "Plots/Fig5b_p2_mixed_Markers_IP8", 
        wd = 3600, ht = 7000)

## density plot
require(Nebulosa)
keyMarkers <- c("YESCAR", "KLRG1", "TBX21", "MKI67",
                "IL7R", "CCR7", "LRRN3", "PLAC8",
                "TCF7", "ADT_CD45RA", "FOXP3", "GATA3",
                "NEAT1", "CCR4", "GZMB", "PRF1", "FAS")

for (gx in keyMarkers){
  p5b_p2 <- plot_density(IP8, gx)
  saveFig(p5b_p2, prefix = paste0("Plots/Density_Plots/Fig5b_p3_density_markers_", gx),
          wd = 1800, ht = 1500)
}


##============================-- Panel 5c Exp and Per Precursors Tracing ----
## part 1: overall distribution
metaDfExp <- FetchData(IP8, vars = c("patient", "annoCoarse", "inExp")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  mutate(group = "Exp_pre") %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patient", "annoCoarse", "inPer")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  mutate(group = "Per_pre") %>%
  select(-inPer)

metaDfOverall <- FetchData(IP8, vars = c("patient", "annoCoarse", "inPer")) %>%
  filter(!is.na(inPer)) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  mutate(group = "Overall")

metaAll <- bind_rows(metaDfExp, metaDfPer)
metaAll <- bind_rows(metaAll, metaDfOverall) %>%
  mutate(group = factor(group, levels = c("Overall", "Exp_pre", "Per_pre")))


p5c_p1 <- ggplot(metaAll, aes(x=group, y=Proportion, fill=annoCoarse)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colset8$annoIP8) +
  theme_bw() + 
  labs(x="Precursor Stages") +
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
  guides(fill=guide_legend(override.aes = list(size=6), nrow = 1))
saveFig(p5c_p1, 
        prefix = "Plots/Fig5c_p1_precursor_annoCoarse_distrib_overall", 
        wd = 2000, ht = 1200)


## part 2: per patient comparisons
metaDfExp <- FetchData(IP8, vars = c("patientShort", "annoCoarse", "inExp")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, patientShort) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  mutate(group = "Exp_pre") %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patientShort", "annoCoarse", "inPer")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, patientShort) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  mutate(group = "Per_pre") %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5c_p2 <- ggboxplot(metaAll, x="annoCoarse", y="Proportion", add = "jitter", facet.by = "group", 
                    fill = "annoCoarse", scales = "free_y", nrow=1) +
  stat_compare_means(ref.group = "Naive", 
                     label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colset8$annoIP8) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black",  angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p5c_p2, 
        prefix = "Plots/Fig5c_p2_precursors_Pheno_Boxplot_annoCoarse", 
        wd = 3000, ht = 2000)

## part 3: overall distribution: split by clonal vs singlet
metaDfExp <- FetchData(IP8, vars = c("patient", "annoCoarse", "inExp", "Exp_Clonal")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, Exp_Clonal) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  rename(isClonal = Exp_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Exp_pre_", isClonal)) %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patient", "annoCoarse", "inPer", "Per_Clonal")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, Per_Clonal) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  rename(isClonal = Per_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Per_pre_", isClonal)) %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5c_p3 <- ggplot(metaAll, aes(x=group, y=Proportion, fill=annoCoarse)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colset8$annoIP8) +
  theme_bw() + 
  labs(x="Precursor Stages") +
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
  guides(fill=guide_legend(override.aes = list(size=6), nrow = 1))
saveFig(p5c_p3, 
        prefix = "Plots/Fig5c_p3_precursor_annoCoarse_distrib_overall_byClonality", 
        wd = 2000, ht = 1200)


## part 4: per patient comparisons
metaDfExp <- FetchData(IP8, vars = c("patientShort", "annoCoarse", "inExp", "Exp_Clonal")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, Exp_Clonal, patientShort) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  rename(isClonal = Exp_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Exp_pre_", isClonal)) %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patientShort", "annoCoarse", "inPer", "Per_Clonal")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, Per_Clonal, patientShort) %>%
  count(annoCoarse) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  rename(isClonal = Per_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Per_pre_", isClonal)) %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5c_p4 <- ggboxplot(metaAll, x="annoCoarse", y="Proportion", add = "jitter", facet.by = "group", 
                    fill = "annoCoarse", scales = "free_y", nrow=1) +
  stat_compare_means(ref.group = "Naive", 
                     label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colset8$annoIP8) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black", angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p5c_p4, 
        prefix = "Plots/Fig5c_p4_precursors_Pheno_Boxplot_annoCoarse_byClonality", 
        wd = 4000, ht = 2000)


##============================-- Panel 5d CAR+ vs CAR- populations ----
dfx <- FetchData(IP8, vars = c("rna_YESCAR", "ADT_CAR-TotalSeqC"))
IP8$CARgroupRNA <- ifelse(dfx$rna_YESCAR > 0, "rPos", "rNeg")
IP8$CARgroupADT <- ifelse(dfx$`ADT_CAR-TotalSeqC` > 0.5, "aPos", "aNeg")
IP8$CARgroup <- paste0(IP8$CARgroupRNA, "_", IP8$CARgroupADT)
colCARgroup <- rev(pal_npg()(4))

## part 1: overall distribution
metaDfExp <- FetchData(IP8, vars = c("patient", "CARgroup", "inExp")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  mutate(group = "Exp_pre") %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patient", "CARgroup", "inPer")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  mutate(group = "Per_pre") %>%
  select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5d_p1 <- ggplot(metaAll, aes(x=group, y=Proportion, fill=CARgroup)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colCARgroup) +
  theme_bw() + 
  labs(x="Precursor Stages") +
  coord_flip() +
  scale_y_continuous(position = "right") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth=1, color = "black")) +
  guides(fill=guide_legend(override.aes = list(size=5), nrow = 1))
saveFig(p5d_p1, 
        prefix = "Plots/Fig5d_p1_precursor_CARgroup_distrib_overall", 
        wd = 2000, ht = 1200)


## part 2: per patient comparisons
metaDfExp <- FetchData(IP8, vars = c("patientShort", "CARgroup", "inExp")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, patientShort) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  mutate(group = "Exp_pre") %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patientShort", "CARgroup", "inPer")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, patientShort) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  mutate(group = "Per_pre") %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5d_p2 <- ggboxplot(metaAll, x="CARgroup", y="Proportion", add = "jitter", facet.by = "group", 
                    fill = "CARgroup", scales = "free_y", nrow=1) +
  stat_compare_means(ref.group = "rNeg_aNeg", 
                     label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colCARgroup) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black",  angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p5d_p2, 
        prefix = "Plots/Fig5d_p2_precursors_Pheno_Boxplot_CARgroup", 
        wd = 3000, ht = 2000)

## part 3: overall distribution: split by clonal vs singlet
metaDfExp <- FetchData(IP8, vars = c("patient", "CARgroup", "inExp", "Exp_Clonal")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, Exp_Clonal) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  rename(isClonal = Exp_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Exp_pre_", isClonal)) %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patient", "CARgroup", "inPer", "Per_Clonal")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, Per_Clonal) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  rename(isClonal = Per_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Per_pre_", isClonal)) %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5d_p3 <- ggplot(metaAll, aes(x=group, y=Proportion, fill=CARgroup)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colCARgroup) +
  theme_bw() + 
  labs(x="Precursor Stages") +
  coord_flip() +
  scale_y_continuous(position = "right") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth=1, color = "black")) +
  guides(fill=guide_legend(override.aes = list(size=5), nrow = 1))
saveFig(p5d_p3, 
        prefix = "Plots/Fig5d_p3_precursor_CARgroup_distrib_overall_byClonality", 
        wd = 2000, ht = 1200)


## part 4: per patient comparisons
metaDfExp <- FetchData(IP8, vars = c("patientShort", "CARgroup", "inExp", "Exp_Clonal")) %>%
  filter(!is.na(inExp)) %>%
  group_by(inExp, Exp_Clonal, patientShort) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inExp=="Yes") %>%
  rename(isClonal = Exp_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Exp_pre_", isClonal)) %>%
  dplyr::select(-inExp)

metaDfPer <- FetchData(IP8, vars = c("patientShort", "CARgroup", "inPer", "Per_Clonal")) %>%
  filter(!is.na(inPer)) %>%
  group_by(inPer, Per_Clonal, patientShort) %>%
  count(CARgroup) %>%
  mutate(Proportion = n/sum(n)) %>%
  ungroup() %>%
  filter(inPer=="Yes") %>%
  rename(isClonal = Per_Clonal) %>%
  mutate(isClonal = plyr::mapvalues(isClonal, from=c("Yes", "No"), to = c("Clonal", "Singlets"))) %>%
  mutate(group = paste0("Per_pre_", isClonal)) %>%
  dplyr::select(-inPer)

metaAll <- bind_rows(metaDfExp, metaDfPer)


p5d_p4 <- ggboxplot(metaAll, x="CARgroup", y="Proportion", add = "jitter", facet.by = "group", 
                    fill = "CARgroup", scales = "free_y", nrow=1) +
  stat_compare_means(ref.group = "rNeg_aNeg", 
                     label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = colCARgroup) +
  theme_prism() +
  theme(axis.text.x = element_text(size=12, color = "black", angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12)) 

saveFig(p5d_p4, 
        prefix = "Plots/Fig5d_p4_precursors_Pheno_Boxplot_CARgroup_byClonality", 
        wd = 4000, ht = 2000)


##============================-- Panel 5e Precursors DEGs ----
## part 1: DEG panel: all Exp vs all Per; Volcano Plots ---
IP8$preGroup <- ifelse(IP8$inExp=="Yes", "Exp_Pre", 
                       ifelse(IP8$inExp=="No" & IP8$inPer=="Yes", "Per_Pre", "Others"))
preDEG <- FindMarkers(IP8, ident.1 = "Exp_Pre", ident.2 = "Per_Pre", group.by = "preGroup", logfc.threshold = 0.01)
saveRDS(preDEG, "RDS/DEG_inIP_Exp_vs_Per_Precursors_May02.rds")

preDf <- preDEG %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val_adj<0.05)

## set parameters
preDEG <- readRDS("RDS/DEG_inIP_Exp_vs_Per_Precursors_May02.rds")
require(ggrepel)
volcano_col <- c(unname(colset8$stageCoarse[2:3]), 'grey80')
Exp_genes <- c("TIGIT", "GZMB", "GZMK", "GZMH", "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
               "FCGR3A", "GNLY", "JUN", "JUND", "CCL5", "NKG7","IFNG", "LAG3", "FOS",
               "KLRD1", "KLRG1", "CD81")
Per_genes <- c("CD52", "YESCAR", "SELL", "CD7", "IL7R", "IL4R","S100A4", "S100A6", "NEAT1",
               "LRRN3", "CD27", "S1PR1", "SLAMF1", "ISG20", "IFITM1", "IFITM2", "IFNAR1", "IFNAR2",
               "IRF7", "CD28", "ZAP70", "IL10RA", "GATA3")
gene_highlight <- c(Exp_genes, Per_genes)

## preprocessing
res <- preDEG %>%
  rownames_to_column(var = "gene") %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(logp = ifelse(logp>100, 100, logp)) %>%
  mutate(avg_log2FC = clipMtx(avg_log2FC, 1.2)) %>%
  mutate(category = case_when(
    logp > 3 & avg_log2FC > 0.25 ~ "Up in Exp Pre",
    logp > 3 & avg_log2FC < -0.25 ~ "Up in Per Pre",
    TRUE ~ "Non-DEGs"
  )) %>%
  mutate(gene_highlight = ifelse(gene %in% gene_highlight, gene, "")) %>%
  mutate(gene_highlight = ifelse(category=="Non-DEGs", "", gene_highlight)) %>%
  mutate(category = factor(category, levels = c("Up in Exp Pre", "Up in Per Pre", "Non-DEGs")))

res_Highlight <- subset(res, gene_highlight!="") 

## plot
p5e_p1 <- ggplot(res, aes(x=avg_log2FC, y=logp, color=category, label=gene_highlight))+
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
saveFig(p5e_p1, 
        prefix = "Plots/Fig5e_p1_volcanoPlot_DEG_Exp_vs_Per_Precursors_inIP", 
        wd = 3000, ht = 2000)


## part 2: DEG panel: select DEGs in Tile Heatmap ---
Exp_genes <- c("TIGIT", "GZMB", "GZMK", "GZMH", "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
               "FCGR3A", "GNLY", "JUN", "JUND", "CCL5", "NKG7","LAG3", "FOS",
               "KLRD1", "KLRG1", "CD81", "IFNG")
Per_genes <- c("CD52", "YESCAR", "SELL", "CD7", "IL7R", "IL4R","S100A4", "S100A6", "NEAT1",
               "LRRN3", "CD27", "S1PR1", "SLAMF1", "ISG20", "IFITM1", "IFITM2", "IFNAR1", "IFNAR2",
               "IRF7", "CD28", "ZAP70", "IL10RA", "GATA3")
gene_highlight <- c(Exp_genes, Per_genes)


sigSet <- list(
  Eff_Func = c("GNLY", "NKG7", "IFNG", "TNF", "GZMB", "GZMH", "GZMK", "KLRG1", "KLRD1", "HLA-DRA", "CCL5"),
  Exh = c("TIGIT", "LAG3"),
  AP1 = c("JUN", "FOS", "JUND"),
  Stem = c("SELL", "IL7R", "CD7", "CD27", "S1PR1"),
  Type2 = c("IL4R", "GATA3", "CXCR4"),
  IFN1 = c("IFITM1", "IFITM2", "IFNAR1", "IFNAR2"),
  others = c("NEAT1", "S100A4", "S100A6", "CD28", "ZAP70", "IL10RA", "CXCR3")
)


geneSelect <- unname(unlist(sigSet))

gsDf <- do.call("rbind", lapply(names(sigSet), function(sx){
  ggx <- sigSet[[sx]]
  data.frame(Gene=ggx, gLab = rep(sx, length(ggx)))
}))

IP8 <- ScaleData(IP8,features = geneSelect, 
                 vars.to.regress = c("percent.mt"))

sigDf <- FetchData(IP8, vars = c(geneSelect, "preGroup"), slot = "scale.data")

sigAvg <- sigDf %>%
  filter(!is.na(preGroup)) %>%
  mutate(preGroup = gsub("Others", "IP_Only", preGroup)) %>%
  pivot_longer(cols = !c("preGroup"), names_to = "Gene", values_to = "Expr") %>%
  group_by(Gene, preGroup) %>%
  summarise(AvgExpr = mean(Expr, na.rm=T)) %>%
  left_join(gsDf, by = "Gene") %>%
  mutate(Gene = factor(Gene, levels = geneSelect)) %>%
  mutate(gLab = factor(gLab, levels = names(sigSet))) %>%
  mutate(preGroup = factor(preGroup, levels = rev(c("Exp_Pre","Per_Pre", "IP_Only")))) %>%
  ungroup() 

p5e_p2 <- ggplot(sigAvg, aes(x=Gene, y=preGroup, fill = AvgExpr)) +
  geom_tile(color = "black", lwd = 1, linetype = 1) +
  theme_minimal() +
  labs(fill = "Normalized\nExpression")+
  facet_grid2(.~gLab, scales = "free_x", space = "free_x") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15, angle=60, hjust=1, 
                                   vjust=1, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black", face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=15, face = "bold")) +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", mid=0) +
  guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
saveFig(p5e_p2, 
        prefix = "Plots/Fig5e_p2_heatmap_DEG_Exp_vs_Per_Precursors_inIP", 
        wd = 3000, ht = 1300)

## part 3: GSEA panel: select GSEA summary dotplot ---
preGSEA <- getGSEA(preDEG)
saveRDS(preGSEA, "RDS/GSEA_inIP_Exp_vs_Per_Precursors_May02.rds")

preGSEA <- readRDS("RDS/GSEA_inIP_Exp_vs_Per_Precursors_May02.rds")
up_terms <- c("HALLMARK_MYC_TARGETS_V1","REACTOME_PD_1_SIGNALING", 
                  "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CHEMOTAXIS",
                  "GOBP_RESPONSE_TO_INTERLEUKIN_1", 
                  "REACTOME_INTERFERON_GAMMA_SIGNALING",
                  "GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
                  "REACTOME_TCR_SIGNALING")
down_terms <- c("GSE28726_NAIVE_VS_ACTIVATED_CD4_TCELL_UP",
                "GSE22886_NAIVE_TCELL_VS_NKCELL_UP", 
                "GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP",
                "GSE22886_NAIVE_CD8_TCELL_VS_NKCELL_UP")
select_terms <- c(up_terms, rev(down_terms))

selectDf <- preGSEA %>%
  filter(ID %in% select_terms) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels=rev(ID))) %>%
  mutate(logq = -log10(qvalue)) %>%
  mutate(enrichment = ifelse(NES>0, "Exp_Pre", "Per_Pre")) %>%
  mutate(enrichment = factor(enrichment, levels=c("Exp_Pre", "Per_Pre")))


p5e_p3 <- ggplot(selectDf, aes(x=ID, y=NES)) +
  geom_segment(aes(y = 0,x = ID, yend = NES, xend = ID), 
               color = "black", linewidth=1) +
  geom_point(aes(color=enrichment), size=6) +
  theme_bw() +
  coord_flip() +
  labs(x="", y="NES") +
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = "top")+
  scale_color_manual(values = unname(colset8$stageCoarse[c(2:3)])) +
  guides(color = guide_legend(title = "",override.aes = list(size=6), nrow = 1)) 
saveFig(p5e_p3, 
        prefix = "Plots/Fig5e_p3_GSEA_Summary_Exp_vs_Per_Precursors_inIP", 
        wd = 2800, ht = 2000)


##============================-- Panel 5f Signature Scores ----
IP8$preGroup <- ifelse(IP8$inExp=="Yes", "Exp_Pre", 
                       ifelse(IP8$inExp=="No" & IP8$inPer=="Yes", "Per_Pre", "IP_Only"))

### part 1: msigdb signature ----
require(UCell)
selectList <- c("GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP", 
                "GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_DN", 
             "GSE9650_NAIVE_VS_EFF_CD8_TCELL_DN",
             "GSE9650_NAIVE_VS_EFF_CD8_TCELL_UP",
             "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP",
             "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN",
             "GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
             "REACTOME_INTERFERON_GAMMA_SIGNALING",
             "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
             "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING",
             "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
             "REACTOME_DNA_REPLICATION",
             "REACTOME_PD_1_SIGNALING")
sigList <- GetMsigdbSet(selectList)
IP8 <- AddModuleScore_UCell(IP8, features = sigList, name = "")

metaDf <- FetchData(IP8, vars = c(names(sigList), "preGroup", "annoCoarse")) %>%
  filter(!is.na(preGroup)) %>%
  mutate(preGroup = factor(preGroup, levels = c("Exp_Pre", "Per_Pre", "IP_Only")))

for (gx in names(sigList)){
  p5f_p1 <- ggviolin(metaDf, x="preGroup", y=gx, facet.by = "annoCoarse", add = "boxplot",
                 add.params = list(fill="white"), fill = "preGroup", 
                 nrow = 1) +
    stat_compare_means(label = "p.signif", 
                       comparisons = list(c("Exp_Pre", "Per_Pre"), c("Per_Pre", "IP_Only"), c("Exp_Pre", "IP_Only")),
                       size=3) +
    labs(y="Signature Score", title=gx, x="time") +
    scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
    scale_fill_manual(values = unname(colset8$stageCoarse[c(2,3,1)])) +
    theme_prism() +
    theme(axis.text.x = element_text(size=12, color = "black", 
                                     angle=60, hjust=1, vjust=1),
          axis.title.x = element_text(size=18, face = "bold"),
          axis.text.y = element_text(size=12, color = "black"),
          axis.title.y = element_text(size=18, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          strip.text = element_text(size=12),
          plot.title = element_text(size=12, hjust=0.5, color="blue", face="bold")) 
  saveFig(p5f_p1, 
          paste0("Plots/SigVlnPlot/Fig5f_Signature_Plot_preGroup_", gx,"_Split"),
          wd = 3000, ht = 1200)
  
  p5f_p2 <- ggviolin(metaDf, x="preGroup", y=gx,add = "boxplot",
                  add.params = list(fill="white"), fill = "preGroup", 
                  nrow = 1) +
    stat_compare_means(label = "p.signif", 
                       comparisons = list(c("Exp_Pre", "Per_Pre"), c("Per_Pre", "IP_Only"), c("Exp_Pre", "IP_Only")),
                       size=3) +
    labs(y="Signature Score", title=gx, x="time") +
    scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
    scale_fill_manual(values = unname(colset8$stageCoarse[c(2,3,1)])) +
    theme_prism() +
    theme(axis.text.x = element_text(size=12, color = "black", 
                                     angle=60, hjust=1, vjust=1),
          axis.title.x = element_text(size=18, face = "bold"),
          axis.text.y = element_text(size=12, color = "black"),
          axis.title.y = element_text(size=18, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          strip.text = element_text(size=12),
          plot.title = element_text(size=12, hjust=0.5, color="blue", face="bold")) 
  saveFig(p5f_p2, 
          paste0("Plots/SigVlnPlot/Fig5f_Signature_Plot_preGroup_", gx, "_Overall"),
          wd = 2000, ht = 2000)
}

##============================-- Panel 5g Custom Signatures ----
require(UCell)
preDEG <- readRDS("RDS/DEG_inIP_Exp_vs_Per_Precursors_May02.rds")
require(msigdbr)
plotGSEA <- function(xx, gID, Ntext=1200, pTitle = NA){
  require(ggpp)
  pathName <- xx@result[gID, "Description"]
  NES <- round(xx@result[gID, "NES"],2)
  qVal <- format(round(xx@result[gID, "qvalue"],3), nsmall=3)
  pVal <- format(round(xx@result[gID, "p.adjust"],3), nsmall=3)
  subText <- paste0("NES=", NES, "\nFDR q-value=", qVal, "\np.adjust= ", pVal)
  
  if (!is.na(pTitle)){
    pathName = pTitle
  }
  
  pRaw <- gseaplot2(xx, geneSetID = gID, 
                    title = NA,
                    subplots = 1:2,
                    rel_heights = c(1.2,0.3), base_size = 10, pvalue_table =F)
  
  pRaw[[1]] <- pRaw[[1]] +
    labs(title = pathName) +
    annotate(geom="text_npc",
             npcx = c("right"),
             npcy = c("top"),
             label=subText,
             color="black", hjust="inward", size=4) +
    theme(plot.title = element_text(size=15, hjust=0.5, face="bold"),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))
  
  pRaw[[2]] <- pRaw[[2]] +
    theme(axis.text.x = element_text(size=12))
  
  return(pRaw)
}


#### 2.1 EMRA vs SCM in JW Human Paper
library(readxl)
library(enrichplot)
dfFile <- "Data/JW_Human_Atlas_DEG_Overall.xlsx"
gps <- excel_sheets(dfFile)
JWList <- unlist(lapply(gps, function(gpx) {
  dfx <- read_excel(path = dfFile, sheet = gpx) 
  dfxClean <- dfx %>%
    filter(padj<0.05) %>%
    filter(abs(log2FoldChange) > 1.5) %>%
    rename(gene = geneSymbol)
  dfxList <- split(dfxClean, f=dfxClean$subset_comparison)
  dfxList
}), recursive = F)

sigSetUP <- lapply(JWList, function(dfx){
  dfx %>% 
    filter(log2FoldChange>0) %>%
    arrange(desc(log2FoldChange)) %>%
    pull(gene) %>% 
    GetGenelist() %>%
    unique()
})
names(sigSetUP) <- paste0(names(sigSetUP), "_UP")

sigSetDN <- lapply(JWList, function(dfx){
  dfx %>% 
    filter(log2FoldChange<0) %>%
    arrange(log2FoldChange) %>%
    pull(gene) %>% 
    GetGenelist() %>%
    unique()
})
names(sigSetDN) <- paste0(names(sigSetDN), "_DN")

sigSet <- c(sigSetUP, sigSetDN)

m_df <- do.call("rbind", lapply(names(sigSet), function(sigx){
  data.frame(gs_name = rep(sigx, length(sigSet[[sigx]])), 
             entrez_gene = as.integer(sigSet[[sigx]])) %>%
    as_tibble()
}))

## part 1: EM EX vs EFF ----
gvec <- GetGenelist(preDEG)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, 
                pAdjustMethod = "fdr", minGSSize = 5, maxGSSize = Inf)

glist<- c("scmR3+_vs_PD1+CD39+_UP", "scmR3+_vs_PD1+CD39+_DN", 
          "naive2_vs_PD1+CD39+_UP", "naive2_vs_PD1+CD39+_DN", 
          "scmR3+_vs_effectorMemory1_UP", "scmR3+_vs_effectorMemory1_DN", 
          "effectorMemory1_vs_centralMemory_UP", "effectorMemory1_vs_centralMemory_DN")

ID_alias = glist

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p5g_p1 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p5g_p1, paste0("Plots/GSEA/Fig5g_p1_GSEA_Exp_vs_Per_Precursor_", pathx),
          wd = 1800, ht = 1200)
}

