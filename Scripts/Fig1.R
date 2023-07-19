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

##============================-- Panel 1c CAR expression ----
## load all T cell Raw Object; memory intensive
AllT <- readRDS("RDS/Temporal_AllT_b123_CD4_CD4_Pos_Apr26_2023.rds")

p1c <- VlnPlot(AllT, group.by = "Final_class", features = "YESCAR", pt.size = 0, 
               split.by = "celltype") +
  stat_compare_means(method = "wilcox.test", label="p.signif") +
  labs(y="Normalized Expression", x="", title="CAR Transgene") +
  theme(axis.text.y = element_text(size=15,color = 'black'),
        axis.text.x = element_text(size=20, angle=0, hjust = 0.5),
        axis.title.y = element_text(size=20,face = 'bold'),
        plot.title = element_text(size=24),
        legend.text = element_text(size=20)) +
  scale_fill_manual(values = c("#845EC2", "#FF9671")) +
  guides(fill = guide_legend(override.aes = list(size=8)))

saveFig(p1c, prefix = "Plots/Fig1c_VlnPlot_CAR_Expr_ALLT", 
        wd = 2000, ht = 1500)



##============================-- Panel 1d Percentage of CAR-T cells ----
dfx <- read.csv("Data/WBC_Percent.csv") %>%
  mutate(CAR_Count = CAR_abundance * WBC_Count * 1000) %>%
  mutate(CAR_perc = CAR_abundance*100) %>%
  mutate(time = as.numeric(gsub("d", "", date))) %>%
  mutate(patientShort = plyr::mapvalues(patient,
                                        from = c("P031", "P033", "P034", "P039", "P048", "P049", "P116"),
                                        to = c("P04", "P01", "P03", "P05", "P06", "P02", "P07"))) %>%
  mutate(patientShort = factor(patientShort, levels=mixedsort(unique(patientShort))))

p1d <- ggplot(dfx, aes(x=time, y=CAR_perc, group=patientShort, color=patientShort)) +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  theme_prism() +
  scale_color_prism() +
  labs(x="Day", y="% CAR") +
  scale_y_continuous(limits = c(0,80)) +
  theme(
    axis.text = element_text(size=10),
    legend.text = element_text(size=10),
    legend.position = "top")

saveFig(p1d, "Plots/Fig1d_CAR_Percentage", wd = 1200, ht = 1200)

### supplementary 2a: CD4 and CD8 breakdown


##============================-- Panel 1e Absolute Count of CAR-T cells ----
dfx <- read.csv("Data/WBC_Percent.csv") %>%
  mutate(CAR_Count = CAR_abundance * WBC_Count * 1000) %>%
  mutate(CAR_perc = CAR_abundance*100) %>%
  mutate(time = as.numeric(gsub("d", "", date))) %>%
  mutate(patientShort = plyr::mapvalues(patient,
                                        from = c("P031", "P033", "P034", "P039", "P048", "P049", "P116"),
                                        to = c("P04", "P01", "P03", "P05", "P06", "P02", "P07"))) %>%
  mutate(patientShort = factor(patientShort, levels=mixedsort(unique(patientShort))))

p1d <- ggplot(dfx, aes(x=time, y=CAR_Count, group=patientShort, color=patientShort)) +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  theme_prism() +
  scale_color_prism() +
  labs(x="Day", y="CAR Count") +
  theme(
    axis.text = element_text(size=10),
    legend.text = element_text(size=10),
    legend.position = "top")

saveFig(p1d, "Plots/Fig1e_CAR_Count", wd = 1200, ht = 1200)


##============================-- Panel 1f&g Overlap Coefficients of CAR-T cells ----
## load all TCR data ----
repCRsCAR <- readRDS("RDS/AggRepList_b123_PI_CAR_CR_wMeta_Patient_stageFine.rds")
repCRsEndo <- readRDS("RDS/AggRepList_b123_PI_Endo_CR_wMeta_Patient_stageFine.rds")

## get overlap ----
overlapCAR <- repOverlap(repCRsCAR$data, .method = "morisita")
overlapEndo <- repOverlap(repCRsEndo$data, .method = "morisita")

## bar plot ----
dfovCAR <- overlapCAR %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = !Sample1, names_to = "Sample2", values_to = "Overlap_Coeff") %>%
  filter(Overlap_Coeff >0) %>%
  mutate(s1Pt = str_split(Sample1, "_", simplify = T)[,2]) %>%
  mutate(s2Pt = str_split(Sample2, "_", simplify = T)[,2]) %>%
  mutate(s1Stage = str_split(Sample1, "_", simplify = T)[,3]) %>%
  mutate(s2Stage = str_split(Sample2, "_", simplify = T)[,3]) %>%
  filter(s1Pt == s2Pt) %>%
  mutate(group = unlist(Map(function(x,y){
    paste0(sort(c(x,y)), collapse = "_")
  }, s1Stage, s2Stage))) %>%
  select(s1Pt, group, Overlap_Coeff) %>%
  distinct() %>%
  mutate(celltype = "CART")

dfovEndo <- overlapEndo %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = !Sample1, names_to = "Sample2", values_to = "Overlap_Coeff") %>%
  filter(Overlap_Coeff >0) %>%
  mutate(s1Pt = str_split(Sample1, "_", simplify = T)[,2]) %>%
  mutate(s2Pt = str_split(Sample2, "_", simplify = T)[,2]) %>%
  mutate(s1Stage = str_split(Sample1, "_", simplify = T)[,3]) %>%
  mutate(s2Stage = str_split(Sample2, "_", simplify = T)[,3]) %>%
  filter(s1Pt == s2Pt) %>%
  mutate(group = unlist(Map(function(x,y){
    paste0(sort(c(x,y)), collapse = "_")
  }, s1Stage, s2Stage))) %>%
  select(s1Pt, group, Overlap_Coeff) %>%
  distinct() %>%
  mutate(celltype = "EndoT")

dfx <- rbind(dfovCAR, dfovEndo)


px <- ggplot(dfx, aes(x=group, y=Overlap_Coeff, fill=group)) +
  stat_summary(geom = "col", fun = mean, width = 0.7) +
  stat_summary(geom = "errorbar",
               fun = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.3) +
  theme_prism() +
  facet_wrap2(~celltype, ncol=2, scales = "free_y") +
  stat_compare_means(comparisons = list(c("S1_S2a", "S1_S2b"),
                                        c("S1_S2b", "S2a_S2b"),
                                        c("S1_S2a", "S2a_S2b")), label = "p.signif",
                     step.increase = c(0,0.04,0.05), tip.length = 0.02, method = "t.test") +
  scale_fill_npg() +
  labs(x="", y="Clonotype Overlap Coefficient") +
  theme(
    strip.text = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text = element_text(size=10),
    legend.text = element_text(size=10)
  )
saveFig(px, "Plots/Fig1fg_Barplot_TCR_clonotype_overlap_coefficient",
        wd = 2000, ht = 1200)


##============================-- Panel 1h Overlap Coefficients of CD8 and CD4 CAR-T cells ----
## load all TCR data ----
repCRsCAR4 <- readRDS("RDS/AggRepList_b123_PI_CAR4_Only_CR_wMeta_Patient_stageFine.rds")
repCRsCAR8 <- readRDS("RDS/AggRepList_b123_PI_CAR8_Only_CR_wMeta_Patient_stageFine.rds")

## get overlap ----
overlapCAR4 <- repOverlap(repCRsCAR4$data, .method = "morisita")
overlapCAR8 <- repOverlap(repCRsCAR8$data, .method = "morisita")

## bar plot ----
dfovCAR4 <- overlapCAR4 %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = !Sample1, names_to = "Sample2", values_to = "Overlap_Coeff") %>%
  filter(Overlap_Coeff >0) %>%
  mutate(Sample1 = gsub("CAR4_", "", Sample1)) %>%
  mutate(Sample2 = gsub("CAR4_", "", Sample2)) %>%
  mutate(s1Pt = str_split(Sample1, "-", simplify = T)[,1]) %>%
  mutate(s2Pt = str_split(Sample2, "-", simplify = T)[,1]) %>%
  mutate(s1Stage = str_split(Sample1, "-", simplify = T)[,2]) %>%
  mutate(s2Stage = str_split(Sample2, "-", simplify = T)[,2]) %>%
  filter(s1Pt == s2Pt) %>%
  mutate(group = unlist(Map(function(x,y){
    paste0(sort(c(x,y)), collapse = "_")
  }, s1Stage, s2Stage))) %>%
  select(s1Pt, group, Overlap_Coeff) %>%
  distinct() %>%
  mutate(celltype = "CAR4")

dfovCAR8 <- overlapCAR8 %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = !Sample1, names_to = "Sample2", values_to = "Overlap_Coeff") %>%
  filter(Overlap_Coeff >0) %>%
  mutate(Sample1 = gsub("CAR8_", "", Sample1)) %>%
  mutate(Sample2 = gsub("CAR8_", "", Sample2)) %>%
  mutate(s1Pt = str_split(Sample1, "-", simplify = T)[,1]) %>%
  mutate(s2Pt = str_split(Sample2, "-", simplify = T)[,1]) %>%
  mutate(s1Stage = str_split(Sample1, "-", simplify = T)[,2]) %>%
  mutate(s2Stage = str_split(Sample2, "-", simplify = T)[,2]) %>%
  filter(s1Pt == s2Pt) %>%
  mutate(group = unlist(Map(function(x,y){
    paste0(sort(c(x,y)), collapse = "_")
  }, s1Stage, s2Stage))) %>%
  select(s1Pt, group, Overlap_Coeff) %>%
  distinct() %>%
  mutate(celltype = "CAR8")

dfx <- rbind(dfovCAR4, dfovCAR8)


px <- ggplot(dfx, aes(x=group, y=Overlap_Coeff, fill=group)) +
  stat_summary(geom = "col", fun = mean, width = 0.7) +
  stat_summary(geom = "errorbar",
               fun = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.3) +
  theme_prism() +
  facet_wrap2(~celltype, ncol=2, scales = "free_y") +
  stat_compare_means(comparisons = list(c("Exp_Per1", "Exp_Per2"),
                                        c("Exp_Per2", "Per1_Per2"),
                                        c("Exp_Per1", "Per1_Per2")), 
                     label = "p.signif",
                     step.increase = c(0,0.04,0.05), 
                     tip.length = 0.02, method = "t.test") +
  scale_fill_npg() +
  labs(x="", y="Clonotype Overlap Coefficient") +
  theme(
    strip.text = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text = element_text(size=10),
    legend.text = element_text(size=10)
  )
saveFig(px, "Plots/Fig1h_Barplot_TCR_clonotype_overlap_coefficient_byCelltype",
        wd = 2000, ht = 1200)

##============================-- Panel 1I Clonality Overlap Map ----
## part 1: overlap matrix ----
CARTCR <- readRDS("RDS/TCR_table_all_IP_and_PI_CAR_Freq_Table_May20_2023.rds")
EndoTCR <- readRDS("RDS/TCR_table_all_IP_and_PI_Endo_Freq_Table_May20_2023.rds")
TCRList <- list(Endo=EndoTCR, CAR=CARTCR)

dfxList <- lapply(names(TCRList), function(gpx){
  dfFC <- TCRList[[gpx]] %>%
    select(patientShort, raw_clonotype_id, Per, Exp) %>%
    group_by(patientShort) %>%
    mutate(PerTotal = sum(Per), PerProp = Per/PerTotal) %>%
    mutate(ExpTotal = sum(Exp), ExpProp = Exp/ExpTotal) %>%
    ungroup() %>%
    mutate(celltype = gpx) %>%
    mutate(Exp_Clonal = case_when(
      Exp==1 ~ "Exp_Singlet",
      Exp>1 & ExpProp < 1e-3 ~ "Exp_Small",
      Exp>1 & ExpProp > 1e-3 ~ "Exp_Large",
      TRUE ~ "Not_Found"
    )) %>%
    mutate(Per_Clonal = case_when(
      Per==1 ~ "Per_Singlet",
      Per>1 & PerProp < 1e-3 ~ "Per_Small",
      Per>1 & PerProp > 1e-3 ~ "Per_Large",
      TRUE ~ "Not_Found"
    )) 
  dfFC
})
names(dfxList) <- names(TCRList)

### version 1: 3 level picture
dfx <- do.call("rbind", lapply(names(dfxList), function(gpx){
  dfxList[[gpx]] %>%
    filter(Exp_Clonal != "Not_Found" & Per_Clonal != "Not_Found") %>%
    mutate(Exp_Clonal = factor(Exp_Clonal, levels = c("Exp_Large", "Exp_Small", "Exp_Singlet"))) %>%
    mutate(Per_Clonal = factor(Per_Clonal, levels = rev(c("Per_Large", "Per_Small", "Per_Singlet")))) %>%
    count(Exp_Clonal, Per_Clonal, .drop = F) %>%
    mutate(celltype = gpx) %>%
    mutate(Prop = round(n/sum(n),3))
}))

### version 2: 2 level picture
dfx <- do.call("rbind", lapply(names(dfxList), function(gpx){
  dfxList[[gpx]] %>%
    filter(Exp_Clonal != "Not_Found" & Per_Clonal != "Not_Found") %>%
    mutate(Exp_Clonal = plyr::mapvalues(Exp_Clonal, from = c("Exp_Large", "Exp_Small"), 
                                        to = c("Exp_Expanded", "Exp_Expanded"))) %>%
    mutate(Per_Clonal = plyr::mapvalues(Per_Clonal, from = c("Per_Large", "Per_Small"), 
                                        to = c("Per_Expanded", "Per_Expanded"))) %>%
    mutate(Exp_Clonal = factor(Exp_Clonal, levels = c("Exp_Expanded", "Exp_Singlet"))) %>%
    mutate(Per_Clonal = factor(Per_Clonal, levels = rev(c("Per_Expanded", "Per_Singlet")))) %>%
    count(Exp_Clonal, Per_Clonal, .drop = F) %>%
    mutate(celltype = gpx) %>%
    mutate(Prop = round(n/sum(n),3))
}))

p1i_p1CAR <- ggplot(subset(dfx, celltype=="CAR"), aes(x=Exp_Clonal, y=Per_Clonal, fill=n)) +
  geom_tile(color = "black", linewidth = 1, linetype = 1) +
  geom_text(aes(label=n), size=5) +
  scale_fill_material("light-green") +
  theme_minimal() +
  labs(fill = "Number of Clones")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=45, 
                                   hjust=0.9, vjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text = element_text(size=15, face = "bold"))
saveFig(p1i_p1CAR, "Plots/Fig1I_p1CAR_TCR_cloneSize_overlap_acrossStage_2Lvl_Count",
        wd = 1200, ht = 1400)

p1i_p1Endo <- ggplot(subset(dfx, celltype=="Endo"), aes(x=Exp_Clonal, y=Per_Clonal, fill=n)) +
  geom_tile(color = "black", linewidth = 1, linetype = 1) +
  geom_text(aes(label=n), size=5) +
  scale_fill_material("orange") +
  theme_minimal() +
  labs(fill = "Number of Clones")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=45, 
                                   hjust=0.9, vjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text = element_text(size=15, face = "bold"))
saveFig(p1i_p1Endo, "Plots/Fig1I_p1Endo_TCR_cloneSize_overlap_acrossStage_2Lvl_Count",
        wd = 1200, ht = 1400)

p1i_p1CAR <- ggplot(subset(dfx, celltype=="CAR"), aes(x=Exp_Clonal, y=Per_Clonal, fill=Prop)) +
  geom_tile(color = "black", linewidth = 1, linetype = 1) +
  geom_text(aes(label=Prop), size=5) +
  scale_fill_material("light-green") +
  theme_minimal() +
  labs(fill = "Proportion of Clones")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=45, 
                                   hjust=0.9, vjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text = element_text(size=15, face = "bold"))
saveFig(p1i_p1CAR, "Plots/Fig1I_p1CAR_TCR_cloneSize_overlap_acrossStage_2Lvl_Prop",
        wd = 1200, ht = 1400)

p1i_p1Endo <- ggplot(subset(dfx, celltype=="Endo"), aes(x=Exp_Clonal, y=Per_Clonal, fill=Prop)) +
  geom_tile(color = "black", linewidth = 1, linetype = 1) +
  geom_text(aes(label=Prop), size=5) +
  scale_fill_material("orange") +
  theme_minimal() +
  labs(fill = "Proportion of Clones")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=45, 
                                   hjust=0.9, vjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text = element_text(size=15, face = "bold"))
saveFig(p1i_p1Endo, "Plots/Fig1I_p1Endo_TCR_cloneSize_overlap_acrossStage_2Lvl_Prop",
        wd = 1200, ht = 1400)

## part 2: Frequency Change ----
CARTCR <- readRDS("RDS/TCR_table_all_IP_and_PI_CAR_Freq_Table_May20_2023.rds")
EndoTCR <- readRDS("RDS/TCR_table_all_IP_and_PI_Endo_Freq_Table_May20_2023.rds")
TCRList <- list(Endo=EndoTCR, CAR=CARTCR)

dfxList <- lapply(names(TCRList), function(gpx){
  dfFC <- TCRList[[gpx]] %>%
    select(patientShort, raw_clonotype_id, Per, Exp) %>%
    group_by(patientShort) %>%
    mutate(PerTotal = sum(Per), PerProp = Per/PerTotal) %>%
    mutate(ExpTotal = sum(Exp), ExpProp = Exp/ExpTotal) %>%
    ungroup() %>%
    filter(PerProp > 1e-3 | ExpProp > 1e-3) %>%
    mutate(log2FC = log1p(Per) - log1p(Exp)) %>%
    mutate(celltype = gpx)

  dfFC
})
names(dfxList) <- names(TCRList)
dfFC <- do.call("rbind", dfxList)
gpCols <- stepped()


p1i_p2 <- ggplot(dfFC, aes(x=celltype, y=log2FC, fill=celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_violin(trim = T, color = "grey50") +
  theme_prism() +
  labs(y = "Clone_Size_Log2FC", fill="",
       title = "All Large TCR Clones (prop>1e-3)") +
  scale_fill_manual(values = c(CAR = "#03A64A", Endo = "#F2811D")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10),
        plot.title = element_text(size=12),
        legend.position = "none")
saveFig(p1i_p2, "Plots/Fig1I_p2_TCR_CloneSize_Log2FoldChange_Exp_to_Per",
        wd = 1500, ht = 1200)
