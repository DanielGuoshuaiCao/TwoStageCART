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

### Major message: robust exhaustion signature but no epigenetic imprinting

##============================-- Panel 6a EX Signature ----
## Section 1: signature Plots ----
require(UCell)
selectList <- c("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN",
                "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP")
msigList <- GetMsigdbSet(selectList)

### JW mouse Tex subsets
JWdf <- read.csv("Data/JWAllClusterDEG.csv")
JWList <- split(JWdf, f=JWdf$cluster)
sigSet <- lapply(JWList, function(dfx){
  dfx %>% 
    arrange(desc(avg_log2FC)) %>%
    head(500) %>%
    pull(gene) %>% 
    unique() %>% 
    HM_Convert(HtoM=F, simple = T) %>% 
    na.omit()
})
names(sigSet) <- paste0("JW_", names(sigSet))
names(sigSet) <- gsub("-", "_", names(sigSet))
names(sigSet) <- gsub(" ", "", names(sigSet))

### JW human atlas
dfAtlas <- read.csv("Data/JWHumanAtlas.csv")
gList1 <- lapply(c(0:6), function(nthresh){
  dfAtlas %>%
    filter(log2FoldChange>1.5 & padj<0.05) %>%
    count(geneSymbol) %>%
    filter(n>nthresh) %>%
    pull(geneSymbol) %>%
    unique()
})
names(gList1) <- paste0("CoreSigIn", 1:7)

gList2 <- lapply(split(dfAtlas, f=dfAtlas$subset_comparison), function(dfx){
  dfx %>%
    filter(log2FoldChange>1.5 & padj<0.05) %>%
    pull(geneSymbol) %>%
    unique()
})
xx <- names(gList2)
x1 <- gsub("\\+", "pos", xx)
x2 <- gsub("\\-", "neg", x1)
names(gList2) <- x2

gListAtlas <- c(gList1, gList2)

### combine
sigList <- c(msigList, sigSet, gListAtlas)

CAR8 <- AddModuleScore_UCell(CAR8, features = sigList, name = "")

metaDf <- FetchData(CAR8, vars = c(names(sigList), "stageCoarse", "annoCoarse")) %>%
  mutate(annoCoarse = factor(annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                  "EM_Eff-like", "CM", "Prolif"))) %>%
  mutate(stageCoarse = factor(stageCoarse, levels=c("Exp", "Per")))

## part 1: msigdb pathways; VlnPlots ----
for (gx in c(names(msigList), "JW_Exh_Term", "JW_Pre_Exh", "CoreSigIn1")){
  metaDfx <- metaDf %>%
    select(stageCoarse, annoCoarse, all_of(gx)) %>%
    rename(SigScore = all_of(gx)) %>%
    filter(SigScore > (mean(SigScore) - 3*IQR(SigScore)))
  
  meanDfx <- metaDfx %>%
    group_by(annoCoarse) %>%
    summarize(meanScore = mean(SigScore)) %>%
    ungroup() %>%
    arrange(desc(meanScore))
  
  metaDfx <- metaDfx %>%
    mutate(annoCoarse = factor(annoCoarse, levels = meanDfx$annoCoarse))
  
  p6a_p1 <- ggviolin(metaDfx, x="stageCoarse", y="SigScore", facet.by = "annoCoarse", add = "boxplot",
                     add.params = list(fill="white"), fill = "stageCoarse", 
                     nrow = 1, trim=T) +
    stat_compare_means(label = "p.signif", 
                       comparisons = list(c("Exp", "Per")),
                       size=3) +
    labs(y="Signature Score", title=gx, x="time") +
    scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
    scale_fill_manual(values = colset8$stageCoarse[c(2,3)]) +
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
          paste0("Plots/SigVlnPlot2/Fig6a_Ex_Signature_Plot_annoCoarse_", gx,"_Split"),
          wd = 3000, ht = 1200)
  
  p6a_p2 <- ggviolin(metaDfx, x="annoCoarse", y="SigScore",add = "boxplot",
                     add.params = list(fill="white"), fill = "annoCoarse", 
                     nrow = 1, trim=T) +
    stat_compare_means(label = "p.signif", 
                       ref.group = "EFF_Exh-like",
                       size=3) +
    labs(y="Signature Score", title=gx, x="time") +
    scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
    scale_fill_manual(values = colset8$annoCoarse) +
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
          paste0("Plots/SigVlnPlot2/Fig6a_EX_Signature_Plot_annoCoarse_", gx, "_Overall"),
          wd = 2000, ht = 2000)
}

## part 2: JW anno; heatmap; exclude prolif ----
sigAvg <- metaDf %>%
  select(-stageCoarse) %>%
  pivot_longer(cols = !annoCoarse, names_to = "Signature", values_to = "SigScore") %>%
  filter(annoCoarse!="Prolif") %>%
  filter(Signature %in% c("JW_Exh_Term", "JW_Pre_Exh", "JW_TransI",
                          "JW_TransII", "JW_TransCTL")) %>%
  group_by(Signature, annoCoarse) %>%
  summarise(AvgScore = median(SigScore)) %>%
  mutate(AvgScore = scale(AvgScore)) %>%
  mutate(annoCoarse = factor(annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                    "EM_Eff-like", "CM")))


p6a_p2 <- ggplot(sigAvg, aes(x=annoCoarse, y=Signature, fill = AvgScore)) +
  geom_tile(color = "black", lwd = 1, linetype = 1) +
  theme_minimal() +
  labs(fill = "Normalized\nScore")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15, angle=60, hjust=1, vjust=1),
        axis.text.y = element_text(size=12),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=10)) +
  scale_fill_continuous_diverging(palette = "Blue-Red3", mid=0) +
  guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
saveFig(p6a_p2, 
        prefix = "Plots/Fig6a_p2_Heatmap_JW_SigScore", 
        wd = 2000, ht = 1200)


## part 3: UMAP highlighting Tex subsets ----
metaDf <- FetchData(CAR8, vars = c("UMAP_1", "UMAP_2", "annoCoarse", "stageCoarse"))
submetaDf <- metaDf %>%
  filter(grepl("Exh", annoCoarse))

p6a_p3 <- ggplot() +
  # geom_point(data = metaDf, aes(x=UMAP_1, y=UMAP_2, color = annoCoarse), 
  #            size=1.5, alpha=0.1) +
  geom_point(data = metaDf,  aes(x=UMAP_1, y=UMAP_2), size=1.5, color="grey90") +
  stat_density_2d(data = submetaDf, aes(x=UMAP_1, y=UMAP_2, color=annoCoarse), 
                  contour_var = "ndensity", linewidth=1) +
  xlab("UMAP1")+ylab("UMAP2")+
  theme_pubr() +
  theme(axis.text = element_text(size=20,color = 'black'),
        axis.title = element_text(size=25,face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values = colset8$annoCoarse) +
  guides(color=guide_legend(override.aes = list(size=8), nrow = 1, keywidth = unit(1.5, "cm")))
saveFig(p6a_p3, 
        prefix = "Plots/Fig6a_p3_DensityPlot_on_UMAP_Exh_Only", 
        wd = 3600, ht = 3200)

##============================-- Panel 6b EX Label Transfer ----
## part 2: barplot of JW label transfer
metaDf <- FetchData(CAR8, vars = c("annoCoarse", "projTILHuman")) %>%
  rename(anno=annoCoarse, label = projTILHuman)

labCols <- pal_npg()(9)[c(4,2,6,5,1,9,3)]
names(labCols) <- c("CD8.NaiveLike", "CD8.CM", "CD8.TPEX", 
                    "CD8.EM", "CD8.TEMRA", "CD8.TEX", "CD8.MAIT")

dfx <- metaDf %>%
  drop_na() %>%
  group_by(anno) %>%
  dplyr::count(label) %>%
  mutate(Proportion = n/sum(n), Freq=n, Total=sum(n)) %>%
  ungroup() %>%
  mutate(label = factor(label, levels = c("CD8.NaiveLike", "CD8.CM", "CD8.TPEX", 
                                          "CD8.EM", "CD8.TEMRA", "CD8.TEX", "CD8.MAIT"))) %>%
  mutate(anno = factor(anno, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                        "EM_Eff-like", "CM", "Prolif")))

p6b <- ggplot(dfx, aes(x=anno, y=Proportion, fill=label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  theme_bw() + 
  labs(x="CD8 subsets") +
  coord_flip() +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values = labCols) +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.y = element_text(size=18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth=1, color = "black")) +
  guides(fill=guide_legend(override.aes = list(size=5), nrow = 2))
saveFig(p6b, 
        prefix = "Plots/Fig6b_Label_Transfer_distrib_overall", 
        wd = 2000, ht = 1200)


##============================-- Panel 6c EX gene density plot ----
require(Nebulosa)
row.names(CAR8@assays$ADT@data) <- gsub("-TotalSeqC", "", row.names(CAR8@assays$ADT@data))
row.names(CAR8@assays$ADT@scale.data) <- gsub("-TotalSeqC", "", row.names(CAR8@assays$ADT@scale.data))

keyMarkers <- c("YESCAR", "KLRG1", "TBX21", 
                "MKI67","IL7R", "CCR7", 
                "TCF7", "ADT_CD45RA", "GZMB", "PRF1",
                "TOX", "TIGIT", "NR4A2", "PDCD1","CTLA4",
                "HAVCR2", "LAG3", "ENTPD1", "IRF4")

keyMarkers <- c("CD27", "B3GAT1", "ADT_CD57", "ADT_CD39", "ADT_CD27", "ADT_CCR7", "ADT_CD127", 
                "ADT_TIM-3", "ADT_LAG-3", "ADT_PD-1", "ADT_TIGIT", "ADT_CD103", "IL7R", "LEF1", 
                "CX3CR1", "NKG7", "GNLY", "IFNG")
for (gx in keyMarkers){
  p6c <- plot_density(CAR8, gx)
  saveFig(p6c, prefix = paste0("Plots/Density_Plots/Fig6c_PI_CAR8_density_markers_", gx),
          wd = 1800, ht = 1500)
}


##============================-- Panel 6d EX DEGs ----
EMDEG <- FindMarkers(CAR8, ident.1 = "EM_Exh-like", ident.2 = "EM_Eff-like", group.by = "annoCoarse", logfc.threshold = 0.01)
EFFDEG <- FindMarkers(CAR8, ident.1 = "EFF_Exh-like", ident.2 = "EFF", group.by = "annoCoarse", logfc.threshold = 0.01)
EMFDEG <- FindMarkers(CAR8, ident.1 = "EM_Exh-like", ident.2 = "EFF_Exh-like", group.by = "annoCoarse", logfc.threshold = 0.01)
degRes <- list(EMex = EMDEG, EFFex = EFFDEG, EMvEFFex = EMFDEG)
saveRDS(degRes, "RDS/DEG_CAR8_EX_May06.rds")


## part 1: EM_EX vs EM ----
require(ggrepel)
EMDEG <- degRes$EMex
volcano_col <- c("#A66D3C", "#B7B545", 'grey80')
EX_genes <- c("CXCR4", "NFKBIA", "DUSP4", "NR4A2", "BCL3", "RUNX3", "CD69", "TOX", "CTLA4",
               "RELB", "PRDM1", "CD38", "HLA-DRB5", "JUNB", "HIF1A", "BHLHE40", "HLA-DRA")
EFF_genes <- c("S100A11", "LTB", "CX3CR1", "KLRD1", "IL7R", "STAT1", "CXCR3", "GNLY", "PRF1", 
               "GZMB", "S1PR1", "GZMM", "KLRB1", "TBX21", "FCGR3A", "IL2RG", "ID2", "KIR2DL3", "KLRG1", 
               "IFITM1", 'IFITM3', "ID2", "CXCR6")
gene_highlight <- c(EX_genes, EFF_genes)

## preprocessing
res <- EMDEG %>%
  rownames_to_column(var = "gene") %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(logp = ifelse(logp>100, 100, logp)) %>%
  mutate(avg_log2FC = clipMtx(avg_log2FC, 1.2)) %>%
  mutate(category = case_when(
    logp > 3 & avg_log2FC > 0.25 ~ "Up in EM_Exh-like",
    logp > 3 & avg_log2FC < -0.25 ~ "Up in EM_Eff-like",
    TRUE ~ "Non-DEGs"
  )) %>%
  mutate(gene_highlight = ifelse(gene %in% gene_highlight, gene, "")) %>%
  mutate(gene_highlight = ifelse(category=="Non-DEGs", "", gene_highlight)) %>%
  mutate(category = factor(category, levels = c("Up in EM_Exh-like", "Up in EM_Eff-like", "Non-DEGs")))

res_Highlight <- subset(res, gene_highlight!="") 

## plot
p6d_p1 <- ggplot(res, aes(x=avg_log2FC, y=logp, color=category, label=gene_highlight))+
  geom_point(size=1) +
  geom_point(data = res_Highlight, size=2, color="black") +
  geom_text_repel(color="black", min.segment.length = 0, max.overlaps = Inf, size=3, force = 1.5) +
  geom_vline(xintercept = 0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_vline(xintercept = -0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_hline(yintercept = 3, linetype="dashed", color="black", linewidth=0.5) +
  theme_pubr()+
  labs(x=bquote(~Log[2]~ 'Fold Change'), y=bquote('-'~Log[10]~P['adj']))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_blank())+
  scale_color_manual(values = volcano_col)+
  scale_y_continuous(expand = expansion(add=30))+
  guides(color = guide_legend(override.aes = list(size=6))) 
saveFig(p6d_p1, 
        prefix = "Plots/Fig6d_p1_volcanoPlot_DEG_EM_EX_vs_EFF_inCAR8", 
        wd = 3000, ht = 2000)


## part 2: EFF_EX vs EFF ----
require(ggrepel)
EFFDEG <- degRes$EFFex
volcano_col <- c( "#5489A3", "#D94436", 'grey80')
EX_genes <- c("TIGIT", "LAG3", "RELB", "CD7", "NFKBIA", "KLRD1", "FCGR3A", "TOX", "CTLA4",
              "HAVCR2", "CD27", "CCL3", "ID2", "BACH2", "BHLHE40", "CXCR6", "EOMES", "IFNG",
              "PDCD1", "CREM", "DUSP5", "BCL3")
EFF_genes <- c("ASCL2", "KLF3", "GZMH", "TGFB1", "LAT", "IFITM1", "LY6E", "PRF1", "S100A10",
               "KLF13", "CD81", "S100A6", "GZMA", "GZMM")
gene_highlight <- c(EX_genes, EFF_genes)

## preprocessing
res <- EFFDEG %>%
  rownames_to_column(var = "gene") %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(logp = ifelse(logp>100, 100, logp)) %>%
  mutate(avg_log2FC = clipMtx(avg_log2FC, 1.2)) %>%
  mutate(category = case_when(
    logp > 3 & avg_log2FC > 0.25 ~ "Up in EFF_Exh-like",
    logp > 3 & avg_log2FC < -0.25 ~ "Up in EFF",
    TRUE ~ "Non-DEGs"
  )) %>%
  mutate(gene_highlight = ifelse(gene %in% gene_highlight, gene, "")) %>%
  mutate(gene_highlight = ifelse(category=="Non-DEGs", "", gene_highlight)) %>%
  mutate(category = factor(category, levels = c("Up in EFF_Exh-like", "Up in EFF", "Non-DEGs")))

res_Highlight <- subset(res, gene_highlight!="") 

## plot
p6d_p2 <- ggplot(res, aes(x=avg_log2FC, y=logp, color=category, label=gene_highlight))+
  geom_point(size=1) +
  geom_point(data = res_Highlight, size=2, color="black") +
  geom_text_repel(color="black", min.segment.length = 0, max.overlaps = Inf, size=3, force = 1.5) +
  geom_vline(xintercept = 0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_vline(xintercept = -0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_hline(yintercept = 3, linetype="dashed", color="black", linewidth=0.5) +
  theme_pubr()+
  labs(x=bquote(~Log[2]~ 'Fold Change'), y=bquote('-'~Log[10]~P['adj']))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_blank())+
  scale_color_manual(values = volcano_col)+
  scale_y_continuous(expand = expansion(add=30))+
  guides(color = guide_legend(override.aes = list(size=6))) 
saveFig(p6d_p2, 
        prefix = "Plots/Fig6d_p2_volcanoPlot_DEG_EFF_EX_vs_EFF_inCAR8", 
        wd = 3000, ht = 2000)


## part 3: EM_EX vs EFF_EX ----
require(ggrepel)
EMFDEG <- degRes$EMvEFFex
volcano_col <- c("#A66D3C", "#5489A3", 'grey80')
EM_genes <- c("GZMK", "CXCR4","DUSP4", "RGS1", "CD69", "TCF7", "TNFAIP3", "SELL",
              "CD74", "NFKBIA", "JUNB", "TGFB1", "NR4A2", "LDHB", "ISG20", "CD44", 
              "KLF13", "KLF3", "S1PR4", "RUNX2", "LEF1", "ZEB1")
EFF_genes <- c("KLRD1", "FCGR3A", "CCL4", "CX3CR1", "GZMB", "TIGIT", "LAG3", "KIR2DL3",
               "S1PR5", "HAVCR2", "KLRG1", 'TBX21', "ID2", "KLRB1", "IL2RG", "ZEB2", "NKG7",
               "PRF1", "STAT1", "CCL3", "TOX", "IRF2", "GZMH", "BCL2", "KLF6", "PDCD1", "BATF",
               "GNLY", "EOMES", "IRF7", "IFNG", "IKZF3", "KLF2")
gene_highlight <- c(EM_genes, EFF_genes)

## preprocessing
res <- EMFDEG %>%
  rownames_to_column(var = "gene") %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(logp = ifelse(logp>100, 100, logp)) %>%
  mutate(avg_log2FC = clipMtx(avg_log2FC, 1.2)) %>%
  mutate(category = case_when(
    logp > 3 & avg_log2FC > 0.25 ~ "Up in EM_Exh-like",
    logp > 3 & avg_log2FC < -0.25 ~ "Up in EFF_Exh-like",
    TRUE ~ "Non-DEGs"
  )) %>%
  mutate(gene_highlight = ifelse(gene %in% gene_highlight, gene, "")) %>%
  mutate(gene_highlight = ifelse(category=="Non-DEGs", "", gene_highlight)) %>%
  mutate(category = factor(category, levels = c("Up in EM_Exh-like", "Up in EFF_Exh-like", "Non-DEGs")))

res_Highlight <- subset(res, gene_highlight!="") 

## plot
p6d_p3 <- ggplot(res, aes(x=avg_log2FC, y=logp, color=category, label=gene_highlight))+
  geom_point(size=1) +
  geom_point(data = res_Highlight, size=2, color="black") +
  geom_text_repel(color="black", min.segment.length = 0, max.overlaps = Inf, size=3, force = 1.5) +
  geom_vline(xintercept = 0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_vline(xintercept = -0.25, linetype="dashed", color="black", linewidth=0.5) +
  geom_hline(yintercept = 3, linetype="dashed", color="black", linewidth=0.5) +
  theme_pubr()+
  labs(x=bquote(~Log[2]~ 'Fold Change'), y=bquote('-'~Log[10]~P['adj']))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_blank())+
  scale_color_manual(values = volcano_col)+
  scale_y_continuous(expand = expansion(add=30))+
  guides(color = guide_legend(override.aes = list(size=6))) 
saveFig(p6d_p3, 
        prefix = "Plots/Fig6d_p3_volcanoPlot_DEG_EM_EX_vs_EFF_EX_inCAR8", 
        wd = 3000, ht = 2000)


##============================-- Panel 6e EX DEG pathway enrichment ----
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

## preload JW
JWdf <- read.csv("Data/JWAllClusterDEG.csv")
JWList <- split(JWdf, f=JWdf$cluster)
sigSet <- lapply(JWList, function(dfx){
  dfx %>% 
    arrange(desc(avg_log2FC)) %>%
    head(500) %>%
    pull(gene) %>% 
    unique() %>% 
    HM_Convert(HtoM=F, simple = T) %>% 
    na.omit() %>%
    GetGenelist()
})
names(sigSet) <- paste0("JW_", gsub(" ", "_", names(sigSet)))
m_df_JW <- do.call("rbind", lapply(names(sigSet), function(sigx){
  data.frame(gs_name = rep(sigx, length(sigSet[[sigx]])), 
             entrez_gene = as.integer(sigSet[[sigx]])) %>%
    as_tibble()
}))

m_df <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  filter(grepl("EXH", gs_name))


## part 1: EM EX vs EFF ----
gvec <- GetGenelist(degRes$EMex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

glist<- c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP",
          "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN",
          "GSE41867_DAY6_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP")

ID_alias = c("EX vs EFF UP (GSE9650)", 
             "MEM vs EX UP (GSE9650)",
             "EFF vs EX UP (GSE41867)")

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6e_p1 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6e_p1, paste0("Plots/GSEA/Fig6e_p1_GSEA_EM_EX_vs_EFF_", pathx),
          wd = 1800, ht = 1200)
}


## part 2: EFF EX vs EFF ----
gvec <- GetGenelist(degRes$EFFex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

glist<- c("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN")

ID_alias = c("EX vs EFF UP (GSE9650)")

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6e_p2 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6e_p2, paste0("Plots/GSEA/Fig6e_p2_GSEA_EFF_EX_vs_EFF_", pathx),
          wd = 1800, ht = 1200)
}

## part 3: EM EX vs EFF EX ----
gvec <- GetGenelist(degRes$EMvEFFex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

glist<- c("GSE41867_DAY6_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP")

ID_alias = c("EFF vs EX UP (GSE41867)")

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6e_p3 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6e_p3, paste0("Plots/GSEA/Fig6e_p3_GSEA_EM_EX_vs_EFF_EX_", pathx),
          wd = 1800, ht = 1200)
}


##============================-- Panel 6f EX DEG pathway enrichment: JW Subsets ----
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

## preload JW
JWList <- readRDS("RDS/Ref_JW_CAR8_deg_List_pairwise_essential_May07_2023.rds")
sigSetUP <- lapply(JWList, function(dfx){
  dfx %>% 
    rownames_to_column(var = "gene") %>%
    filter(avg_log2FC>0) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene) %>% 
    GetGenelist() %>%
    unique()
})
names(sigSetUP) <- paste0(names(sigSetUP), "_UP")

sigSetDN <- lapply(JWList, function(dfx){
  dfx %>% 
    rownames_to_column(var = "gene") %>%
    filter(avg_log2FC<0) %>%
    arrange(avg_log2FC) %>%
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
gvec <- GetGenelist(degRes$EMex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, 
                pAdjustMethod = "fdr", minGSSize = 5, maxGSSize = Inf)

glist<- c("Trans I_vs_Exh-Prog_UP", "Trans I_vs_Exh-Prog_DN",
          "Exh-Pre_vs_Trans I_UP", "Exh-Pre_vs_Trans I_UP")

ID_alias = glist

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6f_p1 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6f_p1, paste0("Plots/GSEA/Fig6f_p1_GSEA_EM_EX_vs_EFF_", pathx),
          wd = 1800, ht = 1200)
}


## part 2: EFF EX vs EFF ----
gvec <- GetGenelist(degRes$EFFex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

dfx <- pathway@result %>%
  mutate(group = gsub("_UP", "", ID)) %>%
  mutate(group = gsub("_DN", "", group)) %>%
  group_by(group) %>%
  summarise(meanNES = mean(abs(NES)))

glist<- c("Trans CTL_vs_Exh-Term_UP", "Trans CTL_vs_Exh-Term_DN",
          "Eff_vs_Exh-Term_UP", "Eff_vs_Exh-Term_UP",
          "CTL_vs_Exh-Int_DN", "CTL_vs_Exh-Int_UP",
          "CTL_vs_Exh-Pre_UP")

ID_alias = glist

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6f_p2 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6f_p2, paste0("Plots/GSEA/Fig6f_p2_GSEA_EFF_EX_vs_EFF_", pathx),
          wd = 1800, ht = 1200)
}

## part 3: EM EX vs EFF EX ----
gvec <- GetGenelist(degRes$EMvEFFex)
pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

dfx <- pathway@result %>%
  mutate(group = gsub("_UP", "", ID)) %>%
  mutate(group = gsub("_DN", "", group)) %>%
  group_by(group) %>%
  summarise(meanNES = mean(abs(NES)),
            Nterm = n())

glist<- c("Exh-Term_vs_Exh-Prog_UP", "Exh-Term_vs_Exh-Prog_DN",
          "CTL_vs_Exh-KLR_UP", "CTL_vs_Exh-KLR_UP")

ID_alias = glist

for (idx in 1:length(glist)){
  pathx <- glist[idx]
  pathAx <- ID_alias[idx]
  
  p6f_p3 <- plotGSEA(pathway, pathx, pTitle = pathAx)
  saveFig(p6f_p3, paste0("Plots/GSEA/Fig6f_p3_GSEA_EM_EX_vs_EFF_EX_", pathx),
          wd = 1800, ht = 1200)
}


##============================-- Panel 6g ATAC UMAP ----
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



