##============================-- 1. Set up Environment ----
## 1.1 load packages ----
rm(list=ls())

pkgs <- c("Seurat", "tidyverse", "pals", "immunarch",
          "ggpubr", "org.Hs.eg.db", "ggsci", "clusterProfiler", 
          "pheatmap", "devtools", "R.utils", "ggh4x", "ggalluvial", 
          "ggpubr", "ggprism", "colorspace")

lapply(pkgs, require, character.only = TRUE)

## 1.2 load auxiliary script ----
source("Scripts/Utility/RNA_Utility.R")

## 1.4 set system colors ----
colset8 <- readRDS("RDS/color_set_CD8.rds")

## 1.5 load CAR4 and CAR8 seurat objects
CAR8 <- readRDS("RDS/CAR8_for_manuscript_Apr262023.rds")


##============================-- Panel 4a Tile Map Stage-specific Gene Expression ----
sigSetS2 <- list(
  Inflammation = c("NFKBIA", "NFKBIE","NFKB1","PDE4D", "DUSP4", "HMGA1", "PIM3", "ZNF580",
                   "CNN2", "CUL3", "NR4A2", "BCL3", "CXCR4", "HIF1A", "ATF4",
                   "HLA-DRA", "HLA-DRB1","CD74", "RHOH", "CD69", "DUSP5", "RELB", "LDHB",
                   "CD44", "RUNX3", "GATA3", "PRDM1", "CCL5", "CD83"),
  TIL = c("PELI1", "RGS1", "PER1", "CTLA4"),
  Immu_Inhib=c("CD300A", "CISH", "GPR171", "KIR2DL3", "IRF2", "ADRB2"),
  IFN_I = c("IRF7", "STAT1", "TRIM22", "IFITM1", "IFITM2"),
  T_Eff = c("CX3CR1", "LAIR2", "IL18RAP", "PNP", "SH2D2A", "KLRD1", "KLRB1", "FCGR3A",
            "PRF1", "GZMB", "KIR3DL1", "ZEB2", "GZMM", "NKG7"),
  T_Migra = c("RGS3", "SEMA4C", "S1PR5", "S1PR1"),
  Others = c("ABI3", "SAMD9L", "SAMD9", "LBH","IGF2R", "TGFBR1", "GIMAP7", "GIMAP5")
)

### 2.2.2 patient level plot ----
geneSelect <- unname(unlist(sigSetS2))

gsDf <- do.call("rbind", lapply(names(sigSetS2), function(sx){
  ggx <- sigSetS2[[sx]]
  data.frame(Gene=ggx, gLab = rep(sx, length(ggx)))
}))

CAR8 <- ScaleData(CAR8,features = geneSelect, 
                  vars.to.regress = c("percent.mt", "CC.Difference"))

sigDf <- FetchData(CAR8, vars = c(geneSelect, "patientShort", "stageCoarse", "annoCoarse"), slot = "scale.data")

sigAvg <- sigDf %>%
  pivot_longer(cols = !c("annoCoarse", "patientShort", "stageCoarse"), names_to = "Gene", values_to = "Expr") %>%
  rename(patient = patientShort) %>%
  group_by(Gene, annoCoarse, patient, stageCoarse) %>%
  summarise(AvgExpr = mean(Expr, na.rm=T)) %>%
  left_join(gsDf, by = "Gene") %>%
  mutate(AvgExpr = clipVec(AvgExpr, 1, lower = F)) %>%
  mutate(AvgExpr = clipVec(AvgExpr, -1, lower = T)) %>%
  mutate(Gene = factor(Gene, levels = geneSelect)) %>%
  mutate(gLab = factor(gLab, levels = names(sigSetS2))) %>%
  mutate(annoCoarse = factor(annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                    "EM_Eff-like", "CM", "Prolif"))) %>%
  mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient)))) %>%
  ungroup()

fancy_strips <- strip_nested(
  background_y = elem_list_rect(fill= c(colset8$stageCoarse[2:3], rep(colset8$annoCoarse, 2)),
                                color = rep("white", 10)),
  text_y = elem_list_text(face="bold",
                          size=10,
                          angle = 0,
                          color= "white"),
  text_x = elem_list_text(size=12, face="bold")
)


p4a_p1 <- ggplot(sigAvg, aes(y=patient, x=Gene, fill = AvgExpr)) +
  geom_tile(color = "black", lwd = 0.2, linetype = 1) +
  theme_minimal() +
  labs(fill = "Normalized\nScore")+
  facet_nested(stageCoarse + annoCoarse ~ gLab, scales = "free", space = "free",
               strip = fancy_strips, switch = "y") +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", mid=0) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=60, hjust=1, 
                                   vjust=1, colour = "black"),
        axis.text.y = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
saveFig(p4a_p1, "Plots/Fig4a_p1_Heatmap_Exp_vs_Per_DEG_patient_Lvl", 
        wd = 4000, ht = 1800)


### 2.2.3 stage level plot ----
geneSelect <- unname(unlist(sigSetS2))

gsDf <- do.call("rbind", lapply(names(sigSetS2), function(sx){
  ggx <- sigSetS2[[sx]]
  data.frame(Gene=ggx, gLab = rep(sx, length(ggx)))
}))

CAR8 <- ScaleData(CAR8,features = geneSelect, 
                  vars.to.regress = c("percent.mt", "CC.Difference"))

sigDf <- FetchData(CAR8, vars = c(geneSelect, "stageCoarse", "annoCoarse"), slot = "scale.data")

sigAvg <- sigDf %>%
  pivot_longer(cols = !c("annoCoarse", "stageCoarse"), names_to = "Gene", values_to = "Expr") %>%
  group_by(Gene, annoCoarse, stageCoarse) %>%
  summarise(AvgExpr = mean(Expr, na.rm=T)) %>%
  left_join(gsDf, by = "Gene") %>%
  mutate(AvgExpr = clipVec(AvgExpr, 0.7, lower = F)) %>%
  mutate(AvgExpr = clipVec(AvgExpr, -0.7, lower = T)) %>%
  mutate(Gene = factor(Gene, levels = geneSelect)) %>%
  mutate(gLab = factor(gLab, levels = names(sigSetS2))) %>%
  mutate(annoCoarse = factor(annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                    "EM_Eff-like", "CM", "Prolif"))) %>%
  ungroup()

fancy_strips <- strip_nested(
  background_y = elem_list_rect(fill= c(colset8$stageCoarse[2:3]),
                                color = "white"),
  text_y = elem_list_text(face="bold",
                          size=10,
                          angle = 0,
                          color= "white"),
  text_x = elem_list_text(size=12, face="bold")
)


p4a_p2 <- ggplot(sigAvg, aes(y=annoCoarse, x=Gene, fill = AvgExpr)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  theme_minimal() +
  labs(fill = "Normalized\nScore")+
  facet_grid2(stageCoarse ~ gLab, scales = "free", space = "free",
              strip = fancy_strips, switch = "y") +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", mid=0) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=60, hjust=1, 
                                   vjust=1, colour = "black"),
        axis.text.y = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.placement = "outside") +
  guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
saveFig(p4a_p2, "Plots/Fig4a_p2_Heatmap_Exp_vs_Per_DEG_stage_Lvl", 
        wd = 3000, ht = 1500)



##============================-- Panel 4b Summary of Stage-specific GSEA ----
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

# mList <- readRDS("RDS/SingleCell_DEG_S1_vs_S2_Seurat_List_Apr262023.rds")
# mStageSC2 <- mList$Overall
mStageSC <- readRDS("RDS/SingleCell_DEG_S1_vs_S2_Raw_May112023.rds")

msigCats <- c("C2", "H", "C5", "C7")
mdfList <- lapply(msigCats, function(catx){
  msigdbr(species = "Homo sapiens", category = catx) %>% 
    dplyr::select(gs_name, entrez_gene)
})
names(mdfList) <- msigCats

## part 1: specific GSEA ----
### Single cell version
gvec <- GetGenelist(mStageSC)

gseaListSC <- lapply(mdfList, function(m_dfx){
  GSEA(gvec, TERM2GENE = m_dfx, pvalueCutoff = 0.5, pAdjustMethod = "fdr")
})


selectListSC <- list(C2 = c("PHONG_TNF_TARGETS_UP", "WP_IL1_SIGNALING_PATHWAY", "TIAN_TNF_SIGNALING_VIA_NFKB",
                    "ALCALA_APOPTOSIS", "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", "HECKER_IFNB1_TARGETS",
                    "MOSERLE_IFNA_RESPONSE", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "BOSCO_TH1_CYTOTOXIC_MODULE",
                    "BIOCARTA_CTL_PATHWAY", "PID_CD8_TCR_PATHWAY"),
             H = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
             C5 = c("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
                    "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", "GOBP_RESPONSE_TO_INTERFERON_BETA", 
                    "GOCC_T_CELL_RECEPTOR_COMPLEX", "GOBP_CELL_KILLING", "GOBP_RESPONSE_TO_TYPE_I_INTERFERON"))

for (catx in names(selectList)){
  glist <- selectListSC[[catx]]
  pathwayx <- gseaListSC[[catx]]
  for (idx in 1:length(glist)){
    pathx <- glist[idx]
    pathAx <- pathx
    
    p4b_p1 <- plotGSEA(pathwayx, pathx, pTitle = pathAx)
    saveFig(p4b_p1, paste0("Plots/GSEA/Fig4b_p1_GSEA_SC_Exp_vs_Per_", pathx),
            wd = 1800, ht = 1200)
  }
}

### pseudoBulk Version
gvec <- GetGenelist(mStagePsudo)

gseaListBulk <- lapply(mdfList, function(m_dfx){
  GSEA(gvec, TERM2GENE = m_dfx, pvalueCutoff = 0.5, pAdjustMethod = "fdr")
})


selectListBulk <- list(C2 = c("REACTOME_CELL_CYCLE", "BENPORATH_CYCLING_GENES", "LEE_EARLY_T_LYMPHOCYTE_UP",
                          "BENPORATH_PROLIFERATION", "GOLDRATH_ANTIGEN_RESPONSE", "ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP",
                          "KAUFFMANN_DNA_REPAIR_GENES","REACTOME_DNA_REPAIR", "PHONG_TNF_RESPONSE_VIA_P38_PARTIAL"),
                   H = c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
                   C5 = c("GOBP_CELL_CYCLE", "GOBP_DNA_REPLICATION", "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY", 
                          "GOBP_POSITIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"))

for (catx in names(selectList)){
  glist <- selectListBulk[[catx]]
  pathwayx <- gseaListBulk[[catx]]
  
  for (idx in 1:length(glist)){
    pathx <- glist[idx]
    pathAx <- pathx
    
    p4b_p1 <- plotGSEA(pathwayx, pathx, pTitle = pathAx)
    saveFig(p4b_p1, paste0("Plots/GSEA/Fig4b_p1_GSEA_PseudoBulk_Exp_vs_Per_", pathx),
            wd = 1800, ht = 1200)
  }
}

## part 2 Overall Summary ----
### Single Cell Version
stageGSEA <-  do.call("rbind", lapply(gseaListSC, function(sobj) {
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  sobj@result
}))

select_terms <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                  "PHONG_TNF_TARGETS_UP", 
                  "WP_IL1_SIGNALING_PATHWAY", 
                  "ALCALA_APOPTOSIS", 
                  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", 
                  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
                  "BIOCARTA_CTL_PATHWAY", 
                  "PID_CD8_TCR_PATHWAY",
                  "HALLMARK_G2M_CHECKPOINT", 
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", 
                  "GOBP_CELL_KILLING", 
                  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON")

selectDf <- stageGSEA %>%
  filter(ID %in% select_terms) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels=rev(ID))) %>%
  mutate(logq = -log10(qvalue)) %>%
  mutate(enrichment = ifelse(NES>0, "Exp", "Per")) %>%
  mutate(enrichment = factor(enrichment, levels=c("Exp", "Per")))


p4b_p2 <- ggplot(selectDf, aes(x=ID, y=NES)) +
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
saveFig(p4b_p2, 
        prefix = "Plots/Fig4b_p2_GSEA_Summary_Exp_vs_Per_SingleCell", 
        wd = 2800, ht = 2000)

### Pseudo-Bulk Version
stageGSEA <-  do.call("rbind", lapply(gseaListBulk, function(sobj) {
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  sobj@result
}))

select_terms <- c("REACTOME_CELL_CYCLE", 
                  "LEE_EARLY_T_LYMPHOCYTE_UP",
                  "GOLDRATH_ANTIGEN_RESPONSE", 
                  "ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP",
                  "KAUFFMANN_DNA_REPAIR_GENES",
                  "REACTOME_DNA_REPAIR", 
                  "PHONG_TNF_RESPONSE_VIA_P38_PARTIAL",
                  "HALLMARK_G2M_CHECKPOINT", 
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  "GOBP_CELL_CYCLE", 
                  "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                  "GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_DN",
                  "GOBP_POSITIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS")

selectDf <- stageGSEA %>%
  filter(ID %in% select_terms) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels=rev(ID))) %>%
  mutate(logq = -log10(qvalue)) %>%
  mutate(enrichment = ifelse(NES>0, "Exp", "Per")) %>%
  mutate(enrichment = factor(enrichment, levels=c("Exp", "Per")))


p4b_p2 <- ggplot(selectDf, aes(x=ID, y=NES)) +
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
saveFig(p4b_p2, 
        prefix = "Plots/Fig4b_p2_GSEA_Summary_Exp_vs_Per_PseudoBulk", 
        wd = 2800, ht = 2000)







##============================-- Panel 4c Top Regulons ----
## load data and define functions ----
RegDataSet <- readRDS("RDS/Fig4_RegDataSet_May082023.rds")
list2env(RegDataSet, .GlobalEnv)

plotRegulonNetwork <- function(TFx, degDf, regList, TFtable, vUp=0.2, labSize=5, toRepel=F){
  require(ggraph)
  require(igraph)
  require(scales)
  require(tidygraph)
  
  ## names of all human TFs
  TFlist <- TFtable$Symbol
  
  ## regulon genes
  reglist <- unlist(regList, recursive = F)
  names(reglist) <- gsub("\\.", "_", names(reglist))
  tfx <- strsplit(TFx, "_")[[1]][2]
  regVec <- unique(c(reglist[[TFx]],tfx))

  ## filter degDF to keep only high DEGs
  subdegDf <- degDf %>% 
    filter(gene %in% regVec)  %>%
    mutate(isTF = ifelse(gene ==tfx, "aTF", 
                         ifelse(gene %in% TFlist, "bTF", "nonTF"))) %>%
    arrange(isTF, desc(avg_log2FC)) %>%
    slice_head(n=11)
    
  if (nrow(subdegDf)<2){
    return(NA)
  } else {
    TFtargets <-  subdegDf %>% 
      filter(isTF!="aTF") %>%
      pull(gene)
    
    lfcmax <- max(abs(subdegDf$avg_log2FC))
    TFfullRange <- c(-lfcmax,lfcmax)
    
    ggdf <- data.frame(from=rep(tfx, length(TFtargets)), to = TFtargets)
    allGenes <- unique(c(ggdf$from, ggdf$to))
    
    vvdf <- degDf %>%
      filter(gene %in% allGenes) %>%
      rename(node=gene)
    
    if (!all(allGenes %in% vvdf$node)){
      dfmat <- do.call("rbind", lapply(setdiff(allGenes, vvdf$node), function(genex){
        data.frame(node = genex, p_val = 1, avg_log2FC = 0, pct.1 = 0, pct.2 = 0, p_val_adj = 1)
      }))
      vvdf <- bind_rows(vvdf, dfmat)
    }
    
    gg <- graph_from_data_frame(ggdf, directed = T, vertices = vvdf)
    gg <- as_tbl_graph(gg)
    gg <- gg %>% mutate(degree=centrality_degree())
    glayout <- create_layout(gg, layout = "stress")
    
    p2d_temp <- ggraph(glayout)+
      geom_edge_link(color="grey40",alpha=0.9, width=0.2)+
      geom_node_point(aes(fill=avg_log2FC), size=8, shape=21, color="black", stroke=1)+
      geom_node_text(aes(label=name), size=labSize, color="black", repel = toRepel, nudge_y=vUp, 
                     hjust="inward") +
      labs(color=NULL, size=NULL) +
      theme(legend.title = element_text(size=20),
            legend.text = element_text(size=15),
            legend.position = "top",
            panel.background = element_rect(fill = "transparent", color="transparent")) +
      scale_fill_gradient2(low=colset8$stageCoarse["Per"], 
                           mid = "#F7F7F7", 
                           high = colset8$stageCoarse["Exp"],
                           limits=TFfullRange) +
      guides(size="none", color="none", fill=guide_colorbar(barwidth = unit(3, "inches"),
                                                            title.vjust = 0.9,
                                                            frame.colour = "black"))
    return(p2d_temp)
  }
}
getTopRegulons <- function(modelObj, topN){
  ## load packages
  require(xgboost)
  require(caret)
  require(Ckmeans.1d.dp)
  require(SHAPforxgboost)
  
  ##-- 3.0.1 soft load relevant data for a BP
  bst <- modelObj$model$finalModel
  cell_trn <- modelObj$trainDf
  cols <- setdiff(colnames(cell_trn), c("Anno"))
  cell_trn_mtx <- data.matrix(cell_trn[, cols])
  
  ##-- 3.0.2 get SHAP values
  shap_values <- shap.values(xgb_model = bst, X_train = cell_trn_mtx)
  shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = cell_trn_mtx)
  shapvDf <- as.data.frame(shap_values$shap_score)
  
  ##-- 3.0.3 get prediction direction
  cell_trn_df <- as.data.frame(cell_trn_mtx)
  shap_force <- unlist(lapply(colnames(shapvDf), function(tfx){
    cor(shapvDf[[tfx]], cell_trn_df[[tfx]])
  }))
  
  ##-- 3.0.4 summarize mean shap df
  shap_mean <- data.frame(regulon=colnames(shapvDf), 
                          meanSHAP=shap_values$mean_shap_score[colnames(shapvDf)],
                          meanForce=shap_force)
  shap_mean$dir <- ifelse(!is.na(shap_mean$meanForce),
                          ifelse(shap_mean$meanForce>0, "Exp", "Per"),
                          NA)
  
  ##-- 3.0.5 prepare Df for plotting
  dfx <- shap_mean %>%
    drop_na() %>%
    group_by(dir) %>%
    arrange(desc(meanSHAP)) %>%
    slice_head(n=topN) %>%
    ungroup()
  
  return(dfx)
}

## part 1: Top Regulon Networks ----
reg8R <- getTopRegulons(xgblist, 30)
reg8RFilt <- reg8R %>%
  group_by(dir) %>%
  filter(grepl(unique(dir), regulon)) %>%
  arrange(desc(meanSHAP)) %>%
  slice_head(n=10) %>%
  ungroup()
  

tfSelect <- reg8RFilt$regulon
for (tx in tfSelect){
  p4c_p1 <- plotRegulonNetwork(TFx = tx, degDf = degDf, regList = regList, TFtable = TFtable)
  if (!any(is.na(p4c_p1))){
    saveFig(pobj = p4c_p1, 
            prefix = paste0("Plots/RegulonNet/Fig4c_p1_NetPlot_ExpvPer_CAR8_", tx), 
            wd = 2400, ht = 1500) 
  }
}

## part 2: Top Regulon Barplot summary ----
p4c_p2 <- ggplot(reg8RFilt, aes(x=meanSHAP, y=reorder(regulon, meanSHAP), fill=dir)) +
  geom_bar(stat = "identity") +
  theme_pubr() +
  scale_fill_manual(values = colset8$stageCoarse[c(2:3)]) +
  labs(x="Mean(|SHAP|)", y="", fill="") +
  facet_wrap(~dir, nrow = 1, scales = "free") +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15),
        strip.text = element_text(size=15, hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = "cm"))
saveFig(pobj = p4c_p2, 
        prefix = "Plots/Fig4c_p2_SHAP_Barplot_ExpvPer_CAR8", 
        wd = 2800, ht = 1500) 

## part 3: Top Regulon VlnPlot ----
metaDf <- FetchData(CAR8, vars = c("stageCoarse", "annoCoarse")) %>%
  rownames_to_column(var = "cell")
annoDf <- annoDf %>%
  rownames_to_column(var = "cell")
annoDfMeta <- annoDf %>%
  left_join(metaDf, by = "cell")


for (gx in reg8RFilt$regulon){
  metaDfx <- annoDfMeta %>%
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
  
  p4c_p3 <- ggviolin(metaDfx, x="stageCoarse", y="SigScore", facet.by = "annoCoarse", add = "boxplot",
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
  saveFig(p4c_p3, 
          paste0("Plots/SigVlnPlot2/Fig4c_p3_Regulon_Signature_Plot_annoCoarse_", gx,"_Split"),
          wd = 3000, ht = 1200)
}

## part 4: Top Regulon VlnPlot: per patient ----
metaDf <- FetchData(CAR8, vars = c("stageCoarse", "annoCoarse", "patientShort")) %>%
  rownames_to_column(var = "cell")
annoDf <- annoDf %>%
  rownames_to_column(var = "cell")
annoDfMeta <- annoDf %>%
  left_join(metaDf, by = "cell")



for (gx in reg8RFilt$regulon){
  metaDfx <- annoDfMeta %>%
    select(stageCoarse, patientShort, all_of(gx)) %>%
    rename(SigScore = all_of(gx)) %>%
    filter(SigScore > (mean(SigScore) - 3*IQR(SigScore))) %>%
    mutate(patientShort = factor(patientShort, levels = gtools::mixedsort(unique(patientShort))))
  
  p4c_p3 <- ggviolin(metaDfx, x="stageCoarse", y="SigScore", facet.by = "patientShort", add = "boxplot",
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
  saveFig(p4c_p3, 
          paste0("Plots/SigVlnPlot2/Fig4c_p3_Regulon_Signature_Plot_Patient_", gx,"_Split"),
          wd = 3000, ht = 1200)
}




## part 5: Plot top regulon number of genes and output table----
dfx <- do.call("rbind", lapply(reg8RFilt$regulon, function(TFx){
  TFlist <- TFtable$Symbol
  reglist <- unlist(regList, recursive = F)
  names(reglist) <- gsub("\\.", "_", names(reglist))
  tfx <- strsplit(TFx, "_")[[1]][2]
  regVec <- unique(c(reglist[[TFx]],tfx))
  
  fullDf <- data.frame(gene = regVec, regulon = TFx)
  
  subdegDf <- mStageSC %>% 
    filter(gene %in% regVec)  %>%
    full_join(fullDf, by = c("gene")) %>%
    mutate(isTF = ifelse(gene ==tfx, "core_TF", 
                         ifelse(gene %in% TFlist, "target_TF", "nonTF"))) %>%
    arrange(isTF, desc(avg_log2FC)) %>%
    mutate(regulon = TFx) %>%
    mutate(regSize = length(regVec)) %>%
    mutate(isDEG = ifelse(!is.na(avg_log2FC) & p_val_adj<0.05 & abs(avg_log2FC)>0.25, "DEG", "non-DEG"))
  
  subdegDf
}))
write.csv(dfx, "Plots/Fig4_SupTable_Top10_Regulon_Genes.csv", quote = F, row.names = F)

## extended version 
reg8R <- getTopRegulons(xgblist, 50)
reg8RFilt <- reg8R %>%
  group_by(dir) %>%
  filter(grepl(unique(dir), regulon)) %>%
  arrange(desc(meanSHAP)) %>%
  slice_head(n=30) %>%
  ungroup()
dfx <- do.call("rbind", lapply(reg8R$regulon, function(TFx){
  TFlist <- TFtable$Symbol
  reglist <- unlist(regList, recursive = F)
  names(reglist) <- gsub("\\.", "_", names(reglist))
  tfx <- strsplit(TFx, "_")[[1]][2]
  regVec <- unique(c(reglist[[TFx]],tfx))
  
  fullDf <- data.frame(gene = regVec, regulon = TFx)
  
  subdegDf <- mStageSC %>% 
    filter(gene %in% regVec)  %>%
    full_join(fullDf, by = c("gene")) %>%
    mutate(isTF = ifelse(gene ==tfx, "core_TF", 
                         ifelse(gene %in% TFlist, "target_TF", "nonTF"))) %>%
    arrange(isTF, desc(avg_log2FC)) %>%
    mutate(regulon = TFx) %>%
    mutate(regSize = length(regVec)) %>%
    mutate(isDEG = ifelse(!is.na(avg_log2FC) & p_val_adj<0.05 & abs(avg_log2FC)>0.25, "DEG", "non-DEG"))
  
  subdegDf
}))
write.csv(dfx, "Plots/Fig4_SupTable_Top30_Regulon_Genes.csv", quote = F, row.names = F)

reg8RFiltLvl <- reg8RFilt %>% 
  mutate(meanSHAP = ifelse(dir=="Exp", meanSHAP, -1 * meanSHAP)) %>%
  arrange(desc(meanSHAP))

dfx_dense <- dfx %>%
  select(regulon, regSize) %>%
  distinct() %>%
  mutate(regulon = factor(regulon, levels = reg8RFiltLvl$regulon))

p4c_p5 <- ggplot(dfx_dense, aes(x=regulon, y=regSize)) + 
  geom_bar(stat = "identity", fill="blue") +
  geom_text(aes(label=regSize), vjust=-0.25) +
  theme_prism() +
  scale_y_log10() +
  labs(y="Size of Regulon") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12, angle=60, hjust=1, 
                                   vjust=1, colour = "black"),
        axis.text.y = element_text(size=12, face="bold"))
saveFig(p4c_p5, 
        paste0("Plots/Fig4c_p5_Regulon_Size_barplot"),
        wd = 2000, ht = 1500)


## part 6: Top Regulon tilePlot: regulons ----
metaDf <- FetchData(CAR8, vars = c("stageCoarse", "annoCoarse")) %>%
  rownames_to_column(var = "cell")
annoDf <- annoDf %>%
  rownames_to_column(var = "cell")
annoDfMeta <- annoDf %>%
  left_join(metaDf, by = "cell")

reg8RFiltLvl <- reg8RFilt %>% 
  mutate(meanSHAP = ifelse(dir=="Exp", meanSHAP, -1 * meanSHAP)) %>%
  arrange(desc(meanSHAP))

metaDfx <- annoDfMeta %>%
  select(stageCoarse, annoCoarse, all_of(reg8RFilt$regulon)) %>%
  pivot_longer(cols = !c(annoCoarse, stageCoarse), names_to = "regulon", values_to = "SigScore") %>%
  group_by(regulon, annoCoarse, stageCoarse) %>%
  summarise(meanScore = mean(SigScore, na.rm=T)) %>%
  ungroup() %>%
  group_by(regulon) %>%
  mutate(meanScore = scale(meanScore)) %>%
  ungroup() %>%
  mutate(regulon = factor(regulon, levels = reg8RFiltLvl$regulon)) %>%
  mutate(annoCoarse = factor(annoCoarse, levels = c("EFF", "EFF_Exh-like", "EM_Exh-like", 
                                                    "EM_Eff-like", "CM", "Prolif")))

fancy_strips <- strip_nested(
  background_y = elem_list_rect(fill= c(colset8$stageCoarse[2:3]),
                                color = "white"),
  text_y = elem_list_text(face="bold",
                          size=10,
                          angle = 0,
                          color= "white"))


p4c_p6 <- ggplot(metaDfx, aes(y=annoCoarse, x=regulon, fill = meanScore)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  theme_minimal() +
  labs(fill = "Normalized\nScore")+
  facet_grid2(stageCoarse ~ ., scales = "free", space = "free",
              strip = fancy_strips, switch = "y") +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", mid=0) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle=60, hjust=1, 
                                   vjust=1, colour = "black"),
        axis.text.y = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.placement = "outside") +
  guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
saveFig(p4c_p6, 
        paste0("Plots/Fig4c_p6_Regulon_TileMap_annoCoarse"),
        wd = 2000, ht = 1500)

# ## part 7: Top Regulon tilePlot: genes ----
# dfx <- read.csv("Plots/Fig4_SupTable_Top10_Regulon_Genes.csv") %>%
#   filter(isDEG=="DEG") %>%
#   filter(isTF!="nonTF")
# 
# Exp_genes <- dfx %>%
#   filter(avg_log2FC>0) %>%
#   arrange(desc(avg_log2FC)) %>%
#   pull(gene)
# Per_genes <- dfx %>%
#   filter(avg_log2FC<0) %>%
#   arrange(desc(avg_log2FC)) %>%
#   pull(gene)
# gene_highlight <- c(Exp_genes, Per_genes)
# 
# 
# sigSet <- list(
#   Exp_Up = c("GNLY", "NKG7", "IFNG", "TNF", "GZMB", "GZMH", "GZMK", "KLRG1", "KLRD1", "HLA-DRA", "CCL5"),
#   Per_Up = c("TIGIT", "LAG3"),
#   AP1 = c("JUN", "FOS", "JUND"),
#   Stem = c("SELL", "IL7R", "CD7", "CD27", "S1PR1"),
#   Type2 = c("IL4R", "GATA3", "CXCR4"),
#   IFN1 = c("IFITM1", "IFITM2", "IFNAR1", "IFNAR2"),
#   others = c("NEAT1", "S100A4", "S100A6", "CD28", "ZAP70", "IL10RA", "CXCR3")
# )
# 
# 
# geneSelect <- unname(unlist(sigSet))
# 
# gsDf <- do.call("rbind", lapply(names(sigSet), function(sx){
#   ggx <- sigSet[[sx]]
#   data.frame(Gene=ggx, gLab = rep(sx, length(ggx)))
# }))
# 
# IP8 <- ScaleData(IP8,features = geneSelect, 
#                  vars.to.regress = c("percent.mt"))
# 
# sigDf <- FetchData(IP8, vars = c(geneSelect, "preGroup"), slot = "scale.data")
# 
# sigAvg <- sigDf %>%
#   filter(!is.na(preGroup)) %>%
#   mutate(preGroup = gsub("Others", "IP_Only", preGroup)) %>%
#   pivot_longer(cols = !c("preGroup"), names_to = "Gene", values_to = "Expr") %>%
#   group_by(Gene, preGroup) %>%
#   summarise(AvgExpr = mean(Expr, na.rm=T)) %>%
#   left_join(gsDf, by = "Gene") %>%
#   mutate(Gene = factor(Gene, levels = geneSelect)) %>%
#   mutate(gLab = factor(gLab, levels = names(sigSet))) %>%
#   mutate(preGroup = factor(preGroup, levels = rev(c("Exp_Pre","Per_Pre", "IP_Only")))) %>%
#   ungroup() 
# 
# p5e_p2 <- ggplot(sigAvg, aes(x=Gene, y=preGroup, fill = AvgExpr)) +
#   geom_tile(color = "black", lwd = 1, linetype = 1) +
#   theme_minimal() +
#   labs(fill = "Normalized\nExpression")+
#   facet_grid2(.~gLab, scales = "free_x", space = "free_x") +
#   theme(panel.grid = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(size=15, angle=60, hjust=1, 
#                                    vjust=1, colour = "black"),
#         axis.text.y = element_text(size=12, colour = "black", face="bold"),
#         legend.position = "top",
#         legend.text = element_text(size=12),
#         legend.title = element_text(size=12),
#         strip.text = element_text(size=15, face = "bold")) +
#   colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", mid=0) +
#   guides(fill = guide_colorbar(title.vjust = 0.9, barwidth = unit(8,"cm")))
# saveFig(p5e_p2, 
#         prefix = "Plots/Fig5e_p2_heatmap_DEG_Exp_vs_Per_Precursors_inIP", 
#         wd = 3000, ht = 1300)


regSelect <- c("CAR8R_S2_FOSL2", "CAR8R_S1_FLI1", "CAR8R_S2_JAZF1", "CAR8R_S2_ELF4", 
               "CAR4R_S2_STAT1", "CAR4R_S2_IRF7")
gDfList <- regDfList[regSelect]
fSelectList <- lapply(regDfList[regSelect], function(dfx){
  dfx %>%
    mutate(pctDiff = abs(pct.1-pct.2)) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    select(gene, sourceReg) %>%
    slice_head(n=5)
})

fDf <- do.call("rbind", fSelectList) %>%
  remove_rownames()

fDfToAdd <- data.frame(gene = c("IRF7", "OAS2", "STAT1", "ID2", "CD247", "NKG7", "CISH"),
                       sourceReg = rep(c("CAR4R_S2_IRF7", "CAR4R_S2_STAT1", 
                                         "CAR8R_S2_JAZF1", "CAR8R_S2_ELF4"), c(1,2,2,2)))
fDf_unique <- bind_rows(fDf_unique, fDfToAdd) %>%
  distinct()

fselect <- unique(fDf_unique$gene)

ALLR <- merge(CAR8R, CAR4R)
ALLR <- NormalizeData(ALLR)

require(colorspace)
subCAR <- subset(ALLR, stageCoarse != "Late")
subCAR$celltype <- ifelse(subCAR$annoCoarse %in% unique(CAR4R$annoCoarse), "CD4-CART", "CD8-CART")
subCAR$groupSuperLong <- paste0(subCAR$celltype, "_", subCAR$patientShort, "_", subCAR$stageCoarse)
exprMtx <- AverageExpression(subCAR, features = fselect, group.by = "groupSuperLong", assays = "RNA",
                             slot="data")$RNA

dfx <- as.data.frame(t(scale(t(exprMtx)))) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = !gene, names_to = "group", values_to = "Expr") %>%
  mutate(celltype = str_split(group, "_", simplify = T)[,1]) %>%
  mutate(patient = str_split(group, "_", simplify = T)[,2]) %>%
  mutate(stage = str_split(group, "_", simplify = T)[,3]) %>%
  mutate(gene = factor(gene, levels = rev(fselect))) %>%
  right_join(fDf_unique, by = "gene")

cellCol <- colSet$CD48Col

fancy_strips <- strip_nested(
  # horizontal strips
  background_x = elem_list_rect(fill= c(cellCol[1:2], stageCol[c(2:3)],stageCol[c(2:3)]),
                                #color=  c(stageCol[c(2:3)], cellCol[1:2], cellCol[1:2])),
                                color = rep("transparent", 6)),
  text_x = elem_list_text(face=rep("bold",6), 
                          size=rep(12,6), 
                          angle = rep(0, 6), 
                          color=c(rep("black", 2), rep("white",4))),
  text_y = elem_list_text(face = rep("bold", 5),
                          size = rep(12, 5),
                          angle = rep(0,5),
                          color = rep("black",5),
                          hjust = rep(0, 5)),
  size = "variable"
)


ps4b_v2 <- ggplot(dfx, aes(x=patient, y=gene, fill = Expr)) +
  geom_tile(color = "black", lwd = 1, linetype = 1) +
  facet_nested(sourceReg ~ celltype + stage, strip = fancy_strips, scales = "free_y", 
               space = "free_y") +
  theme_minimal() +
  labs(fill = "Normalized\nExpression")+
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
saveFig(pobj = ps4b_v2, 
        prefix = "Plots/Temporal/DraftOne/Fig5g_Heatmap_Regulon_Target_Gene_v2", 
        wd = 2400, ht = 2400)



