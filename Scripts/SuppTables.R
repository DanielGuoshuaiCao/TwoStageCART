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

##============================-- 2. Tbl1 Patient-Cell ----
### 2.1 post infusion
CAR8 <- readRDS("RDS/CAR8_for_manuscript_Apr262023.rds")
CAR4 <- readRDS("RDS/CAR4_for_manuscript_Apr262023.rds")

CAR4pc <- FetchData(CAR4, vars = c("patientShort", "time")) %>%
  group_by(patientShort, time) %>%
  count() %>%
  rename(PBMC_CD4_CART=n)

CAR8pc <- FetchData(CAR8, vars = c("patientShort", "time")) %>%
  group_by(patientShort, time) %>%
  count() %>%
  rename(PBMC_CD8_CART=n)

PItbl <- full_join(CAR4pc, CAR8pc, by = c("patientShort", "time"))
write.csv(PItbl, "Tables/SupTable1_Patient_Cell_Distribution_PostInfusion.csv", 
          row.names = F, quote = F)


### 2.2 IP 
IP_all <- readRDS("RDS/IP_b1b2_AxiCel_processed_merged_Jan02_2023.rds")
IP_all$patientShort <- gsub("CR", "P0", IP_all$patientShort)
IPpc <- FetchData(IP_all, vars = c("patientShort","Final_class")) %>%
  filter(Final_class %in% c("CD4", "CD8")) %>%
  group_by(patientShort, Final_class) %>%
  count() %>%
  pivot_wider(names_from = "Final_class", values_from = "n") %>%
  mutate(time="IP") %>%
  rename(CD4_CART = CD4, CD8_CART = CD8)

### 2.3 combine
PCtbl <- read.csv("Tables/SupTable1_Patient_Cell_Distribution_PostInfusion.csv") %>%
  rename(CD4_CART = PBMC_CD4_CART, CD8_CART = PBMC_CD8_CART) %>%
  bind_rows(IPpc) %>%
  arrange(patientShort, time)
write.csv(PCtbl, "Tables/SupTable1_Patient_Cell_Distribution_Overall_afterQC.csv", 
          row.names = F, quote = F)


##============================-- 3. Tbl2-3 Stage-specific DEG & GSEA ----
mStageSCRaw <- readRDS("RDS/SingleCell_DEG_S1_vs_S2_Raw_May112023.rds")
mStageSCfilt <- mStageSCRaw %>%
  filter(p_val_adj<0.05 & abs(avg_log2FC>0.25)) %>%
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC))
write.csv(mStageSCfilt, "Tables/SupTable2_Stage-Specific_DEG_SingleCellMethod.csv", row.names = F, quote = F)


msigCats <- c("C2", "H", "C5", "C7")
mdfList <- lapply(msigCats, function(catx){
  msigdbr(species = "Homo sapiens", category = catx) %>% 
    dplyr::select(gs_name, entrez_gene)
})
names(mdfList) <- msigCats

gvec <- GetGenelist(mStageSCRaw)

gseaListSC <- lapply(mdfList, function(m_dfx){
  GSEA(gvec, TERM2GENE = m_dfx, pvalueCutoff = 0.5, pAdjustMethod = "fdr")
})

lapply(names(gseaListSC), function(cx){
  sobj <- gseaListSC[[cx]]
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  mdfx <- sobj@result %>% filter(p.adjust<0.05)
  write.csv(mdfx, paste0("Tables/SupTable3_Stage-Sepcific_GSEA_Msigdb_Cat_", cx, ".csv"),
            row.names = F, quote = F)
})

### double check NES
stageGSEA <-  do.call("rbind", lapply(gseaListSC, function(sobj) {
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  sobj@result %>% filter(p.adjust<0.05)
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

selectDf %>% select(ID, NES)


##============================-- 3. Tbl4-5 Pre-specific DEG & GSEA ----
preDEG <- readRDS("RDS/DEG_inIP_Exp_vs_Per_Precursors_May02.rds")
preDEGfilt <- preDEG %>%
  filter(p_val_adj<0.05 & abs(avg_log2FC>0.25)) %>%
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC))
write.csv(preDEGfilt, "Tables/SupTable4_Precursor-Specific_DEG_SingleCellMethod.csv", row.names = F, quote = F)


preGSEA <- readRDS("RDS/GSEA_inIP_Exp_vs_Per_Precursors_May02.rds") %>%
  rownames_to_column(var = "fullID") %>%
  mutate(Cat = str_split(fullID, "\\.", simplify = T)[,1])


gseaListPre <- lapply(split(preGSEA, f=preGSEA$Cat), function(dfx){
  dfx %>% select(-fullID, -Cat)
})

lapply(c("H", "C2", "C5", "C7"), function(cx){
  mdfx <- gseaListPre[[cx]]
  write.csv(mdfx, paste0("Tables/SupTable5_Precursor-Sepcific_GSEA_Msigdb_Cat_", cx, ".csv"),
            row.names = F, quote = F)
})

### double check NES
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
  mutate(enrichment = ifelse(NES>0, "Exp", "Per")) %>%
  mutate(enrichment = factor(enrichment, levels=c("Exp", "Per")))

selectDf %>% select(ID, NES, p.adjust, qvalue, Cat)

##============================-- 4. Tbl6-7 Ex-specific DEG & GSEA ----
## EM-Ex vs EM
degRes <- readRDS("RDS/DEG_CAR8_EX_May06.rds")
EMDEG <- degRes$EMex
EMDEGfilt <- EMDEG %>%
  filter(p_val_adj<0.05 & abs(avg_log2FC>0.25))
write.csv(EMDEGfilt, "Tables/SupTable6_EM-Ex_vs_EM_DEG.csv")

##============================-- 5. Tbl7-8 Ex-specific DEG & GSEA ----
mList <- readRDS("RDS/PseudoBulk_DEG_S1_vs_S2_Libra_List_Apr262023.rds")
mStagePseudoRaw <- mList$All
mStagePseudofilt <- mStagePseudoRaw %>%
  filter(p_val_adj<0.05 & abs(avg_logFC > 1.5))
write.csv(mStagePseudofilt, "Tables/SupTable7_Stage-Specific_DEG_PseudoBulk.csv", row.names = F, quote = F)

pseudlBulkGSEA <- getGSEA(mStagePseudoRaw)

msigCats <- c("C2", "H", "C5", "C7")
mdfList <- lapply(msigCats, function(catx){
  msigdbr(species = "Homo sapiens", category = catx) %>% 
    dplyr::select(gs_name, entrez_gene)
})
names(mdfList) <- msigCats

gvec <- GetGenelist(mStagePseudoRaw)

gseaListPseudo <- lapply(mdfList, function(m_dfx){
  GSEA(gvec, TERM2GENE = m_dfx, pvalueCutoff = 0.5, pAdjustMethod = "fdr")
})

lapply(names(gseaListPseudo), function(cx){
  sobj <- gseaListPseudo[[cx]]
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  mdfx <- sobj@result %>% filter(p.adjust<0.05)
  write.csv(mdfx, paste0("Tables/SupTable8_Stage-Sepcific_GSEA_Msigdb_Cat_", cx, "_Pseudo.csv"),
            row.names = F, quote = F)
})

### double check NES
stageGSEA <-  do.call("rbind", lapply(gseaListPseudo, function(sobj) {
  sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
  sobj@result %>% filter(p.adjust<0.1)
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
  filter(ID %in% select_terms)
selectDf %>% select(ID, NES)


