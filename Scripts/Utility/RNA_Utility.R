## take both head and tail top n
headAndTail <- function(vecx,n1,n2=n1){
  if (n1 > floor(length(vecx)/2)){
    return(vecx)
  } else {
    return(c(head(vecx, n1), tail(vecx, n2))) 
  }
}

RelabelPatient <- function(sobj){
  sobj$patientShort <- plyr::mapvalues(sobj$patient, from = paste0("P", c("033", "028",
                                                                          "049", "052",
                                                                          "034", "032")),
                                       to = c("CR1", "NR1", 
                                              "CR2", "NR2",
                                              "CR3", "NR3")) 
  return(sobj)
}


## convert human symbol to mouse symbol
HM_Convert <- function(genes, HtoM=T, simple=T){
  # require pre-compiled conversion table
  HM_Table <- readRDS("RDS/HM_Conversion_Table.rds")
  
  # strip any "." in genes
  genes.clean <- unlist(lapply(genes, function(str){
    strsplit(str, "\\.")[[1]][1]
  }))
  
  if (HtoM){
    refvec <- HM_Table$Mouse
    names(refvec) <- HM_Table$Human
  } else {
    refvec <- HM_Table$Human
    names(refvec) <- HM_Table$Mouse
  }
  
  # convert
  gout <- refvec[genes.clean]
  names(gout) <- genes
  
  # simplify
  if(simple){
    gout <- unname(gout)
    gout <- gout[!is.na(gout)] 
  }
  
  return(gout)
}


## Generic function to split tradeseq result df into lineage list of dfs
breakTradeDf <- function(resDf, Nlin){
  LinList <- lapply(c(1:Nlin), function(idx){
    dfx <- resDf[ , grepl(as.character(idx), colnames(resDf))]
    xx <- colnames(dfx)
    colnames(dfx) <- unlist(lapply(xx, function(ss) {
      s1 <- strsplit(ss, "_")[[1]][1]
      s2 <- gsub(paste0("lineage", idx), "", s1)
      s2
      }))
    dfx$padj <- p.adjust(dfx$pvalue, method = "BH")
    dfx <- dfx[dfx$padj<0.05,]
    if (any(grepl("log", colnames(dfx)))){
      dfx <- dfx[order(dfx$logFC, decreasing = T),] 
    } else {
      dfx <- dfx[order(dfx$waldStat, decreasing = T),] 
    }
    dfx <- na.omit(dfx)
    dfx$gene <- row.names(dfx)
    dfx
  }) 
  names(LinList) <- paste0("Lineage", 1:Nlin)
  return(LinList)
}

## clip values of a matrix
clipMtx <- function(mtx, cap){
  mtx[mtx > cap] <- cap
  mtx[mtx < -cap] <- -cap
  return(mtx)
}

## clip values of a matrix
clipVec <- function(vec, cap, val=cap ,lower=T){
  if (lower) {
    vec[vec<=cap] <- val
  } else {
    vec[vec>=cap] <- val
  }
  return(vec)
}


## search a gene
getGene <- function(str, sobj=TIiso){
  ss <- row.names(sobj)
  return(sort(ss[grepl(str, ss)]))
}


## regList to Network
regToNetwork <- function(regL){
  TFs <- names(regL)
  netDf <- do.call("rbind",lapply(TFs, function(tfx){
    targets <- regL[[tfx]]
    data.frame(from=rep(tfx, length(targets)), to=targets)
  }))
  row.names(netDf) <- NULL
  return(netDf)
}

## transfer Labels
TransferLabels <- function(fromObj, toObj, cols){
  from_meta <- FetchData(fromObj, cols)
  to_meta <- FetchData(toObj, "orig.ident")
  combine_meta <- base::merge(to_meta, from_meta, by=0, all.x=T)
  row.names(combine_meta) <- combine_meta$Row.names
  combine_meta$Row.names <- NULL
  toObj <- AddMetaData(toObj, combine_meta)
  return(toObj)
}

## trim totalseq names
TrimTotalSeqNames <- function(sobj){
  row.names(sobj@assays$ADT@counts) <- gsub("-TotalSeqC", "", row.names(sobj@assays$ADT@counts))
  row.names(sobj@assays$ADT@data) <- gsub("-TotalSeqC", "", row.names(sobj@assays$ADT@data))
  row.names(sobj@assays$ADT@scale.data) <- gsub("-TotalSeqC", "", row.names(sobj@assays$ADT@scale.data))
  return(sobj)
}

## conver large sparse matrix to dense matrix
asDenseMatrix <- function(norm_counts){
  
  if (ncol(norm_counts)>10000){
    xx <- as.matrix(norm_counts[, c(1:10000)])
    NN <- ceiling(ncol(norm_counts)/10000)
    for (ii in 2:NN){
      nstart <- (ii-1)*10000+1
      nstop <- min(ncol(norm_counts), ii*10000)
      xx <- cbind(xx,as.matrix(norm_counts[, nstart:nstop]))
    }
  } else {
    xx <- as.matrix(norm_counts)
  }
  return(xx)
}



##============================-- Quick Plot Functions --=======================================##

### Dimplot with default ident
DP <- function(sobj, ident=NULL, red="umap", ...) {
  if (is.null(ident)){
    DimPlot(sobj, label = T, repel = T, reduction = red, ...) + NoLegend() 
  } else {
    DimPlot(sobj, label = T, repel = T, group.by = ident, reduction = red, ...) + NoLegend() 
  }
}

## binary plot
BinaryDP <- function(sobj, identX, metacol){
  sobj$binary <- ifelse(sobj[[metacol]][,1]==identX, "yes", "no")
  DimPlot(sobj, group.by = "binary")
}

## quick violin plot
VP <- function(sobj, gx, ident=NULL, xtSize=25, gSize=30){
  if (is.null(ident)){
    VlnPlot(sobj, features = gx, pt.size = 0.01) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = xtSize),
            plot.title = element_text(size=gSize)) +
      NoLegend()
  } else {
    VlnPlot(sobj, features = gx, pt.size = 0.01, group.by = ident) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = xtSize),
            plot.title = element_text(size=gSize)) +
      NoLegend()
  }
}


## SigScore Plot
# SigVlnPlot <- function(sobj,sigx, toPrint=T){
#   require(rstatix)
#   
#   ftable <- FetchData(sobj, vars = c("annoExtraFine2", "disease", sigx))
#   colnames(ftable)[c(1,3)] <- c("Anno_Fine","SigScore")
#   stat.test <- ftable %>% 
#     group_by(Anno_Fine) %>%
#     wilcox_test(SigScore ~ disease, p.adjust.method="fdr") %>%
#     add_significance() %>%
#     add_xy_position(x="Anno_Fine", dodge=0.8)
#   
#   ymargin = diff(range(ftable$SigScore)) * 0.3
#   ymax = max(ftable$SigScore) + ymargin
#   
#   
#   psigx <- ggviolin(ftable, x="Anno_Fine", y="SigScore", fill="disease") +
#     geom_boxplot(width=0.15, position = position_dodge(width = 0.8), outlier.colour = "transparent", fill="white",
#                  aes(color=disease), show.legend = F) +
#     stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0,
#                        bracket.nudge.y = ymargin*1/3, label.size = 12) +
#     labs(y=sigx, x="", title = sigx) +
#     theme_pubr() +
#     ylim(0,ymax) +
#     theme(axis.text.y = element_text(size=12,color = 'black'),
#           axis.text.x = element_text(size=15,color = 'black', 
#                                      angle = 45, hjust=1),
#           axis.title.y = element_text(size=18, margin = margin(r=10)),
#           plot.title = element_text(size=20, hjust = 0.5, margin = margin(b=15), face = "bold")) +
#     scale_fill_manual(values = disease_cols) +
#     scale_color_manual(values = c("black", "black")) +
#     NoLegend()
#   
#   if (toPrint){
#     png(paste0("Plots/Fig4/Signature_Plot_", sigx, ".png"), width = 3600, height = 1500, res = 300)
#     print(psigx)
#     dev.off()  
#   } else {
#     print(psigx)
#   }
# }

## quickly visualize dynamic change of a particular gene
Bulk_Dynamics <- function(gx, rna_expr, saveToDrive=F){
  rx <- subset(rna_expr, gene==gx)
  lbls <- colnames(rna_expr)[-43]
  expr <- as.numeric(rx[1,1:42])
  gps <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][1]))
  days <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][2]))
  reps <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][3]))
  dat <- data.frame(Expr = expr, Time=days, Group=gps, Rep=reps, Cat=paste0(gps, "_", days))
  dat <- do.call("rbind", lapply(split(dat, f=dat$Cat), function(ddf){
    expr <- mean(ddf$Expr)
    sdex <- sd(ddf$Expr)
    dfx <- data.frame(Expr=expr, sd=sdex, Group=unique(ddf$Group), Time=unique(ddf$Time))
    dfx
  }))
  dat <- subset(dat, Group!="DLN")
  diseases <- c("Cl13", "Arm","Tumor")
  dat_sub <- dat[rep("Naive_D0",3),]
  dat <- dat[dat$Group!="Naive",]
  dat_sub$Group <- diseases
  row.names(dat_sub) <- paste0(dat_sub$Group, "_", dat_sub$Time)
  dat <- rbind(dat, dat_sub)
  dat$Time <- factor(dat$Time, levels = c("D0", "D4", "D8", "D21", "D90", "D300"))
  dat$Group <- factor(dat$Group, levels = diseases)
  delta <- max((max(dat$Expr) - min(dat$Expr))*0.3, 1)

  p <- ggplot(data=dat, aes(x=Time, y=Expr, color=Group, group=Group))+
    geom_errorbar(aes(ymin=Expr-sd, ymax=Expr+sd), width=0.1, size=1.5, show.legend = F)+
    geom_line(size=1.5)+ geom_point()+theme_pubr()+
    xlab("")+ylab("Normalized Expression")+
    labs(color=gx)+
    theme(axis.text.y = element_text(size=20,color = 'black'),
          axis.text.x = element_text(size=20,color = 'black'),
          axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          legend.key.width = unit(1.5, "cm"))+
    guides(color = guide_legend(override.aes = list(size=2)))+
    scale_color_manual(values = pal_npg()(3)[c(1,3,2)]) +
    scale_y_continuous(expand = expansion(add = c(delta,delta)))
  
  if(saveToDrive){
    png(paste0("Plots/Fig1/Gene_expr_",gx,".png"), width = 2400, height = 1800, res = 300)
    print(p)
    dev.off()  
  } else {
    print(p)
  }
}


## plot gene signature or genes
SignatureDynamics <- function(sDf, saveToDrive=T){
  smtx <- sDf
  sname <- colnames(sDf)[1]
  lbls <- row.names(smtx)
  expr <- as.numeric(smtx[,1])
  gps <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][1]))
  days <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][2]))
  reps <- unlist(lapply(lbls, function(str) strsplit(str, "_")[[1]][3]))
  dat <- data.frame(Expr = expr, Time=days, Group=gps, Rep=reps, Cat=paste0(gps, "_", days))
  dat <- do.call("rbind", lapply(split(dat, f=dat$Cat), function(ddf){
    expr <- mean(ddf$Expr)
    sdex <- sd(ddf$Expr)
    dfx <- data.frame(Expr=expr, sd=sdex, Group=unique(ddf$Group), Time=unique(ddf$Time))
    dfx
  }))
  
  diseases <- c("Cl13", "Arm","Tumor")
  dat_sub <- dat[rep("Naive_D0",3),]
  dat <- dat[dat$Group!="Naive",]
  dat_sub$Group <- diseases
  row.names(dat_sub) <- paste0(dat_sub$Group, "_", dat_sub$Time)
  dat <- rbind(dat, dat_sub)
  dat$Time <- factor(dat$Time, levels = c("D0", "D4", "D8", "D21", "D90", "D300"))
  dat$Group <- factor(dat$Group, levels = diseases)
  
  p <- ggplot(data=dat, aes(x=Time, y=Expr, color=Group, group=Group))+
    geom_errorbar(aes(ymin=Expr-sd, ymax=Expr+sd), width=0.1, size=1.5, show.legend = F)+
    geom_line(size=1.5)+ geom_point()+theme_pubr()+
    xlab("")+ylab("AUC Score")+
    labs(color="")+
    theme(axis.text.y = element_text(size=20,color = 'black'),
          axis.text.x = element_text(size=20,color = 'black'),
          axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          legend.text = element_text(size=20),
          legend.title = element_blank(),
          legend.key.width = unit(1.5, "cm"))+
    guides(color = guide_legend(override.aes = list(size=2)))+
    scale_color_manual(values = pal_npg()(3)[c(1,3,2)]) +
    scale_y_continuous(expand = expansion(add = c(0.1, 0.1)))
  
  if(saveToDrive){
    png(paste0("Plots/Fig1/",sname,".png"), width = 2400, height = 1800, res = 300)
    print(p)
    dev.off()  
  } else {
    print(p)
  }
}


TrajUmap <- function(curveList, sobj, segLenVec, colors=brewer.pal(10, "Set3"), widthCurve = 1){
  curveList <- lapply(curveList, function(cx){
    curve1 <- as.data.frame(cx$s) %>% arrange(desc(UMAP_1))
    curve1
  })
  
  if (length(segLenVec)>1){
    if (length(segLenVec)!=length(curveList)){
      stop("wrong length of segLenVec")
    } else{
      curveList <- Map(function(cx, sx){
        cx <- cx[1:round(nrow(cx)*sx),]
      }, curveList, segLenVec)
    }
  } else {
    curveList <- lapply(curveList, function(cx) {
      cx <- cx[1:round(nrow(cx)*segLenVec),]
    })
  }
  
  pp <- DimPlot(sobj, reduction = "umap", group.by = "annoFine",
                pt.size = 1, label=F, raster = F, shuffle = T)+
    xlab("UMAP1")+ylab("UMAP2")+
    theme_pubr() +
    theme(axis.text = element_text(size=20,color = 'black'),
          axis.title = element_text(size=25,face = 'bold'),
          legend.text = element_text(size=24),
          plot.title = element_blank(),
          legend.key.width = unit(1.5,"inches")) +
    scale_color_manual(values = anno_cols)  +
    guides(color=guide_legend(override.aes = list(size=8)))
  
  for (ii in 1:length(curveList)){
    curvex <- curveList[[ii]]
    colorx <- colors[ii]
    pp <- pp + geom_line(data = curvex, aes(x=UMAP_1, y=UMAP_2), 
                         color=colorx, size = widthCurve,
                         arrow=arrow(ends = "first", length = unit(0.2, "inches"), 
                                     type="closed",angle = 15))
  }
  return(pp)
}


TrajUmapExtra <- function(curveList, sobj, segLenVec, colors=brewer.pal(10, "Set3"), widthCurve = 1){
  curveList <- lapply(curveList, function(cx){
    curve1 <- as.data.frame(cx$s) %>% arrange(desc(UMAP_1))
    curve1
  })
  
  if (length(segLenVec)>1){
    if (length(segLenVec)!=length(curveList)){
      stop("wrong length of segLenVec")
    } else{
      curveList <- Map(function(cx, sx){
        cx <- cx[1:round(nrow(cx)*sx),]
      }, curveList, segLenVec)
    }
  } else {
    curveList <- lapply(curveList, function(cx) {
      cx <- cx[1:round(nrow(cx)*segLenVec),]
    })
  }
  
  pp <- DimPlot(sobj, reduction = "umap", group.by = "annoExtraFine2",
                pt.size = 1.5, label=F, raster = F, shuffle = T)+
    xlab("UMAP1")+ylab("UMAP2")+
    theme(axis.text = element_text(size=15,color = 'black'),
          axis.title = element_text(size=18,face = 'bold'),
          legend.text = element_text(size=15),
          plot.title = element_blank(),
          legend.position = "top",
          legend.key.width = unit(1.3,"cm"))+
    scale_color_manual(values = anno_cols)  +
    guides(color=guide_legend(override.aes = list(size=6), 
                              nrow = 4))
  
  for (ii in 1:length(curveList)){
    curvex <- curveList[[ii]]
    colorx <- colors[ii]
    pp <- pp + geom_line(data = curvex, aes(x=UMAP_1, y=UMAP_2), 
                         color=colorx, size = widthCurve,
                         arrow=arrow(ends = "first", length = unit(0.2, "inches"), 
                                     type="closed",angle = 15))
  }
  return(pp)
}

## split violin plot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                 newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



## Trajectory Gene Dynamics Plot
TrajDynamics <- function(gx, curveIdx, sobj, curveObj, norm=F, branchList=NULL, colpal=TrajCols ){
  curveIdx <- sort(curveIdx)
  branchTime <- list("1_2" = 1.218, "3_4"=2.885, "2_4"=2.885, "2_3" = 6.245,
                     "1_3" = 1.218, "1_4"=1.218)
  TrajColTemp <- colpal[curveIdx]
  cellTime <- as.data.frame(slingPseudotime(curveObj))
  gtDf <- do.call("rbind", lapply(curveIdx, function(ii){
    dfxTemp <- cellTime[, ii,drop=F]
    colnames(dfxTemp) <- "Pseudotime"
    cellSelect <- row.names(dfxTemp)
    dfy <- FetchData(sobj, vars = c(gx))
    colnames(dfy) <- "NormExp"
    dfy[] <- dfy[cellSelect, ]
    dfout <- cbind(dfy, dfxTemp)
    dfout <- na.omit(dfout)
    xx <- dfout$Pseudotime
    
    if (norm){
      dfout$Pseudotime <- (xx - min(xx))/(max(xx) - min(xx)) 
    }
    
    dfout$curve <- paste0("Trajectory_",ii)
    dfout
  }))
  
  pp <- ggplot(gtDf, aes(x=Pseudotime, y=NormExp, color=curve)) +
    geom_smooth(method = "loess", size=1) +
    theme_pubr() +
    scale_color_manual(values = TrajColTemp)
  
  
  if (!is.null(branchList)){
    for (ii in 1:length(branchList)){
      branchT <- branchTime[[branchList[ii]]]
      pp <- pp + geom_vline(xintercept = branchT, linetype="dashed", size=0.6) +
        geom_vline(xintercept = branchT-1, linetype="dashed", size=0.6)
    }
  }
  
  
  if (norm){
    pp <- pp + labs(x="Normalized Pseudotime",
                    y="Normalized Expression", 
                    title = gx)
  } else {
    pp <- pp + labs(x="Pseudotime", 
                    y="Normalized Expression",
                    title = gx)
  }
  
  pp <- pp + 
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          plot.title = element_text(size=25, face="bold", hjust = 0.5),
          legend.position = "none") 
  
  return(pp)
}


## Trajectory Gene Dynamics Plot
TrajDynamics_byDisease <- function(gx, curveIdx, sobj, curveObj, norm=F, branchList=NULL, colpal=TrajCols){
  curveIdx <- sort(curveIdx)
  branchTime <- list("1_2" = 1.218, "3_4"=2.885, "2_4"=2.885, "2_3" = 6.245,
                     "1_3" = 1.218, "1_4"=1.218)
  TrajColTemp <- colpal[curveIdx]
  cellTime <- as.data.frame(slingPseudotime(curveObj))
  gtDf <- do.call("rbind", lapply(curveIdx, function(ii){
    dfxTemp <- cellTime[, ii,drop=F]
    colnames(dfxTemp) <- "Pseudotime"
    cellSelect <- row.names(dfxTemp)
    dfy <- FetchData(sobj, vars = c(gx, "disease"))
    colnames(dfy) <- c("NormExp", "disease")
    dfy[] <- dfy[cellSelect, ]
    dfout <- cbind(dfy, dfxTemp)
    dfout <- na.omit(dfout)
    xx <- dfout$Pseudotime
    
    if (norm){
      dfout$Pseudotime <- (xx - min(xx))/(max(xx) - min(xx)) 
    }
    
    dfout$curve <- paste0("Trajectory_",ii)
    dfout
  }))
  
  pp <- ggplot(gtDf, aes(x=Pseudotime, y=NormExp, color=disease)) +
    geom_smooth(method = "loess", size=1) +
    theme_pubr() +
    scale_color_manual(values = disease_cols) +
    facet_wrap(~curve, nrow = 2)
  
  
  if (!is.null(branchList)){
    for (ii in 1:length(branchList)){
      branchT <- branchTime[[branchList[ii]]]
      pp <- pp + geom_vline(xintercept = branchT, linetype="dashed", size=0.6) +
        geom_vline(xintercept = branchT-1, linetype="dashed", size=0.6)
    }
  }
  
  
  if (norm){
    pp <- pp + labs(x="Normalized Pseudotime",
                    y="Normalized Expression", 
                    title = gx)
  } else {
    pp <- pp + labs(x="Pseudotime", 
                    y="Normalized Expression",
                    title = gx)
  }
  
  pp <- pp + 
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          plot.title = element_text(size=25, face="bold", hjust = 0.5),
          legend.position = "none",
          strip.text.x = element_text(face = "bold", size=12),
          strip.background = element_rect(color="transparent", fill = "transparent")) 
  
  return(pp)
}

## Trajectory Gene Dynamics Plot with multiple genes
TrajDynamicsGrid <- function(gvec, curveIdx, sobj, curveObj, norm=F, Ncols=2, colpal=TrajCols){
  curveIdx <- sort(curveIdx)
  TrajColTemp <- colpal[curveIdx]
  cellTime <- as.data.frame(slingPseudotime(curveObj))
  norm <- T
  
  gtDf <- do.call("rbind", lapply(gvec, function(gx){
    dfx <-  do.call("rbind", lapply(curveIdx, function(ii){
      dfxTemp <- cellTime[, ii,drop=F]
      colnames(dfxTemp) <- "Pseudotime"
      cellSelect <- row.names(dfxTemp)
      dfy <- FetchData(sobj, vars = c(gx))
      colnames(dfy) <- "NormExp"
      dfy[] <- dfy[cellSelect, ]
      dfout <- cbind(dfy, dfxTemp)
      dfout <- na.omit(dfout)
      xx <- dfout$Pseudotime
      
      if (norm){
        dfout$Pseudotime <- (xx - min(xx))/(max(xx) - min(xx)) 
      }
      
      dfout$curve <- paste0("Trajectory_",ii)
      dfout$gene <- gx
      dfout
    }))
    dfx
  }))
  
  gtDf$gene <- factor(gtDf$gene, levels = gvec)
  
  
  pp <- ggplot(gtDf, aes(x=Pseudotime, y=NormExp, color=curve)) +
    geom_smooth(method = "loess", size=1) +
    theme_pubr() +
    scale_color_manual(values = TrajColTemp) +
    facet_wrap(~gene, ncol = Ncols, scales = 'free_y')
  
  if (norm){
    pp <- pp + labs(x="Normalized Pseudotime",
                    y="Normalized Expression")
  } else {
    pp <- pp + labs(x="Pseudotime", 
                    y="Normalized Expression")
  }
  
  pp <- pp + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.background = element_rect(fill="transparent", color="transparent"),
          strip.text = element_text(size=12),
          plot.title = element_blank(),
          panel.grid.major.x = element_line(color="grey90"),
          legend.position = "none") 
  return(pp)
}

## Bulk heatmap



##============================-- Preprocess --=======================================##
PreprocessDefault <- function(sobj, nVgenes=5000, minDistUmap=0.3, nNeighborUmap=30, 
                              cellCycle=F, harmVar=c("patient"), Human=T, normADT=F, nDims=30, nRes=0.3){
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(sobj, nfeatures = nVgenes)
  if (cellCycle){
    if (Human){
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    } else {
      s.genes <- HM_Convert(cc.genes$s.genes, simple = T)
      g2m.genes <- HM_Convert(cc.genes$g2m.genes, simple = T) 
    }
    sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
    sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
    sobj <- ScaleData(sobj, vars.to.regress = c("percent.mt", "CC.Difference"))
  } else {
    sobj <- ScaleData(sobj, vars.to.regress = c("percent.mt"))
  }
  
  
  
  sobj <- RunPCA(sobj, npcs = 100, verbose = FALSE)
  ElbowPlot(sobj, ndims=100)
  require(harmony)
  
  if (!is.na(harmVar)){
    sobj <- RunHarmony(sobj, group.by.vars = harmVar)
    sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:nDims,
                    min.dist = minDistUmap, n.neighbors = nNeighborUmap)
    sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:nDims)
    sobj <- FindClusters(sobj, resolution = nRes)
  } else {
    sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:nDims,
                    min.dist = minDistUmap, n.neighbors = nNeighborUmap) 
    sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:nDims)
    sobj <- FindClusters(sobj, resolution = nRes)
  }
  
  ## re normalize ADT
  if (normADT){
    sobj <- NormalizeData(sobj, assay = "ADT", normalization.method = "CLR")
    sobj <- ScaleData(sobj, assay = "ADT") 
  }
  
  return(sobj)
}



##============================-- Geneset Scores --=======================================##
# submsig <- msigdbr() %>%
#   filter(gs_cat %in% c("H", "C2", "C3", "C5", "C7", "C8"))
# saveRDS(submsig, "RDS/msigdbr_subset_HomoSapiens.rds")

## extract gene list of a particular term
GetMsigdbSet <- function(term, vague=F){
  submsig <- readRDS("RDS/msigdbr_subset_HomoSapiens.rds")
  require(msigdbr)
  if (vague){
    gsetDf <- submsig[grepl(term, submsig$gs_name), c("gs_name", "gene_symbol")]
  } else {
    gsetDf <- submsig[submsig$gs_name %in% term, c("gs_name", "gene_symbol")]
  }
  glist <- list()
  for (nx in unique(gsetDf$gs_name)) {
    gsetx <- gsetDf %>%
      filter(gs_name==nx) %>%
      pull(gene_symbol)
    gsetx <- unique(gsetx)
    glist[[nx]] <- gsetx
  }
  return(glist)
}

## extract gene list of a particular gene
GetGeneMsigdbSet <- function(gene){
  submsig <- readRDS("RDS/msigdbr_subset_HomoSapiens.rds")
  require(msigdbr)
  gsetDf <- submsig[submsig$gene_symbol %in% gene, c("gs_name", "gene_symbol")]
  glist <- list()
  for (nx in unique(gsetDf$gs_name)) {
    gsetx <- gsetDf %>%
      filter(gs_name==nx) %>%
      pull(gene_symbol)
    gsetx <- unique(gsetx)
    glist[[nx]] <- gsetx
  }
  return(glist)
}





ScoreMsigdbSet <- function(sobj, term, vague=T){
  require(msigdbr)
  if (!exists("submsig")){
    submsig <- readRDS("RDS/msigdbr_subset_HomoSapiens.rds") 
  }
  if (vague){
    gsetDf <- submsig %>%
      filter(grepl(term, gs_name)) %>%
      dplyr::select(gs_name, gene_symbol) 
  } else {
    gsetDf <- submsig %>%
      filter(gs_name %in% term) %>%
      dplyr::select(gs_name, gene_symbol)
  }
  
  
  for (nx in unique(gsetDf$gs_name)) {
    gsetx <- gsetDf %>%
      filter(gs_name==nx) %>%
      pull(gene_symbol)
    sobj <- AddModuleScore(sobj, features = list(gsetx), name = nx)
  }
  return(sobj)
}


CalGeneSetSignature <- function(geneset, CR, Sname, nc=7, maxrankcell=0.05, raw=F){
  require("GSEABase")
  require("AUCell")
  xx <- geneset
  memset <- GeneSet(xx, setName=Sname)
  Cells_AUC <- AUCell_calcAUC(memset, CR, nCores = nc, 
                              aucMaxRank = ceiling(maxrankcell*nrow(CR)))
  if (raw) {
    return(Cells_AUC)
  }
  mtx_score <- as.data.frame(t(getAUC(Cells_AUC)))
  return(mtx_score)
}

ScoreSeurat <- function(sobj, regulon_list, ncore=7, mthresh=0.05, returnRaw=F, prefix=""){
  require("AUCell")
  require("GSEABase")
  counts <- sobj@assays$RNA@data
  ranking <- AUCell_buildRankings(asDenseMatrix(counts), plotStats = F)
  TFs <- names(regulon_list)
  Score_list <- lapply(TFs, function(tfx){
    gsx <- regulon_list[[tfx]]
    CalGeneSetSignature(geneset = gsx, CR=ranking, Sname = tfx, nc=ncore, maxrankcell = mthresh)
  })
  score_df <- as.data.frame(do.call("cbind", Score_list))
  colnames(score_df) <- paste0(prefix, TFs)

  if (returnRaw){
    return(score_df)
  }
  
  sobj <- AddMetaData(sobj, score_df)
  return(sobj)
}

getRawScoreSeurat <- function(sobj, regulon_list, tfx, ncore=7, mthresh=0.2){
  require("AUCell")
  require("GSEABase")
  counts <- sobj@assays$RNA@counts
  ranking <- AUCell_buildRankings(asDenseMatrix(counts), plotStats = F)
  gsx <- regulon_list[[tfx]]
  scoreobj <- CalGeneSetSignature(geneset = gsx, CR=ranking, Sname = tfx, nc=ncore, 
                                  maxrankcell = mthresh, raw=T)
  return(scoreobj)
}

getBulkScoreSeurat <- function(mtx, signatureVec, Sname, ncore=3, mthresh = 0.7, raw=F) {
  require("AUCell")
  require("GSEABase")
  ranking <- AUCell_buildRankings(mtx, plotStats = F)
  memset <- GeneSet(signatureVec, setName=Sname)
  Cells_AUC <- AUCell_calcAUC(memset, ranking, nCores = ncore, 
                              aucMaxRank = ceiling(mthresh*nrow(ranking)))
  if (raw) {
    return(Cells_AUC)
  }
  mtx_score <- as.data.frame(t(getAUC(Cells_AUC)))
  return(mtx_score)
  
}


##============================-- Pathway Enrichment --=======================================##
## convert differential 
GetGenelist <- function(vecOrDf){
  if (is.data.frame(vecOrDf)){
    ddf <- vecOrDf
    xx <- colnames(ddf)[grep("log",colnames(ddf))]
    gvec <- ddf[[xx]]
    if (!("gene" %in% colnames(ddf))) ddf$gene <- row.names(ddf)
  } else {
    gvec <- vecOrDf
    ddf <- data.frame(gene=gvec, stringsAsFactors = F)
  }
  gdf <- biomaRt::select(org.Hs.eg.db, ddf$gene, "ENTREZID", "SYMBOL")
  gdf <- gdf[!base::duplicated(gdf$SYMBOL),]
  names(gvec) <- as.character(gdf$ENTREZID)
  gvec <- sort(gvec, decreasing = T)
  gvec <- gvec[!is.na(names(gvec))]
  if (is.data.frame(vecOrDf)){
    return(gvec)
  } else {
    return(names(gvec) )
  }
}


## 1.5 Wrapper function to do enrichGO from clusterprofiler with preset pval, fdr, ont
myEnrichGO <- function(genes, pthresh=0.05){
  BP <- enrichGO(genes, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "BP", readable = T)
  CC <- enrichGO(genes, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "CC", readable = T)
  MF <- enrichGO(genes, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "MF", readable = T)
  pathway <- list(BP, CC, MF)
  names(pathway) <- c("BP", "CC", "MF")
  return(pathway)
}

## 1.5.2 Wrapper function to do enrichGO from clusterprofiler with preset pval, fdr, ont
myEnrichGO_msigdb <- function(genes, pthresh=0.05){
  require(msigdbr)
  reslist <- list()
  mobj <- readRDS("RDS/msigdbr_subset_HomoSapiens.rds")
  catlbl <- unique(mobj$gs_cat)
  for (cat in catlbl){
    m_df <- mobj %>% 
      filter(gs_cat==cat) %>%
      dplyr::select(gs_name, entrez_gene)
    pathway <- enricher(genes, TERM2GENE = m_df, pvalueCutoff = pthresh, pAdjustMethod = "fdr")
    reslist[[cat]] <- pathway
  }
  return(reslist)
}


## 1.6.1 Wrapper function to do GO GSEA from clusterprofiler
myEnrichGSEA <- function(gvec, pthresh=0.05){
  BP <- gseGO(gvec, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "BP", minGSSize = 1, maxGSSize = 20000)
  CC <- gseGO(gvec, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "CC")
  MF <- gseGO(gvec, OrgDb = org.Hs.eg.db, pvalueCutoff = pthresh, pAdjustMethod = "fdr", ont = "MF")
  pathway <- list(BP, CC, MF)
  names(pathway) <- c("BP", "CC", "MF")
  return(pathway)
}


## 1.6.2 Wrapper function to do MESH GSEA from clusterprofiler
myEnrichMesh <- function(gvec, pthresh=0.05){
  require(meshes)
  reslist <- list()
  catlbl <- c("A", "C", "D", "E", "G", "L", "N")
  for (cat in catlbl){
    pathway <- gseMeSH(gvec, MeSHDb = "MeSH.Hsa.eg.db", database = "gene2pubmed", category = cat,
                       pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 5000)
    reslist[[cat]] <- pathway
  }
  return(reslist)
}

## 1.5.2 Wrapper function to do enrichGO from clusterprofiler with preset pval, fdr, ont
myEnrichGSEA_msigdb <- function(gvec, pthresh=0.05){
  require(msigdbr)
  reslist <- list()
  mobj <- readRDS("RDS/msigdbr_subset_HomoSapiens.rds")
  catlbl <- unique(mobj$gs_cat)
  for (cat in catlbl){
    m_df <- mobj %>% 
      filter(gs_cat==cat) %>%
      dplyr::select(gs_name, entrez_gene)
    pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = pthresh, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 5000)
    reslist[[cat]] <- pathway
  }
  return(reslist)
}


## 1.6 Extract genes from the GO result
ExtractGenes <- function(dfx, idx){
  xx <- dfx[idx, "geneID"]
  xx <- strsplit(xx, "\\/")[[1]]
  dfy <- select(org.Hs.eg.db, xx, "SYMBOL", "ENTREZID")
  gx <- dfy$SYMBOL
  return(gx)
}


##============================-- Differential Gene Expression (Bulk) --=======================================##
## function to quickly get the result of a particular contrast
getDEG <- function(ddsobj,c1,c2, pmethod="fdr") {
  res <- results(ddsobj, contrast = c("Condition", c1, c2), pAdjustMethod = pmethod)
  resDf <- res %>%
    as.data.frame() %>%
    drop_na() %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::rename(log2FC=log2FoldChange) %>%
    filter(padj < 0.05)
  return(resDf)
}

##  Double Differential Gene Expression
getDoubleDGE <- function(CvA, TvA, maxlfc=5, thresh_lfc=1.5, thresh_padj=0.05, geneHL=NULL){
  colnames(CvA) <- paste0(colnames(CvA), "_CvA")
  colnames(TvA) <- paste0(colnames(TvA), "_TvA")
  CvA$gene <- row.names(CvA)
  TvA$gene <- row.names(TvA)
  dfx <- merge(CvA, TvA, by="gene", all=F)
  dfx$log2FC_CvA[dfx$log2FC_CvA > maxlfc] <- maxlfc
  dfx$log2FC_TvA[dfx$log2FC_TvA > maxlfc] <- maxlfc
  dfx$log2FC_CvA[dfx$log2FC_CvA < -maxlfc] <- -maxlfc
  dfx$log2FC_TvA[dfx$log2FC_TvA < -maxlfc] <- -maxlfc
  dfx$group <- "Non-DEGs"
  dfx$maxPadj <- pmax(dfx$padj_CvA, dfx$padj_TvA)
  dfx$maxLogMean <- log(pmax(dfx$baseMean_CvA, dfx$baseMean_TvA))
  dfx$ExprLvl <- ifelse(dfx$maxLogMean>5.2, "high", "low")
  dfx$group[abs(dfx$log2FC_CvA)> thresh_lfc  &  abs(dfx$log2FC_TvA)< thresh_lfc & dfx$maxPadj < thresh_padj] <- "Cl13 Specific"
  dfx$group[abs(dfx$log2FC_CvA)< thresh_lfc  &  abs(dfx$log2FC_TvA)> thresh_lfc & dfx$maxPadj < thresh_padj] <- "B16F10-GP Specific"
  dfx$group[dfx$log2FC_CvA > thresh_lfc  &  dfx$log2FC_TvA > thresh_lfc & dfx$maxPadj < thresh_padj] <- "Shared Up"
  dfx$group[dfx$log2FC_CvA < -thresh_lfc  &  dfx$log2FC_TvA < -thresh_lfc & dfx$maxPadj < thresh_padj] <- "Shared Down"
  dfx$group <- factor(dfx$group, levels = c("Cl13 Specific", "B16F10-GP Specific", "Shared Up", "Shared Down", "Non-DEGs"))
  dfx$log2FC <- ifelse(dfx$group=="Cl13 Specific", dfx$log2FC_CvA, 
                           ifelse(dfx$group =="B16F10-GP Specific", dfx$log2FC_TvA, 
                                  rowMeans(dfx[,c('log2FC_TvA', 'log2FC_TvA')], na.rm=TRUE)))
  
  
  
  if (!is.null(geneHL)){
    gene_highlight <- geneHL
    is_highlight <- (dfx$gene %in% gene_highlight) & (dfx$group != "Non-DEGs")
    dfx$gene_highlight <- ifelse(is_highlight, dfx$gene, '')
  }
  
  return(dfx)
}

##  Double Differential Gene Expression: Seurat Version
getDoubleDGE_Seurat <- function(CdPD1, TdPD1, maxlfc=3, thresh_lfc=0.2, thresh_padj=0.05, geneHL=NULL){
  colnames(CdPD1) <- c("p_val", "log2FC", "pct.1", "pct.2", "padj")
  colnames(TdPD1) <- c("p_val", "log2FC", "pct.1", "pct.2", "padj")
  colnames(CdPD1) <- paste0(colnames(CdPD1), "_CdPD1")
  colnames(TdPD1) <- paste0(colnames(TdPD1), "_TdPD1")
  CdPD1$gene <- row.names(CdPD1)
  TdPD1$gene <- row.names(TdPD1)
  dfx <- merge(CdPD1, TdPD1, by="gene", all=F)
  dfx$log2FC_CdPD1[dfx$log2FC_CdPD1 > maxlfc] <- maxlfc
  dfx$log2FC_TdPD1[dfx$log2FC_TdPD1 > maxlfc] <- maxlfc
  dfx$log2FC_CdPD1[dfx$log2FC_CdPD1 < -maxlfc] <- -maxlfc
  dfx$log2FC_TdPD1[dfx$log2FC_TdPD1 < -maxlfc] <- -maxlfc
  dfx$group <- "Non-DEGs"
  dfx$maxPadj <- pmax(dfx$padj_CdPD1, dfx$padj_TdPD1)
  dfx$group[abs(dfx$log2FC_CdPD1)> thresh_lfc  &  abs(dfx$log2FC_TdPD1)< thresh_lfc & dfx$maxPadj < thresh_padj] <- "Cl13 Specific"
  dfx$group[abs(dfx$log2FC_CdPD1)< thresh_lfc  &  abs(dfx$log2FC_TdPD1)> thresh_lfc & dfx$maxPadj < thresh_padj] <- "B16F10-GP Specific"
  dfx$group[dfx$log2FC_CdPD1 > thresh_lfc  &  dfx$log2FC_TdPD1 > thresh_lfc & dfx$maxPadj < thresh_padj] <- "Shared Up"
  dfx$group[dfx$log2FC_CdPD1 < -thresh_lfc  &  dfx$log2FC_TdPD1 < -thresh_lfc & dfx$maxPadj < thresh_padj] <- "Shared Down"
  dfx$group <- factor(dfx$group, levels = c("Cl13 Specific", "B16F10-GP Specific", "Shared Up", "Shared Down", "Non-DEGs"))
  dfx$log2FC <- ifelse(dfx$group=="Cl13 Specific", dfx$log2FC_CdPD1, 
                       ifelse(dfx$group =="B16F10-GP Specific", dfx$log2FC_TdPD1, 
                              rowMeans(dfx[,c('log2FC_TdPD1', 'log2FC_TdPD1')], na.rm=TRUE)))
  
  
  
  if (!is.null(geneHL)){
    gene_highlight <- geneHL
    is_highlight <- (dfx$gene %in% gene_highlight) & (dfx$group != "Non-DEGs")
    dfx$gene_highlight <- ifelse(is_highlight, dfx$gene, '')
  }
  
  return(dfx)
}



## Inspect doubleDEG to identity genes to highlight
inspectDoubleDEG <- function(ddegDfx) {
  geneSpecList <- list()
  gseaTempList <- list()
  
  for (gpx in c("Cl13 Specific", "B16F10-GP Specific", "Shared Up", "Shared Down")){
    dfspec <- subset(ddegDfx, group==gpx) %>% 
      dplyr::select(gene, log2FC, maxLogMean, ExprLvl) %>% 
      arrange(desc(log2FC))
    gseaTemp <- getGSEA(dfspec) %>% 
      filter(p.adjust < 0.05) %>% 
      arrange(desc(NES))
    geneSpecList[[gpx]] <- dfspec
    gseaTempList[[gpx]] <- gseaTemp
  }
  
  tempResList <- list(specGenes = geneSpecList, gseaTemp=gseaTempList)
}


## get gsea results
getGO <- function(degDfx){
  gseaRes <- myEnrichGO_msigdb(GetGenelist(degDfx))
  resDf <- do.call("rbind", lapply(gseaRes, function(sobj) {
    sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
    sobj@result
  }))
}



## get gsea results
getGSEA <- function(degDfx){
  gseaRes <- myEnrichGSEA_msigdb(GetGenelist(degDfx))
  resDf <- do.call("rbind", lapply(gseaRes, function(sobj) {
    sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
    sobj@result
  }))
}

## get gsea results
getMesh <- function(degDfx){
  gseaRes <- myEnrichMesh(GetGenelist(degDfx))
  resDf <- do.call("rbind", lapply(gseaRes, function(sobj) {
    sobj <- setReadable(sobj, 'org.Hs.eg.db', keyType = "ENTREZID")
    sobj@result
  }))
}




##============================-- ggplot utility  --=======================================##
saveFig <- function(pobj, prefix, wd, ht, drawPDF=T){
  png(filename = paste0(prefix, ".png"), width = wd, height = ht, res = 300, units = "px")
  print(pobj)
  dev.off()
  
  if (drawPDF){
    ggsave(plot = pobj, filename = paste0(prefix, ".pdf"), width = wd, height = ht,
           units = "px", device = "pdf") 
  }
}