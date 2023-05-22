################################################################################
################################################################################
################################################################################
################################################################################
##################################### SETUP ####################################

setwd("/efs/projects/tea-biscuit-retro/")
# set working directory to where your directory is

################################ LOAD PACKAGES #################################

library(tidyverse)
library(data.table)
library(PCAtools)
library(dplyr)
library(matrixStats)
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)
library(sva)
library(ashr)
library(cowplot)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(pheatmap)
library(scales)

################################## LOAD REFS ###################################

# This file contains the TE annotation, and the mapping of gene IDs to symbols
load("refs/refs.Rdata")

################################## LOAD DATA ###################################

# read frontal cortex samples
fc.counts.tx <- readRDS("data/results.BRAIN-frontal_cortex/counts.gene.rds")
fc.counts.rtx <- readRDS("data/results.BRAIN-frontal_cortex/counts.retro.rds")
fc.samples <- readRDS("data/results.BRAIN-frontal_cortex/samples.rds")

# read hippocampus samples
hippo.counts.tx <- readRDS("data/results.BRAIN-hippocampus/counts.gene.rds")
hippo.counts.rtx <- readRDS("data/results.BRAIN-hippocampus/counts.retro.rds")
hippo.samples <- readRDS("data/results.BRAIN-hippocampus/samples.rds")

# sanity check that all the TEs in the count data are in the annotation
stopifnot(all(rownames(hippo.counts.rtx) == retro.annot.v2$locus))
# sanity check that the transcriptoeme and retrotranscriptome have the same samples
stopifnot(all(names(fc.counts.tx) == names(fc.counts.rtx)))

# sanity check that all the TEs in the count data are in the annotation
stopifnot(all(rownames(fc.counts.rtx) == retro.annot.v2$locus))
# sanity check that the transcriptoeme and retrotranscriptome have the same samples
stopifnot(all(names(fc.counts.tx) == names(fc.counts.rtx)))

############################### COMBINE METADATA  ##############################

# extract common metadata
common.cols <- intersect(colnames(fc.samples), colnames(hippo.samples))

# combine sample sheets by common metadata
brain.samples <- rbind(fc.samples[common.cols], hippo.samples[common.cols])
# set rownames
rownames(brain.samples) <- brain.samples$SAMPID

################################ COMBINE SAMPLES ###############################

# combine .tx and .rtx counts for all lymphoma samples + sanity check
brain.rtx <- cbind(fc.counts.rtx, hippo.counts.rtx)
brain.tx <- cbind(fc.counts.tx, hippo.counts.tx)

brain.comb <- rbind(brain.tx, brain.rtx)
stopifnot(all(names(brain.tx) == names(brain.rtx)))
stopifnot(all(names(brain.tx) == rownames(brain.samples)))

brain.herv <- brain.comb[rownames(retro.annot.v2[retro.annot.v2$Class == "HERV",]),]
brain.l1 <- brain.comb[rownames(retro.annot.v2[retro.annot.v2$Class == "L1",]),]

############################### BATCH CORRECTION ###############################

ischemic.time = ifelse(brain.samples$SMTSISCH <= 300, 1, 
                          ifelse(brain.samples$SMTSISCH > 300 & brain.samples$SMTSISCH <= 600, 2,
                                 ifelse(brain.samples$SMTSISCH > 600 & brain.samples$SMTSISCH <= 900, 3,
                                        ifelse(brain.samples$SMTSISCH > 900 & brain.samples$SMTSISCH <= 1200, 4, 5)))) 


tissue = sapply(as.character(brain.samples$SMTSD),
            switch, "Brain - Frontal Cortex (BA9)" = 1,
            "Brain - Hippocampus" = 2,
            USE.NAMES = F)

brain.comb.corrected = ComBat_seq(counts = as.matrix(brain.comb),
                            batch = ischemic.time,
                            group = tissue,
                            full_mod = TRUE)

save(brain.comb.corrected, file="brain.comb.correct.Rds")

############################# SUBSET HERVs and L1s #############################

brain.comb.corrected <- as.data.frame(brain.comb.corrected)

brain.cor.tx <- as.data.frame(brain.comb.corrected[rownames(brain.tx),])
brain.cor.rtx <- as.data.frame(brain.comb.corrected[rownames(brain.rtx),])
brain.cor.herv <- as.data.frame(brain.cor.rtx[rownames(retro.annot.v2[retro.annot.v2$Class == "HERV",]),])
brain.cor.l1 <- as.data.frame(brain.cor.rtx[rownames(retro.annot.v2[retro.annot.v2$Class == "L1",]),])

########################### FILTER ORIGINAL DATA ###############################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(brain.comb) * 0.015)


brain.filt.tx <- brain.tx[rowSums(brain.tx > cutoff.count) > cutoff.samp, ]
brain.filt.rtx <- brain.rtx[rowSums(brain.rtx > cutoff.count) > cutoff.samp, ]
brain.filt.comb <- rbind(brain.filt.tx, brain.filt.rtx)
brain.filt.herv <- brain.herv[rowSums(brain.herv > cutoff.count) > cutoff.samp, ]
brain.filt.l1 <- brain.l1[rowSums(brain.l1 > cutoff.count) >cutoff.samp, ]

save(brain.filt.tx, brain.filt.rtx, brain.filt.comb, brain.filt.herv, brain.filt.l1,
     brain.samples, file="r_outputs/brain.filt.Rdata")

############################ FILTER CORRECTED DATA ############################# 

cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(brain.comb.corrected) * 0.015)


brain.filt.cor.tx <- brain.cor.tx[rowSums(brain.cor.tx > cutoff.count) > cutoff.samp, ]
brain.filt.cor.rtx <- brain.cor.rtx[rowSums(brain.cor.rtx > cutoff.count) > cutoff.samp, ]
brain.filt.cor.comb <- rbind(brain.filt.cor.tx, brain.filt.cor.rtx)
brain.filt.cor.herv <- brain.cor.herv[rowSums(brain.cor.herv > cutoff.count) > cutoff.samp, ]
brain.filt.cor.l1 <- brain.cor.l1[rowSums(brain.cor.l1 > cutoff.count) >cutoff.samp, ]

save(brain.filt.cor.tx, brain.filt.cor.rtx, brain.filt.cor.comb, brain.filt.cor.herv, brain.filt.cor.l1,
     brain.samples, file="r_outputs/brain.filt.cor.Rdata")

############################# COUNT READS PER TYPE ############################# 

te_percent <- function(herv.df, rtx.df, comb.df, metadata, metadata.col) {
  
  herv.reads <- as.data.frame(colSums(herv.df))
  te.reads <- as.data.frame(colSums(rtx.df))
  all.reads <- as.data.frame(colSums(comb.df))
  
  colnames(herv.reads) <- c("reads")
  colnames(te.reads) <- c("reads")
  colnames(all.reads) <- c("reads")
  
  herv.reads$sample  <- rownames(herv.reads)
  te.reads$sample <- rownames(te.reads)
  all.reads$sample <- rownames(all.reads)
  
  herv.reads$type <- metadata[[metadata.col]]
  te.reads$type <- metadata[[metadata.col]]
  all.reads$type <- metadata[[metadata.col]]
  
  stopifnot(all(all.reads$sample == te.reads$sample))
  
  te.reads$proportion <- te.reads$reads/all.reads$reads*100
  herv.reads$proportion <- herv.reads$reads/all.reads$reads*100
  
  output <- list(herv.reads = herv.reads,
                 te.reads = te.reads)
  
  return(output)
}

###################### COUNT READS PER AGIRRE TISSUE TYPE ###################### 

brain.te.percent <-
  te_percent(brain.filt.cor.herv, brain.filt.cor.rtx, brain.filt.cor.comb,
             brain.samples, "SMTSD")


brain.te.percent$herv.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

brain.te.percent$te.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

pdf("plots/brain_te_percent.pdf", height=3, width=3)
brain.te.percent$te.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of TE Fragments") +
  ylim(0, 0.75) +
  scale_fill_manual(values = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1],
                               "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7])) + 
  scale_x_discrete(labels=c("Brain - Frontal Cortex (BA9)" = "FC", 
                            "Brain - Hippocampus" = "HC")) +
  theme(aspect.ratio = 1) 
dev.off()

pdf("plots/brain_herv_percent.pdf", height=3, width=3)
brain.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of HERV Fragments") +
  ylim(0, 0.75) +
  scale_fill_manual(values = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1],
                               "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7])) + 
  scale_x_discrete(labels=c("Brain - Frontal Cortex (BA9)" = "FC", 
                            "Brain - Hippocampus" = "HC")) +
  theme(aspect.ratio = 1)
dev.off()

################################# FUNCTION SCREE PLOT #################################

pca_standard <- function(tform, metadata, var) {
  
  removeVar <- var
  pca.obj <- PCAtools::pca(assay(tform), 
                           metadata=metadata, 
                           removeVar=removeVar)
  
  cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
              removeVar*100, length(pca.obj$xvars)))
  
  varline <- 50
  varline.x <- min(which(cumsum(pca.obj$variance) >= varline))
  
  horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
  elbow <- PCAtools::findElbowPoint(pca.obj$variance)
  
  screeplot <-PCAtools::screeplot(pca.obj,
                                  axisLabSize = 6,
                                  components = getComponents(pca.obj, 1:30),
                                  title=paste("Retrotranscriptome SCREE",
                                              metadata$cancer_type[1],
                                              sep=" "),
                                  hline=varline, vline=c(varline.x, horn$n, elbow)
  ) +
    geom_label(aes(x=varline.x+1, y=50, 
                   label = paste0(varline, '% var'), vjust = -1)) +
    geom_label(aes(x = horn$n + 1, y = 50,
                   label = 'Horn\'s', vjust = -1)) +
    geom_label(aes(x = elbow + 1, y = 50,
                   label = 'Elbow method', vjust = -1))
  
  
  cat(sprintf('%d PCs for Elbow method\n', elbow))
  cat(sprintf('%d PCs for Horn method\n', horn$n))
  cat(sprintf('%d PCs needed to explain %d percent of variation\n', 
              varline.x, varline))
  
  print(screeplot)
  
  return(pca.obj)
}


####################### DIFFERENTIAL EXPRESSION ORIGINAL ####################### 

brain.gh.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.comb,
                                              colData = brain.samples,
                                              design = ~ SMTSD)

brain.gh.dds <- DESeq2::DESeq(brain.gh.dds, parallel=T)
brain.gh.tform <- DESeq2::varianceStabilizingTransformation(brain.gh.dds, 
                                                           blind=FALSE)

brain.gh.pca.obj <-
  pca_standard(tform = brain.gh.tform, 
               metadata = brain.samples, 
               var = 0.1)

all(rownames(brain.gh.pca.obj$loadings) %in% rownames(gene_table))
rownames(brain.gh.pca.obj$loadings) <- 
  gene_table[rownames(brain.gh.pca.obj$loadings), 'display']

pdf("plots/brain_pca_genes_hervs_uncorrected.pdf", height=5, width=7)
biplot(brain.gh.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "SMTSD",
       colkey = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1], 
                  "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7]),
       legendPosition = "right")  +
  theme(aspect.ratio = 1) +
  theme_cowplot()
dev.off()

save(brain.gh.dds, brain.gh.pca.obj, file = "r_outputs/brain.gh.dds.pca.Rdata")

####################### DIFFERENTIAL EXPRESSION CORRECTED ###################### 

brain.gh.cor.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.cor.comb,
                                               colData = brain.samples,
                                               design = ~ SMTSD)

brain.gh.cor.dds <- DESeq2::DESeq(brain.gh.cor.dds, parallel=T)
brain.gh.cor.tform <- DESeq2::varianceStabilizingTransformation(brain.gh.cor.dds, 
                                                            blind=FALSE)

brain.gh.cor.pca.obj <-
  pca_standard(tform = brain.gh.cor.tform, 
               metadata = brain.samples, 
               var = 0.1)

all(rownames(brain.gh.cor.pca.obj$loadings) %in% rownames(gene_table))
rownames(brain.gh.cor.pca.obj$loadings) <- 
  gene_table[rownames(brain.gh.cor.pca.obj$loadings), 'display']

pdf("plots/brain_pca_genes_hervs_corrected.pdf", height=5, width=7)
biplot(brain.gh.cor.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "SMTSD",
       colkey = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1], 
                  "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7]),
       legendPosition = "right")  +
  theme_cowplot()
dev.off()

save(brain.gh.cor.dds, brain.gh.cor.pca.obj, file = "r_outputs/brain.gh.cor.dds.pca.Rdata")


#################### DIFFERENTIAL EXPRESSION CORRECTED HERV ####################

brain.herv.cor.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.cor.herv,
                                                   colData = brain.samples,
                                                   design = ~ SMTSD)

brain.herv.cor.dds <- DESeq2::DESeq(brain.herv.cor.dds, parallel=T)
brain.herv.cor.tform <- DESeq2::varianceStabilizingTransformation(brain.herv.cor.dds, 
                                                                blind=FALSE)

brain.herv.cor.pca.obj <-
  pca_standard(tform = brain.herv.cor.tform, 
               metadata = brain.samples, 
               var = 0.1)

pdf("plots/brain_pca_hervs_corrected.pdf", height=5, width=7)
biplot(brain.herv.cor.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "SMTSD",
       colkey = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1], 
                  "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7]),
       legendPosition = "right")  +
  theme_cowplot()
dev.off()

save(brain.herv.cor.dds, brain.herv.cor.pca.obj, file = "r_outputs/brain.herv.cor.dds.pca.Rdata")

################################ SET THRESHOLDS ################################

lfc.cutoff <- 1.5
pval=0.001 # p value threshold

############################# TOP GENES / HERVS ################################

brain.res <- DESeq2::results(brain.gh.cor.dds, contrast=c("SMTSD", "Brain - Hippocampus", 
                                                          "Brain - Frontal Cortex (BA9)"), alpha=pval)

brain.res$display <- gene_table[rownames(brain.res),]$display
brain.res$class <- gene_table[rownames(brain.res),]$gene_type

sig.gh <- subset(brain.res, padj < pval & abs(log2FoldChange) > lfc.cutoff)
sig.gh <- sig.gh[order(sig.gh$padj),]


############################### TOP HERVs ONLY #################################

brain.res.herv <- DESeq2::results(brain.herv.cor.dds, contrast=c("SMTSD", "Brain - Hippocampus", 
                                         "Brain - Frontal Cortex (BA9)"), alpha=pval)

brain.res.herv$display <- gene_table[rownames(brain.res.herv),]$display
brain.res.herv$class <- gene_table[rownames(brain.res.herv),]$gene_type

sig.hervs <- subset(brain.res.herv, padj < pval & abs(log2FoldChange) > lfc.cutoff)
sig.hervs <- sig.hervs[order(sig.hervs$padj),]

################################ LIST TOP HERVS ################################ 

upvars <- rownames(subset(sig.hervs, log2FoldChange>0))
downvars <- rownames(subset(sig.hervs, log2FoldChange<0))

################################# HEATMAPS #####################################

################################### SETUP ######################################

cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(brain.herv.cor.dds)[,c("SMTS","SMTSD", "SMGEBTCH")])
df <- subset(df, select = -c(1))

# Create colors for each group
annoCol <- c(pal_jco("default", alpha = 0.8)(7)[1], pal_jco("default", alpha = 0.8)(7)[7])

names(annoCol) <- unique(brain.samples$SMTSD)
annoCol <- list(SMTSD = annoCol)


######################### UPREGULATED IN ALL GROUPS ############################

annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
annoRow <- annoRow[rownames(sig.hervs),]
annoRow <- subset(annoRow, select = -c(2))

pdf("plots/brain_sig_hervs_heatmap.pdf", height=10, width=10)
pheatmap(assay(brain.herv.cor.tform)[rownames(sig.hervs),], 
         main="Differentially Expressed HERVs",
         cluster_rows=TRUE,
         show_rownames=FALSE,
         show_colnames = FALSE,
         color = cols,
         scale="row",
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[rownames(sig.hervs),]$display,
         cluster_cols=TRUE, 
         treeheight_row=0,
         annotation_col=df,
         annotation_row=annoRow,
         annotation_colors = annoCol)
dev.off()

############################# VOLCANO DZ VS LZ #################################

pdf("plots/brain_herv_volcano.pdf", height=7, width=7)
EnhancedVolcano(sig.hervs,
                lab = rownames(sig.hervs),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Frontal Cortex vs Hippocampus')

dev.off()

########################## PLOT INDIVIDUAL HERVs ###############################

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="SMTSD", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=SMTSD, y=count, fill=SMTSD))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cell Type") +
    ylab("Counts") +
    scale_x_discrete(labels=c("Brain - Frontal Cortex (BA9)" = "FC", 
                              "Brain - Hippocampus" = "HC")) +
    scale_fill_manual(values = c("Brain - Frontal Cortex (BA9)" = pal_jco("default", alpha = 0.8)(7)[1],
                                 "Brain - Hippocampus" = pal_jco("default", alpha = 0.8)(7)[7])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene) + 
    scale_y_log10(labels = label_comma()) +
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    stat_compare_means(comparisons = list(c("Brain - Frontal Cortex (BA9)", "Brain - Hippocampus")),
                       method = "t.test", 
                       label = "p.signif",
                       hide.ns = FALSE)
}


pdf("plots/brain_HERVH_12q13.2b.pdf", height=4, width=3)
plot.counts(brain.herv.cor.dds, "HERVH_12q13.2b")
dev.off()

pdf("plots/brain_ERVLB4_20q13.12a.pdf", height=4, width=3)
plot.counts(brain.herv.cor.dds, "ERVLB4_20q13.12a")
dev.off()

############################### FAMILY LEVEL UP ############################### 

upreg.hervs.df <- as.data.frame(sig.hervs)
upreg.hervs.df$upregin <- ifelse(upreg.hervs.df$log2FoldChange > 0, "Hippocampus", "Frontal Cortex")
upreg.hervs.df$family <- retro.annot.v2$Family[match(upreg.hervs.df$display, 
                                                     retro.annot.v2$Locus)]

upreg.families <-
  upreg.hervs.df %>% dplyr::count(family, upregin, sort = TRUE) 

pdf("plots/brain_family_level.pdf", height=8, width=6)
ggplot(upreg.families, aes(fill=reorder(family, -n), y=upregin, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot.v2$Family),
                    labels = unique(retro.annot.v2$Family)) + 
  coord_flip() +
  theme_cowplot() +  
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  guides(fill=guide_legend(title="HERV Family")) +
  ylab(NULL) +
  xlab("Number of HERV Loci") + 
  theme(legend.position = c("right"),
        plot.margin = margin(10, 10, 10, 40),
        axis.line=element_blank()) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

################################### PATHWAYS ###################################

pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")

################################ FGSEA FUNCTION ################################

make.fsgsea <- function(pathway, fgsea.res, clust_name, pathway_name) {
  
  fgsea.res$SYMBOL <- gene_table[rownames(fgsea.res),]$display
  
  fgsea.res <- as.data.frame(fgsea.res) %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))
  
  fgsea.ranks <- deframe(fgsea.res)
  
  fgsea.out <- fgsea(pathways=pathway, 
                     stats=fgsea.ranks, 
                     nPermSimple = 10000,
                     eps=0)
  return(fgsea.out)
  
  assign(paste0(clust_name, ".", pathway_name, ".fgsea.out"), fgsea.out, envir = .GlobalEnv )
}

################################ HALLMARK FGSEA ################################

fsgsea.hallmarks <- list(
  "Hippocampus" = make.fsgsea(pathways.hallmark, brain.res,"hallmark"))

fgseaResTidy <- fsgsea.hallmarks[["Hippocampus"]] %>%
  as_tibble() %>%
  arrange(desc(NES))

rownames(fgseaResTidy) <- fgseaResTidy$pathway

pdf("plots/brain_hallmark.pdf", height=10, width=7)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hippocampus") + 
  theme_cowplot()
dev.off()

# Get longform df
fsgsea.hallmarks.summary <- rbindlist(fsgsea.hallmarks, idcol = "index")

# Recode df
fsgsea.hallmarks.summary$pathway <- gsub("HALLMARK_","",fsgsea.hallmarks.summary$pathway)
fsgsea.hallmarks.summary$pathway <- gsub("_"," ",fsgsea.hallmarks.summary$pathway)

# Bubble plot
pdf("plots/brain_hallmark_bubble.pdf", height=8.5, width=6)
ggplot(fsgsea.hallmarks.summary, aes(x = index, 
                                     y = pathway, 
                                     size = -log(padj), 
                                     color = NES)) +
  geom_point() +
  scale_size(name = "-log (P value)", range = c(1, 10)) + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  scale_colour_gradientn(colors = viridis_pal()(10)) +
  xlab("Hippocampus") +
  ylab("Hallmark Pathway") 
dev.off()
