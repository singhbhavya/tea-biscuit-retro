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
       colkey = c("Brain - Frontal Cortex (BA9)" = "orange", 
                  "Brain - Hippocampus" = "lightblue"),
       legendPosition = "right")  +
  theme_cowplot()

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

biplot(brain.gh.cor.pca.obj, 
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
       colkey = c("Brain - Frontal Cortex (BA9)" = "orange", 
                  "Brain - Hippocampus" = "lightblue"),
       legendPosition = "right")  +
  theme_cowplot()

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

biplot(brain.herv.cor.pca.obj, 
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
       colkey = c("Brain - Frontal Cortex (BA9)" = "orange", 
                  "Brain - Hippocampus" = "lightblue"),
       legendPosition = "right")  +
  theme_cowplot()

save(brain.herv.cor.dds, brain.herv.cor.pca.obj, file = "r_outputs/brain.herv.cor.dds.pca.Rdata")

################################ SET THRESHOLDS ################################

lfc.cutoff <- 1.5
pval=0.001 # p value threshold

