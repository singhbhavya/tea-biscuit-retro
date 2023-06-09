# Welcome to Day 2: Bulk Retrotranscriptomics


Our goals for today are the following:
1. Combine data from two different datasets (GTEx Frontal Cortex and GTEx Hippocampus)
2. Conduct batch correction on the combined data
3. Filter out lowly expressed genes and HERVs
4. Use PCA to compare batch-corrected filtered data to the original filtered data
5. Compare PCA sample clustering when based on genes, versus when based on only HERVs
6. Calculate % TE and % HERV in the two tissue types
7. Conduct differential expression (DE) analysis to identify DE HERVs.
8. Create heatmaps and volcano plots to identify DE HERVs
9. Plot individual HERV counts and conduct t-tests to compare means between tissue types
10. Identify family-level HERV abundance per tissue type
11. Gene-set enrichment analysis with Hallmark pathways

## Let's begin!

First, we will set our working directories and load some packages. You should set your working directory to wherever you have downloaded or symlinked the "tea-biscuit-retro" directory. To figure out where your directory was downloaded, you can go into the directory and type "pwd".Second, we will be loading some packages required for today's tutorials. All packages should already be downloaded in R, so you shouldn't have any issues here. If you do, let me know!

```
# Set working directory
setwd("/efs/users/[YOURUSERNAME]/tea-biscuit-user/")

# Load packages
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
```

## Loading References

```
# This file contains the Telescope TE annotation, and the mapping of gene IDs to symbols
load("refs/refs.Rdata")

# View the first few rows of this file: 
head(retro.annot.v2)

```

## Loading Data

Now, we will be loading the frontal cortex and hippocampus sample matrices.

```
# Read frontal cortex samples
fc.counts.tx <- readRDS("data/results.BRAIN-frontal_cortex/counts.gene.rds")
fc.counts.rtx <- readRDS("data/results.BRAIN-frontal_cortex/counts.retro.rds")
fc.samples <- readRDS("data/results.BRAIN-frontal_cortex/samples.rds")

# Read hippocampus samples
hippo.counts.tx <- readRDS("data/results.BRAIN-hippocampus/counts.gene.rds")
hippo.counts.rtx <- readRDS("data/results.BRAIN-hippocampus/counts.retro.rds")
hippo.samples <- readRDS("data/results.BRAIN-hippocampus/samples.rds")

# Look at the counts
#head(fc.counts.tx)

# Look at the sample sheet
#head(fc.samples)
```

Now that the data is all loaded, we want to combine the metadata.

```
# Extract common metadata columns
common.cols <- intersect(colnames(fc.samples), colnames(hippo.samples))

print(common.cols)

# Combine sample sheets by common metadata
brain.samples <- rbind(fc.samples[common.cols], hippo.samples[common.cols])
# Set rownames
rownames(brain.samples) <- brain.samples$SAMPID
```

And now, combine the samples:

```
# Combine .tx and .rtx counts for frontal cortex and hippocampus respectively 
brain.rtx <- cbind(fc.counts.rtx, hippo.counts.rtx)
brain.tx <- cbind(fc.counts.tx, hippo.counts.tx)

# Combine the .tx and .rtx counts together 
brain.comb <- rbind(brain.tx, brain.rtx)

# Sanity check. All the samples .tx, .rtx, and sample sheets should match@
stopifnot(all(names(brain.tx) == names(brain.rtx)))
stopifnot(all(names(brain.tx) == rownames(brain.samples)))

# Extract the HERV and LINE1 elements from all samples
brain.herv <- brain.comb[rownames(retro.annot.v2[retro.annot.v2$Class == "HERV",]),]
brain.l1 <- brain.comb[rownames(retro.annot.v2[retro.annot.v2$Class == "L1",]),]

head(brain.herv)
```

## Batch Correction

GTEx consortium datasets are known to have batch issues, such as ischemic time, the actual batch that the samples were processed in, gender, etc. For the sake of keeping this brief, today we will be batch correcting for just ischemic time. SVA, which is the batch correction tool we are using, only takes discrete variables. Since ischemic time is a continuous variable, we will be dividing it into chunks.

```
max(brain.samples$SMTSISCH) 
```

The max ischemic.time is 1360, meaning we will need 5 chunks.

```
ischemic.time = ifelse(brain.samples$SMTSISCH <= 300, 1, 
                          ifelse(brain.samples$SMTSISCH > 300 & 
                                   brain.samples$SMTSISCH <= 600, 2,
                                 ifelse(brain.samples$SMTSISCH > 600 &
                                          brain.samples$SMTSISCH <= 900, 3,
                                        ifelse(brain.samples$SMTSISCH > 900 &
                                                 brain.samples$SMTSISCH <= 1200, 4, 5)))) 


print(ischemic.time)
```

Our variable of interest is the tissue (frontal cortex vs. hippocampus), so we will be accounting for that in our batch correction as well.

```
tissue = sapply(as.character(brain.samples$SMTSD),
            switch, "Brain - Frontal Cortex (BA9)" = 1,
            "Brain - Hippocampus" = 2,
            USE.NAMES = F)
```

Now time for the actual batch correction! You may run the following code. It will take a little bit of time!

```
brain.comb.corrected = ComBat_seq(counts = as.matrix(brain.comb),
                            batch = ischemic.time,
                            group = tissue,
                            full_mod = TRUE)

# If you encountered issues in running the above, you can magically load the output below: 
load("r_outputs/brain.comb.correct.Rdata")
```

## Set Up HERV and L1 Matrices

```
# Convert the corrected, batch corrected combined matrix to a dataframe
brain.comb.corrected <- as.data.frame(brain.comb.corrected)

# Re-create subsets of the corrected data, including the .tx, .rtx, .herv., and L1 matrices
brain.cor.tx <- as.data.frame(brain.comb.corrected[rownames(brain.tx),])
brain.cor.rtx <- as.data.frame(brain.comb.corrected[rownames(brain.rtx),])
brain.cor.herv <- as.data.frame(brain.cor.rtx[rownames(retro.annot.v2[retro.annot.v2$Class == "HERV",]),])
brain.cor.l1 <- as.data.frame(brain.cor.rtx[rownames(retro.annot.v2[retro.annot.v2$Class == "L1",]),])
```

## Filter Original and Batch-Corrected data

First, let's filter the original data:

```
# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(brain.comb) * 0.015)


brain.filt.tx <- brain.tx[rowSums(brain.tx > cutoff.count) > cutoff.samp, ]
brain.filt.rtx <- brain.rtx[rowSums(brain.rtx > cutoff.count) > cutoff.samp, ]
brain.filt.comb <- rbind(brain.filt.tx, brain.filt.rtx)
brain.filt.herv <- brain.herv[rowSums(brain.herv > cutoff.count) > cutoff.samp, ]
brain.filt.l1 <- brain.l1[rowSums(brain.l1 > cutoff.count) >cutoff.samp, ]

# If you were unable to run the above, you can load the following data
load("r_outputs/brain.filt.Rdata")
```

And now the batch-corrected data:

```
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(brain.comb.corrected) * 0.015)

brain.filt.cor.tx <- brain.cor.tx[rowSums(brain.cor.tx > cutoff.count) > cutoff.samp, ]
brain.filt.cor.rtx <- brain.cor.rtx[rowSums(brain.cor.rtx > cutoff.count) > cutoff.samp, ]
brain.filt.cor.comb <- rbind(brain.filt.cor.tx, brain.filt.cor.rtx)
brain.filt.cor.herv <- brain.cor.herv[rowSums(brain.cor.herv > cutoff.count) > cutoff.samp, ]
brain.filt.cor.l1 <- brain.cor.l1[rowSums(brain.cor.l1 > cutoff.count) >cutoff.samp, ]


# If you were unable to run the above, you can load the following data
load("r_outputs/brain.filt.cor.Rdata")
```

**Question: Is there a difference in the number of genes, HERVs, and LINEs preserved after filtering in the original and corrected data?**

## Principal Component Analysis (PCA)

We will now be visualizing the corrected and uncorrected data based on a gene-and-HERV-based PCA. However, since PCA requires transformed data matrices as input, we need to run DESeq2 to generate those matrices. We won't be conducting differential expression analysis quite yet, but we will run the main program.

The following is the function we will use to create PCA objects:

```
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

```

Now for DESeq2 on our original data:

```

# Note that we are using the filtered and combined data matrix here. This means that the matrix contains both transcriptomic and retrotranscriptomic data, from both frontal cortex and hippocampus. 

# Generating the DESeq dataset:
brain.gh.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.comb,
                                              colData = brain.samples,
                                              design = ~ SMTSD)

# Run DESeq:
brain.gh.dds <- DESeq2::DESeq(brain.gh.dds, parallel=T)

# Create variance-transformed matrix (for PCA)
brain.gh.tform <- DESeq2::varianceStabilizingTransformation(brain.gh.dds, 
                                                           blind=FALSE)
```

With the tform matrix, we can now run PCA on the original data:

```
# PCA on the original data, using the function from earlier 
brain.gh.pca.obj <-
  pca_standard(tform = brain.gh.tform, 
               metadata = brain.samples, 
               var = 0.1)

# Rename the gene IDS to names, based on the gene_table
all(rownames(brain.gh.pca.obj$loadings) %in% rownames(gene_table))
rownames(brain.gh.pca.obj$loadings) <- 
  gene_table[rownames(brain.gh.pca.obj$loadings), 'display']

# If you were unable to run any of the above, you can load the output data here:
load("r_outputs/brain.gh.dds.pca.Rdata")
```

Now we can plot the PCA!

```
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
```

Let's do DESeq2 on our corrected data:

```
# Generating the DESeq dataset:
brain.gh.cor.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.cor.comb,
                                               colData = brain.samples,
                                               design = ~ SMTSD)

# Run DESeq
brain.gh.cor.dds <- DESeq2::DESeq(brain.gh.cor.dds, parallel=T)

# Create variance-transformed matrix (for PCA)

brain.gh.cor.tform <- DESeq2::varianceStabilizingTransformation(brain.gh.cor.dds, 
                                                            blind=FALSE)
```

PCA on corrected data:

```
# PCA on the original data, using the function from earlier 
brain.gh.cor.pca.obj <-
  pca_standard(tform = brain.gh.cor.tform, 
               metadata = brain.samples, 
               var = 0.1)

# Rename the gene IDS to names, based on the gene_table
all(rownames(brain.gh.cor.pca.obj$loadings) %in% rownames(gene_table))
rownames(brain.gh.cor.pca.obj$loadings) <- 
  gene_table[rownames(brain.gh.cor.pca.obj$loadings), 'display']


# If you were unable to run any of the above, you can load the output data here:
load("r_outputs/brain.gh.cor.dds.pca.Rdata")
```

Now we can plot the PCA of our corrected data!

```
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
  ```
  
  **Question: Do the PCA plots look different for the corrected and uncorrected data? Why? Why not?**
  
  ## Gene-based vs. HERV-based PCA: Do the two tell us different things?
  
  We will be creating the same PCA plots as before, but instead of using the combined matrices, we will use ONLY HERV counts. Let's see if things look any different!
  
  ```
  # Generating the DESeq dataset:
brain.herv.cor.dds <- DESeq2::DESeqDataSetFromMatrix(countData = brain.filt.cor.herv,
                                                   colData = brain.samples,
                                                   design = ~ SMTSD)

# Run DESeq:
brain.herv.cor.dds <- DESeq2::DESeq(brain.herv.cor.dds, parallel=T)

# Create variance-transformed matrix (for PCA)
brain.herv.cor.tform <- DESeq2::varianceStabilizingTransformation(brain.herv.cor.dds, 
                                                                blind=FALSE)

# PCA based on batch-corrected HERVs: 
brain.herv.cor.pca.obj <-
  pca_standard(tform = brain.herv.cor.tform, 
               metadata = brain.samples, 
               var = 0.1)

# If you were unable to run any of the above, you can load the output data here:
load("r_outputs/brain.herv.cor.dds.pca.Rdata")
```

Now let's plot the PCA!

```
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
```

**Question: Do things look any different compared to the gene-based PCA? Why? Why not?**

## What % of HERVs and TEs are in our tissue types?

Function to count HERV and TE reads per tissue-type:

```
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
```

What % of the total reads are made up by TEs and HERVs?

```

# Use function to count reads
brain.te.percent <-
  te_percent(brain.filt.cor.herv, brain.filt.cor.rtx, brain.filt.cor.comb,
             brain.samples, "SMTSD")

# Calculate herv percent
brain.te.percent$herv.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

# Calculate te percent
brain.te.percent$te.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))
```

Let's plot our results!

```
# Plot TE percent
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

# Plot HERV percent 
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
```

**Question: Why do you think the TE and/or HERV percents would be different?**

## Extract Differential Expression Analysis Results

Here, we will be using the DESEq2 objects from earlier to identify and visualize differentially-expressed HERVs.

First, let's set our p-value and log2foldchange thresholds:

```
lfc.cutoff <- 1.5
pval=0.001 
```

Now, let's extract our significantly differentially-expressed genes and HERVs:

```
# Extract DE results, comparing Hippocampus to the Frontal Cortex. In this case, the hippocampus is the numerator, meaning that any genes with a POSITIVE log2fold change have higher expression in the hippocampus, and lower expression in the frontal cortex.

brain.res <- DESeq2::results(brain.gh.cor.dds, contrast=c("SMTSD", "Brain - Hippocampus", 
                                                          "Brain - Frontal Cortex (BA9)"),
                             alpha=pval)

# Add a column with gene names, because gene IDs are confusing
brain.res$display <- gene_table[rownames(brain.res),]$display

# Add a column with gene type. This will allow us to distinguish the HERVs.
brain.res$class <- gene_table[rownames(brain.res),]$gene_type

# Extract significant genes 
sig.gh <- subset(brain.res, padj < pval & abs(log2FoldChange) > lfc.cutoff)
# Order results by p-value
sig.gh <- sig.gh[order(sig.gh$padj),]
head(sig.gh)
```

Using our HERV-only DESEq object, let's extract just the HERVs:

```
# Extract DE results, comparing Hippocampus to the Frontal Cortex. In this case, the hippocampus is the numerator, meaning that any genes with a POSITIVE log2fold change have higher expression in the hippocampus, and lower expression in the frontal cortex.

brain.res.herv <- DESeq2::results(brain.herv.cor.dds, contrast=c("SMTSD", "Brain - Hippocampus", 
                                         "Brain - Frontal Cortex (BA9)"), alpha=pval)

brain.res.herv$display <- gene_table[rownames(brain.res.herv),]$display
brain.res.herv$class <- gene_table[rownames(brain.res.herv),]$gene_type

# Extract significant genes 
sig.hervs <- subset(brain.res.herv, padj < pval & abs(log2FoldChange) > lfc.cutoff)
# Order results by p-value
sig.hervs <- sig.hervs[order(sig.hervs$padj),]
head(sig.hervs)
```

List top HERVs:

```
upvars <- rownames(subset(sig.hervs, log2FoldChange>0)) # HERVs upregulated in hippocampus
downvars <- rownames(subset(sig.hervs, log2FoldChange<0)) # HERVs upregulated in frontal cortex
```

## Visualize DE Results: Heatmaps, Volcano Plots, Oh My!

Let's set up some colors for the heatmap:

```
# Heatmap colors
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(brain.herv.cor.dds)[,c("SMTS","SMTSD", "SMGEBTCH")])
df <- subset(df, select = -c(1))

# Create colors for each group
annoCol <- c(pal_jco("default", alpha = 0.8)(7)[1], pal_jco("default", alpha = 0.8)(7)[7])
 
# Create annotation column
names(annoCol) <- unique(brain.samples$SMTSD)
annoCol <- list(SMTSD = annoCol)
```

Okay, time to rock and roll.

```
# Set up row annotations
annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
annoRow <- annoRow[rownames(sig.hervs),]
annoRow <- subset(annoRow, select = -c(2))

# Make the heatmap:
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
```
  
**Question: What do you see? Why are some samples not clustering with their tissue-type?**

We can also visualize HERVs with the highest log2fold changes and p-values using volcano plots:

```
EnhancedVolcano(sig.hervs,
                lab = rownames(sig.hervs),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Frontal Cortex vs Hippocampus')
```

## Visualize Individual HERVs of Interest

Here is a custom function I created to plot HERV counts:

```
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
```

Let's pick out a few random HERVs from the volcano plot to visualize:

```
plot.counts(brain.herv.cor.dds, "HERVH_12q13.2b")

plot.counts(brain.herv.cor.dds, "ERVLB4_20q13.12a")
```

**Activity: Play around with the matrices and try to find more HERVs of interest. Do any of these HERVs show up in the literature?**

## A Broader Look: What about Family-Level?

Different HERV families can have different broad functions. Sometimes, it can be useful to take a look at the bigger picture.

```

# Count the number of HERV loci, and tally which families they belong to.
# Pull out significant HERVs
upreg.hervs.df <- as.data.frame(sig.hervs)
# HERVs with log2fold change > 0 are upregulated in the Hippocampus, while others are upregulated in the FC.
upreg.hervs.df$upregin <- ifelse(upreg.hervs.df$log2FoldChange > 0, "Hippocampus", "Frontal Cortex")
# Which families do these HERVs belong to?
upreg.hervs.df$family <- retro.annot.v2$Family[match(upreg.hervs.df$display, 
                                                     retro.annot.v2$Locus)]
# Count number of loci per family
upreg.families <-
  upreg.hervs.df %>% dplyr::count(family, upregin, sort = TRUE) 
```

Let's plot it!

```
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
  ```
  
  **Question: Which families are more prominent in one tissue type vs the other? Why could this be?**
  
  ## A brief introduction to gene-set enrichment analysis
  
  Import pathways:
  
  ```
  pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")
  head(pathways.hallmark)
  ```
  
  Here is a function I wrote to quickly do a fast gene-set enrichment analysis.
  
  ```
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
```

Let's run the FGSEA for the hallmark pathways!

```
# Run the function 
fsgsea.hallmarks <- list(
  "Hippocampus" = make.fsgsea(pathways.hallmark, brain.res,"hallmark"))

# Tidy the data
fgseaResTidy <- fsgsea.hallmarks[["Hippocampus"]] %>%
  as_tibble() %>%
  arrange(desc(NES))

# Let's take a look
fgseaResTidy

# Set rownames
rownames(fgseaResTidy) <- fgseaResTidy$pathway
```

This is a classic way to plot GSEA:

```
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hippocampus") + 
  theme_cowplot()
```

**Question: Do the upregulated and downregulated pathways make sense?**

If you'd like a slightly fancier way to plot the same, you can also do bubble plots!

```
# Get longform df
fsgsea.hallmarks.summary <- rbindlist(fsgsea.hallmarks, idcol = "index")

# Recode df
fsgsea.hallmarks.summary$pathway <- gsub("HALLMARK_","",fsgsea.hallmarks.summary$pathway)
fsgsea.hallmarks.summary$pathway <- gsub("_"," ",fsgsea.hallmarks.summary$pathway)

# Bubble plot
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

```

**Question: Do you think there is a benefit to plotting this way compared to the previous way?**

