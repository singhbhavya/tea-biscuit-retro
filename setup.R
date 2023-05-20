################################################################################
################################################################################
################################################################################
################################################################################
##################################### SETUP ####################################

setwd("/efs/projects/tea-biscuit-retro/")
# set working directory to where your directory is

################################ LOAD PACKAGES #################################

library(tidyverse)
library(sva)

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

############################### BATCH CORRECTION ###############################

ischemic.time = ifelse(brain.samples$SMTSISCH <= 300, 1, 
                          ifelse(brain.samples$SMTSISCH > 300 & brain.samples$SMTSISCH <= 600, 2,
                                 ifelse(brain.samples$SMTSISCH > 600 & brain.samples$SMTSISCH <= 900, 3,
                                        ifelse(brain.samples$SMTSISCH > 900 & brain.samples$SMTSISCH <= 1200, 4, 5)))) 


tissue = sapply(as.character(brain.samples$SMTSD),
            switch, "Brain - Frontal Cortex (BA9)" = 1,
            "Brain - Hippocampus" = 2,
            USE.NAMES = F)

BL.counts.comb.corrected = ComBat_seq(counts = as.matrix(brain.comb),
                            batch = ischemic.time,
                            group = tissue,
                            full_mod = TRUE)

