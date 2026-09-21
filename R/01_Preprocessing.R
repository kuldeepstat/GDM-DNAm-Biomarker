# ==============================================================================
# 01_Preprocessing.R
#
# Purpose:
# Quality control, normalization, probe filtering, and preprocessing of
# Illumina DNA methylation data used in the GDM biomarker study.
#
# Main steps:
#   1. Load required R packages
#   2. Read sample sheet and raw IDAT files
#   3. Perform sample-level quality control using detection P-values
#   4. Normalize methylation data
#   5. Assess group and batch structure using MDS
#   6. Perform probe-level quality control
#   7. Assess post-filtering data structure
#   8. Extract beta values and M-values
#
# Note:
# The analytical commands in this script correspond to those used in the
# original study analysis.
# ==============================================================================


# ==============================================================================
# 1. Load required packages
# ==============================================================================

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(methylationArrayAnalysis)
library(maxprobes)


# ==============================================================================
# 2. Set data directory and read raw methylation data
# ==============================================================================

# Set up a path to the data directory
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")

# List the files
list.files(dataDirectory, recursive = TRUE)

# Obtain EPIC array annotation
ann450k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Read in the sample sheet for the experiment
targets <- read.metharray.sheet(
  dataDirectory,
  pattern = "SampleSheet4.csv"
)

targets

# Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets = targets)
rgSet

# Give the samples descriptive names
targets$ID <- paste(targets$Sample_Name)
sampleNames(rgSet) <- targets$ID
rgSet


# ==============================================================================
# 3. Sample-level quality control
# ==============================================================================

# Calculate the detection P-values
detP <- detectionP(rgSet)
head(detP)


# ------------------------------------------------------------------------------
# 3.1 Examine mean detection P-values across samples
# ------------------------------------------------------------------------------

pal <- brewer.pal(8, "Dark2")

par(mfrow = c(1, 2))

barplot(
  colMeans(detP),
  col = pal[factor(targets$Sample_Group)],
  las = 2,
  cex.names = 0.01,
  ylim = c(0, 1),
  ylab = "Mean detection p-values"
)

abline(h = 0.05, col = "red")

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  fill = pal,
  bg = "white",
  cex = 0.7
)


barplot(
  colMeans(detP[, 1:258]),
  col = pal[factor(targets$Sample_Group)],
  las = 2,
  cex.names = 0.5,
  ylim = c(0, 0.026),
  ylab = "Mean detection p-values"
)

abline(h = -log10(0.05), col = "red")

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  fill = pal,
  bg = "white",
  cex = 0.7
)


# ------------------------------------------------------------------------------
# 3.2 Generate QC report
# ------------------------------------------------------------------------------

qcReport(
  rgSet,
  sampNames = targets$ID,
  sampGroups = targets$Sample_Group,
  pdf = "qcReport.pdf"
)


# ------------------------------------------------------------------------------
# 3.3 Remove poor-quality samples
# ------------------------------------------------------------------------------

keep <- colMeans(detP) < 0.05

rgSet <- rgSet[, keep]
rgSet

targets <- targets[keep, ]
targets[, 1:5]

# Remove poor-quality samples from detection P-value table
detP <- detP[, keep]
dim(detP)


# ==============================================================================
# 4. Normalize methylation data
# ==============================================================================

# Normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)

# Create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)


# ------------------------------------------------------------------------------
# 4.1 Visualize data before and after normalization
# ------------------------------------------------------------------------------

par(mfrow = c(1, 2))

densityPlot(
  rgSet,
  sampGroups = targets$Sample_Group,
  main = "Raw",
  legend = FALSE
)

legend(
  "top",
  legend = levels(factor(targets$Sample_Group)),
  text.col = brewer.pal(7, "Dark2")
)


densityPlot(
  getBeta(mSetSq),
  sampGroups = targets$Sample_Group,
  main = "Normalized",
  legend = FALSE
)

legend(
  "top",
  legend = levels(factor(targets$Sample_Group)),
  text.col = brewer.pal(8, "Dark2")
)


# ==============================================================================
# 5. Pre-filtering assessment using multidimensional scaling
# ==============================================================================


# ------------------------------------------------------------------------------
# 5.1 Group-wise MDS
# ------------------------------------------------------------------------------

par(mfrow = c(1, 2))

plotMDS(
  getM(mSetSq),
  top = 700000,
  gene.selection = "common",
  main = "Group wise",
  pch = 16,
  col = pal[factor(targets$Sample_Group)]
)

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  text.col = pal,
  bg = "white",
  cex = 0.7
)


# ------------------------------------------------------------------------------
# 5.2 Define batch variable
# ------------------------------------------------------------------------------

batch <- substr(targets$Sample_Name, 11, 13)

batch[batch == "B1"] <- "One"
batch[batch == "B2"] <- "One"
batch[batch == "B3"] <- "One"
batch[batch == "B4"] <- "One"
batch[batch == "B5"] <- "One"
batch[batch == "B6"] <- "One"
batch[batch == "B7"] <- "One"

batch[batch == "B8"] <- "Two"
batch[batch == "B9"] <- "Two"
batch[batch == "B10"] <- "Two"
batch[batch == "B11"] <- "Two"
batch[batch == "B12"] <- "Two"


# ------------------------------------------------------------------------------
# 5.3 Batch-wise MDS
# ------------------------------------------------------------------------------

plotMDS(
  getM(mSetSq),
  top = 100000,
  gene.selection = "common",
  main = "Batch wise",
  pch = 16,
  col = as.factor(batch)
)

legend(
  "topleft",
  legend = levels(as.factor(batch)),
  text.col = pal,
  bg = "white",
  cex = 0.7
)


# ------------------------------------------------------------------------------
# 5.4 Examine higher MDS dimensions
# ------------------------------------------------------------------------------

par(mfrow = c(1, 3))

plotMDS(
  getM(mSetSq),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Group)],
  dim = c(1, 3)
)

legend(
  "top",
  legend = levels(factor(targets$Sample_Group)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


plotMDS(
  getM(mSetSq),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Group)],
  dim = c(2, 3)
)

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


plotMDS(
  getM(mSetSq),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Group)],
  dim = c(3, 4)
)

legend(
  "topright",
  legend = levels(factor(targets$Sample_Group)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


# ==============================================================================
# 6. Probe-level quality control
# ==============================================================================


# ------------------------------------------------------------------------------
# 6.1 Match probe order between methylation and detection P-value objects
# ------------------------------------------------------------------------------

detP <- detP[
  match(featureNames(mSetSq), rownames(detP)),
]


# ------------------------------------------------------------------------------
# 6.2 Remove probes failing detection P-value criteria
# ------------------------------------------------------------------------------

# Remove probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)

table(keep)

mSetSqFlt <- mSetSq[keep, ]
mSetSqFlt


# ------------------------------------------------------------------------------
# 6.3 Detection P-value diagnostic plot
# ------------------------------------------------------------------------------

barplot(
  rowMeans(detP[1:25000, ]),
  col = pal[factor(targets$Sample_Group)],
  las = 2,
  cex.names = 0.01,
  ylim = c(0, 1),
  ylab = "Mean detection p-values",
  xlab = "CpGs"
)

abline(h = 0.05, col = "red")


# ------------------------------------------------------------------------------
# 6.4 Remove probes located on sex chromosomes
# ------------------------------------------------------------------------------

# If data include males and females, remove probes on sex chromosomes
keep <- !(
  featureNames(mSetSqFlt) %in%
    ann450k$Name[ann450k$chr %in% c("chrX", "chrY")]
)

table(keep)

mSetSqFlt <- mSetSqFlt[keep, ]


# ------------------------------------------------------------------------------
# 6.5 Remove probes containing SNPs at the CpG site
# ------------------------------------------------------------------------------

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt


# ------------------------------------------------------------------------------
# 6.6 Remove cross-reactive probes
# ------------------------------------------------------------------------------

xReactiveProbes <- xreactive_probes(array_type = "EPIC")

keep <- !(
  featureNames(mSetSqFlt) %in% xReactiveProbes
)

table(keep)

mSetSqFlt <- mSetSqFlt[keep, ]
mSetSqFlt


# ==============================================================================
# 7. Post-filtering quality assessment
# ==============================================================================


# ------------------------------------------------------------------------------
# 7.1 Group-wise and batch-wise MDS
# ------------------------------------------------------------------------------

par(mfrow = c(1, 2))

plotMDS(
  getM(mSetSqFlt),
  top = 700000,
  gene.selection = "common",
  main = "Group wise",
  pch = 16,
  col = pal[factor(targets$Sample_Group)],
  cex = 0.8
)

legend(
  "topright",
  legend = levels(factor(targets$Sample_Group)),
  text.col = pal,
  cex = 0.65,
  bg = "white"
)


plotMDS(
  getM(mSetSqFlt),
  top = 700000,
  gene.selection = "common",
  main = "Batch wise",
  pch = 16,
  col = as.factor(Datam$Batch),
  cex = 0.8
)

legend(
  "topright",
  legend = levels(as.factor(Datam$Batch)),
  text.col = pal,
  cex = 0.65,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 7.2 Examine higher MDS dimensions
# ------------------------------------------------------------------------------

par(mfrow = c(1, 3))

plotMDS(
  getM(mSetSqFlt),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Source)],
  dim = c(1, 3)
)

legend(
  "right",
  legend = levels(factor(targets$Sample_Source)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


plotMDS(
  getM(mSetSqFlt),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Source)],
  dim = c(2, 3)
)

legend(
  "topright",
  legend = levels(factor(targets$Sample_Source)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


plotMDS(
  getM(mSetSqFlt),
  top = 1000,
  gene.selection = "common",
  col = pal[factor(targets$Sample_Source)],
  dim = c(3, 4)
)

legend(
  "right",
  legend = levels(factor(targets$Sample_Source)),
  text.col = pal,
  cex = 0.7,
  bg = "white"
)


# ==============================================================================
# 8. Extract M-values and beta values
# ==============================================================================

# Calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[, 1:5])

# Calculate beta values
bVals <- getBeta(mSetSqFlt)
head(bVals[, 1:5])


# ==============================================================================
# 9. Visualize final beta-value and M-value distributions
# ==============================================================================

par(mfrow = c(1, 2))

densityPlot(
  bVals,
  sampGroups = targets$Sample_Group,
  main = "Beta values",
  legend = FALSE,
  xlab = "Beta values"
)

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  text.col = brewer.pal(8, "Dark2")
)


densityPlot(
  mVals,
  sampGroups = targets$Sample_Group,
  main = "M-values",
  legend = FALSE,
  xlab = "M values"
)

legend(
  "topleft",
  legend = levels(factor(targets$Sample_Group)),
  text.col = brewer.pal(8, "Dark2")
)
