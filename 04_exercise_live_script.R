## Quality Control - Practical

# In the course materials we performed QC and filtering of 2 samples from each
# of the sample groups. For this challenge we would like you to perform QC and 
# filtering on all of the samples from the Caron data set. 

# In this script we have provided you with much of the code for the analysis.
# There are a few exercises where you will need to write your own code.
# Please refer to the accompanying page on the website for further insights and
# suggestions.

## Initial set up of resources

# Load packages 
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)

# Load the samplesheet 
samplesheet <- read_tsv("Data/sample_sheet.tsv")

# Set up parallelisation 
bp.params <- MulticoreParam(workers = 7)


## -- Exercise 1 -- ############################################################
# Load the data from Cell Ranger 

# In order to load the CellRanger data for all of the samples, you will first
# need to create a named vector of the paths to the filtered count matrix
# folders called `list_of_files` and then use this in the `read10xCounts`
# command.

samples <- samplesheet$Sample
list_of_files <- str_c("Data/CellRanger_Outputs/", 
                       samples, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samples
list_of_files

sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)
sce

################################################################################

# Check Single Cell Experiment object - does it contain the samples you are
# expecting?

colData(sce) %>%
    as.data.frame() %>% 
    select(Sample) %>% 
    distinct()
sce

# Add sample metadata to the droplet annotation

sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), samplesheet, by="Sample", sort=FALSE)
rownames(colData(sce)) <- sce$Barcode

# Check that this is what we want
colData(sce)

# Filter out undetected 
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes,]

## -- Exercise 2 -- ############################################################

# What proportion of the genes have been detected in at least 1 cell?

detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)


################################################################################

## Add the chromosome information to the feature (gene) annotation

## Add per cell QC metrics to the `colData` table

# we will use `addPerCellQC` to do this, as well as the metrics for all genes
# we want to add metrics that pertain just to the mitochondrial genes.

is.mito <- rowData(sce)[grepl('^MT-', rowData(sce)$Symbol),]$ID

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

## -- Exercise 3 -- ############################################################

# Use the `scater` function `plotColData` to generate plots showing the
# distributions of the total number of UMIs, the number of genes detected and
# percentage of UMIs aligned to mitochondrial genes across all cells for each
# sample.

# This time we will use `facet_wrap` to split the plot according to the Sample
# Group. We will need to include the sample group in the plot data using the 
# `other_fields` argument in `plotColData` (see the help page for details - 
# `?plotColData`) in order that we can use it in the `facet_wrap` command. 

# The code for plotting the total number of UMIs is shown below. You will also
# need to plot the the number of genes detected and percentage of UMIs aligned
# to mitochondrial.

# Total UMI plot
plotColData(sce, x="SampleName", y="sum", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Total count")

# Detected genes 
plotColData(sce, x="SampleName", y="detected", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Detected features")

# Mitochondrial content 
plotColData(sce, x="SampleName", y="subsets_Mito_percent", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") +
   ggtitle("Mito percent")



################################################################################

## Identification of low-quality cells with adaptive thresholds

## -- Exercise 4 -- ############################################################

# Use the scater function `quickPerCellQC` to assess cell quality based
# on the three metrics. Name the object generated `cell_qc_results`.  

# When running the command, consider the distribution plots above and decide
# whether to use the `batch` option and if so, at what level it should be
# applied.

batch.cell_qc_results <- quickPerCellQC(colData(sce), 
                                        sub.fields = TRUE,
                                        batch=sce$Sample)

# How many cells will be removed from the data set?

qc_tbl <- batch.cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(where(is.logical), sum))

sum(qc_tbl$discard)

################################################################################
 
# Add the filtering results to the colData of the SCE object


colData(sce) <- cbind(colData(sce), batch.cell_qc_results)

## Visualize the results of the filtering

# Library size filtering 

plotColData(sce, 
            x="SampleName", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))

# Detected genes filtering 

plotColData(sce, 
            x="SampleName", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))

# MT content filtering

plotColData(sce, 
        x="SampleName", 
        y="subsets_Mito_percent",
        other_fields="SampleGroup", 
        colour_by = "high_subsets_Mito_percent") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


## Filter out the poor quality cells

sce <- sce[, !sce$discard]

# Save the final object

saveRDS(sce, "filtered_cells.rds")

# Load object back in for next time

sce <- readRDS("filtered_cells.rds")

