# load packages 
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(ggvenn)

# load samplesheet 
samplesheet <- read_tsv("Data/sample_sheet.tsv")

# set up parallelisation 
bp.params <- MulticoreParam(workers = 7)

# load a single sample 
sample.path <- "Data/CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
sce.sing <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)

## The single cell experiment (sce) object

# The counts matrix
dim(counts(sce.sing))

counts(sce.sing)[1:10, 1:10]

# feature metadata 

# droplet metadata 
colData(sce.sing)

## Properties of scRNAseq data

# How many genes have been detected in at least 1 cell?
sum(rowSums(counts(sce.sing)) > 0)

# genes_per_cell 
genesPerCell <- colSums(counts(sce.sing) > 0)
plot(density(genesPerCell), main="", xlab="Genes per cell")

# expression_v_detected 
plot(rowSums(counts(sce.sing)) / rowSums(counts(sce.sing) > 0),
     rowMeans(counts(sce.sing) > 0),
     log = "x",
     xlab="Mean UMIs per cell",
     ylab="proportion of cells expressing the gene"
)

# top_20_genes 
rel_expression <- t( t(counts(sce.sing)) / colSums(counts(sce.sing))) * 100
rownames(rel_expression) <- rowData(sce.sing)$Symbol
most_expressed <- sort(rowSums( rel_expression ), decreasing = T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)

## Quality Control - filtering poor quality cells

# Loading multiple samples 
samples <- samplesheet$Sample[c(1,5,7,9)]
list_of_files <- str_c("Data/CellRanger_Outputs/", 
                       samples, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samples

sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)

# Modify the droplet annotation 
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), samplesheet, by="Sample", sort=FALSE)
rownames(colData(sce)) <- sce$Barcode

# Remove Undetected genes 
detected_genes <- rowSums(counts(sce)) > 0

sce <- sce[detected_genes,]

# Annotate genes 
#ah <- AnnotationHub()
#ens.hs.107<- query(ah, c("Homo sapiens", "EnsDb", 107))[[1]] 

#genes <- rowData(sce)$ID
#gene_annot <- AnnotationDbi::select(ens.hs.107, 
#                                    keys = genes,
#                                    keytype = "GENEID",
#                                    columns = c("GENEID", "SEQNAME")) %>%
#    set_names(c("ID", "Chromosome"))
#rowData(sce) <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE)
  
## Add per cell QC metrics 

#is.mito <- which(rowData(sce)$Chromosome=="MT")

is.mito <- rowData(sce)[grepl('^MT-', rowData(sce)$Symbol),]$ID

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

## QC metric distribution

# Total counts 
plotColData(sce, x="SampleName", y="sum") + 
    scale_y_log10() + 
    ggtitle("Total count")

# Detected genes 
plotColData(sce, x="SampleName", y="detected") + 
    scale_y_log10() + 
    ggtitle("Detected features")


# Mitochondrial content 
plotColData(sce, x="SampleName", y="subsets_Mito_percent") + 
    ggtitle("Mito percent")

## Identification of low-quality cells with adaptive thresholds

# Libary size filtering 
low_lib_size <- isOutlier(sce$sum, log=TRUE, type="lower")

# Detected genes 
low_n_features <- isOutlier(sce$detected, log=TRUE, type="lower")

colData(sce)$low_n_features <- low_n_features
plotColData(sce, x="SampleName", y="detected", colour_by = "low_n_features") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


# Mitochondrial gene content 

high_Mito_percent <- isOutlier(sce$subsets_Mito_percent, type="higher")

colData(sce)$high_Mito_percent <- high_Mito_percent
plotColData(sce,
            x="SampleName",
            y="subsets_Mito_percent",
            colour_by = "high_Mito_percent") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


# Summary table 
tibble(low_lib_size, low_n_features, high_Mito_percent) %>%
  mutate(discard = low_lib_size | low_n_features | high_Mito_percent) %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(where(is.logical), sum))


## All three filters at once - `quickPerCellQC`

cell_qc_results <- quickPerCellQC(colData(sce), sub.fields = TRUE)

cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(where(is.logical), sum))


# Separate thresholds for each sample

batch.cell_qc_results <- quickPerCellQC(colData(sce), 
                                         sub.fields = TRUE,
                                         batch=sce$Sample)

batch.cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(where(is.logical), sum))

# compare_thresholds 
all.thresholds <- tibble(`SampleName`="All",
       `Library Size`=attr(cell_qc_results$low_lib_size, "thresholds")[1],
       `Genes detected`=attr(cell_qc_results$low_n_features, "thresholds")[1],
       `Mitochondrial UMIs`=attr(cell_qc_results$high_subsets_Mito_percent, "thresholds")[2])


tibble(`Sample`=names(attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,]),
       `Library Size`=attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,],
       `Genes detected`=attr(batch.cell_qc_results$low_n_features, "thresholds")[1,],
       `Mitochondrial UMIs`=attr(batch.cell_qc_results$high_subsets_Mito_percent, "thresholds")[2,]) %>% 
    left_join(samplesheet) %>% 
    select(SampleName, `Library Size`, `Genes detected`, `Mitochondrial UMIs`) %>% 
    bind_rows(all.thresholds) %>% 
    mutate(across(where(is.numeric), round, digits=2))


# replace_filters_in_sce

sce$low_lib_size <- batch.cell_qc_results$low_lib_size
sce$low_n_features <- batch.cell_qc_results$low_n_features
sce$high_Mito_percent <- batch.cell_qc_results$high_subsets_Mito_percent
sce$discard <- batch.cell_qc_results$discard


# plot_library_size_batch_filters 
plotColData(sce, x="SampleName", y="sum", colour_by = "low_lib_size") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))


# plot_detected_genes_batch_filters 
plotColData(sce, x="SampleName", y="detected", colour_by = "low_n_features") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


# plot_MT_content_batch_filters 
plotColData(sce, 
        x="Sample", 
        y="subsets_Mito_percent",
        colour_by = "high_Mito_percent") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))

# Finally, remove the poor quality cells
sce.filtered <- sce[, !sce$discard]

## MT content v library_size 

plotColData(sce, 
            x="sum", 
            y="subsets_Mito_percent", 
            other_fields="SampleName",
            colour_by="discard") +
    facet_wrap(~SampleName, ncol=5, scale="free_x")



