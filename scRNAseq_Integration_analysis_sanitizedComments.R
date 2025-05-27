# original analysis from Dec 2024 - not used for the current version of the manuscript. See x_0092-analysis-paper-figs_2025-05-27.Rmd for current version of analysis. 

# setup

# #install packages
# install.packages("tidyverse")
# install.packages('Seurat')
# install.packages("here")
# install.packages("patchwork")
# #Seurat does not require, but makes use of, packages developed by other labs that can substantially enhance speed and performance. These include presto (Korunsky/Raychaudhari labs), BPCells (Greenleaf Lab), and glmGamPoi (Huber Lab). We recommend users install these along with Seurat:
# 
# setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
# install.packages(c("BPCells", "presto", "glmGamPoi"))
#load packages
library(tidyverse)
library(Seurat)
library(here)
library(patchwork)
library(BPCells)
library(presto)
library(glmGamPoi)


# # Install the remotes package
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# install.packages('Signac')
# remotes::install_github("satijalab/seurat-data", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

# Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("scater", "scran", "DropletUtils"))

# Load libraries
library(scater)
library(scran)
library(DropletUtils)

# install.packages("RColorBrewer")
# install.packages("ggnewscale")

library(RColorBrewer)
library(ggnewscale)

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################


# import data
# this uses the list of samples given to cellranger to specify the location of input files, 
## and then loops through those files, reading them in and assigning the group metadata variable to each. 
## it also adds the percent mitochondrial reads to the metadata. 

# Import sample names from sample sheet
samp_names <- read_lines(here("..", "cellranger_out", "sample.list"))

# Loop through sample names to create Seurat objects and assign groups
samp_list <- map(samp_names, ~ {
  data_dir <- here("..", "cellranger_out", .x, "outs", "raw_feature_bc_matrix")
  
  # Verify required files
  stopifnot(all(c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz") %in% list.files(data_dir)))
  
  # Read data and create Seurat object
  dataobj <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(
    counts = dataobj,
    project = .x,
    min.cells = 3,
    min.features = 200
  )
  
  # Assign group
  seurat_obj$group <- .x
  
  # Add percentage of mitochondrial genes (genes starting with "DR311")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^DR311")
  
  seurat_obj
})
names(samp_list) <- samp_names


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


# this looks messy but is necessary. 
# merging with a loop in a similar fashion to reading in cellranger counts above - while vastly cleaner-looking here - can cause infuriating and nearly uninterpretable errors in the integration steps. 
# This is because of how Seurat handles successive merge operations. The metadata get garbled. It was worth avoiding, even though it means something that looks as silly as the following code. 

# Pull individual seurat objects out of the list and assign to individual objects
LV_5FBS <- samp_list[["Lv_5FBS"]]
LV_10FBS <- samp_list[["Lv_10FBS"]]
LV_15FBS <- samp_list[["Lv_15FBS"]]

Lv_2hpf <- samp_list[["Lv-2hpf"]]
Lv_3hpf <- samp_list[["Lv-3hpf"]]
Lv_4hpf <- samp_list[["Lv-4hpf"]]
Lv_5hpf <- samp_list[["Lv-5hpf"]]
Lv_6hpf <- samp_list[["Lv-6hpf"]]
Lv_7hpf <- samp_list[["Lv-7hpf"]]
Lv_8hpf <- samp_list[["Lv-8hpf"]]
Lv_9hpf <- samp_list[["Lv-9hpf"]]
Lv_10hpf <- samp_list[["Lv-10hpf"]]
Lv_11hpf <- samp_list[["Lv-11hpf"]]
Lv_12hpf <- samp_list[["Lv-12hpf"]]
Lv_13hpf <- samp_list[["Lv-13hpf"]]
Lv_14hpf <- samp_list[["Lv-14hpf"]]
Lv_15hpf <- samp_list[["Lv-15hpf"]]
Lv_16hpf <- samp_list[["Lv-16hpf"]]
Lv_18hpf <- samp_list[["Lv-18hpf"]]
Lv_20hpf <- samp_list[["Lv-20hpf"]]
Lv_24hpf <- samp_list[["Lv-24hpf"]]

# Assign experiment identifiers
LV_5FBS[["samplesource"]] <- "cultures"
LV_10FBS[["samplesource"]] <- "cultures"
LV_15FBS[["samplesource"]] <- "cultures"

Lv_2hpf[["samplesource"]] <- "atlas"
Lv_3hpf[["samplesource"]] <- "atlas"
Lv_4hpf[["samplesource"]] <- "atlas"
Lv_5hpf[["samplesource"]] <- "atlas"
Lv_6hpf[["samplesource"]] <- "atlas"
Lv_7hpf[["samplesource"]] <- "atlas"
Lv_8hpf[["samplesource"]] <- "atlas"
Lv_9hpf[["samplesource"]] <- "atlas"
Lv_10hpf[["samplesource"]] <- "atlas"
Lv_11hpf[["samplesource"]] <- "atlas"
Lv_12hpf[["samplesource"]] <- "atlas"
Lv_13hpf[["samplesource"]] <- "atlas"
Lv_14hpf[["samplesource"]] <- "atlas"
Lv_15hpf[["samplesource"]] <- "atlas"
Lv_16hpf[["samplesource"]] <- "atlas"
Lv_18hpf[["samplesource"]] <- "atlas"
Lv_20hpf[["samplesource"]] <- "atlas"
Lv_24hpf[["samplesource"]] <- "atlas"

# Merge Seurat objects for cultures and atlas. Explicitly specify all samples instead of looping through them. 
LVcultures <- merge(x = LV_5FBS, y = c(LV_10FBS, LV_15FBS))
LV_2_to_24hr <- merge(x = Lv_2hpf, y = c(
  Lv_3hpf, Lv_4hpf, Lv_5hpf, Lv_6hpf, Lv_7hpf, Lv_8hpf, 
  Lv_9hpf, Lv_10hpf, Lv_11hpf, Lv_12hpf, Lv_13hpf, Lv_14hpf,
  Lv_15hpf, Lv_16hpf, Lv_18hpf, Lv_20hpf, Lv_24hpf
))

#saving the raw seurat objects in case I need to backtrack. 
save(LV_2_to_24hr,  file = here("RawSeuratObjects","Lvar_Atlas_2_to_24_raw.Rda"))
save(LVcultures,  file = here("RawSeuratObjects", "Lvar_cellcultures_raw.Rda"))

#filtering and joining (joining now apepars necessary to avoid weird behavior in integration steps)

LV_2_to_24hr_filtered <- subset(
  LV_2_to_24hr,
  subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
    nCount_RNA > 600 & nCount_RNA < 10000 &
    percent.mt < 30
)
LV_2_to_24hr_filtered<- JoinLayers(LV_2_to_24hr_filtered)

LVcultures_filtered <- subset(
  LVcultures,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
    nCount_RNA > 400 & nCount_RNA < 10000 &
    percent.mt < 30
)
LVcultures_filtered <- JoinLayers(LVcultures_filtered)


##################################################################################
##################################################################################

#SCTransform here replaces eplaces NormalizeData(), ScaleData(), and FindVariableFeatures() from the standard workflow
#SCTransform does not assume equal amounts of RNA per cell in normalization. 
#This assumption is probably not valid for the integrated analysis, or we'd use the standard workflow. We go this route to be conservative.

LV_2_to_24hr_sct <- SCTransform(LV_2_to_24hr_filtered, method = "glmGamPoi", verbose = FALSE)
LV_2_to_24hr_sct <- LV_2_to_24hr_sct %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:30) #dimension choice based on analysis done for sc atlas and on elbow plot of these data (not included in this minimal reproduction of the analysis)

LVcultures_sct <- SCTransform(LVcultures_filtered, method = "glmGamPoi", verbose = FALSE)
LVcultures_sct <- LVcultures_sct %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:30) #as above

#probably dont need to make this explicit in retrospect
DefaultAssay(LV_2_to_24hr_sct) <- "SCT"
DefaultAssay(LVcultures_sct) <- "SCT"

# Step 2: Find integration features
features <- SelectIntegrationFeatures(object.list = list(LV_2_to_24hr_sct, LVcultures_sct), nfeatures = 3000)

# Step 3: Prepare SCT Integration
integrationlist <- PrepSCTIntegration(object.list = list(LV_2_to_24hr_sct, LVcultures_sct), anchor.features = features)
save(integrationlist, file = here("integrationlist.Rda"))

#integration
#using RPCA here which is supposedly more efficient in finding anchors vs CCA (default, which hangs when I've tried it - or else just takes way too long)
LVanchors <- FindIntegrationAnchors(object.list = integrationlist, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
LVint <- IntegrateData(anchorset = LVanchors, normalization.method = "SCT", dims = 1:30)
LVint <- FindNeighbors(LVint, dims = 1:50, verbose = FALSE) %>% FindClusters(resolution = 2, verbose = FALSE) #dims determined by elbow plot. resolution makes a similar number of clusters as we observe cell types in LV scRNA atlas, which seems appropriate for our understanding of the biology. 
LVint <- RunPCA(LVint, verbose = FALSE)
LVint <- RunUMAP(LVint, reduction = "pca", dims = 1:30)
save(LVint, file = "LVint.Rda") #This file is available from Wray lab on request. 

#visualizations
p1 <- DimPlot(LVint, reduction = "umap", group.by = "samplesource")
p1 #Fig 4A
ggsave(filename = "UMAP-source.png", path = "IntegratedUMAPs")
p2 <- DimPlot(LVint, reduction = "umap", group.by = "orig.ident", label = TRUE,
              repel = TRUE)
p2 #unused for manuscript
ggsave(filename = "umap-orig-ident.png", path = "IntegratedUMAPs")
p1 + p2 #unused for manuscript

#plot samples

# Ensure orig.ident is a factor meaningful order
LVint$orig.ident <- factor(LVint$orig.ident, levels = c(
  "Lv-2hpf", "Lv-3hpf", "Lv-4hpf", "Lv-5hpf", "Lv-6hpf", 
  "Lv-7hpf", "Lv-8hpf", "Lv-9hpf", "Lv-10hpf", "Lv-11hpf",
  "Lv-12hpf", "Lv-13hpf", "Lv-14hpf", "Lv-15hpf", "Lv-16hpf",
  "Lv-18hpf", "Lv-20hpf", "Lv-24hpf", "Lv_5FBS", "Lv_10FBS", "Lv_15FBS"
))

# Plot UMAP showing `orig.ident` with correct order
DimPlot(LVint, reduction = "umap", group.by = "orig.ident", label = TRUE, ) #unused for manuscript
ggsave(filename = "UMAP-sample.png", path = "IntegratedUMAPs")


#plot clusters
p3 <- DimPlot(LVint, reduction = "umap", group.by = "seurat_clusters") #unused for manuscript
p3
ggsave(filename = "umap-seurat-clusters.png", path = "IntegratedUMAPs")

#split by sample source and plot clusters - easier to see that there is indeed overlap given the big change in cell-type abundance from embryos to cultures. 
#fig 4D
p4 <- DimPlot(LVint, reduction = "umap", group.by = "seurat_clusters", split.by = "samplesource", label = TRUE, label.size = 3) +
  theme(
    legend.position = "bottom", 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 10), 
    legend.key.size = unit(3, "lines"), 
    legend.box = "horizontal", 
    legend.spacing.x = unit(0.1, "cm") 
  ) + #legend was kind of a pain
  guides(
    color = guide_legend(
      nrow = 3, 
      byrow = TRUE 
    )
  )
p4 #fig 4D
ggsave(filename = "splitclusterUMAP.png", path = "IntegratedUMAPs")
ggsave(filename = "splitclust2l.png", width = 11, height = 6, units = "in")

#marker genes  - select favorites
# LOC121428528 - pks1 - in cluster markers
#using this in supplement with modifications. modifying below for better visibility of sample source.
VlnPlot(LVint, features = "LOC121428528")
ggsave(filename = "vln-pks1.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121428528", label = T, )  + labs(title = "LOC121428528 - pks1")#pks1
ggsave(filename = "umap-pks1.png", path = "IntegratedUMAPs/")


####### SUPP FIG 6A (PKS1) #####

# Fetch UMAP coordinates, gene expression, and metadata for pks1
plot_data <- FetchData(LVint, vars = c("umap_1", "umap_2", "LOC121428528", "samplesource"))

# Fetch UMAP coordinates, expression, and metadata & Normalize expression values for each group
plot_data <- plot_data %>%
  group_by(samplesource) %>%
  mutate(norm_expression = scales::rescale(LOC121428528)) %>% # Normalize to [0, 1] per group for better color scaling. mention in text. 
  ungroup()

plot_data <- plot_data %>%
  mutate(LOC121428528 = pmin(LOC121428528, 20)) #impose max cutoff to keep rare high count observations from blowing out color scaling. mention in text.


# Adjust Brewer palettes to start further up the scale for better visibility
greens_adjusted <- colorRampPalette(brewer.pal(9, "Greens")[3:9])(100) # Remove very light greens
reds_adjusted <- colorRampPalette(brewer.pal(9, "Reds")[3:9])(100)    # Remove very light reds

# Plot with adjusted gradients and dual color scale
ggplot(plot_data, aes(x = umap_1, y = umap_2)) +
  # Plot atlas with adjusted red scale 
  geom_point(
    data = subset(plot_data, samplesource == "atlas"),
    aes(color = LOC121428528),
    size = 0.01
  ) +
  scale_color_gradientn(colors = reds_adjusted, name = "Expression\n(Atlas)") +
  new_scale_color() +
  # Plot cultures with adjusted green scale 
  geom_point(
    data = subset(plot_data, samplesource == "cultures"),
    aes(color = LOC121428528),
    size = 0.01
  ) +
  scale_color_gradientn(colors = greens_adjusted, name = "Expression\n(Cultures)") +
theme_minimal() +
  labs(title = "LOC121428528 - Pks1") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

ggsave(filename = "umap-pks1_v2.png", path = "IntegratedUMAPs/") ####### SUPP FIG 6A (PKS1) #####

# LOC447794 - alx-1 
VlnPlot(LVint, features = "LOC447794") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-alx1.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC447794", min.cutoff = 'q1', label = T)  + labs(title = "LOC447794 - alx1")#alx1  ##unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-alx1.png", path = "IntegratedUMAPs/")

# LOC446151 - Brachyury
VlnPlot(LVint, features = "LOC446151")  #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-brachyury.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC446151", min.cutoff = 'q1', label = T) + labs(title = "LOC446151 - brachyuri") #Brachyury  #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-brachyury.png", path = "IntegratedUMAPs/")

#Seawi LOC121431753
VlnPlot(LVint, features = "LOC121431753")  #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-seawi.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121431753", label = T) + labs(title = "LOC121431753 - seawi") #seawi  #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-seawi.png", path = "IntegratedUMAPs/")

#Vasa LOC121424246
VlnPlot(LVint, features = "LOC121424246") #unused in manuscript but probably interesting for those exploring these data
FeaturePlot(object = LVint, features = "LOC121424246", label = T) + labs(title = "LOC121424246 - vasa")#vasa #unused in manuscript but probably interesting for those exploring these data

#Nanos2 LOC121408828
VlnPlot(LVint, features = "LOC121408828")#unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-nanos2.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121408828", label = T) + labs(title = "LOC121408828 - nanos") #nanos#unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-nanos2.png", path = "IntegratedUMAPs/")

#tcf - LOC121432141
VlnPlot(LVint, features = "LOC121431848")#unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-tcf.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121431848", label = T) + labs(title = "LOC121431848 - tcf")#tcf#unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-tcf.png", path = "IntegratedUMAPs/")

#endo16 LOC121406173 --- alter this to make it easier to see sample source and include as supp fig
VlnPlot(LVint, features = "LOC121406173")
ggsave(filename = "vln-endo16.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121406173", label = T, max.cutoff = 5) + labs(title = "LOC121406173 - endo16")
ggsave(filename = "umap-endo16.png", path = "IntegratedUMAPs/")

#alt endo16 plot - two separate color scales to show sample source (atlas vs cultures)
####### SUPP FIG 6B (ENDO 16) #####

# Fetch UMAP coordinates, gene expression, and metadata for endo16
plot_data <- FetchData(LVint, vars = c("umap_1", "umap_2", "LOC121406173", "samplesource"))

# Normalize expression values for each group
plot_data <- plot_data %>%
  group_by(samplesource) %>%
  mutate(norm_expression = scales::rescale(LOC121406173)) %>% # Normalize to [0, 1] per group (creates consistent color scaling). mention in text
  ungroup()

plot_data <- plot_data %>%
  mutate(LOC121406173 = pmin(LOC121406173, 10)) #impose max count cutoff per cell of 10 (above 10 just registers as the maximum). this prevents a few very strongly expressing cells from blowing out the scaling. mention in text


# Adjust Brewer palettes to start further up the scale for better visibility
greens_adjusted <- colorRampPalette(brewer.pal(9, "Greens")[3:9])(100) # Remove very light greens
reds_adjusted <- colorRampPalette(brewer.pal(9, "Reds")[3:9])(100)    # Remove very light reds

# Plot with adjusted dual gradients
ggplot(plot_data, aes(x = umap_1, y = umap_2)) +
  # Plot cultures with adjusted green scale 
  geom_point(
    data = subset(plot_data, samplesource == "cultures"),
    aes(color = LOC121406173),
    size = 0.01
  ) +
  scale_color_gradientn(colors = greens_adjusted, name = "Expression\n(Cultures)") +
  new_scale_color() +
  # Plot atlas with adjusted red scale
  geom_point(
    data = subset(plot_data, samplesource == "atlas"),
    aes(color = LOC121406173),
    size = 0.01
  ) +
  scale_color_gradientn(colors = reds_adjusted, name = "Expression\n(Atlas)") +
  theme_minimal() +
  labs(title = "LOC121406173 - Endo16") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

ggsave(filename = "umap-endo16_v2.png", path = "IntegratedUMAPs/") ####### SUPP FIG 6B (ENDO 16) #####

#foxa1 LOC121419318
VlnPlot(LVint, features = "LOC121419318") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-foxa1.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121419318", label = T) + labs(title = "LOC121419318 - foxa1") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-foxa1.png", path = "IntegratedUMAPs/")


# tektin-1-like LOC121408954
VlnPlot(LVint, features = "LOC121408954") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-tektin-1-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121408954", min.cutoff = 'q1', label = T) + labs(title = "LOC121408954 - tektin-1-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-tektin-1-like.png", path = "IntegratedUMAPs/")

# cilia-_and_flagella-associated_protein_251-like LOC121411688
VlnPlot(LVint, features = "LOC121411688") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-cilia-_and_flagella-associated_protein_251-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121411688",  label = T) + labs(title = "LOC121411688 - cilia-_and_flagella-associated_protein_251-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-cilia-_and_flagella-associated_protein_251-like.png", path = "IntegratedUMAPs/")

#LOC121405827	- G1/S-specific_cyclin-D1-like
VlnPlot(LVint, features = "LOC121405827") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-G1.S-specific_cyclin-D1-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121405827", min.cutoff = 'q1', label = T) + labs(title = "LOC121405827	- G1/S-specific_cyclin-D1-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-G1.S-specific_cyclin-D1-like.png", path = "IntegratedUMAPs/")

#LOC121417436 - proliferation-associated_protein_2G4-like
VlnPlot(LVint, features = "LOC121417436") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-proliferation-associated_protein_2G4-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121417436", min.cutoff = 'q1', label = T) + labs(title = "LOC121417436 - proliferation-associated_protein_2G4-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-proliferation-associated_protein_2G4-like.png", path = "IntegratedUMAPs/")


#LOC121428657	-	hepatocyte_nuclear_factor_4-alpha-like (FOXA1)
VlnPlot(LVint, features = "LOC121428657") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-hnf4-foxa1.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121428657", min.cutoff = 'q1', label = T) + labs(title = "LOC121428657	-	hepatocyte_nuclear_factor_4-alpha-like (FOXA1)") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-hnf4-foxa1.png", path = "IntegratedUMAPs/")

#LOC121419318	-	hepatocyte_nuclear_factor_3-beta-like (FOXA2)
VlnPlot(LVint, features = "LOC121419318") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-hnf3bl-foxa2.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121419318", min.cutoff = 'q1', label = T) + labs(title = "LOC121419318	-	hepatocyte_nuclear_factor_3-beta-like (FOXA2)") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-hnf3bl-foxa2.png", path = "IntegratedUMAPs/")

# LOC121430133 -	transcription_factor_Sox-2-like
VlnPlot(LVint, features = "LOC121430133") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-sox2-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121430133", min.cutoff = 'q1', label = T) + labs(title = "LOC121430133 -	transcription_factor_Sox-2-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "umap-sox2-like.png", path = "IntegratedUMAPs/")

# LOC121418887	-	homeobox_protein_six1-like
VlnPlot(LVint, features = "LOC121418887") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-homeobox_protein_six1-like.png", path = "IntegratedUMAPs/")
FeaturePlot(object = LVint, features = "LOC121418887", min.cutoff = 'q1', label = T) + labs(title = "LOC121418887 	-	homeobox_protein_six1-like") #unused in manuscript but probably interesting for those exploring these data
ggsave(filename = "vln-homeobox_protein_six1-like.png", path = "IntegratedUMAPs/")

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


######## Manuscript figure 4C #############
#which clusters share cells??
metatable_orig.ident <- table(LVint@meta.data$integrated_snn_res.2, LVint@meta.data$orig.ident, LVint@meta.data$samplesource) %>% as_tibble(.name_repair = c("minimal"))
colnames(metatable_orig.ident) <- c("cluster", "sample", "sample_source", "counts")
metatable <- table(LVint@meta.data$integrated_snn_res.2, LVint@meta.data$samplesource) %>% as_tibble(.name_repair = c("minimal"))
colnames(metatable) <- c("cluster", "sample_source", "counts")

# Filter clusters with counts >= 5 in both conditions
filtered_metatable <- metatable %>%
  group_by(cluster) %>%
  filter(all(counts[sample_source == "atlas"] >= 5, counts[sample_source == "cultures"] >= 5)) #defining shared contributions to clusters as contributions of >= 5 cells. This might be a touch too conservative. 

library(forcats) # For fct_reorder

# Order clusters by counts in cultures - asthetic decision
ordered_clusters <- filtered_metatable %>%
  filter(sample_source == "cultures") %>%
  arrange(desc(counts)) %>%
  pull(cluster) %>%
  unique()

# Apply the ordering to the cluster factor
filtered_metatable <- filtered_metatable %>%
  mutate(cluster = factor(cluster, levels = ordered_clusters))

# make plot
ggplot(filtered_metatable, aes(x = counts, y = cluster, color = sample_source)) + #plot counts by cluster, color by culture vs atlas
  geom_line(aes(group = cluster), color = "black", size = 0.5) + #lines connecting the dots, behind dots
  geom_point(size = 3) + #dots
  scale_x_log10() + #log scaling for counts
  scale_color_manual(
    values = c("atlas" = "#F37971", "cultures" = "#0FBCC2") #match the colors in FIG 4A 
  ) +
  theme_minimal()  + #minor theme adjustments. Labeling is mostly done in Adobe illustrator 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    legendL
  ) 
ggsave("cellclustdiff.png", width = 3, height = 4, units = "in") #save - plot will get cropped but this makes correct size for compound figure. 


###### VENN DIAGRAM FIG 4B ######
#this is not scaled. venn diagrams are awful but here I think it nicely gets the conceptual point across. 
# Install and load ggvenn 
if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")
library(ggvenn)
# Filter clusters with counts < 5 as having no counts
filtered_data <- metatable %>%
  mutate(counts = ifelse(counts < 5, 0, counts)) %>%
  group_by(cluster) %>%
  summarize(
    atlas = any(sample_source == "atlas" & counts > 0),
    cultures = any(sample_source == "cultures" & counts > 0)
  )
# Create a list of sets for Venn diagram
venn_sets <- list(
  Atlas = filtered_data$cluster[filtered_data$atlas],
  Cultures = filtered_data$cluster[filtered_data$cultures]
)
# Plot the Venn diagram
ggvenn::ggvenn(
  venn_sets,
  fill_color = c("#F37971", "#0FBCC2"), #match colors to plot 4A
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 5
) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave("clustervenn.png")

# alternative plot for Fig 4C (not used)
# alt line-dot plot order (increasing absolute difference top to bottom)

# Prepare the data
filtered_metatable <- metatable %>%
  mutate(counts = ifelse(counts < 5, 0, counts)) %>% # Treat counts < 5 as 0
  pivot_wider(names_from = sample_source, values_from = counts, values_fill = 0) %>% # Wide format
  filter(atlas >= 5 & cultures >= 5) %>% # Include only clusters where both atlas and cultures >= 5
  mutate(abs_diff = abs(atlas - cultures)) %>% # Calculate absolute difference
  arrange(abs_diff) %>% # Order by increasing absolute difference
  mutate(cluster = fct_reorder(cluster, -abs_diff)) # Reorder clusters so smallest difference is at the top

# Reshape data to long format for plotting
long_data <- filtered_metatable %>%
  pivot_longer(cols = c(atlas, cultures), names_to = "sample_source", values_to = "counts")

# Create the plot
ggplot(long_data, aes(x = counts, y = cluster, color = sample_source)) +
  geom_line(aes(group = cluster), color = "black", size = 0.5) + # Plot lines first with line width in points
  geom_point(size = 5) + # Plot points with point size in points
  scale_x_log10() + # Logarithmic X scale
  scale_color_manual(
    values = c("atlas" = "#F37971", "cultures" = "#0FBCC2"),
    labels = c("Atlas", "Cultures")
  ) +
  theme_minimal() +
  labs(
    title = "Cluster Counts (Both Atlas and Cultures >= 5)",
    x = "Log10(Counts)",
    y = "Cluster",
    color = "Sample Source"
  ) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "bottom"
  )
ggsave("cellclustdiff_altorder.png", width = 3, height = 4, units = "in")

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################



sessioninfo::session_info(pkgs = "loaded", include_base = T, info = "all", to_file = "session_info.txt")  
  
