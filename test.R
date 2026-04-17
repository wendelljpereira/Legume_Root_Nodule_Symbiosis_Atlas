#setwd("/Users/wendellpereira/Library/CloudStorage/Dropbox-UF/Wendell Pereira/Root_nodule_symbiosis_atlas")

require(Seurat)

## SATURN

saturn <- readRDS("app_ready_integration/saturn/seurat_build.rds")
saturn
DimPlot(saturn,
        group.by = "saturn_species", 
        reduction = "umap")

saturn_C <- readRDS("app_ready_integration/saturn/clustered_dataset.rds")
saturn_C
DimPlot(saturn_C,
        group.by = "saturn_species", 
        reduction = "umap")

scCustomize::FeaturePlot_scCustom(saturn,
                                  features = "medicago::MtrunA17Chr1g0197491", 
                                  pt.size = 1.5)

saturn <- NormalizeData(saturn)
saturn <- FindVariableFeatures(saturn,
                               selection.method = "vst",
                               nfeatures = 2000)

all.genes <- rownames(saturn)
saturn <- ScaleData(saturn, features = all.genes)

scCustomize::FeaturePlot_scCustom(saturn,
                                  features = "medicago::MtrunA17Chr1g0197491", 
                                  pt.size = 1.5)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#"species", "sample_name", "time_point", "saturn_label", "saturn_ref_label", "saturn_species_code","saturn_labels", "saturn_labels2", "saturn_ref_labels", "saturn_species",

all_vars <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "new_ident")

unique(saturn@meta.data[, all_vars[1] ] )

scCustomize::FeaturePlot_scCustom(saturn,
            features = "medicago::MtrunA17Chr1g0197491", 
            pt.size = 1.5)

saturn <- NormalizeData(saturn)
saturn <- FindVariableFeatures(saturn,
                               selection.method = "vst",
                               nfeatures = 2000)

all.genes <- rownames(saturn)
saturn <- ScaleData(saturn, features = all.genes)

scCustomize::FeaturePlot_scCustom(saturn,
                                  features = "medicago::MtrunA17Chr1g0197491", 
                                  pt.size = 1.5)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


## CAMEX

camex_C <- readRDS("camex/clustered_dataset.rds")
camex_S <- readRDS("camex/seurat_build.rds")

readRDS("camex/clustered_dataset.rds")

DimPlot(camex_C)
DimPlot(camex_S, split.by = "species")


