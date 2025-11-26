library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sceasy)
library(scCustomize)

# Read full Cd3+ object from Caushi et. al. 2021
ser <- readRDS("./mydata/tcells/Robjs/ser.integrate.rds")
DimPlot(ser, group.by = "seurat_clusters", label = T)
MAIT = subset(ser, seurat_clusters == "6")
dim(MAIT)
saveRDS(MAIT, file = "./mydata/tcells/Robjs/MAIT.rds")

# Read the MAIT object
MAIT = readRDS("./mydata/tcells/Robjs/MAIT.rds")
unique(MAIT$CellType)
DefaultAssay(MAIT) <- 'integrated'

# Cluster MAIT cells
MAIT <- FindClusters(MAIT, resolution = 0.3, cluster.name = "new.clusters")
DimPlot(MAIT, group.by = "new.clusters")

DefaultAssay(MAIT) <- "RNA"
Idents(MAIT) <- MAIT$new.clusters
MAIT.markers <- FindAllMarkers(MAIT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
MAIT.markers %>% 
 group_by(cluster) %>% 
top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(MAIT.markers, "./mydata/tcells/csvs/mait.markers.res0.3.csv")

MAIT <- ScaleData(MAIT)
DimPlot(MAIT, reduction = "umap", group.by = c("patient_id"))
DimPlot(MAIT, reduction = "umap", group.by = c("new.clusters"))

Idents(MAIT) <- MAIT$new.clusters
DimPlot(MAIT)
MAIT <- RenameIdents(MAIT,
                     "0" = "Cytotoxic MAIT 1",
                     "1" =  "Naïve MAIT 1", 
                     "2" = "MAIT 17", 
                     "3" =  "Activated MAIT 1", 
                     "4" =  "TRM", 
                     "5" = "Activated tissue resident")
MAIT$new.cluster.names <- Idents(MAIT)

# Figure 1A
DimPlot(MAIT, label = T)

# Figure 1B
DotPlot(MAIT, features = unique(top10$gene), assay = "RNA", cols = c("lightgray", "brown"), 
        scale.by = "size", dot.scale = 5.0, dot.min = 0.01) +
  coord_flip() +
  # scale_colour_viridis_c(direction = -1)+
  theme_bw() + scale_color_gradient(low = "whitesmoke", high = "brown", trans = "exp") + 
  theme(plot.background=element_blank(),
        panel.grid = element_line(size = 0.1),
        legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 8, family = "Helvetica"), 
        legend.text = element_text(colour = "black", size = 6, family = "Helvetica"),
        legend.spacing = unit(0, "pt"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(colour = "black", size = 8, family = "Helvetica", angle = 0), # element_blank(), # 
        axis.text.x = element_text(colour = "black", size = 8, family = "Helvetica", angle = 45, vjust = 1, hjust = 1),
        plot.title=element_blank())

# Fixed the NA values for NR patient
table(MAIT$response, useNA = "ifany")
MAIT$response[is.na(MAIT$response)] <- "NR"
table(MAIT$response, useNA = "ifany")

# Save the MAIT object 
saveRDS(MAIT, file = "./mydata/tcells/Robjs/MAIT_newclusters.rds")
mait <- readRDS("./mydata/tcells/Robjs/MAIT_newclusters.rds")

# Subset the cytotoxic cluster from MAIT object
cytotoxic <- subset(mait, subset = new.clusters == "0")
DimPlot(cytotoxic, group.by = "new.clusters")
# cytotoxic <- FindClusters(cytotoxic, resolution = 0.3, cluster.name = "subclusters")
# cytotoxic <- RunUMAP(cytotoxic, dims = 1:10)
DimPlot(cytotoxic, group.by = "subclusters")
DimPlot(cytotoxic, group.by = "response")
DimPlot(cytotoxic, group.by = "tissue")

DefaultAssay(cytotoxic) <- "RNA"
Idents(cytotoxic) <- cytotoxic$response

# Subset the tumor cells from cytotoxic cluster
cytotoxic_tumor <- subset(cytotoxic, subset = tissue == "tumor")
DimPlot(cytotoxic_tumor, group.by = "tissue")
DefaultAssay(cytotoxic_tumor) <- "RNA"
Idents(cytotoxic_tumor) <- cytotoxic_tumor$response
response.markers <- FindMarkers(cytotoxic_tumor, ident.1 = "R", ident.2 = "NR")

# Supp Figure 1A
EnhancedVolcano(response.markers,
                lab = rownames(response.markers), 
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, labSize = 2.5,
                pointSize = 0.3, FCcutoff = 0.5, colAlpha = 0.8) + #xlim(-5, 5) +
  theme_classic(base_size = 8) + 
  theme(text = element_text(size = 8, colour = "black"), 
        plot.margin = margin(0,0,0,0),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7, color = "black"), 
        axis.line = element_line(size = 0.4),
        legend.position = "none", legend.key.size = unit(0.2,"cm"),
        legend.justification = "center",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-5,-10,-10))


# Ontology analysis for response markers from tumor cells in cytotoxic cluster
marker_data = response.markers[response.markers$avg_log2FC > 0.0 & response.markers$p_val_adj < 0.01,]
organism_dbs <- c("org.Mm.eg.db", "org.Hs.eg.db")
orgamisms = c("Mus musculus", "Homo sapiens")

gene.df <- bitr(rownames(marker_data), fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = organism_dbs[2])

egoRes <- enrichGO(gene        =  gene.df$ENTREZID,
                   OrgDb         = organism_dbs[2],
                   keyType       = 'ENTREZID',
                   ont           = c("BP"),
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05, 
                   readable = T )

goRes <- clusterProfiler::simplify(egoRes)
summaryEgoRes <- as.data.frame(egoRes)

dotplot(egoRes, showCategory = 20)
egoRes@result$minuslogpval  = -1 * log10(egoRes@result$p.adjust)
egoRes@result$minuslogpval = as.numeric(egoRes@result$minuslogpval)
barplot(egoRes, showCategory = 15, x = "minuslogpval", label_format = 100) & theme_classic(base_size = 6) & theme(plot.margin = margin(0,0,0,0), axis.text = element_text(color= "black", size = 6), plot.title = element_text(size = 6), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt")) & scale_fill_gradient(low = "#5b5b5b", high = "#5b5b5b") & xlab("-log adjusted p-value") & ylab("Term") & ggtitle("Gene ontology terms enriched for genes upregulated on DMB")

VlnPlot(mait, features = str_split(summaryEgoRes[1,"geneID"], "/")[[1]], 
        group.by = "response", pt.size = 0.0)

marker_data[str_split(summaryEgoRes[1,"geneID"], "/")[[1]],]

# R vs NR analysis on all MAIT cells
table(mait$response)
Idents(mait) <- mait$response
allmait.response.markers <- FindMarkers(mait, ident.1 = "R", ident.2 = "NR", 
                                logfc.threshold = 0.0, 
                                min.pct = 0.25)
library(EnhancedVolcano)
EnhancedVolcano(allmait.response.markers,
                lab = rownames(allmait.response.markers), 
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, labSize = 2.5,
                pointSize = 0.3, FCcutoff = 0.5, colAlpha = 0.8) + xlim(-2, 2) +
  theme_classic(base_size = 8) + 
  theme(text = element_text(size = 8, colour = "black"), 
        plot.margin = margin(0,0,0,0),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7, color = "black"), 
        axis.line = element_line(size = 0.4),
        legend.position = "none", legend.key.size = unit(0.2,"cm"),
        legend.justification = "center",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-5,-10,-10))

# Read file for TCR sequencing data from Caushi et. al. 2021 

vdj = read.csv('./mydata/tcells/csvs/MAIT_vdj_information.txt') 
vdj$barcode.barcode1 = vdj$barcode.barcode1 %>%  
  strsplit(split = '\t') %>%  
  sapply(function(x) x[1]) 
colnames(vdj)[1] = 'barcode' 

table(vdj$chain)
colnames(vdj)
vdj_a = vdj %>% filter(chain == "TRA")
vdj_b = vdj %>% filter(chain == "TRB")

allconventional = vdj %>% filter(v_gene=='TRAV1-2' & (j_gene=='TRAJ33' | j_gene=='TRAJ12' | j_gene=='TRAJ20'))
sum(colnames(mait) %in% allconventional$barcode)
mait = AddMetaData(mait, colnames(mait) %in% allconventional$barcode, col.name = 'allconventional')
table(mait$allconventional)

J33 = vdj %>% filter(v_gene=='TRAV1-2' & j_gene=='TRAJ33')
mait = AddMetaData(mait, colnames(mait) %in% J33$barcode, col.name = 'J33') 

J20 = vdj %>% filter(v_gene=='TRAV1-2' & j_gene=='TRAJ20')
mait = AddMetaData(mait, colnames(mait) %in% J20$barcode, col.name = 'J20') 

J12 = vdj %>% filter(v_gene=='TRAV1-2' & j_gene=='TRAJ12')
mait = AddMetaData(mait, colnames(mait) %in% J12$barcode, col.name = 'J12') 

# Read file for TCR sequencing from Caushi et. al. 2021 

trab.wide = readRDS("./mydata/tcells/Robjs/trab.wide.public.rds") 

sum(trab.wide$barcode %in% colnames(mait)) 
meta = mait@meta.data 
rownames(trab.wide) <- trab.wide$barcode
meta = meta %>% left_join(trab.wide, by = "barcode") 
rownames(meta) = meta$barcode 
colnames(meta)
mait = AddMetaData(mait, metadata = meta) 

top10_clono = table(mait$clonotype) %>% sort(decreasing = T) %>% head(10) 
top10_clono = sapply(mait$clonotype, function(x) { 
  if (x %in% names(top10_clono)) { 
    return(x) 
  } else { 
    return(NA) 
  } 
}) 
mait = AddMetaData(mait, top10_clono, col.name = 'top10_clonotypes') 
DimPlot(mait, group.by = 'top10_clonotypes', order = T) 
mait$j_gene <- ifelse(mait$J12 == TRUE, "J12",
                      ifelse(mait$J20 == TRUE, "J20",
                             ifelse(mait$J33 == TRUE, "J33", "Non-conventional")))

# Save python object for more plots
as.anndata(x = mait, file_path = "./mydata/tcells/h5ad/", file_name = "MAIT_final.h5ad")


