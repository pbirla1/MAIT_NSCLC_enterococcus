library(EnhancedVolcano)

markers_all = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo_deg_all.csv", row.names = 1)
text_size = 14
point_size = 2.0
EnhancedVolcano(markers_all,
                lab = markers_all$gene_name,
                x = 'log2FoldChange', 
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~- Log[10]~ 'adjusted p-value'),
                title = 'Riboflavin vs. EF + Ribo (|LFC| > 1 & adj. p-value < 0.05)',
                titleLabSize = text_size,
                subtitle = NULL,
                caption = NULL,
                legendPosition = "none",
                pCutoff = 0.05,
                drawConnectors = T,
                axisLabSize = text_size, 
                labSize = text_size/3,
                pointSize = point_size,
                FCcutoff = 1.0,
                colAlpha = 0.8) + 
  # xlim(-4.5,3.5) +
  xlim(-21, 3) +
  ylim(1, 176) +
  theme_classic() +
  theme(plot.margin = margin(0,0,0,0),
        text = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size, color = "black"), 
        axis.line = element_line(size = 0.4), 
        legend.position = "none", 
        legend.justification = "center", 
        legend.margin=margin(0,0,0,0),  
        legend.box.margin=margin(0,0,0,0))

##############################

markers_all = read.csv("./mydata/tcells/files_from_c/EfaecalisvsEcoli_deg_all.csv", row.names = 1)
text_size = 14
point_size = 2.0
EnhancedVolcano(markers_all,
                lab = markers_all$gene_name,
                x = 'log2FoldChange', 
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~- Log[10]~ 'adjusted p-value'),
                title = 'EF vs. E coli (|LFC| > 1 & adj. p-value < 0.05)',
                titleLabSize = text_size,
                subtitle = NULL,
                caption = NULL,
                legendPosition = "none",
                pCutoff = 0.05,
                drawConnectors = T,
                axisLabSize = text_size, 
                labSize = text_size/3,
                pointSize = point_size,
                FCcutoff = 1.0,
                colAlpha = 0.8) + 
  xlim(-3,2.0) +
  ylim(1, 155) +
  theme_classic() +
  theme(plot.margin = margin(0,0,0,0),
        text = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size, color = "black"), 
        axis.line = element_line(size = 0.4), 
        legend.position = "none", 
        legend.justification = "center", 
        legend.margin=margin(0,0,0,0),  
        legend.box.margin=margin(0,0,0,0))

##############################

markers_all = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo_deg_all.csv", row.names = 1)
go_resuts = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_GOenrich_significant.csv")
gene_list1 = strsplit(go_resuts[go_resuts$Description == "protein targeting to ER", "geneName"], split = "/")[[1]]
gene_list2 = strsplit(go_resuts[go_resuts$Description == "ribonucleoprotein complex biogenesis", "geneName"], split = "/")[[1]]
kegg_resuts = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_KEGGenrich_significant.csv")
gene_list3 = strsplit(kegg_resuts[kegg_resuts$Description == "Endocytosis", "geneName"], split = "/")[[1]]
gene_list4 = strsplit(kegg_resuts[kegg_resuts$Description == "Protein processing in endoplasmic reticulum", "geneName"], split = "/")[[1]]

highlight_genes <- unique(c(gene_list1, gene_list2, gene_list3, gene_list4))
markers_all <- markers_all[order(markers_all$gene_name %in% highlight_genes), , drop = FALSE]

keyvals <- ifelse(markers_all$gene_name %in% gene_list1, "blue",
                  ifelse(markers_all$gene_name %in% gene_list2, "red",
                         ifelse(markers_all$gene_name %in% gene_list3, "green",
                                ifelse(markers_all$gene_name %in% gene_list4, "yellow",
                                       "lightgrey"))))
keyvals[is.na(keyvals)] <- "lightgrey"

names(keyvals)[keyvals == 'blue'] <- 'GO: protein targeting to ER'
names(keyvals)[keyvals == 'red'] <- 'GO: ribonucleoprotein complex biogenesis'
names(keyvals)[keyvals == 'green'] <- 'KEGG: protein targeting to ER'
names(keyvals)[keyvals == 'yellow'] <- 'KEGG: ribonucleoprotein complex biogenesis'
names(keyvals)[keyvals == 'lightgrey'] <- "Others"
table(keyvals)

text_size = 14
point_size = 2.0

EnhancedVolcano(markers_all,
    lab = markers_all$gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = unique(c(gene_list1, gene_list2, gene_list3, gene_list4)),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~- Log[10]~ 'adjusted p-value'),
    title = 'Riboflavin vs. EF + Ribo',
    subtitle = NULL,
    caption = NULL,
    pCutoff = 0.01,
    FCcutoff = 1.0,
    pointSize = point_size,
    labSize = text_size/3,
    shape = c(18, 15, 16, 17),
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'left',
    legendLabSize = text_size,
    legendIconSize = point_size,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    arrowheads = FALSE) +
  labs(color = "Pathway") + 
  xlim(-3.0,3.0) +
  ylim(1, 60) +
  theme_classic() +
  theme(plot.margin = margin(0,0,0,0),
        text = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size, color = "black"), 
        axis.line = element_line(size = 0.4), 
        legend.position = "right", 
        legend.justification = "center", 
        legend.margin=margin(0,0,0,0),  
        legend.box.margin=margin(0,0,0,0))

########################################

markers_all <- read.csv(
  "./mydata/tcells/files_from_c/riboflavinvsEf_ribo_deg_all.csv",
  row.names = 1
)
go_resuts   <- read.csv(
  "./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_GOenrich_significant.csv"
)
kegg_resuts <- read.csv(
  "./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_KEGGenrich_significant.csv"
)

gene_list1 <- strsplit(
  go_resuts[go_resuts$Description == "protein targeting to ER", "geneName"],
  split = "/"
)[[1]]

gene_list2 <- strsplit(
  go_resuts[go_resuts$Description == "ribonucleoprotein complex biogenesis", "geneName"],
  split = "/"
)[[1]]

gene_list3 <- strsplit(
  kegg_resuts[kegg_resuts$Description == "Endocytosis", "geneName"],
  split = "/"
)[[1]]

gene_list4 <- strsplit(
  kegg_resuts[kegg_resuts$Description == "Protein processing in endoplasmic reticulum", "geneName"],
  split = "/"
)[[1]]

highlight_genes <- unique(c(gene_list1, gene_list2, gene_list3, gene_list4))

sig_highlight_genes <- with(
  markers_all,
  gene_name[(gene_name %in% highlight_genes) & (padj <= 0.05) & (abs(log2FoldChange) > log2(1.5))]
)

markers_all <- markers_all[
  order(markers_all$gene_name %in% sig_highlight_genes),
  ,
  drop = FALSE
]

gene_list1 <- gene_list1[gene_list1 %in% sig_highlight_genes]
gene_list2 <- gene_list2[gene_list2 %in% sig_highlight_genes]
gene_list3 <- gene_list3[gene_list3 %in% sig_highlight_genes]
gene_list4 <- gene_list4[gene_list4 %in% sig_highlight_genes]

keyvals <- ifelse(
  markers_all$gene_name %in% gene_list1, "yellow",
  ifelse(
    markers_all$gene_name %in% gene_list2, "red",
    ifelse(
      markers_all$gene_name %in% gene_list3, "green",
      ifelse(
        markers_all$gene_name %in% gene_list4, "blue",
        "lightgrey"
      )
    )
  )
)

keyvals[is.na(keyvals)] <- "lightgrey"
names(keyvals)[keyvals == "yellow"]      <- "GO: protein targeting to ER"
names(keyvals)[keyvals == "red"]       <- "GO: ribonucleoprotein complex biogenesis"
names(keyvals)[keyvals == "green"]     <- "KEGG: Endocytosis"
names(keyvals)[keyvals == "blue"]    <- "KEGG: Protein processing in ER"
names(keyvals)[keyvals == "lightgrey"] <- "Others"

table(keyvals)

text_size  <- 14
point_size <- 2.0

EnhancedVolcano(
  markers_all,
  lab            = markers_all$gene_name,
  x              = "log2FoldChange",
  y              = "padj",
  selectLab      = sig_highlight_genes,                                      
  xlab           = bquote(~Log[2]~ "fold change"),
  ylab           = bquote(~-Log[10]~ "adjusted p-value"),
  title          = 'Riboflavin vs. EF + Ribo (|LFC| > 0.585 & adj. p-value < 0.05)',
  subtitle       = NULL,
  caption        = NULL,
  pCutoff        = 0.05,
  FCcutoff       = log2(1.5),                                                         
  pointSize      = point_size,
  labSize        = text_size / 3,
  shape          = c(18, 15, 16, 17),
  colCustom      = keyvals,
  colAlpha       = 1,
  legendPosition = "left",
  legendLabSize  = text_size,
  legendIconSize = point_size,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors   = "black",
  arrowheads      = FALSE, 
  max.overlaps = 20
) +
  labs(color = "Pathway") +
  xlim(-3, 3) +
  ylim(1, 60) +
  theme_classic() +
  theme(
    plot.margin       = margin(0, 0, 0, 0),
    text              = element_text(size = text_size),
    axis.title        = element_text(size = text_size),
    axis.text         = element_text(size = text_size, colour = "black"),
    legend.text       = element_text(size = text_size, colour = "black"),
    axis.line         = element_line(size = 0.4),
    legend.position   = "right",
    legend.justification = "center",
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

##############################

markers_all = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo_deg_all.csv", row.names = 1)
go_resuts = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_GOenrich_significant.csv")
gene_list1 = strsplit(go_resuts[go_resuts$Description == "protein targeting to ER", "geneName"], split = "/")[[1]]
gene_list2 = strsplit(go_resuts[go_resuts$Description == "ribonucleoprotein complex biogenesis", "geneName"], split = "/")[[1]]
kegg_resuts = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo.down_KEGGenrich_significant.csv")
gene_list3 = strsplit(kegg_resuts[kegg_resuts$Description == "Endocytosis", "geneName"], split = "/")[[1]]
gene_list4 = strsplit(kegg_resuts[kegg_resuts$Description == "Protein processing in endoplasmic reticulum", "geneName"], split = "/")[[1]]

highlight_genes <- unique(c(gene_list1, gene_list2, gene_list3, gene_list4))
markers_all <- markers_all[order(markers_all$gene_name %in% highlight_genes), , drop = FALSE]

keyvals <- ifelse(markers_all$gene_name %in% gene_list1, "blue",
                  ifelse(markers_all$gene_name %in% gene_list2, "red",
                         ifelse(markers_all$gene_name %in% gene_list3, "green",
                                ifelse(markers_all$gene_name %in% gene_list4, "yellow",
                                       "lightgrey"))))
keyvals[is.na(keyvals)] <- "lightgrey"

names(keyvals)[keyvals == 'blue'] <- 'GO: protein targeting to ER'
names(keyvals)[keyvals == 'red'] <- 'GO: ribonucleoprotein complex biogenesis'
names(keyvals)[keyvals == 'green'] <- 'KEGG: protein targeting to ER'
names(keyvals)[keyvals == 'yellow'] <- 'KEGG: ribonucleoprotein complex biogenesis'
names(keyvals)[keyvals == 'lightgrey'] <- "Others"
table(keyvals)

text_size = 14
point_size = 2.0

EnhancedVolcano(markers_all,
                lab = markers_all$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = unique(c(gene_list1, gene_list2, gene_list3, gene_list4)),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~- Log[10]~ 'adjusted p-value'),
                title = 'Riboflavin vs. EF + Ribo',
                subtitle = NULL,
                caption = NULL,
                pCutoff = 0.01,
                FCcutoff = 1.0,
                pointSize = point_size,
                labSize = text_size/3,
                shape = c(18, 15, 16, 17),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = text_size,
                legendIconSize = point_size,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                arrowheads = FALSE) +
  labs(color = "Pathway") + 
  xlim(-3.0,3.0) +
  ylim(1, 60) +
  theme_classic() +
  theme(plot.margin = margin(0,0,0,0),
        text = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size, color = "black"), 
        axis.line = element_line(size = 0.4), 
        legend.position = "right", 
        legend.justification = "center", 
        legend.margin=margin(0,0,0,0),  
        legend.box.margin=margin(0,0,0,0))

########################################

markers_all = read.csv("./mydata/tcells/files_from_c/riboflavinvsEf_ribo_deg_all.csv", row.names = 1)
text_size = 14
point_size = 2.0
EnhancedVolcano(markers_all,
                lab = markers_all$gene_name,
                selectLab = unique(c("DEPP1", "FICD", "MT1G")),
                x = 'log2FoldChange', 
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~- Log[10]~ 'adjusted p-value'),
                title = 'Riboflavin vs. EF + Ribo (|LFC| > 1 & adj. p-value < 0.05)',
                titleLabSize = text_size,
                subtitle = NULL,
                caption = NULL,
                legendPosition = "none",
                pCutoff = 0.05,
                drawConnectors = T,
                axisLabSize = text_size, 
                labSize = text_size/3,
                pointSize = point_size,
                FCcutoff = 1.0,
                colAlpha = 0.8) + 
  xlim(-4.5,3.5) +
  # xlim(-21, 3) +
  ylim(1, 176) +
  theme_classic() +
  theme(plot.margin = margin(0,0,0,0),
        text = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size, color = "black"), 
        axis.line = element_line(size = 0.4), 
        legend.position = "none", 
        legend.justification = "center", 
        legend.margin=margin(0,0,0,0),  
        legend.box.margin=margin(0,0,0,0))
