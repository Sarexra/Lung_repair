if (!require("devtools")) 
  install.packages("devtools")

library(readxl)
library(tidyverse)
library(dplyr)


## 1. Exvivo tidal ventilation delayes alveolar epithelial wound closure ####
### 1.1. Wound healing BEAS-2B + BALF ####
data <- read_excel("Wound_healing_BEAS.xlsx",
                   sheet = 1)
time=7
time_zero<- subset(data, data$Hora == 0)
time_zero<-time_zero[rep(1:nrow(time_zero),time),]
data <-data %>%
  mutate(Porcentaje=data$Area/time_zero$Area)
data$Porcentaje <- data$Area/time_zero$Area 
data<-data %>% 
  group_by(Condicion,Hora, BALF) %>% 
  mutate(mean=mean(Porcentaje), sd=sd(Porcentaje))

data <-data %>% 
  mutate(group = paste(Condicion, BALF, Replica, sep="_"))
data <- data %>% 
  filter( group != "Control_Control_1")
ag.data <- data %>% 
  group_by(Hora, Condicion)%>% 
  summarise(media=mean(Porcentaje), sd=sd(Porcentaje), se=sd/sqrt(length(Porcentaje)), upp=se+media, down=media-se)

ag.data %>%
  filter(Condicion != "Control") %>% 
  ggplot(aes(x=Hora, y=media, col=Condicion))+
  geom_line(size=1.25, alpha=0.8)+
  geom_point(size=2)+
  geom_ribbon(aes(ymin=down, ymax=upp, fill=Condicion), alpha=0.5,  color = NA)+
  scale_y_continuous(labels= scales::percent_format(scale= 100), limits=c(0,1))+
  labs(title="BEAS-2B wound healing",x="Time(h)", y = "Percentage", color = "Condition",   
       fill  = "Condition")+
  scale_color_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                     labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  scale_fill_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                    labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  theme_minimal()+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 19), 
    legend.position = "right",
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 18),
    axis.line = element_line(color = "grey50", linewidth = 1))

#### Statistic of wound healing BEAS-2B ####
library(emmeans)
modelo <- aov(Porcentaje ~ Hora + Condicion + Error(Replica/Hora),
              data = data)

emm <- emmeans(modelo, ~ Condicion | Hora)
pairs(emm, adjust = "tukey") 


### 1.2. RNAseq####
## Sheet2 -> BEAS  
library(readODS)
samples_BEAS <- read_ods("Samples_BEAS.ods", sheet=1)
samples_BEAS$condition[samples_BEAS$condition == "UA"] <- "VM"
samples_BEAS$condition[samples_BEAS$condition == "UC"] <- "CPAP"
samples_BEAS$patient <- c("VM1", "VM2", "VM3", "CPAP1", "CPAP2", "CPAP3")
 
### 1.2.1. Count matrix ####
library(AnnotationHub)
library(ensembldb)
library(tximport)
library(DESeq2)
ah_BEAS <- AnnotationHub()
edb_BEAS <- ah_BEAS[["AH73986"]]

txi_BEAS<- readRDS("transcripts_BEAS_github_rds")
View(txi_BEAS$counts)


### 1.2.2. DESeq2 ####
#BiocManager::install("DESeq2")
samples_BEAS$condition <- factor(samples_BEAS$condition, levels = c("CPAP","VM"))
ddsTxi_BEAS <- DESeqDataSetFromTximport(txi_BEAS, colData=samples_BEAS, design= ~ condition)  
dim(ddsTxi_BEAS) #34528

keep_BEAS<-rowSums(counts(ddsTxi_BEAS))>=ncol(ddsTxi_BEAS)
ddsTxi_BEAS<-ddsTxi_BEAS[keep_BEAS,]
dim(ddsTxi_BEAS)#14217

dds_BEAS<-DESeq(ddsTxi_BEAS)
genenames_BEAS<-ensembldb::select(edb_BEAS, keys=rownames(dds_BEAS), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds_BEAS)<-genenames_BEAS$SYMBOL


### 1.2.3. Differential expression ####
res_BEAS<-results(dds_BEAS)
resOrdered_BEAS <- res_BEAS[order(res_BEAS$log2FoldChange),] 
resOrdered_BEAS
summary(res_BEAS, alpha=0.05) 
sig.genes_BEAS<-resOrdered_BEAS[!is.na(resOrdered_BEAS$padj) & resOrdered_BEAS$padj<0.05,]
#154 genes con ED, 103-up. 51-down
sig.genes_BEAS <- as.data.frame(sig.genes_BEAS)
sig.genes_BEAS$genes <- rownames(sig.genes_BEAS)
resOrdered_BEAS <- as.data.frame(resOrdered_BEAS) 
resultsNames(dds_BEAS) # "condition_VM_vs_CPAP"


### 1.2.4. Volcano Plot ####
library(ggrepel)
library(ggplot2)
my_gene_BEAS <- c("IL6", "EGFR", "PPIA", "TGFB2", "KRT7", "FN1", "CASP3", "FOSL2", "ICAM1", "PPP2CA", "LRP1", "CANX")

sig.genes_BEAS$significance[sig.genes_BEAS$log2FoldChange > 1 & sig.genes_BEAS$padj < 0.05] <- "Up"
sig.genes_BEAS$significance[sig.genes_BEAS$log2FoldChange < -1 & sig.genes_BEAS$padj < 0.05] <- "Down"

my_gene_BEAS <- resOrdered_BEAS[rownames(resOrdered_BEAS) %in% my_gene_BEAS,]

ggplot(as.data.frame(res_BEAS), aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(col = padj<0.05), alpha = 0.8)+
  xlim(-11,11)+
  ylim(0,13)+
  geom_point(data = my_gene_BEAS, color = "black", size = 1.8)+
  geom_label_repel(data=my_gene_BEAS, aes(label=rownames(my_gene_BEAS)),
                   force = 2, 
                   nudge_y = 1,
                   size=6)+
  coord_cartesian(ylim = c(0,8), xlim = c(-5, 5), expand=FALSE)+
  labs(title = "BEAS-2B MV vs CPAP", 
       x = "log2 Fold Change", 
       y = "-log10(padj)") +
  scale_color_manual(
    values = c("FALSE" = "#40425c", "TRUE" = "orange", "NA" = none),
    labels = c("FALSE" = "p > 0.05", "TRUE" = "p < 0.05"),
    name = "Padj",
    na.translate = FALSE  
  )+
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 15), 
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14))


### 1.2.5. Heatmap ####
library(ComplexHeatmap)
library(circlize)

mat_BEAS <- counts(dds_BEAS, normalized = TRUE)[rownames(sig.genes_BEAS), ]
mat_BEAS <- (mat_BEAS - rowMeans(mat_BEAS)) / rowSds(mat_BEAS)
select_BEAS <- order(rowMeans(mat_BEAS), decreasing = TRUE)
mat_BEAS <- mat_BEAS[select_BEAS, ]
df_BEAS <- as.data.frame(colData(dds_BEAS)[,c("condition")])
rownames(df_BEAS) <- colnames(mat_BEAS)
colnames(df_BEAS) <- "Condition"

df_BEAS$Condition <- as.character(df_BEAS$Condition)
df_BEAS$Condition[df_BEAS$Condition == "VM"]   <- "MV"
df_BEAS$Condition[df_BEAS$Condition == "CPAP"] <- "CPAP"
df_BEAS$Condition <- factor(df_BEAS$Condition, levels = c("CPAP", "MV"))

Condition_colors <- c("MV" = "orange", "CPAP" = "#343e68")
top_annot_BEAS <- HeatmapAnnotation(
  Condition = df_BEAS$Condition,
  col = list(Condition = Condition_colors),
  show_annotation_name = TRUE)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "orangered"))

Heatmap(mat_BEAS,
        name = "Z-score",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp    = gpar(fontsize = 5),
        row_dend_side   = "left",
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title_gp  = gpar(fontsize = 6),
          labels_gp = gpar(fontsize = 6)),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = Condition_colors, col  = NA, alpha=0.6),
                                                            labels = c("MV", "CPAP"), 
                                                            labels_gp = gpar(col = "black", fontsize = 9)),
                                           height = unit(0.4, "cm")),
        column_km = 2,
        column_title = NULL,
        width           = unit(8, "cm"),
        height          = unit(25, "cm"),)


### 1.2.6. Enrichment analysis ####
library(clusterProfiler)
dds_BEAS<-DESeqDataSetFromTximport(txi_BEAS, colData = samples_BEAS, design = ~condition) 
dim(dds_BEAS)#34528     
keep_BEAS<-rowSums(counts(dds_BEAS))>=ncol(dds_BEAS)
dds_BEAS<-dds_BEAS[keep_BEAS,]
dim(dds_BEAS)#14517
dds_BEAS<-DESeq(dds_BEAS)
res_BEAS <- results(dds_BEAS)
original_gene_list <- res_BEAS$stat 
names(original_gene_list) <- rownames(res_BEAS)
gene_list_BEAS<-na.omit(original_gene_list)
gene_list_BEAS = sort(gene_list_BEAS, decreasing = TRUE)

### 1.2.7. GSE ####
library("org.Hs.eg.db", character.only = TRUE)
set.seed(1000) 
gse_BEAS <- gseGO(geneList=gene_list_BEAS, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 5, 
                  maxGSSize = 8000,
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")
dim(gse_BEAS) #83
View(gse_BEAS@result)

res_anon_BEAS <- res_BEAS
genenames_BEAS<-ensembldb::select(edb_BEAS, keys=rownames(res_anon_BEAS), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(res_anon_BEAS) <- genenames_BEAS$SYMBOL

### 1.2.8. Network gen-gen ####
library(ensembldb)
library(stringr)
library(biomaRt)
string <- read.table("9606.protein.links.v11.5.txt", header = T)
string$protein1 <- str_replace(string$protein1, "^9606\\.", "")
string$protein2 <- str_replace(string$protein2, "^9606\\.", "")
listEnsembl()
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "ensembl_peptide_id", "ensembl_transcript_id"),
                 filters="ensembl_peptide_id",
                 values = c(unique(string$protein1),unique(string$protein2)),
                 mart = mart)

string$gene1 <- results$ensembl_gene_id[match(string$protein1, results$ensembl_peptide_id)]
string$gene2 <- results$ensembl_gene_id[match(string$protein2, results$ensembl_peptide_id)]
string$gene1_symbol <- results$hgnc_symbol[match(string$protein1, results$ensembl_peptide_id)]
string$gene2_symbol <- results$hgnc_symbol[match(string$protein2, results$ensembl_peptide_id)]

### 1.2.9. BP - GO:0006952 defense response(BEAS)---- 
gse_1 <-gse_BEAS@result
gse_1 <- subset(gse_1, gse_1$ID=="GO:0006952") #IL6
gse_1$core_enrichment
res_GO_1 <- strsplit(gse_1$core_enrichment, "/")[[1]]

genenames_BEAS<-ensembldb::select(edb_BEAS, keys=res_GO_1, keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
res_GO_1<-as.vector(genenames_BEAS$SYMBOL  )
res_GO_1 <- as.data.frame(res_anon_BEAS[rownames(res_anon_BEAS) %in% res_GO_1,])

edges_1 <- string[string$combined_score >= 900,]
edges_1 <- edges_1[edges_1$gene1_symbol %in% rownames(res_GO_1)[!is.na(res_GO_1$padj)&res_GO_1$pvalue<0.05],]
edges_1 <- edges_1[edges_1$gene2_symbol %in% rownames(res_GO_1)[!is.na(res_GO_1$padj)&res_GO_1$pvalue<0.05 ],]
edges_1 <- edges_1[!is.na(edges_1$gene1_symbol),]
edges_1 <- edges_1[!is.na(edges_1$gene2_symbol),c("gene1", "gene2","combined_score")]

table(is.na(edges_1$gene1))
table(is.na(edges_1$gene1_symbol))
table(is.na(edges_1$gene2))
table(is.na(edges_1$gene2_symbol))

nodes_1 <- data.frame(nodes_1=unique(c(edges_1$gene1, edges_1$gene2)),
                      symbol = results$hgnc_symbol[match(unique(c(edges_1$gene1, edges_1$gene2)), results$ensembl_gene_id)])
nodes_1$LFC <- ifelse(nodes_1$symbol %in% rownames(res_GO_1)[res_GO_1$log2FoldChange>0],"UP", "DOWN")
nodes_1$FC <- res_GO_1$log2FoldChange[match(nodes_1$symbol, rownames(res_GO_1))]
nodes_1$FC.st <- (nodes_1$FC-mean(nodes_1$FC))/sd(nodes_1$FC)
nodes_1$sig <- ifelse(nodes_1$symbol %in% rownames(res_anon_BEAS)[res_BEAS$padj<0.05],"Sig", "NO_Sig" )

table(is.na(nodes_1))
library(tidygraph)
library(critcolors)

routes_tidy_1 <- tbl_graph(nodes=nodes_1,
                           edges=edges_1,
                           directed = T)
ggraph(routes_tidy_1, layout = "kk") + 
  geom_edge_link(aes(alpha=combined_score), width=0.2) + 
  geom_node_point(aes(col = sig, shape=LFC), size = 6, alpha = 0.7) +
  geom_node_label(aes(label = symbol),  size=3, nudge_y=0.2,
                  family = "sans") +
  theme_graph(base_family = "sans") +
  scale_colour_critcolors(discrete = TRUE)


### 1.2.10. BP - GO:0009653 anatomical structure morphogenesis (BEAS)---- 
gse_2 <-gse_BEAS@result
gse_2 <- subset(gse_2, gse_2$ID=="GO:0009653") #anatomical structure morphogenesis #EGFR #TGFB2 #FN1 #CASP3 #WDR1

gse_2$core_enrichment
res_GO_2 <- strsplit(gse_2$core_enrichment, "/")[[1]]

genenames_BEAS<-ensembldb::select(edb_BEAS, keys=res_GO_2, keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
res_GO_2<-as.vector(genenames_BEAS$SYMBOL  )
res_GO_2 <- as.data.frame(res_anon_BEAS[rownames(res_anon_BEAS) %in% res_GO_2,])

edges_2 <- string[string$combined_score >= 600,]
edges_2 <- edges_2[edges_2$gene1_symbol %in% rownames(res_GO_2)[!is.na(res_GO_2$padj)&res_GO_2$pvalue<0.05],]
edges_2 <- edges_2[edges_2$gene2_symbol %in% rownames(res_GO_2)[!is.na(res_GO_2$padj)&res_GO_2$pvalue<0.05 ],]
edges_2 <- edges_2[!is.na(edges_2$gene1_symbol),]
edges_2 <- edges_2[!is.na(edges_2$gene2_symbol),c("gene1", "gene2","combined_score")]

table(is.na(edges_2$gene1))
table(is.na(edges_2$gene1_symbol))
table(is.na(edges_2$gene2))
table(is.na(edges_2$gene2_symbol))

nodes_2 <- data.frame(nodes=unique(c(edges_2$gene1, edges_2$gene2)),
                      symbol = results$hgnc_symbol[match(unique(c(edges_2$gene1, edges_2$gene2)), results$ensembl_gene_id)])
nodes_2$LFC <- ifelse(nodes_2$symbol %in% rownames(res_GO_2)[res_GO_2$log2FoldChange>0],"UP", "DOWN")
nodes_2$FC <- res_GO_2$log2FoldChange[match(nodes_2$symbol, rownames(res_GO_2))]
nodes_2$FC.st <- (nodes_2$FC-mean(nodes_2$FC))/sd(nodes_2$FC)
nodes_2$sig <- ifelse(nodes_2$symbol %in% rownames(res_anon_BEAS)[res_BEAS$padj<0.05],"Sig", "NO_Sig" )

table(is.na(nodes_2))

routes_tidy_2 <- tbl_graph(nodes=nodes_2,
                           edges=edges_2,
                           directed = T)
ggraph(routes_tidy_2, layout = "kk") + 
  geom_edge_link(aes(alpha=combined_score), width=0.2) + 
  geom_node_point(aes(col = sig, shape=LFC), size = 6, alpha = 0.7) +
  geom_node_label(aes(label = symbol),  size=3, nudge_y=0.2, family = "sans") +
  theme_graph(base_family = "sans") +
  scale_colour_critcolors(discrete = TRUE)



# 2. Exvivo tidal ventilation promotes a fibroproliferative response####
### 2.1. Wound healing MRC5 + BALF####
data <- read_excel("Wound_healing_MRC5.xlsx",
                   sheet = 1)
time_MRC5=9 
time_zero_MRC5<- subset(data, data$Hora == 0)
time_zero_MRC5<-time_zero_MRC5[rep(1:nrow(time_zero_MRC5),time_MRC5),]
data <-data %>%
  mutate(Porcentaje=data$Area/time_zero_MRC5$Area)
data$Porcentaje <- data$Area/time_zero_MRC5$Area
data<-data %>% 
  group_by(Condicion,Hora, BALF) %>% 
  mutate(mean=mean(Porcentaje), sd=sd(Porcentaje))
data <-data %>% 
  mutate(group = paste(Condicion, BALF, sep="_"))
data <- data %>% 
  filter(group != "CPAP_4") 
ag.data <- data %>% 
  group_by(Hora, Condicion)%>% 
  summarise(media=mean(Porcentaje), sd=sd(Porcentaje), se=sd/sqrt(length(Porcentaje)), upp=se+media, down=media-se)

ag.data %>%
  ggplot(aes(x=Hora, y=media, col=Condicion))+
  geom_line(size=1.25, alpha=0.8)+
  geom_point(size=2)+
  geom_ribbon(aes(ymin=down, ymax=upp, fill=Condicion), alpha=0.5,  color = NA)+
  scale_y_continuous(labels= scales::percent_format(scale= 100), limits=c(0,1))+
  labs(title="MRC5 wound healing",x="Time(h)", y = "Percentage", color = "Condition",   
       fill  = "Condition")+
  scale_color_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                     labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  scale_fill_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                    labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  theme_minimal()+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 19), #, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 18),
    axis.line = element_line(color = "grey50", linewidth = 1))


#### Statistic of wound healing MRC5 ####
library(emmeans)
modelo <- aov(Porcentaje ~ Hora + Condicion + Error(Replica/Hora),
              data = data)

emm <- emmeans(modelo, ~ Condicion | Hora)
pairs(emm, adjust = "tukey") 


### 2.2. RNAseq####
samples_MRC5 <- read_ods("Samples_MRC5.ods", sheet = 1)
samples_MRC5$condition[samples_MRC5$condition == "Sin CEC"] <- "VM"
samples_MRC5$condition[samples_MRC5$condition == "Con CEC"] <- "CPAP"

### 2.2.1. Count matrix ####

ah_MRC5 <- AnnotationHub()
edb_MRC5 <- ah_MRC5[["AH73986"]]

txi_MRC5<- readRDS("transcripts_MRC5_github.rds")
View(txi_MRC5$counts)


### 2.2.2. DESeq2 ####
#BiocManager::install("DESeq2")
samples_MRC5$condition <- factor(samples_MRC5$condition, levels = c("CPAP","VM"))
ddsTxi_MRC5 <- DESeqDataSetFromTximport(txi_MRC5, colData=samples_MRC5, design= ~condition) 
dim(ddsTxi_MRC5) #34528

keep_MRC5<-rowSums(counts(ddsTxi_MRC5))>=ncol(ddsTxi_MRC5)
ddsTxi_MRC5<-ddsTxi_MRC5[keep_MRC5,]
dim(ddsTxi_MRC5)#14512

dds_MRC5<-DESeq(ddsTxi_MRC5)
genenames_MRC5<-AnnotationDbi::select(edb_MRC5, keys=rownames(dds_MRC5), keytype = "GENEID", columns=c("SYMBOL", "GENEID"))
rownames(dds_MRC5)<-genenames_MRC5$SYMBOL


### 2.2.3. Differential expression ####
res_MRC5<-results(dds_MRC5)
resOrdered_MRC5 <- res_MRC5[order(res_MRC5$log2FoldChange),] 
resOrdered_MRC5
summary(res_MRC5, alpha=0.05)
#583 genes con ED, 279-up. 304-down
sig.genes_MRC5<-resOrdered_MRC5[!is.na(resOrdered_MRC5$padj) & resOrdered_MRC5$padj<0.05,]
sig.genes_MRC5 <- as.data.frame(sig.genes_MRC5)
sig.genes_MRC5$genes <- rownames(sig.genes_MRC5)
resOrdered_MRC5 <- as.data.frame(resOrdered_MRC5)

resultsNames(dds_MRC5) #"condition_VM_vs_CPAP"


### 2.2.4. Volcano Plot ####
my.genes_MRC5 <- c("TGFB1", "IL6", "LIF", "MMP1", "MMP10", "PPIA", "BMP2",
                   "APLP1", "ITGB8", "LRP6", "ITGB5", "ITGA2", "LRP8", "SDC2", "HAVCR2", "F2R", "F3", "ERBB4", "FAS",
                   "IL1RAPL1", "LPAR3", "DDR1", "TSPAN14", "ADORA1", "TNFRSF19")
sig.genes_MRC5$significance[sig.genes_MRC5$log2FoldChange > 1 & sig.genes_MRC5$padj < 0.05] <- "Up"
sig.genes_MRC5$significance[sig.genes_MRC5$log2FoldChange < -1 & sig.genes_MRC5$padj < 0.05] <- "Down"

my.genes_MRC5 <- resOrdered_MRC5[rownames(resOrdered_MRC5) %in% my.genes_MRC5,]

ggplot(as.data.frame(res_MRC5), aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(col = padj<0.05), alpha = 0.8)+
  xlim(-11,11)+
  ylim(0,13)+
  geom_point(data = my.genes_MRC5, color = "black", size = 1.8)+
  geom_label_repel(data=my.genes_MRC5, aes(label=rownames(my.genes_MRC5)),
                   force = 2, 
                   nudge_y = 1,
                   size=6)+
  coord_cartesian(ylim = c(0,8), xlim = c(-5, 5), expand=FALSE)+
  labs(title = "MRC5 MV vs CPAP", 
       x = "log2 Fold Change", 
       y = "-log10(padj)") +
  scale_color_manual(
    values = c("FALSE" = "#40425c", "TRUE" = "orange", "NA" = none),
    labels = c("FALSE" = "p > 0.05", "TRUE" = "p < 0.05"),
    name = "Padj",
    na.translate = FALSE 
  )+
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 15), 
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14))


### 2.2.5. Heatmap ####
library(ComplexHeatmap)
library(circlize)

mat_MRC5 <- counts(dds_MRC5, normalized = TRUE)[rownames(sig.genes_MRC5), ]
mat_MRC5 <- (mat_MRC5 - rowMeans(mat_MRC5)) / rowSds(mat_MRC5)
select_MRC5 <- order(rowMeans(mat_MRC5), decreasing = TRUE)
mat_MRC5 <- mat_MRC5[select_MRC5, ]
mat_MRC5 <- mat_MRC5[, order(samples_MRC5$condition)]
mat_MRC5 <- mat_MRC5[select_MRC5, ]
n_genes_MRC5 <- 154
mat_MRC5 <- mat_MRC5[1:n_genes_MRC5, ]

df_MRC5 <- as.data.frame(colData(dds_MRC5)[,c("condition")])
rownames(df_MRC5) <- colnames(mat_MRC5)
colnames(df_MRC5) <- "Condition"

df_MRC5$Condition <- as.character(df_MRC5$Condition)
df_MRC5$Condition[df_MRC5$Condition == "VM"]   <- "MV"
df_MRC5$Condition[df_MRC5$Condition == "CPAP"] <- "CPAP"
df_MRC5$Condition <- factor(df_MRC5$Condition, levels = c("CPAP", "MV"))

Condition_colors <- c("MV" = "orange", "CPAP" = "#343e68")
top_annot_MRC5 <- HeatmapAnnotation(
  Condition = df_MRC5$Condition,
  col = list(Condition = Condition_colors),
  show_annotation_name = TRUE)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "orangered"))

Heatmap(mat_MRC5,
        name = "Z-score",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp    = gpar(fontsize = 5),
        row_dend_side   = "left",
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title_gp  = gpar(fontsize = 6),
          labels_gp = gpar(fontsize = 6)),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = Condition_colors, col  = NA, alpha=0.6),
                                                            labels = c("MV", "CPAP"), 
                                                            labels_gp = gpar(col = "black", fontsize = 9)),
                                           height = unit(0.4, "cm")),
        column_km = 2,
        column_title = NULL,
        width           = unit(8, "cm"),
        height          = unit(25, "cm"),)



### 2.2.6. Enrichment analysis ####
dds_MRC5<-DESeqDataSetFromTximport(txi_MRC5, colData = samples_MRC5, design = ~condition) 
dim(dds_MRC5)#34528     
keep_MRC5<-rowSums(counts(dds_MRC5))>=ncol(dds_MRC5)
dds_MRC5<-dds_MRC5[keep_MRC5,]
dim(dds_MRC5)#14512
dds_MRC5<-DESeq(dds_MRC5)
res_MRC5 <- results(dds_MRC5)
original_gene_list_MRC5 <- res_MRC5$stat 
names(original_gene_list_MRC5) <- rownames(res_MRC5)
gene_list_MRC5<-na.omit(original_gene_list_MRC5)
gene_list_MRC5 = sort(gene_list_MRC5, decreasing = TRUE)

### 2.2.7. GSE ####
set.seed(1000) 
gse_MRC5 <- gseGO(geneList=gene_list_MRC5, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 20, 
                  maxGSSize = 100,
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")
dim(gse_MRC5) #188
View(gse_MRC5@result)

res_anon_MRC5 <- res_MRC5
genenames_MRC5<-ensembldb::select(edb_MRC5, keys=rownames(res_anon_MRC5), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(res_anon_MRC5) <- genenames_MRC5$SYMBOL


### 2.2.8. Network gen-gen ####
string <- read.table("9606.protein.links.v11.5.txt", header = T)
string$protein1 <- str_replace(string$protein1, "^9606\\.", "")
string$protein2 <- str_replace(string$protein2, "^9606\\.", "")
listEnsembl()
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "ensembl_peptide_id", "ensembl_transcript_id"),
                 filters="ensembl_peptide_id",
                 values = c(unique(string$protein1),unique(string$protein2)),
                 mart = mart)

string$gene1 <- results$ensembl_gene_id[match(string$protein1, results$ensembl_peptide_id)]
string$gene2 <- results$ensembl_gene_id[match(string$protein2, results$ensembl_peptide_id)]
string$gene1_symbol <- results$hgnc_symbol[match(string$protein1, results$ensembl_peptide_id)]
string$gene2_symbol <- results$hgnc_symbol[match(string$protein2, results$ensembl_peptide_id)]


### 2.2.9. MF - GO:0008083 growth factor activity(MRC5)---- 
library(ggraph)
gse_4 <-gse_MRC5@result
gse_4  <- subset(gse_4 , gse_4 $ID=="GO:0008083") 
gse_4 $core_enrichment
res_GO_4  <- strsplit(gse_4$core_enrichment, "/")[[1]]

genenames_MRC5<-ensembldb::select(edb_MRC5, keys=res_GO_4, keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
res_GO_4 <-as.vector(genenames_MRC5$SYMBOL  )
res_GO_4  <- as.data.frame(res_anon_MRC5[rownames(res_anon_MRC5) %in% res_GO_4,])

edges_4 <- string[string$combined_score >= 400,]
edges_4 <- edges_4[edges_4$gene1_symbol %in% rownames(res_GO_4)[!is.na(res_GO_4$padj)&res_GO_4$pvalue<0.05],]
edges_4 <- edges_4[edges_4$gene2_symbol %in% rownames(res_GO_4)[!is.na(res_GO_4$padj)&res_GO_4$pvalue<0.05 ],]
edges_4 <- edges_4[!is.na(edges_4$gene1_symbol),]
edges_4 <- edges_4[!is.na(edges_4$gene2_symbol),c("gene1", "gene2","combined_score")]

table(is.na(edges_4$gene1))
table(is.na(edges_4$gene1_symbol))
table(is.na(edges_4$gene2))
table(is.na(edges_4$gene2_symbol))

nodes_4 <- data.frame(nodes=unique(c(edges_4$gene1, edges_4$gene2)),
                      symbol = results$hgnc_symbol[match(unique(c(edges_4$gene1, edges_4$gene2)), results$ensembl_gene_id)])
nodes_4$LFC <- ifelse(nodes_4$symbol %in% rownames(res_GO_4)[res_GO_4$log2FoldChange>0],"UP", "DOWN")
nodes_4$FC <- res_GO_4$log2FoldChange[match(nodes_4$symbol, rownames(res_GO_4))]
nodes_4$FC.st <- (nodes_4$FC-mean(nodes_4$FC))/sd(nodes_4$FC)
nodes_4$sig <- ifelse(nodes_4$symbol %in% rownames(res_anon_MRC5)[res_MRC5$padj<0.05],"Sig", "NO_Sig" )

table(is.na(nodes_4))

routes_tidy_4 <- tbl_graph(nodes=nodes_4,
                           edges=edges_4,
                           directed = T)
ggraph(routes_tidy_4, layout = "kk") + 
  geom_edge_link(aes(alpha=combined_score), width=0.2) + 
  geom_node_point(aes(col = sig, shape=LFC), size = 6, alpha = 0.7) +
  geom_node_label(aes(label = symbol),  size=3, nudge_y=0.2, family = "sans") +
  theme_graph(base_family = "sans") +
  scale_colour_critcolors(discrete = TRUE)

# 3. Tocilizumab reverts tidal ventilation exvivo effects in epithelial cells####

### 3.1. Wound healing-BEAS-2B + BALF + Tocilizumab ####
data <- read_excel("Wound_healing_BEAS.xlsx",                   
                   sheet = 2)

time_BEAS_toci=7 
time_zero_BEAS_toci<- subset(data, data$Hora == 0)
time_zero_BEAS_toci<-time_zero_BEAS_toci[rep(1:nrow(time_zero_BEAS_toci),time_BEAS_toci),]

data <-data %>%
  mutate(Porcentaje=data$Area/time_zero_BEAS_toci$Area)
data$Porcentaje <- data$Area/time_zero_BEAS_toci$Area
data<-data %>% 
  group_by(Condicion,Hora, BALF) %>% 
  mutate(mean=mean(Porcentaje), sd=sd(Porcentaje))
data <-data %>% 
  mutate(group = paste(Condicion, BALF, sep="_"))
ag.data <- data %>% 
  group_by(Hora, Condicion)%>% 
  summarise(media=mean(Porcentaje), sd=sd(Porcentaje), se=sd/sqrt(length(Porcentaje)), upp=se+media, down=media-se)

ag.data %>%
  ggplot(aes(x=Hora, y=media, col=Condicion))+
  geom_line(size=1.25, alpha=0.8)+
  geom_point(size=2)+
  geom_ribbon(aes(ymin=down, ymax=upp, fill=Condicion), alpha=0.5,  color = NA)+
  scale_y_continuous(labels= scales::percent_format(scale= 100), limits=c(0,1))+
  labs(title="BEAS+Tocilizumab (10ug/mL) wound healing",x="Time(h)", y = "Percentage", color = "Condition",  
       fill  = "Condition")+
  scale_color_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                     labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  scale_fill_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                    labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  theme_minimal()+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 19),
    legend.position = "right",
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 18),
    axis.line = element_line(color = "grey50", linewidth = 1))


#### Statistic of wound healing BEAS + tocilizumab ####
library(emmeans)
modelo <- aov(Porcentaje ~ Hora + Condicion + Error(Replica/Hora),
              data = data)

emm <- emmeans(modelo, ~ Condicion | Hora)
pairs(emm, adjust = "tukey") #p=0.1559

### 3.2. RNAseq ####
samples_BEAS_toci <- read_ods("Samples_BEAS.ods", sheet=2)

samples_BEAS_toci$condition <- as.factor(samples_BEAS_toci$condition)

### 3.2.1. Count matrix ####
ah_BEAS_toci <- AnnotationHub()
edb_BEAS_toci <- ah_BEAS_toci[["AH73986"]]

txi_BEAS_toci<- readRDS("transcripts_BEAS_toci_github.rds") 
View(txi_BEAS_toci$counts)

### 3.2.2. DESeq2 ####
#BiocManager::install("DESeq2")

ddsTxi_BEAS_toci <- DESeqDataSetFromTximport(txi_BEAS_toci, colData=samples_BEAS_toci, design= ~condition) 
dim(ddsTxi_BEAS_toci) #34528

keep_BEAS_toci<-rowSums(counts(ddsTxi_BEAS_toci))>=ncol(ddsTxi_BEAS_toci)
ddsTxi_BEAS_toci<-ddsTxi_BEAS_toci[keep_BEAS_toci,]
dim(ddsTxi_BEAS_toci)#13371

dds_BEAS_toci<-DESeq(ddsTxi_BEAS_toci)
genenames_BEAS_toci<-ensembldb::select(edb_BEAS_toci, keys=rownames(dds_BEAS_toci), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds_BEAS_toci)<-genenames_BEAS_toci$SYMBOL

### 3.2.3. Differential expression ####
res_BEAS_toci<-results(dds_BEAS_toci)
resOrdered_BEAS_toci <- res_BEAS_toci[order(res_BEAS_toci$log2FoldChange),] 
resOrdered_BEAS_toci
summary(res_BEAS_toci, alpha=0.05) 
#BEAS+toci:24 genes con ED, 14-up. 10-down
sig.genes_BEAS_toci<-resOrdered_BEAS_toci[!is.na(resOrdered_BEAS_toci$padj) & resOrdered_BEAS_toci$padj<0.05,]
sig.genes_BEAS_toci <- as.data.frame(sig.genes_BEAS_toci)
sig.genes_BEAS_toci$genes <- rownames(sig.genes_BEAS_toci)
resOrdered_BEAS_toci <- as.data.frame(resOrdered_BEAS_toci)

resultsNames(dds_BEAS_toci) #SinCEC vs ConCEC


### 3.2.4. Volcano Plot #### 
my_gene_BEAS_toci <- c("IL6", "EGFR", "PPIA", "TGFB2", "KRT7", "FN1", "CASP3", "FOSL2", "ICAM1", "PPP2CA", "LRP1", "CANX")

sig.genes_BEAS_toci$significance[sig.genes_BEAS_toci$log2FoldChange > 1 & sig.genes_BEAS_toci$padj < 0.05] <- "Up"
sig.genes_BEAS_toci$significance[sig.genes_BEAS_toci$log2FoldChange < -1 & sig.genes_BEAS_toci$padj < 0.05] <- "Down"

my_gene_BEAS_toci <- resOrdered_BEAS_toci[rownames(resOrdered_BEAS_toci) %in% my_gene_BEAS_toci,]

ggplot(as.data.frame(res_BEAS_toci), aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(col = padj<0.05), alpha = 0.8)+
  xlim(-11,11)+
  ylim(0,13)+
  geom_point(data = my_gene_BEAS_toci, color = "black", size = 1.8)+
  geom_label_repel(data=my_gene_BEAS_toci, aes(label=rownames(my_gene_BEAS_toci)),
                   force = 2, 
                   nudge_y = 1,
                   size=6)+
  coord_cartesian(ylim = c(0,8), xlim = c(-5, 5), expand=FALSE)+
  labs(title = "BEAS-2B MV vs CPAP + Tocilizumab (10ug/mL)", 
       x = "log2 Fold Change", 
       y = "-log10(padj)") +
  scale_color_manual(
    values = c("FALSE" = "#40425c", "TRUE" = "orange", "NA" = none),
    labels = c("FALSE" = "p > 0.05", "TRUE" = "p < 0.05"),
    name = "Padj",
    na.translate = FALSE  
  )+
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 15), #, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14))



# 4. Tocilizumab reverts tidal ventilation exvivo effects in mesenchymal cells####
### 4.1. Wound healing-MRC5 + BALF + Tocilizumab ####

data <- read_excel("Wound_healing_MRC5.xlsx",                   
                   sheet = 2)
time_MRC5_toci=9
time_zero_MRC5_toci<- subset(data, data$Hora == 0)
time_zero_MRC5_toci<-time_zero_MRC5_toci[rep(1:nrow(time_zero_MRC5_toci),time_MRC5_toci),]

data <-data %>%
  mutate(Porcentaje=data$Area/time_zero_MRC5_toci$Area)
data$Porcentaje <- data$Area/time_zero_MRC5_toci$Area

data<-data %>% 
  group_by(Condicion,Hora, BALF) %>% 
  mutate(mean=mean(Porcentaje), sd=sd(Porcentaje))
data <-data %>% 
  mutate(group = paste(Condicion, BALF, sep="_"))
ag.data <- data %>% 
  group_by(Hora, Condicion)%>% 
  summarise(media=mean(Porcentaje), sd=sd(Porcentaje), se=sd/sqrt(length(Porcentaje)), upp=se+media, down=media-se)

ag.data %>%
  ggplot(aes(x=Hora, y=media, col=Condicion))+
  geom_line(size=1.25, alpha=0.8)+
  geom_point(size=2)+
  geom_ribbon(aes(ymin=down, ymax=upp, fill=Condicion), alpha=0.5,  color = NA)+
  scale_y_continuous(labels= scales::percent_format(scale= 100), limits=c(0,1))+
  labs(title="MRC5+Tocilizumab (10ug/mL) wound healing",x="Time(h)", y = "Percentage", color = "Condition",   
       fill  = "Condition")+
  scale_color_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                     labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  scale_fill_manual(values = c("CPAP" = "#2a3d7f", "VM" = "orange"),
                    labels = c("CPAP" = "CPAP", "VM" = "MV"))+
  theme_minimal()+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 19), 
    legend.position = "right",
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 18),
    axis.line = element_line(color = "grey50", linewidth = 1))


#### Statistic of wound healing MRC5 + tocilizumab ####
library(emmeans)
modelo <- aov(Porcentaje ~ Hora + Condicion + Error(Replica/Hora),
              data = data)

emm <- emmeans(modelo, ~ Condicion | Hora)
pairs(emm, adjust = "tukey") #p= 0.3826


### 4.2. RNAseq ####
samples_MRC5_toci <- read_ods("Samples_MRC5.ods", sheet = 2)

samples_MRC5_toci$condition <- as.factor(samples_MRC5_toci$condition)

### 4.2.1. Matriz de conteo ####
ah_MRC5_toci <- AnnotationHub()
edb_MRC5_toci <- ah_MRC5_toci[["AH73986"]]

txi_MRC5_toci<- readRDS("transcripts_MRC5_toci_github.rds") 
View(txi_MRC5_toci$counts)


### 4.2.2. DESeq2 ####
#BiocManager::install("DESeq2")
ddsTxi_MRC5_toci <- DESeqDataSetFromTximport(txi_MRC5_toci, colData=samples_MRC5_toci, design= ~condition)
dim(ddsTxi_MRC5_toci) #34528

keep_MRC5_toci<-rowSums(counts(ddsTxi_MRC5_toci))>=ncol(ddsTxi_MRC5_toci)
ddsTxi_MRC5_toci<-ddsTxi_MRC5_toci[keep_MRC5_toci,]
dim(ddsTxi_MRC5_toci)#13904

dds_MRC5_toci<-DESeq(ddsTxi_MRC5_toci)
genenames_MRC5_toci<-ensembldb::select(edb_MRC5_toci, keys=rownames(dds_MRC5_toci), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds_MRC5_toci)<-genenames_MRC5_toci$SYMBOL

### 4.2.3. Differential expression ####
res_MRC5_toci<-results(dds_MRC5_toci)
resOrdered_MRC5_toci <- res_MRC5_toci[order(res_MRC5_toci$log2FoldChange),] 
resOrdered_MRC5_toci
summary(res_MRC5_toci, alpha=0.05)
#MRC5+toci:218 genes con ED, 126-up. 92-down
sig.genes_MRC5_toci<-resOrdered_MRC5_toci[!is.na(resOrdered_MRC5_toci$padj) & resOrdered_MRC5_toci$padj<0.05,]

sig.genes_MRC5_toci <- as.data.frame(sig.genes_MRC5_toci)
sig.genes_MRC5_toci$genes <- rownames(sig.genes_MRC5_toci)
resOrdered_MRC5_toci <- as.data.frame(resOrdered_MRC5_toci)

resultsNames(dds_MRC5_toci) #SinCEC vs ConCEC



### 4.2.4. Volcano Plot ####
my.genes_MRC5_toci <- c("TGFB1", "IL6", "LIF", "MMP1", "MMP10", "PPIA", "BMP2",
                        "APLP1", "ITGB8", "LRP6", "ITGB5", "ITGA2", "LRP8", "SDC2", "HAVCR2", "F2R", "F3", "ERBB4", "FAS",
                        "IL1RAPL1", "LPAR3", "DDR1", "TSPAN14", "ADORA1", "TNFRSF19")
sig.genes_MRC5_toci$significance[sig.genes_MRC5_toci$log2FoldChange > 1 & sig.genes_MRC5_toci$padj < 0.05] <- "Up"
sig.genes_MRC5_toci$significance[sig.genes_MRC5_toci$log2FoldChange < -1 & sig.genes_MRC5_toci$padj < 0.05] <- "Down"

my.genes_MRC5_toci <- resOrdered_MRC5_toci[rownames(resOrdered_MRC5_toci) %in% my.genes_MRC5_toci,]

ggplot(as.data.frame(res_MRC5_toci), aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(col = padj<0.05), alpha = 0.8)+
  xlim(-11,11)+
  ylim(0,13)+
  geom_point(data = my.genes_MRC5_toci, color = "black", size = 1.8)+
  geom_label_repel(data=my.genes_MRC5_toci, aes(label=rownames(my.genes_MRC5_toci)),
                   force = 2, 
                   nudge_y = 1,
                   size=6)+
  coord_cartesian(ylim = c(0,8), xlim = c(-5, 5), expand=FALSE)+
  labs(title = "MRC5 MV vs CPAP+ Tocilizumab (10ug/mL)", 
       x = "log2 Fold Change", 
       y = "-log10(padj)") +
  scale_color_manual(
    values = c("FALSE" = "#40425c", "TRUE" = "orange", "NA" = none),
    labels = c("FALSE" = "p > 0.05", "TRUE" = "p < 0.05"),
    name = "Padj",
    na.translate = FALSE 
  )+
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")+
  theme(
    aspect.ratio = 1/1.62,
    plot.title = element_text(size = 15), #, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14))


# 5. BALF proteomic ####

expr <- readRDS("proteins_abundance_github.rds")

conditions <- factor(c(rep("VM", 8), rep("CPAP", 4)),
                     levels = c("VM", "CPAP"))
design <- model.matrix(~0 + conditions)
colnames(design) <- levels(conditions)
rownames(design) <- colnames(expr)
all(rownames(design) == colnames(expr))  #TRUE

### Log2 transform ###
expr_log <- log2(expr + 1)

### Filter zero-variance proteins ###
expr_log <- expr_log[apply(expr_log, 1, var, na.rm = TRUE) > 0, ]

### Fit model ###
library(limma)
fit <- lmFit(expr_log, design)

### Contrast ###
contrast_matrix <- makeContrasts(
  VM_vs_CPAP = VM - CPAP,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

#### Results_MV_CPAP ####
results_MV_CPAP <- topTable(fit2,
                            coef = "VM_vs_CPAP",
                            number = Inf,
                            adjust.method = "BH")

signif_MV_CPAP <- results_MV_CPAP %>%
  filter(adj.P.Val < 0.1 & abs(logFC) > 1)

nrow(signif_MV_CPAP)  #24


### 5.1. Heatmap MV vs CPAP ####
#install.packages("pheatmap")
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(circlize)
library(dplyr)

proteins_sig <- (signif_MV_CPAP$ID) 
expr_log_sig <- expr_log[proteins_sig, ] 

expr_scaled <- t(scale(t(expr_log_sig))) 
annotation_col <- data.frame(
  Condition = factor(c(rep("VM", 8), rep("CPAP", 4)))
)
rownames(annotation_col) <- colnames(expr_scaled)
condition <- annotation_col$Condition
names(condition) <- rownames(annotation_col)

#Only > 50% of proteins per sample
filter_na_by_condition <- function(x, cond, threshold = 0.5) {
  tapply(x, cond, function(vals) mean(is.na(vals))) %>% 
    { all(. <= threshold) } #TRUE 
}

genes_to_keep <- apply(expr_scaled, 1, filter_na_by_condition, cond = condition, threshold = 0.5)
table(genes_to_keep) #18

cat(sum(!genes_to_keep), "genes eliminados por exceso de NAs\n") #6 genes eliminated cause NAs excess


genes_to_keep <- as.vector(genes_to_keep)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "orangered"))
cols_keep <- c(
  "VM1", "VM3", "VM4", "VM5", "VM7", "VM8",
  "CPAP1", "CPAP2", "CPAP3", "CPAP4"
)

expr_filtered <- expr_scaled[genes_to_keep, cols_keep]

colnames(expr_filtered) <- c(
  rep("MV", 6),
  rep("CPAP", 4)
)

Condition_colors <- c(
  "CPAP" = "#343e68", 
  "VM" = "orange"  
)  

column_split <- c(rep("MV", 6), rep("CPAP", 4))  

  Heatmap(expr_filtered,
          name = "Z-score",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE, 
          column_split = column_split, 
          show_row_names = TRUE,
          show_column_names = FALSE,
          show_column_dend = FALSE,
          top_annotation = HeatmapAnnotation(
            foo = anno_block(
              gp = gpar(fill = Condition_colors, alpha = 0.6),
              labels = c("CPAP", "MV"),
              labels_gp = gpar(col = "black", fontsize = 9)
            )
          ),
          column_names_gp = gpar(fontsize = 7),
          row_names_gp = gpar(fontsize = 8),
          row_dend_side = "left",
          column_title = NULL,
          width = unit(8, "cm"),
          height = unit(20, "cm"))

## 6. Proteomic BALF ligands-Transcriptomic receptors ####
# DB de CellTalkDB
  
lr_pairs <- readRDS("human_lr_pair.rds")

BALF_ligand <- lr_pairs[lr_pairs$ligand_gene_symbol %in% results_MV_CPAP$ID,]
dim(BALF_ligand) #1003

#And as receptor genes in BEAS and MRC5
interactions_BEAS <- BALF_ligand[BALF_ligand$receptor_gene_symbol %in% rownames(sig.genes_BEAS), 
                                 c("ligand_gene_symbol", "receptor_gene_symbol")]
dim(interactions_BEAS)#47 

colnames(interactions_BEAS) <- c("from", "to")

interactions_MRC5 <- BALF_ligand[BALF_ligand$receptor_gene_symbol %in% rownames(sig.genes_MRC5), 
                                 c("ligand_gene_symbol", "receptor_gene_symbol")]
dim(interactions_MRC5)#55

colnames(interactions_MRC5) <- c("from", "to")


### 6.1. BEAS networks based on receptors ####
networks <- read.table("9606.protein.links.v11.5.txt", header = TRUE)
dim(networks)
head(networks)
length(unique(networks$protein1)) #19385

networks$protein1 <- str_replace(networks$protein1, "^9606\\.", "")
networks$protein2 <- str_replace(networks$protein2, "^9606\\.", "")

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "hgnc_symbol"),
                 filters = "ensembl_peptide_id", values = c(unique(networks$protein1), unique(networks$protein2)),
                 mart = mart)

networks$gene1 <- results$ensembl_gene_id[match(networks$protein1, results$ensembl_peptide_id)]
networks$name1 <- results$hgnc_symbol[match(networks$protein1, results$ensembl_peptide_id)]
networks$gene2 <- results$ensembl_gene_id[match(networks$protein2, results$ensembl_peptide_id)]
networks$name2 <- results$hgnc_symbol[match(networks$protein2, results$ensembl_peptide_id)]


nets_BEAS <- na.omit(networks[networks$name1 %in% interactions_BEAS$to,])
nets_BEAS <- nets_BEAS[nets_BEAS$name2 %in% rownames(sig.genes_BEAS)[which(sig.genes_BEAS$padj < 0.05)],]
nets_BEAS <- nets_BEAS[nets_BEAS$combined_score > 600, c("name1", "name2")]
dim(nets_BEAS) #21
colnames(nets_BEAS) <- c("from", "to")

position <- c(rep("ligand", length(unique(interactions_BEAS$from))), 
              rep("receptor", length(unique(interactions_BEAS$to))),
              rep("intracellular", length(unique(nets_BEAS$to))))
names(position) <- c(unique(interactions_BEAS$from), unique(interactions_BEAS$to), unique(nets_BEAS$to))
position

#### build the network ####

edges_BEAS <- rbind(interactions_BEAS, nets_BEAS)
dim(edges_BEAS) #68

nodes_BEAS <- data.frame(nodes = unique(c(unique(edges_BEAS$from), unique(edges_BEAS$to))))
nodes_BEAS$position <- position[match(nodes_BEAS$nodes, names(position))]

FC <- c(results_MV_CPAP$logFC[results_MV_CPAP$ID %in% nodes_BEAS$nodes[nodes_BEAS$position == "ligand"]],
        sig.genes_BEAS$log2FoldChange[rownames(sig.genes_BEAS) %in% nodes_BEAS$nodes[nodes_BEAS$position == "receptor"]],
        sig.genes_BEAS$log2FoldChange[rownames(sig.genes_BEAS) %in% nodes_BEAS$nodes[nodes_BEAS$position == "intracellular"]])

names(FC) <- c(results_MV_CPAP$ID[results_MV_CPAP$ID %in% nodes_BEAS$nodes[nodes_BEAS$position == "ligand"]],
               rownames(sig.genes_BEAS)[rownames(sig.genes_BEAS) %in% nodes_BEAS$nodes[nodes_BEAS$position == "receptor"]],
               rownames(sig.genes_BEAS)[rownames(sig.genes_BEAS) %in% nodes_BEAS$nodes[nodes_BEAS$position == "intracellular"]])

nodes_BEAS$FC <- FC[match(nodes_BEAS$nodes, names(FC))]

nodes_BEAS$position <- factor(nodes_BEAS$position, levels = c("ligand", "receptor", "intracellular"))

nodes_BEAS$color <- nodes_BEAS$FC
nodes_BEAS$color[nodes_BEAS$position == "ligand"] <- NA

nodes_BEAS$x <- NA
nodes_BEAS$y <- NA

routes_tidy_BEAS <- tbl_graph(nodes = nodes_BEAS, edges = edges_BEAS, directed = TRUE)

routes_tidy_BEAS <- routes_tidy_BEAS %>%
  activate(nodes) %>%
  mutate(
    x = x,
    y = y,
    color = color,
    position = position,
    name = nodes 
  )

nodes_df_BEAS <- nodes_BEAS

for(pos in unique(nodes_df_BEAS$position)){
  idx <- which(nodes_df_BEAS$position == pos)
  n <- length(idx)
  
  if(pos == "receptor"){
    spread <- 3 
  } else {
    spread <- 11  
  }
  nodes_df_BEAS$x[idx] <- seq(-spread, spread, length.out = n)
}

ic <- which(nodes_df_BEAS$position == "intracellular")
nodes_df_BEAS$y[ic] <- -0.3

rc <- which(nodes_df_BEAS$position == "receptor")
nodes_df_BEAS$y[rc] <- 0

lig <- which(nodes_df_BEAS$position == "ligand")
nodes_df_BEAS$y[lig] <- 0.3

plot(nodes_df_BEAS$x, nodes_df_BEAS$y)


ggraph(routes_tidy_BEAS, 
       layout = "manual", 
       x = nodes_df_BEAS$x, 
       y = nodes_df_BEAS$y) +
  geom_edge_diagonal(
    color = "grey30",
    mapping = NULL,
    data = get_edges(),
    position = "identity",
    arrow = arrow(length = unit(1, "mm"), type = "closed"),
    end_cap = circle(6, 'mm'),
    strength = 1,
    flipped = FALSE,
    n = 100,
    lineend = "butt", #round, butt, square
    linejoin = "round", #round, mitre, bevel
    linemitre = 1,
    label_colour = "black",
    label_alpha = 1.5,
    label_parse = FALSE,
    check_overlap = FALSE,
    angle_calc = "rot",
    force_flip = TRUE)+
  geom_node_point(aes(fill = color, shape = position, size = position)) +
  geom_node_text(aes(label = nodes), size = 2) +
  theme_graph() +
  theme(legend.position = "bottom")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  scale_shape_manual(values = c(ligand = 21, receptor = 22, intracellular = 21)) +
  scale_size_manual(values = c(ligand = 16, receptor = 19, intracellular = 15)) 


####BP downstream intracell BEAS ####
library(biomaRt)
intracell_genes_BEAS <- nodes_BEAS$nodes[nodes_BEAS$position == "intracellular" ]
genes_BEAS <- nodes_BEAS$nodes[nodes_BEAS$position %in% c("intracellular", "receptor")]
length(intracell_genes_BEAS)
length(genes_BEAS)
head(intracell_genes_BEAS)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

found_genes_BEAS <- getBM(
  attributes = "hgnc_symbol",
  filters = "hgnc_symbol",
  values = genes_BEAS,
  mart = mart
)
GO_slim_BEAS <- getBM(
  attributes = c(
    "hgnc_symbol",
    "go_id",
    "name_1006",               
    "namespace_1003",        
    "goslim_goa_accession",   
    "goslim_goa_description"   
  ),
  filters = "hgnc_symbol",
  values = found_genes_BEAS,
  mart = mart
)

GO_slim_BP_BEAS <- GO_slim_BEAS[GO_slim_BEAS$namespace_1003 == "biological_process" &
                                  !is.na(GO_slim_BEAS$goslim_goa_description), ]

BP_edges_BEAS <- data.frame(
  from = GO_slim_BP_BEAS$hgnc_symbol,
  to = GO_slim_BP_BEAS$goslim_goa_description
)
BP_counts_BEAS <- table(BP_edges_BEAS$to)
BP_keep_BEAS <- names(BP_counts_BEAS)[BP_counts_BEAS >= 3]
BP_edges_filtered_BEAS <- BP_edges_BEAS[BP_edges_BEAS$to %in% BP_keep_BEAS, ]


gene_fc_BEAS <- nodes_BEAS %>%
  filter(position == "intracellular") %>%
  select(nodes, FC)
gene_fc_BEAS <- nodes_BEAS %>%
  filter(position %in% c("intracellular", "receptor")) %>%
  select(nodes, FC)


GO_gene_BP_BEAS <- BP_edges_filtered_BEAS %>%
  rename(
    gene = from,
    GO = to
  ) %>%
  distinct()


plot_df_BEAS <- GO_gene_BP_BEAS %>%
  left_join(gene_fc_BEAS, by = c("gene" = "nodes")) %>%
  filter(!is.na(FC))

  ggplot(plot_df_BEAS,aes(x = gene,y = GO,fill = FC)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Log2FC"
  ) +
  labs(
    x = "Genes (downstream)",
    y = "Biological Process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)#,
    #panel.grid = element_blank()
  )


### 6.2. BEAS + toci networks based on receptors ####

ligands_used_BEAS   <- nodes_BEAS$nodes[nodes_BEAS$position == "ligand"]
receptors_used_BEAS <- nodes_BEAS$nodes[nodes_BEAS$position == "receptor"]
intracell_used_BEAS <- nodes_BEAS$nodes[nodes_BEAS$position == "intracellular"]
edges_used_BEAS <- edges_BEAS

res_BEAS_toci <- as.data.frame(res_BEAS_toci)
res_BEAS_toci <- mutate(res_BEAS_toci, gen = rownames(res_BEAS_toci))

res_BEAS_toci_filtered <- res_BEAS_toci %>%
  filter(gen %in% c(ligands_used_BEAS,
                    receptors_used_BEAS,
                    intracell_used_BEAS))

nodes_BEAS_toci <- data.frame(
  nodes = c(ligands_used_BEAS,
            receptors_used_BEAS,
            intracell_used_BEAS),
  position = c(
    rep("ligand", length(ligands_used_BEAS)),
    rep("receptor", length(receptors_used_BEAS)),
    rep("intracellular", length(intracell_used_BEAS))
  )
)

nodes_BEAS_toci <- left_join(nodes_BEAS_toci, res_BEAS_toci_filtered, 
                             by = c("nodes" = "gen"))

nodes_BEAS_toci$color <- nodes_BEAS_toci$log2FoldChange

edges_BEAS_toci <- edges_used_BEAS

routes_BEAS_toci <- tbl_graph(nodes = nodes_BEAS_toci, edges = edges_BEAS_toci, directed = TRUE)

nodes_BEAS_toci$x <- nodes_BEAS$x 
nodes_BEAS_toci$y <- nodes_BEAS$y 

length(nodes_BEAS$nodes) #BEAS 63
length(nodes_BEAS_toci$nodes) #BEAS + toci 63


for(pos in unique(nodes_BEAS_toci$position)){
  idx <- which(nodes_BEAS_toci$position == pos)
  n <- length(idx)
  
  if(pos == "receptor"){
    spread <- 3  
  } else {
    spread <- 11 
  }
  nodes_BEAS_toci$x[idx] <- seq(-spread, spread, length.out = n)
}

nodes_BEAS_toci$color[nodes_BEAS_toci$position == "ligand"] <- NA

routes_BEAS_toci <- tbl_graph(nodes = nodes_BEAS_toci, edges = edges_BEAS_toci, directed = TRUE)

routes_BEAS_toci <- routes_BEAS_toci %>%
  activate(nodes) %>%   # activar la tabla de nodos dentro del tbl_graph
  mutate(
    x = x,
    y = y,
    color = color,
    position = position,
    name = nodes  
  )

routes_BEAS_toci %>% activate(nodes) %>% as_tibble() %>% head()


ggraph(routes_BEAS_toci, 
       layout = "manual", 
       x = nodes_df_BEAS$x, 
       y = nodes_df_BEAS$y) +
  geom_edge_diagonal(
    color = "grey30",
    mapping = NULL,
    data = get_edges(),
    position = "identity",
    arrow = arrow(length = unit(1, "mm"), type = "closed"),
    end_cap = circle(6, 'mm'),
    strength = 1,
    flipped = FALSE,
    n = 100,
    lineend = "butt", 
    linejoin = "round", 
    linemitre = 1,
    label_colour = "black",
    label_alpha = 1.5,
    label_parse = FALSE,
    check_overlap = FALSE,
    angle_calc = "rot",
    force_flip = TRUE)+
  geom_node_point(aes(fill = color, shape = position, size = position)) +
  geom_node_text(aes(label = nodes), size = 2) +
  theme_graph() +
  theme(legend.position = "bottom")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  scale_shape_manual(values = c(ligand = 21, receptor = 22, intracellular = 21)) +
  scale_size_manual(values = c(ligand = 16, receptor = 19, intracellular = 15)) 

### 6.3. MRC5 networks based on receptors ####

networks <- read.table("9606.protein.links.v11.5.txt", header = TRUE)
dim(networks)
head(networks)
length(unique(networks$protein1)) #19385

networks$protein1 <- str_replace(networks$protein1, "^9606\\.", "")
networks$protein2 <- str_replace(networks$protein2, "^9606\\.", "")

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "hgnc_symbol"),
                 filters = "ensembl_peptide_id", values = c(unique(networks$protein1), unique(networks$protein2)),
                 mart = mart)

networks$gene1 <- results$ensembl_gene_id[match(networks$protein1, results$ensembl_peptide_id)]
networks$name1 <- results$hgnc_symbol[match(networks$protein1, results$ensembl_peptide_id)]
networks$gene2 <- results$ensembl_gene_id[match(networks$protein2, results$ensembl_peptide_id)]
networks$name2 <- results$hgnc_symbol[match(networks$protein2, results$ensembl_peptide_id)]


nets_MRC5 <- na.omit(networks[networks$name1 %in% interactions_MRC5$to,])
nets_MRC5 <- nets_MRC5[nets_MRC5$name2 %in% rownames(sig.genes_MRC5)[which(sig.genes_MRC5$padj < 0.05)],]
nets_MRC5<- nets_MRC5[nets_MRC5$combined_score > 600, c("name1", "name2")]
dim(nets_MRC5) #92
colnames(nets_MRC5) <- c("from", "to")

position <- c(rep("ligand", length(unique(interactions_MRC5$from))), 
              rep("receptor", length(unique(interactions_MRC5$to))),
              rep("intracellular", length(unique(nets_MRC5$to))))
names(position) <- c(unique(interactions_MRC5$from), unique(interactions_MRC5$to), unique(nets_MRC5$to))
position

#### build the network ####

edges_MRC5 <- rbind(interactions_MRC5, nets_MRC5)
dim(edges_MRC5) #147

nodes_MRC5 <- data.frame(nodes = unique(c(unique(edges_MRC5$from), unique(edges_MRC5$to))))
nodes_MRC5$position <- position[match(nodes_MRC5$nodes, names(position))]

FC <- c(results_MV_CPAP$logFC[results_MV_CPAP$ID %in% nodes_MRC5$nodes[nodes_MRC5$position == "ligand"]],
        sig.genes_MRC5$log2FoldChange[rownames(sig.genes_MRC5) %in% nodes_MRC5$nodes[nodes_MRC5$position == "receptor"]],
        sig.genes_MRC5$log2FoldChange[rownames(sig.genes_MRC5) %in% nodes_MRC5$nodes[nodes_MRC5$position == "intracellular"]])

names(FC) <- c(results_MV_CPAP$ID[results_MV_CPAP$ID %in% nodes_MRC5$nodes[nodes_MRC5$position == "ligand"]],
               rownames(sig.genes_MRC5)[rownames(sig.genes_MRC5) %in% nodes_MRC5$nodes[nodes_MRC5$position == "receptor"]],
               rownames(sig.genes_MRC5)[rownames(sig.genes_MRC5) %in% nodes_MRC5$nodes[nodes_MRC5$position == "intracellular"]])

nodes_MRC5$FC <- FC[match(nodes_MRC5$nodes, names(FC))]

nodes_MRC5$position <- factor(nodes_MRC5$position, levels = c("ligand", "receptor", "intracellular"))

nodes_MRC5$color <- nodes_MRC5$FC
nodes_MRC5$color[nodes_MRC5$position == "ligand"] <- NA

nodes_MRC5$x <- NA
nodes_MRC5$y <- NA

routes_tidy_MRC5 <- tbl_graph(nodes = nodes_MRC5, edges = edges_MRC5, directed = TRUE)

routes_tidy_MRC5 <- routes_tidy_MRC5 %>%
  activate(nodes) %>%
  mutate(
    x = x,
    y = y,
    color = color,
    position = position,
    name = nodes 
  )

nodes_df_MRC5 <- nodes_MRC5

for(pos in unique(nodes_df_MRC5$position)){
  idx <- which(nodes_df_MRC5$position == pos)
  n <- length(idx)
  
  if(pos == "receptor"){
    spread <- 6 
  } else {
    spread <- 11  
  }
  nodes_df_MRC5$x[idx] <- seq(-spread, spread, length.out = n)
}

ic <- which(nodes_df_MRC5$position == "intracellular")
nodes_df_MRC5$y[ic] <- -0.3

rc <- which(nodes_df_MRC5$position == "receptor")
nodes_df_MRC5$y[rc] <- 0

lig <- which(nodes_df_MRC5$position == "ligand")
nodes_df_MRC5$y[lig] <- 0.3

plot(nodes_df_MRC5$x, nodes_df_MRC5$y)


ggraph(routes_tidy_MRC5, 
       layout = "manual", 
       x = nodes_df_MRC5$x, 
       y = nodes_df_MRC5$y) +
  geom_edge_diagonal(
    color = "grey30",
    mapping = NULL,
    data = get_edges(),
    position = "identity",
    arrow = arrow(length = unit(1, "mm"), type = "closed"),
    end_cap = circle(6, 'mm'),
    strength = 1,
    flipped = FALSE,
    n = 100,
    lineend = "butt", 
    linejoin = "round",
    linemitre = 1,
    label_colour = "black",
    label_alpha = 1.5,
    label_parse = FALSE,
    check_overlap = FALSE,
    angle_calc = "rot",
    force_flip = TRUE)+
  geom_node_point(aes(fill = color, shape = position, size = position)) +
  geom_node_text(aes(label = nodes), size = 2) +
  theme_graph() +
  theme(legend.position = "bottom")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  scale_shape_manual(values = c(ligand = 21, receptor = 22, intracellular = 21)) +
  scale_size_manual(values = c(ligand = 16, receptor = 19, intracellular = 15)) 


####BP downstream intracell MRC5 ####
library(biomaRt)
intracell_genes_MRC5 <- nodes_MRC5$nodes[nodes_MRC5$position == "intracellular"]
genes_MRC5 <- nodes_MRC5$nodes[nodes_MRC5$position %in% c("intracellular", "receptor")]
length(genes_MRC5)
length(intracell_genes_MRC5)
head(intracell_genes_MRC5)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

found_genes_MRC5 <- getBM(
  attributes = "hgnc_symbol",
  filters = "hgnc_symbol",
  values = genes_MRC5,
  mart = mart
)
GO_slim_MRC5 <- getBM(
  attributes = c(
    "hgnc_symbol",
    "go_id",
    "name_1006",               
    "namespace_1003",         
    "goslim_goa_accession",    
    "goslim_goa_description"   
  ),
  filters = "hgnc_symbol",
  values = found_genes_MRC5,
  mart = mart
)

GO_slim_BP_MRC5 <- GO_slim_MRC5[GO_slim_MRC5$namespace_1003 == "biological_process" &
                                  !is.na(GO_slim_MRC5$goslim_goa_description), ]

BP_edges_MRC5 <- data.frame(
  from = GO_slim_BP_MRC5$hgnc_symbol,
  to = GO_slim_BP_MRC5$goslim_goa_description
)
BP_counts_MRC5 <- table(BP_edges_MRC5$to)
BP_keep_MRC5 <- names(BP_counts_MRC5)[BP_counts_MRC5 >= 3]
BP_edges_filtered_MRC5 <- BP_edges_MRC5[BP_edges_MRC5$to %in% BP_keep_MRC5, ]


gene_fc_MRC5 <- nodes_MRC5 %>%
  filter(position == "intracellular") %>%
  select(nodes, FC)

gene_fc_MRC5 <- nodes_MRC5 %>%
  filter(position %in% c("intracellular", "receptor")) %>%
  select(nodes, FC)

GO_gene_BP_MRC5 <- BP_edges_filtered_MRC5 %>%
  rename(
    gene = from,
    GO = to
  ) %>%
  distinct()


plot_df_MRC5 <- GO_gene_BP_MRC5 %>%
  left_join(gene_fc_MRC5, by = c("gene" = "nodes")) %>%
  filter(!is.na(FC))


  ggplot(plot_df_MRC5,aes(x = gene,y = GO,fill = FC)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Log2FC"
  ) +
  labs(
    x = "Genes (downstream)",
    y = "Biological Process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



### 6.4. MRC5 + toci networks based on receptors ####

ligands_used_MRC5 <- nodes_MRC5$nodes[nodes_MRC5$position == "ligand"]
receptors_used_MRC5 <- nodes_MRC5$nodes[nodes_MRC5$position == "receptor"]
intracell_used_MRC5 <- nodes_MRC5$nodes[nodes_MRC5$position == "intracellular"]
edges_used_MRC5 <- edges_MRC5

res_MRC5_toci <- as.data.frame(res_MRC5_toci)
res_MRC5_toci <- mutate(res_MRC5_toci, gen = rownames(res_MRC5_toci))

res_MRC5_toci_filtered <- res_MRC5_toci %>%
  filter(gen %in% c(ligands_used_MRC5,
                    receptors_used_MRC5,
                    intracell_used_MRC5))

nodes_MRC5_toci <- data.frame(
  nodes = c(ligands_used_MRC5,
            receptors_used_MRC5,
            intracell_used_MRC5),
  position = c(
    rep("ligand", length(ligands_used_MRC5)),
    rep("receptor", length(receptors_used_MRC5)),
    rep("intracellular", length(intracell_used_MRC5))
  )
)

nodes_MRC5_toci <- left_join(nodes_MRC5_toci, res_MRC5_toci_filtered, 
                             by = c("nodes" = "gen"))

nodes_MRC5_toci$color <- nodes_MRC5_toci$log2FoldChange

edges_MRC5_toci <- edges_used_MRC5

routes_MRC5_toci <- tbl_graph(nodes = nodes_MRC5_toci, edges = edges_MRC5_toci, directed = TRUE)

nodes_MRC5_toci$x <- nodes_MRC5$x 
nodes_MRC5_toci$y <- nodes_MRC5$y 

length(nodes_MRC5$nodes) #MRC5 99
length(nodes_MRC5_toci$nodes) #MRC5 + toci 99


for(pos in unique(nodes_MRC5_toci$position)){
  idx <- which(nodes_MRC5_toci$position == pos)
  n <- length(idx)
  
  if(pos == "receptor"){
    spread <- 3  
  } else {
    spread <- 11  
  }
  nodes_MRC5_toci$x[idx] <- seq(-spread, spread, length.out = n)
}

nodes_MRC5_toci$color[nodes_MRC5_toci$position == "ligand"] <- NA

routes_MRC5_toci <- tbl_graph(nodes = nodes_MRC5_toci, edges = edges_MRC5_toci, directed = TRUE)

routes_MRC5_toci <- routes_MRC5_toci %>%
  activate(nodes) %>%  
  mutate(
    x = x,
    y = y,
    color = color,
    position = position,
    name = nodes  
  )

routes_MRC5_toci %>% activate(nodes) %>% as_tibble() %>% head()


ggraph(routes_MRC5_toci, 
       layout = "manual", 
       x = nodes_df_MRC5$x, 
       y = nodes_df_MRC5$y) +
  geom_edge_diagonal(
    color = "grey30",
    mapping = NULL,
    data = get_edges(),
    position = "identity",
    arrow = arrow(length = unit(1, "mm"), type = "closed"),
    end_cap = circle(6, 'mm'),
    strength = 1,
    flipped = FALSE,
    n = 100,
    lineend = "butt", #round, butt, square
    linejoin = "round", #round, mitre, bevel
    linemitre = 1,
    label_colour = "black",
    label_alpha = 1.5,
    label_parse = FALSE,
    check_overlap = FALSE,
    angle_calc = "rot",
    force_flip = TRUE)+
  geom_node_point(aes(fill = color, shape = position, size = position)) +
  geom_node_text(aes(label = nodes), size = 2) +
  theme_graph() +
  theme(legend.position = "bottom")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  scale_shape_manual(values = c(ligand = 21, receptor = 22, intracellular = 21)) +
  scale_size_manual(values = c(ligand = 16, receptor = 19, intracellular = 15)) 

# 7. Delayed epithelium wound closure during tidal ventilation depends on PPIA####
### 7.1. Interaction analysis ####
samples_BEAS_toci 
samples_BEAS 
samples_BEAS$treatment <- "Without_Toci"
samples_BEAS$ID <- c("CRI1", "CRI2", "CRI3", "CRI4", "CRI5", "CRI6")
samples_BEAS <- samples_BEAS [ ,c("sample_name", "ID", "patient", "condition", "cell_line", "treatment")]
samples_BEAS <-rename(samples_BEAS, Sample_name = sample_name)

samples_interaction <- rbind(samples_BEAS_toci, samples_BEAS)
samples_interaction$patient <- c("VM1t", "VM2t", "VM3t", "CPAP1t", "CPAP2t", "CPAP3t", "VM1", "VM2", "VM3", "CPAP1", "CPAP2", "CPAP3")
samples_interaction$condition[samples_interaction$condition == "Sin CEC"] <- "VM"
samples_interaction$condition[samples_interaction$condition == "Con CEC"] <- "CPAP"
samples_interaction <-rename(samples_interaction, ventilation = condition)
samples_interaction <-rename(samples_interaction, condition = treatment)
samples_interaction$ventilation <- factor(samples_interaction$ventilation,
                                          levels = c("CPAP", "VM"))

samples_interaction$condition <- factor(samples_interaction$condition,
                                        levels = c("Toci", "Without_Toci"))
### 7.1.1. Count matrix ####

ah_interaction <- AnnotationHub()
edb_interaction <- ah_interaction[["AH73986"]]

txi_interaction<- readRDS("transcripts_interaction_analysis_github.rds") 
View(txi_interaction$counts)

### 7.1.2. DESeq2 ####
#BiocManager::install("DESeq2")
ddsTxi_interaction <- DESeqDataSetFromTximport(txi_interaction, colData=samples_interaction, design= ~ventilation + condition + ventilation:condition)
#ventilation: effect of MV vs CPAP; condition: effect of toci vs without toci; ventilation:condiiton: if the ventilation effect depends of the treatment
keep_interaction<-rowSums(counts(ddsTxi_interaction))>=ncol(ddsTxi_interaction)
ddsTxi_interaction<-ddsTxi_interaction[keep_interaction,]

dds_interaction <-DESeq(ddsTxi_interaction)
resultsNames(dds_interaction)
#[1] "Intercept"                           "ventilation_VM_vs_CPAP"              "condition_Without_Toci_vs_Toci"     
#[4] "ventilationVM.conditionWithout_Toci"


resVentilation <- results(dds_interaction, name="ventilation_VM_vs_CPAP")
#Results of the principal effect of the ventilation: MV vs CPAP

resInteraction <- results(dds_interaction, name="ventilationVM.conditionWithout_Toci")
#Results of the interaction between ventilation and treatment(condition), that is if the effect of MV vs CPAP changes with the treatment

sigVentilation <- subset(resVentilation, padj < 0.05)
#Significant genes in ventilation MV vs CPAP with tocilizumab

nonsigInteraction <- subset(resInteraction, padj >= 0.05)
# No significant genes in the interacción, that mean the effect of MV vs CPAP doesn´t change with the treatment

genesExclusiveVentilation <- intersect(rownames(sigVentilation), rownames(nonsigInteraction))
exclusiveVentilation <- sigVentilation[genesExclusiveVentilation, ]
#These genes show a consistent effect of MV vs CPAP, independent of treatment

genes_interaction <- rownames(exclusiveVentilation)

symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes_interaction,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rownames(exclusiveVentilation) <- symbols %>% View()

#"BCL3"    "OAS2"    "RSAD1"   "USF3"    "DCTPP1"  "H2BC4"   "BLOC1S5" "PPIA"    "PPT2"
#set of genes whose expression changes with the type of ventilation (MV vs CPAP), and this effect is 
#robust and independent of Tocilizumab use

# 8. qPCR BEAS + BALF MV vs CPAP with CsA ####
### 8.1. qPCR hPPIA ####

Doc <- read_excel("qPCR_BEAS_CsA.xls", sheet = 1, range = "B8:J89")

Datos_utiles <- Doc[c(1:2,9)]

Datos_utiles$"Cт" <- as.numeric(Datos_utiles$"Cт", na.rm = TRUE)

Datos_utiles <- rename(Datos_utiles, ciclos = "Cт", dilucion = "Sample Name", gen = "Target Name")
Datos_utiles <- na.omit(Datos_utiles)
Datos_utiles$ciclos<- as.numeric(Datos_utiles$ciclos) 

Datos_utiles <-Datos_utiles %>% group_by (gen, dilucion) %>% mutate(mean = mean(ciclos), sd = sd(ciclos))

a = Datos_utiles$ciclos
b = Datos_utiles$mean
Datos_utiles<- cbind(Datos_utiles, dif =abs(a-b))

data <- rbind(Datos_utiles %>% 
                filter(sd < 0.5) %>% 
                summarise(mean = mean(ciclos), sd = sd(ciclos)),
              Datos_utiles %>% 
                filter(sd >=0.5) %>% 
                group_by(dilucion, gen) %>%
                filter(dif < max(dif)) %>% 
                summarise(mean = mean(ciclos), sd = sd(ciclos)))
view(data)

values_hGAPDH <- subset(data, gen == "hGAPDH")
View(values_hGAPDH)

values_hPPIA <- data.frame(filter(data, gen == "hPPIA"))
view(values_hPPIA)

tabla2 <- merge(values_hPPIA, values_hGAPDH, by = "dilucion", suffixes = c("_hPPIA", "_hGAPDH"))
tabla2 <-tabla2 %>%
  mutate(dCt = mean_hPPIA - mean_hGAPDH) %>% 
  mutate(dCt2 = 2^(-dCt))

filter1<-tabla2[tabla2$dilucion=="1 CON CEC+CICLOSPORINA 0.4 uM",]
filter2<-tabla2[tabla2$dilucion=="2 CON CEC+CICLOSPORINA 0.4 uM",]
filter3<-tabla2[tabla2$dilucion=="3 CON CEC+CICLOSPORINA 0.4 uM",]
mean1<-mean(filter1$dCt2)
mean2<-mean(filter2$dCt2)
mean3<-mean(filter3$dCt2)


tabla2$avbasal<-ifelse(tabla2$dilucion %in% c("1 CON CEC+CICLOSPORINA 0.4 uM", "2 CON CEC+CICLOSPORINA 0.4 uM", "3 CON CEC+CICLOSPORINA 0.4 uM"), mean(c(mean1, mean2, mean3)),mean(c(mean1, mean2, mean3)))

tabla2<-tabla2 %>% 
  mutate(avbasal2dct=dCt2/avbasal)

library(ggsignif)
tabla2$condicion <- c("CPAP + CsA", "MV + CsA", "CPAP + CsA", "MV + CsA", "CPAP + CsA", "MV + CsA", "Control", "Control" )
library(ggpubr)
  tabla2 %>% 
  filter(condicion != "Control") %>% 
  ggplot(aes(x = condicion, y = avbasal2dct, color=condicion))+
  geom_boxplot(aes(fill = condicion), alpha = 0.3, width = 0.2)+
  geom_point(size = 2)+
  theme_minimal() +
  labs(y="qPCR relative expression", x="Condition", title = "PPIA in BEAS + BALF with CsA 0.4 uM")+
  scale_color_manual(values = c("CPAP + CsA" = "#2a3d7f", "MV + CsA" = "orange"), labels = c("CPAP", "MV"), name = "Condition")+
  scale_fill_manual(values = c("CPAP + CsA" = "#2a3d7f", "MV + CsA" = "orange"), labels = c("CPAP", "MV"), name = "Condition")+
  #scale_x_discrete(labels = c("Est" = "Static", "Flex" = "Stretch")) +
  theme(
    aspect.ratio = 1/1.32,
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.line = element_line(color = "grey50", linewidth = 1),
    plot.title = element_text(size = 17), #, face = "bold"),
    legend.position = "right"
  )+
  stat_compare_means(method = "t.test", comparisons = list(c("CPAP + CsA", "MV + CsA"))) #p=0.26


### 8.2. qPCR hIL6 ####
Doc <- read_excel("qPCR_BEAS_CsA.xls", sheet = 1, range = "B8:J89")

Datos_utiles <- Doc[c(1:2,9)]

Datos_utiles$"Cт" <- as.numeric(Datos_utiles$"Cт", na.rm = TRUE)

Datos_utiles <- rename(Datos_utiles, ciclos = "Cт", dilucion = "Sample Name", gen = "Target Name")
Datos_utiles <- na.omit(Datos_utiles)
Datos_utiles$ciclos<- as.numeric(Datos_utiles$ciclos) 

Datos_utiles <-Datos_utiles %>% group_by (gen, dilucion) %>% mutate(mean = mean(ciclos), sd = sd(ciclos))

a = Datos_utiles$ciclos
b = Datos_utiles$mean
Datos_utiles<- cbind(Datos_utiles, dif =abs(a-b))

data <- rbind(Datos_utiles %>% 
                filter(sd < 0.5) %>% 
                summarise(mean = mean(ciclos), sd = sd(ciclos)),
              Datos_utiles %>% 
                filter(sd >=0.5) %>% 
                group_by(dilucion, gen) %>%
                filter(dif < max(dif)) %>% 
                summarise(mean = mean(ciclos), sd = sd(ciclos)))
view(data)

values_hGAPDH <- subset(data, gen == "hGAPDH")
View(values_hGAPDH)

values_hIL6 <- data.frame(filter(data, gen == "hIL6"))
view(values_hIL6)

tabla2 <- merge(values_hIL6, values_hGAPDH, by = "dilucion", suffixes = c("_hIL6", "_hGAPDH"))
tabla2 <-tabla2 %>%
  mutate(dCt = mean_hIL6 - mean_hGAPDH) %>% 
  mutate(dCt2 = 2^(-dCt))

filter1<-tabla2[tabla2$dilucion=="1 CON CEC+CICLOSPORINA 0.4 uM",]
filter2<-tabla2[tabla2$dilucion=="2 CON CEC+CICLOSPORINA 0.4 uM",]
filter3<-tabla2[tabla2$dilucion=="3 CON CEC+CICLOSPORINA 0.4 uM",]
mean1<-mean(filter1$dCt2)
mean2<-mean(filter2$dCt2)
mean3<-mean(filter3$dCt2)


tabla2$avbasal<-ifelse(tabla2$dilucion %in% c("1 CON CEC+CICLOSPORINA 0.4 uM", "2 CON CEC+CICLOSPORINA 0.4 uM", "3 CON CEC+CICLOSPORINA 0.4 uM"), mean(c(mean1, mean2, mean3)),mean(c(mean1, mean2, mean3)))

tabla2<-tabla2 %>% 
  mutate(avbasal2dct=dCt2/avbasal)

library(ggsignif)
tabla2$condicion <- c("CPAP + CsA", "MV + CsA", "CPAP + CsA", "MV + CsA", "CPAP + CsA", "MV + CsA", "Control", "Control" )

  tabla2 %>% 
  filter(condicion != "Control") %>% 
  ggplot(aes(x = condicion, y = avbasal2dct, color=condicion))+
  geom_boxplot(aes(fill = condicion), alpha = 0.3, width = 0.2)+
  geom_point(size = 2)+
  theme_minimal() +
  labs(y="qPCR relative expression", x="Condition", title = "IL6 in BEAS + BALF with CsA 0.4 uM")+
  scale_color_manual(values = c("CPAP + CsA" = "#2a3d7f", "MV + CsA" = "orange"), labels = c("CPAP", "MV"), name = "Condition")+
  scale_fill_manual(values = c("CPAP + CsA" = "#2a3d7f", "MV + CsA" = "orange"), labels = c("CPAP", "MV"), name = "Condition")+
  #scale_x_discrete(labels = c("Est" = "Static", "Flex" = "Stretch")) +
  theme(
    aspect.ratio = 1/1.32,
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.line = element_line(color = "grey50", linewidth = 1),
    plot.title = element_text(size = 17), #, face = "bold"),
    legend.position = "right"
  )+
  stat_compare_means(method = "t.test", comparisons = list(c("CPAP + CsA", "MV + CsA"))) #p=0.37

