############################################################################
##############    Male vs female SPLiT-seq all sequencing     ##############
############################################################################

# Libraries ----
library(plyr)         # CRAN v1.8.9
library(dplyr)        # CRAN v1.1.4
library(Seurat)       # CRAN v5.1.0
library(ggplot2)      # CRAN v3.5.1
library(readxl)       # CRAN v1.4.3
library(stringr)      # CRAN v1.5.1
library(readr)        # CRAN v2.1.5
library(SeuratDisk)   # [github::mojaveazure/seurat-disk] v0.0.0.9021
library(DropletUtils) # Bioconductor v1.24.0
library(cowplot)      # CRAN v1.1.3
library(MetBrewer)    # CRAN v0.2.0
library(ape)          # CRAN v5.8
library(gt)           # CRAN v0.11.0

setwd('~/Desktop/data/Split-seq/All_analysis/')

load('./sc_Agamb_new.RData')

# Import matrix ----
L77_2_L88_1_G50 <- read.table(file = "./matrice/L77_2_L88_1_final_G50_MQ0_matrix.txt.gz", header =TRUE, row.names = 1, colClasses =c("character", rep("numeric")))
L77_3_L88_2_G50 <- read.table(file = "./matrice/L77_3_L88_2_final_G50_MQ0_matrix.txt.gz", header =TRUE, row.names = 1, colClasses =c("character", rep("numeric")))
L77_5_L88_3_G50 <- read.table(file = "./matrice/L77_5_L88_3_final_G50_MQ0_matrix.txt.gz", header =TRUE, row.names = 1, colClasses =c("character", rep("numeric")))
L77_6_L88_4_G50 <- read.table(file = "./matrice/L77_6_L88_4_final_G50_MQ0_matrix.txt.gz", header =TRUE, row.names = 1, colClasses =c("character", rep("numeric")))

#Create Seurat objects
#Minimum cutoff (1 gene, 1 cell)
sL77_2_L88_1_G50 <- CreateSeuratObject(counts = L77_2_L88_1_G50, project = "Agamb", min.cells = 1, min.features = 50)
sL77_3_L88_2_G50 <- CreateSeuratObject(counts = L77_3_L88_2_G50, project = "Agamb", min.cells = 1, min.features = 50)
sL77_5_L88_3_G50 <- CreateSeuratObject(counts = L77_5_L88_3_G50, project = "Agamb", min.cells = 1, min.features = 50)
sL77_6_L88_4_G50 <- CreateSeuratObject(counts = L77_6_L88_4_G50, project = "Agamb", min.cells = 1, min.features = 50)

#Assign identity to each library:
sL77_2_L88_1_G50 [["library"]] <- "L77_2_L88_1"
sL77_2_L88_1_G50@meta.data$library <- "L77_2_L88_1"
sL77_3_L88_2_G50 [["library"]] <- "L77_3_L88_2"
sL77_3_L88_2_G50@meta.data$library <- "L77_3_L88_2"
sL77_5_L88_3_G50 [["library"]] <- "L77_5_L88_3"
sL77_5_L88_3_G50@meta.data$library <- "L77_5_L88_3"
sL77_6_L88_4_G50 [["library"]] <- "L77_6_L88_4"
sL77_6_L88_4_G50@meta.data$library <- "L77_6_L88_4"

#Adding ID to the cells and create final Seurat objects:
scAngamb <- merge(sL77_2_L88_1_G50, y = c(sL77_3_L88_2_G50, sL77_5_L88_3_G50, sL77_6_L88_4_G50),
                 add.cell.ids = c("sL77_2_L88_1","sL77_3_L88_2","sL77_5_L88_3","sL77_6_L88_4"),
                 project = "scAngamb",
                 merge.data = F)
scAngamb=JoinLayers(scAngamb)

#Convert into H5ad & 10xCounts:
SaveH5Seurat(scAngamb, filename = 'scAngamb.h5Seurat')
Convert('scAngamb.h5Seurat', dest ='h5ad')
write10xCounts(x = scAngamb@assays$RNA$counts, version="3", path = "./Matrix_scAngamb")
write.csv(scAngamb@meta.data, 'metadata_scAngamb.csv')


# Get BC for conditions ----
BC = read_excel("~/Desktop/data/Split-seq/All_analysis/BC_plate_male_vs_female.xlsx")

# Assign condition to read ----
sex=list()
replicate=list()
sex_rep=list()
for (i in 1:length(names(scAngamb@active.ident))){
  a=which(lapply(as.vector(BC$Barcode), function(x) grepl(x,unlist(str_split(names(scAngamb@active.ident)[i], "_"))[5]))==T)
  s=BC$Sex[a]
  r=BC$Replicate[a]
  sex = append(sex, s)
  replicate = append(replicate, r)
  sex_rep = append(sex_rep, paste0(s, '_', r))
}
sexe=unlist(sex)
replicate=unlist(replicate)
sexe_rep=unlist(sex_rep)
scAngamb=AddMetaData(scAngamb, sexe, col.name = "Sex")
scAngamb=AddMetaData(scAngamb, replicate, col.name = "Replicate")
scAngamb=AddMetaData(scAngamb, sexe_rep, col.name = "Sex replicate")
remove(a, i, r, s, sexe, sexe_rep, replicate, sex_rep, sex)

###########################
### SAVE R.DATA HERE! ----
###########################

### Normalisation and clustering ----

# show count and features by proportion of reads
VlnPlot(scAngamb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# get list genes used
all.genes <- rownames(scAngamb)
length(all.genes) # 9239 genes
# Normalisation
scAngamb.n <- NormalizeData(scAngamb, normalization.method = "LogNormalize", scale.factor = 10000)
# Find genes of interest
scAngamb.n <- FindVariableFeatures(scAngamb.n, selection.method = "vst", nfeatures = 8000)
VariableFeaturePlot(scAngamb.n)
# Scaling data
scAngamb.n <- ScaleData(scAngamb.n, features = all.genes)
# PCA
scAngamb.n <- RunPCA(scAngamb.n, features = VariableFeatures(object = scAngamb.n))
DimPlot(scAngamb.n, reduction = "pca")
# Cluster the dataset
scAngamb.c <- JackStraw(scAngamb.n, num.replicate = 100, dims = 50) #~4h
scAngamb.c <- ScoreJackStraw(scAngamb.c, dims = 1:50)
JackStrawPlot(scAngamb.c, dims = 1:50)
ElbowPlot(scAngamb.c, ndims=50)













# Normalisation
scAngamb.nlow <- NormalizeData(scAngamb, normalization.method = "LogNormalize", scale.factor = 10000)
# Find genes of interest
scAngamb.nlow <- FindVariableFeatures(scAngamb.nlow, selection.method = "vst", nfeatures = 6000)
VariableFeaturePlot(scAngamb.nlow)
# Scaling data
scAngamb.nlow <- ScaleData(scAngamb.nlow, features = all.genes)
# PCA
scAngamb.nlow <- RunPCA(scAngamb.nlow, features = VariableFeatures(object = scAngamb.nlow))
DimPlot(scAngamb.nlow, reduction = "pca")
# Cluster the dataset
scAngamb.clow <- JackStraw(scAngamb.nlow, num.replicate = 50, dims = 50) #~4h
scAngamb.clow <- ScoreJackStraw(scAngamb.clow, dims = 1:50)
JackStrawPlot(scAngamb.clow, dims = 1:50)
ElbowPlot(scAngamb.clow, ndims=50)
#k = 11 resolution 0.8
scAngamb.c_k11.0.8_low <- FindNeighbors(scAngamb.clow, dims = 1:50, k.param = 12)
scAngamb.c_k11.0.8_low <- FindClusters(scAngamb.c_k11.0.8_low, resolution = 0.7)
scAngamb.c_k11.0.8_low=RunUMAP(scAngamb.c_k11.0.8_low, dims = 1:50)
DimPlot(scAngamb.c_k11.0.8_low, label = T)+coord_fixed()

VlnPlot(scAngamb.c_k11.0.8_low, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)















#k = 12
#scAngamb.c_k12 <- FindNeighbors(scAngamb.c, dims = 1:50, k.param = 12)
#scAngamb.c_k12 <- FindClusters(scAngamb.c_k12, resolution = 1)
#scAngamb.c_k12=RunUMAP(scAngamb.c_k12, dims = 1:50)
#k12=DimPlot(scAngamb.c_k12, label = T)+coord_fixed()

#k = 10
#scAngamb.c_k10 <- FindNeighbors(scAngamb.c, dims = 1:50, k.param = 10)
#scAngamb.c_k10 <- FindClusters(scAngamb.c_k10, resolution = 1)
#scAngamb.c_k10=RunUMAP(scAngamb.c_k10, dims = 1:50)
#k10=DimPlot(scAngamb.c_k10, label = T)+coord_fixed()

# k=8
#scAngamb.c_k8 <- FindNeighbors(scAngamb.c, dims = 1:50, k.param = 8)
#scAngamb.c_k8 <- FindClusters(scAngamb.c_k8, resolution = 1)
#scAngamb.c_k8=RunUMAP(scAngamb.c_k8, dims = 1:50)
#k8=DimPlot(scAngamb.c_k8, label = T)+coord_fixed()

# k=6
#scAngamb.c_k6 <- FindNeighbors(scAngamb.c, dims = 1:50, k.param = 6)
#scAngamb.c_k6 <- FindClusters(scAngamb.c_k6, resolution = 1)
#scAngamb.c_k6=RunUMAP(scAngamb.c_k8, dims = 1:50)
#k6=DimPlot(scAngamb.c_k8, label = T)+coord_fixed()

#plot_grid(k6, k8, k10, k12, labels = c('k6', 'k8', 'k10', 'k12'))

#k = 11 resolution 0.8
#scAngamb.c_k11.0.8 <- FindNeighbors(scAngamb.c, dims = 1:50, k.param = 11)
scAngamb.c_k12.0.95_PC40 <- FindNeighbors(scAngamb.c, dims = 1:42, k.param = 12)
scAngamb.c_k12.0.95_PC40 <- FindClusters(scAngamb.c_k12.0.95_PC40, resolution = 0.8)
scAngamb.c_k12.0.95_PC40=RunUMAP(scAngamb.c_k12.0.95_PC40, dims = 1:42)
cl=DimPlot(scAngamb.c_k12.0.95_PC40, label = T)+coord_fixed()+theme(legend.position = 'none')

#cl=DimPlot(scAngamb.c_k11.0.8, label = T)+coord_fixed()

cl_1=DimPlot(scAngamb.c_k12.0.95_PC40, label = T)+coord_fixed()+theme(legend.position = 'none')
cl_sex=DimPlot(scAngamb.c_k12.0.95_PC40, label = F, group.by = 'Sex')+coord_fixed()+labs(title = '')+theme(legend.position = 'right')
cl_comb=plot_grid(cl_1, cl_sex, labels = c('A)', 'B)'), rel_widths = c(1,1.2))
vol_cluster=VlnPlot(scAngamb.c_k12.0.95_PC40, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot_grid(cl_comb, vol_cluster, nrow = 2, rel_heights = c(1,1.2), labels = c('', 'C)'))

###########################
### SAVE R.DATA HERE! ----
###########################

# Plot proportion sex each cluster
b1=ggplot(scAngamb.c_k12.0.95_PC40@meta.data, aes(x=seurat_clusters, fill=Sex)) + 
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
  xlab('Cluster')+
  ylab('% cells')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 

# Plot proportion sex and replicate each cluster
b2=ggplot(scAngamb.c_k12.0.95_PC40@meta.data, aes(x=seurat_clusters, fill=`Sex replicate`)) + 
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
  scale_fill_met_d('Hiroshige')+
  xlab('Cluster')+
  ylab('% cells')+
  theme_minimal()+
  theme(legend.position = 'bottom')

b3=ggplot(scAngamb.c_k12.0.95_PC40@meta.data, aes(x=seurat_clusters, fill=`library`)) +
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
  scale_fill_met_d('Austria')+
  xlab('Cluster')+
  ylab('% cells')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  guides(fill=guide_legend(nrow=2, byrow=TRUE, title = 'Sub-library')) 

plot_grid(b1, b2, b3, nrow = 1, labels = c('A)', 'B)', 'C)'))

### Get markers ----

markers=FindAllMarkers(scAngamb.c_k12.0.95_PC40)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(scAngamb.c_k12.0.95_PC40, features = top10$gene) + NoLegend()

### Test markers ----

# Photoreceptor cells (ARR2)
photo=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP006263')+ggtitle('AGAP006263 \n Arrestin 2', subtitle = 'Photoreceptor cells')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()
# Muscle cells
muscle=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP006186')+ggtitle('AGAP006186 \n Calcium-transporting ATPase sarcoplasmic/endoplasmic reticulum type', subtitle = 'Muscle cells')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()
# Monoamine
mono=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP004476')+ggtitle('AGAP004476 \n Solute carrier family 18 (vesicular monoamine), member 2', subtitle = 'Monoamine cells')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()

plot_grid(cl, photo, muscle, mono, labels = c('A)', 'B)', 'C)', 'D)'))





# Get annotation markers
Agamb_gff=read.gff('./VectorBase-68_AgambiaePEST.gff')
Agamb_gff=Agamb_gff[which(Agamb_gff$type=='protein_coding_gene'|    
                            Agamb_gff$type=='ncRNA_gene'), c(1,3,9)]
Agamb_gff$geneID=sapply(Agamb_gff$attributes, function(x) unlist(str_split(unlist(str_split(x, ';'))[1], '='))[2])
Agamb_gff$description=sapply(Agamb_gff$attributes, function(x) unlist(str_split(unlist(str_split(x, ';'))[2], '='))[2])
Agamb_gff=Agamb_gff[ ,-3]

get_annot <- function(x){
  for (i in 1:nrow(x)){
    x$description[i]=Agamb_gff$description[which(rownames(x)[i]==Agamb_gff$geneID)]
  }
  return(x)
}


markers6 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 6)
markers6=markers6[which(markers6$p_val_adj<0.05),]
m6=get_annot(markers6)
gt(m6[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 6")

markers11 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 11)
markers11=markers11[which(markers11$p_val_adj<0.05),]
m11=get_annot(markers11)
gt(m11[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 11")

markers13 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 13)
markers13=markers13[which(markers13$p_val_adj<0.05),]
m13=get_annot(markers13)
gt(m13[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 13")

markers12 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 12)
markers12=markers12[which(markers12$p_val_adj<0.05),]
m12=get_annot(markers12)
gt(m12[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 12")

markers10 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 10)
markers10=markers10[which(markers10$p_val_adj<0.05),]
m10=get_annot(markers10)
gt(m10[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 10")

markers14 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 14)
markers14=markers14[which(markers14$p_val_adj<0.05),]
m14=get_annot(markers14)
gt(m14[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 14")

markers3 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 3)
markers3=markers3[which(markers3$p_val_adj<0.05),]
m3=get_annot(markers3)
gt(m3[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 3")

markers0 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 0)
markers0=markers0[which(markers0$p_val_adj<0.05),]
m0=get_annot(markers0)
gt(m0[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 0")

markers1 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 1)
markers1=markers1[which(markers1$p_val_adj<0.05),]
m1=get_annot(markers1)
gt(m1[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 0")

markers9 = FindMarkers(scAngamb.c_k12.0.95_PC40, ident.1 = 9)
markers9=markers9[which(markers9$p_val_adj<0.05),]
m9=get_annot(markers9)
gt(m9[,c(5,6)], rownames_to_stub = T)|> tab_header(title = "Markers cluster 9")


# Elav
t=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP000965')+theme(legend.position = 'none')+coord_fixed()
# N-cadherin
t=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP009723')+theme(legend.position = 'none')+coord_fixed()

# timeless 
t=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP008288')+theme(legend.position = 'none')+coord_fixed()
# period
t=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP001856')+theme(legend.position = 'none')+coord_fixed()

# Dhc93AB Dynein heavy chain at 93AB
# JO
t=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP010435')+theme(legend.position = 'none')+coord_fixed()


DotPlot(scAngamb.c_k12.0.95_PC40, features = c('AGAP010435', 'AGAP006263'))

FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP000045')+theme(legend.position = 'none')+coord_fixed()

DotPlot(scAngamb.c_k12.0.95_PC40, features = c('AGAP002560', 'AGAP004596'))
plot_grid(cl, t, ncol=2)

##



### Opsins ----

opsin=c("AGAP001161", 'AGAP001162', "AGAP002443", "AGAP002444", "AGAP002462", "AGAP005356", "AGAP006126", "AGAP007548","AGAP010089", "AGAP012982", "AGAP012985", "AGAP013149")

op_dotplot=DotPlot(scAngamb.c_k12.0.95_PC40, features = opsin, split.by ='Sex') + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))
DoHeatmap(subset(scAngamb.c_k12.0.95_PC40, downsample = 100), features = opsin, size = 3)


op1=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP010089", "AGAP006126"), blend = TRUE)+coord_fixed()
op2=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP010089", "AGAP012982"), blend = TRUE)+coord_fixed()
op3=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP006126", "AGAP012982"), blend = TRUE)+coord_fixed()
op4=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP012982", "AGAP012985"), blend = TRUE)+coord_fixed()
op5=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP012985", 'AGAP013149'), blend = TRUE)+coord_fixed()
op6=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = c("AGAP012982", "AGAP013149"), blend = TRUE)+coord_fixed()

plot_grid(op1, op2, op3, labels=c('Blue-UV', 'Blue-LW', 'UV-LW'), ncol=1)
plot_grid(op4,op5,op6, labels = c('LW-LW', 'LW-LW', 'LW-LW'), ncol=1)
###




### Annotation ----
new.cluster.ids=c('0', '1', '2', '3', 'Monoamine cells', 'Photoreceptor cells', 'Orcokinin+ neurons', 'Muscle cells', 'Photoreceptor cells', '9', 'Fat cells?', '11', 'Auditory neurons', 'Epithelial cells?', 'Pdf neurons')
names(new.cluster.ids)=levels(scAngamb.c_k12.0.95_PC40$seurat_clusters)
scAgamb_annotated=RenameIdents(scAngamb.c_k12.0.95_PC40, new.cluster.ids)
levels(scAgamb_annotated$seurat_clusters)=new.cluster.ids


DimPlot(scAgamb_annotated, label = T, repel=F)+coord_fixed()+theme(legend.position = 'none')


pclust=DimPlot(scAgamb_annotated, label = T, repel=F, label.size = 2.9)+coord_fixed()+theme(legend.position = 'none')
pclust_sex=DimPlot(scAgamb_annotated, label = F, repel=F, group.by = 'Sex')+coord_fixed()+ggtitle('')+theme(legend.position = 'right')
elav=FeaturePlot(scAgamb_annotated, features = 'AGAP000965')+ggtitle('AGAP000965 \n elav', subtitle = 'neurons')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()
cadher=FeaturePlot(scAgamb_annotated, features = 'AGAP009723')+ggtitle('AGAP009723 \n N-cadherin', subtitle = 'glial cells')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()


plot_grid(pclust, pclust_sex, elav, cadher, ncol=2, labels = c('A)', 'B)', 'C)', 'D)'), rel_widths = c(1,1.2,1,1))




opsin=c("AGAP001161", 'AGAP001162', "AGAP006126", "AGAP007548","AGAP010089", "AGAP012982", "AGAP012985", "AGAP013149")

DotPlot(scAgamb_annotated, features = opsin) + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))
DoHeatmap(subset(scAgamb_annotated, downsample = 50), features = opsin, size = 3)

DotPlot(scAgamb_annotated, features = opsin[1:9]) + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))
DotPlot(scAgamb_annotated, features = opsin[10:12]) + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))

op1=FeaturePlot(scAgamb_annotated, features = c("AGAP010089", "AGAP006126"), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
op2=FeaturePlot(scAgamb_annotated, features = c("AGAP010089", "AGAP012982"), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
op3=FeaturePlot(scAgamb_annotated, features = c("AGAP006126", "AGAP012982"), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
op4=FeaturePlot(scAgamb_annotated, features = c("AGAP012982", "AGAP012985"), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
op5=FeaturePlot(scAgamb_annotated, features = c("AGAP012985", 'AGAP013149'), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
op6=FeaturePlot(scAgamb_annotated, features = c("AGAP012982", "AGAP013149"), cells= which(scAgamb_annotated@active.ident=='Photoreceptor cells'), blend = TRUE)+coord_fixed()
plot_grid(op1, op2, op3, labels=c('A)', 'B)', 'C)'), ncol=1)
plot_grid(op4,op5,op6, labels=c('A)', 'B)', 'C)'), ncol=1)
###


cm1=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP010435')+ggtitle('AGAP010435 \n dynein heavy chain 93AB', subtitle = 'Auditory neurons')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()
cm2=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP009594')+ggtitle('AGAP009594 \n dynein assembly factor 1%2C axonemal homolog', subtitle = 'Auditory neurons')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()

cm3=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP005776')+ggtitle('AGAP005776 \n pigment-dispersing hormone', subtitle = 'Pdf neurons')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()

cm4=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP012220')+ggtitle('AGAP012220 \n Orcokinin', subtitle = 'Orcokinin+ neurons')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()

cm5=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP005563')+ggtitle('AGAP005563 \n Facilitated trehalose transporter 1', subtitle = 'Fat cells?')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()
cm6=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP008227')+ggtitle('AGAP008227 \n Trehalose 6-phosphate synthase/phosphatase', subtitle = 'Fat cells?')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()

cm7=FeaturePlot(scAngamb.c_k12.0.95_PC40, features = 'AGAP009399')+ggtitle('AGAP009399 \n nuclear hormone receptor FTZ-F1 beta', subtitle = 'Epithelial cells?')+theme(legend.position = 'none', plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(size = 10))+coord_fixed()




plot_grid(cl, cm1, cm3, cm4, cm5, cm7, labels = c('A)', 'B)', 'C)', 'D)', 'E)', 'F)'))


