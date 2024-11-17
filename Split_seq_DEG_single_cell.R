############################################################################
######################### Male vs female L77 DEG ##########################
############################################################################

# Libraries ----

library(Matrix)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggvenn)

# Import functions Elena Emili et al. 2023 ----
# url: https://www.biorxiv.org/content/10.1101/2023.11.01.565140v1
# Github url: https://github.com/scbe-lab/planarian_cell_type_allometry
source('~/Desktop/data/Split-seq/R_functions_Alberto/clean_sampletable.R')
source('~/Desktop/data/Split-seq/R_functions_Alberto/deseq_pseudobulk.R')
source('~/Desktop/data/Split-seq/R_functions_Alberto/deseq_sc.R')
source('~/Desktop/data/Split-seq/R_functions_Alberto/dynamic_colors_annotation.R')

# set working directory and import clusters ----
setwd('~/Desktop/data/Split-seq/All_analysis/')
load('./DEG_annotated.RData')

#load('./sc_Agamb_new.RData')


# Generate the matrix ----

all.genes <- rownames(scAgamb_annotated)

# Load raw matrix
mat_x <- readMM(file = "./Matrix_scAngamb/matrix.mtx.gz")
colnames(mat_x) <-  read.table("./Matrix_scAngamb/barcodes.tsv.gz")[,1]
rownames(mat_x) <-  read.table("./Matrix_scAngamb/features.tsv.gz")[,1]
#mat_x <- mat_x[rownames(mat_x) %in% all.genes,colnames(mat_x) %in% names(head_Agamb_clust_k14_l$orig.ident)]
mat_x <- mat_x[rownames(mat_x) %in% all.genes,colnames(mat_x) %in% names(scAgamb_annotated$orig.ident)]

# define variables
#identities = factor(head_Agamb_clust_k14_l$seurat_clusters, levels = unique(head_Agamb_clust_k14_l$seurat_clusters))
#conditions = factor(head_Agamb_clust_k14_l$Sex, levels = unique(head_Agamb_clust_k14_l$Sex))
#replicates = factor(head_Agamb_clust_k14_l$Replicate, levels = unique(head_Agamb_clust_k14_l$Replicate))
identities = factor(scAgamb_annotated$seurat_clusters, levels = unique(scAgamb_annotated$seurat_clusters))
conditions = factor(scAgamb_annotated$Sex, levels = unique(scAgamb_annotated$Sex))
replicates = factor(scAgamb_annotated$Replicate, levels = unique(scAgamb_annotated$Replicate))


# retrieve unique instances
idents <- levels(identities); print("found idents")
conds <- levels(conditions); print("found conds")
reps <- levels(replicates); print("found reps")

# Create the matrix
mtx <- matrix(
    0, 
    nrow = length(all.genes), # all same genes as in input data
    ncol = length(idents)*length(conds)*length(reps) # amount of combinations of cluster, experiment, replicate
  )
rownames(mtx) <- all.genes # all same genes as in input data
colnames(mtx) <- apply(expand.grid(idents,conds, reps), 1, paste, collapse="_") # combinations of cluster, experiment, replicate

# create a data frame with information for DGE
sampletable <- data.frame(
  sample = paste0("S",formatC(1:ncol(mtx),width=2,flag="0")),
  id_combined = apply(expand.grid(idents,conds, reps), 1, paste, collapse="_"),
  ctype = "",
  condition = "",
  replicate = ""
)


# begin nested loop
for (ident in idents) {
  for (cond in conds) {
    for (rep in reps) {
      sample <- paste(ident,cond,rep,sep="_") # recreate a sample name
      print(paste0("Starting sample ",sample))
      
      filt = identities == ident & conditions == cond & replicates == rep # combination of conditions: cells from cluster i in experiment j replicate z
      
      ncells_filt = length(which(filt)==TRUE)
      
      print(paste0("Found ",ncells_filt," cells of sample ",sample))
      
      if(ncells_filt == 0) {
        print(paste0("No cells in ", sample, ", skipping"))
        next 
      } 
      else {
        x_i <- mat_x[,filt] # filter the matrix of counts x cells, keep cells from combination filter
      }
      
      if(ncells_filt > 1){ # if we have more than one cell passing the filter
        pseudocounts_sample_i <- rowSums(x_i) # sum the counts of all the same cells, output at a gene level
      } 
      else {
        pseudocounts_sample_i <- x_i # otherwise this column will be the columns of that single cell
      }
      which_sample <- which(colnames(mtx) == sample)
      mtx[,which_sample] <- pseudocounts_sample_i #  dump those values in the corresponding sample column of the output matrix
      sampletable[which_sample,3:5] <- c(ident,cond,rep) # dump the values of sample etc in the corresponding row of sampletable data
      print(paste0("Done sample ",sample)) #  message to screen
    }
  }
}

colnames(mtx) <- sampletable$sample[match(sampletable$id_combined,colnames(mtx))]

# create results list
res <- list(
  matrix = mtx, # output matrix
  sampletable = sampletable # output data frame with information for DGE
)


# Clean sampletable
sex_sampletable <- clean_sampletable(res$sampletable)
sex_matrix <- 
  res$matrix[
    rownames(res$matrix) %in% all.genes,
    colnames(res$matrix) %in% sex_sampletable$sample
  ]

# SINGLE CELL ANALYSIS OVER ALL CELL TYPES

sex_DGE_all <- list()
for(i in unique(sex_sampletable$ctype)){
  print(paste0("starting with cell type ",i))
  sex_DGE_all[[i]] <-
    deseq_pseudobulk(
      count_matrix = sex_matrix,
      samples_info = sex_sampletable[,-2],
      celltype = i,
      filter_by = "pvalue", p_threshold = 0.05,
      contrast_info = c("condition","Male","Female"),
      plot_results = F, min_passing_samples = 2,
      min_counts_per_sample = 1,
      keep_dubious = FALSE
    )
  print(paste0("done cell type ",i))
}





















cluster='Pdf neurons'

ggplot(sex_DGE_all[[cluster]]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[[cluster]]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[[cluster]]$res$pvalue>0.05, NA, rownames(sex_DGE_all[[cluster]]$res))))+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')

library(ape)
library(stringr)
# Get annotation
Agamb_gff=read.gff('./VectorBase-68_AgambiaePEST.gff')
Agamb_gff=Agamb_gff[which(Agamb_gff$type=='protein_coding_gene'| 	
                            Agamb_gff$type=='ncRNA_gene'), c(1,3,9)]
Agamb_gff$geneID=sapply(Agamb_gff$attributes, function(x) unlist(str_split(unlist(str_split(x, ';'))[1], '='))[2])
Agamb_gff$description=sapply(Agamb_gff$attributes, function(x) unlist(str_split(unlist(str_split(x, ';'))[2], '='))[2])
Agamb_gff=Agamb_gff[ ,-3]

get_annot <- function(x){
  for (i in 1:nrow(x)){
    x$Log2FC[i]=sex_DGE_all[[which(names(sex_DGE_all)==DEG$Cluster[i])]]$res$log2FoldChange[which(rownames(sex_DGE_all[[which(names(sex_DGE_all)==DEG$Cluster[i])]]$res)==DEG$DEG[i])]
    x$pval[i]=sex_DGE_all[[which(names(sex_DGE_all)==DEG$Cluster[i])]]$res$pval[which(rownames(sex_DGE_all[[which(names(sex_DGE_all)==DEG$Cluster[i])]]$res)==DEG$DEG[i])]
    x$Description[i]=Agamb_gff$description[which(DEG$DEG[i]==Agamb_gff$geneID)]
  }
  return(x)
}

# Extract DEG per cell type
get_deg_per_cluster <- function(x){
  Cluster=x
  DEG = sex_DGE_all[[which(names(sex_DGE_all)==x)]]$diffgenes
  data.frame(Cluster, DEG)
}

# Get DEG 
DEG=do.call(rbind, lapply(as.list(names(sex_DGE_all)), get_deg_per_cluster))
# Remove NA
DEG=DEG[-which(is.na(DEG$DEG)),]

# Add annotation, Log2FC and p-value
DEG_des=get_annot(DEG)
# Reorder table
DEG_des=DEG_des %>% arrange(Cluster)


library(gt)
colnames(DEG_des)[2]='Gene ID'
gt(DEG_des)






#### Heatmap
library(ggpubr)
library(QuantNorm)
library(MetBrewer)
library(DESeq2)
library(ComplexHeatmap)
library(qtlTools)
library(cowplot)
colnames(DEG_des)[2]='Gene ID'

# Set colours
a=met.brewer('Egypt', n=4, 'discrete')
b=met.brewer('Johnson', n=5, 'discrete')
my_colors = list(
  condition = c(Male = a[1] , Female =a[3]), 
  replicate = c(`1`=b[1], `2`=b[3], `3`=b[4]))

# Choose cluster to generate the figure for
cluster='Photoreceptor cells'

# Take the counts 
signif=rownames(sex_DGE_all[[cluster]]$res)[which(sex_DGE_all[[cluster]]$res$pvalue<0.05)]
mat = assay(sex_DGE_all[[cluster]]$dds)[rownames(sex_DGE_all[[cluster]]$dds) %in% signif,]
mat = mat - rowMeans(mat)

df2 = as.data.frame(cbind(colData(sex_DGE_all[[cluster]]$dds)[,c("replicate"),drop=FALSE], colData(sex_DGE_all[[cluster]]$dds)[,c("condition"),drop=FALSE]))
#colnames(df2)<-c("condition")
rownames(df2)=c('Male 1', 'Female 1', 'Male 3', 'Female 3', 'Male 2', 'Female 2')
# Reorder rows in the table
#mat=mat[match(logorder, rownames(mat)),]

# Get transcript annotation
#annotgenes = read_excel('~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/new_list_genes.xlsx', sheet = 'Feuil1', col_names = F)
#for (i in 1:nrow(mat)){
#  a=which(annotgenes$...1==rownames(mat)[i])
#  rownames(mat)[i]=annotgenes$...6[a]
#}

mat=as.matrix(scale(mat))
# Change sample names
colnames(mat)=rownames(df2)

ha_sampletable <- HeatmapAnnotation(
  df = df2,
  col = my_colors, 
  show_legend = FALSE
)
anno_legend_list = lapply(ha_sampletable@anno_list,
                          function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
prp=Heatmap(mat,
  name = "z-score",
  clustering_method_rows = "ward.D2",
  bottom_annotation = ha_sampletable, 
  column_names_rot = 45, 
  show_heatmap_legend = F,
  column_order = c('Female 1', 'Female 2', 'Female 3', 'Male 1', 'Male 2', 'Male 3')
)
heatmap_legend = color_mapping_legend(prp@matrix_color_mapping, plot = FALSE)

prp_n=draw(prp, heatmap_legend_list = c(list(heatmap_legend), anno_legend_list), padding = unit(c(2, 10, 2, 2), "mm"))

prv=ggplot(sex_DGE_all[[cluster]]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[[cluster]]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[[cluster]]$res$pvalue>0.05, NA, rownames(sex_DGE_all[[cluster]]$res))))+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
DEG_des_f=DEG_des[DEG_des$Cluster==cluster,c(2:5)] %>% arrange(-Log2FC)
#DEG_des_f$Description[3]='C2H2-type domain-containing protein'
prtab=ggtexttable(DEG_des_f, rows = NULL, theme = ttheme(base_size = 9))

test1=plot_grid(grid.grabExpr(draw(prp_n)), prv, rel_widths = c(2,1), nrow = 1, labels = c('A)', 'B)'))
plot_grid(test1, prtab, nrow=2, rel_heights = c(1,1.5), labels = c('', 'C)'))

#













library(Seurat)
DimPlot(scAgamb_annotated)


p0=ggplot(sex_DGE_all[['0']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['0']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['0']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['0']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')

p1=ggplot(sex_DGE_all[['1']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['1']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['1']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['1']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p2=ggplot(sex_DGE_all[['2']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['2']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['2']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['2']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p3=ggplot(sex_DGE_all[['3']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['3']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['3']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['3']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p4=ggplot(sex_DGE_all[['4']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['4']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['4']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['4']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p5=ggplot(sex_DGE_all[['5']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['5']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['5']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['5']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p6=ggplot(sex_DGE_all[['6']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['6']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['6']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['6']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p7=ggplot(sex_DGE_all[['7']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['7']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  #geom_text_repel(aes(label=ifelse(sex_DGE_all[['7']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['7']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p8=ggplot(sex_DGE_all[['8']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['8']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['8']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['8']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p9=ggplot(sex_DGE_all[['9']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['9']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['9']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['9']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p10=ggplot(sex_DGE_all[['10']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['10']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['10']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['10']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p11=ggplot(sex_DGE_all[['11']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['11']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['11']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['11']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p12=ggplot(sex_DGE_all[['12']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['12']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['12']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['12']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p13=ggplot(sex_DGE_all[['13']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['13']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['13']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['13']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p14=ggplot(sex_DGE_all[['14']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['14']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['14']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['14']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')
p15=ggplot(sex_DGE_all[['15']]$res, aes(x = log2FoldChange, y = -log10(pvalue), colour = ifelse(sex_DGE_all[['15']]$res$pvalue>0.05, "DEG", "noDEG"))) + 
  geom_point(size=1.5)+
  geom_text_repel(aes(label=ifelse(sex_DGE_all[['15']]$res$pvalue>0.05, NA, rownames(sex_DGE_all[['15']]$res))), size=3)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  scale_colour_manual(values=c(alpha('grey', .5),alpha('red', .8)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = 'none')

plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 2) #, labels = c('Cluster 0','Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6', 'Cluster 7', 'Cluster 8', 'Cluster 9'), label_size = 10, hjust = -.005)

plot_grid(p10, p11, p12, p13, p14, p15, nrow=2)
