############################################################################
###############        BUSCO plot with tree  + table        ################
############################################################################


### Libraries ----
library(ggtree)
library(readxl)
library(stringr)
library(ggplot2)
library(cowplot)
library(xtable)



### Import informations and tree ----

# Import informations
info=read_excel("~/Desktop/data/Phylogenomics/0_final_list.xlsx", sheet = "busco_plot") 
# 1st column = ID, then columns S, D, F, M
info=info[,c(1, 3:6, 8:11)]
# Import tree 
tree=read.tree('~/Desktop/data/Phylogenomics/Phylogenomics/Fig_trees_phylogenomics/Astral_ultrametric.tre')
# Format tip names
tree$tip.label=sapply(tree$tip.label, function (x) str_replace(x, '_', ' '))
names(tree$tip.label)=tree$tip.label


### Format BUSCO table for plotting ----

# Function get busco table graph form 
get_busco_table2 <- function(table){
  l = length(1:nrow(table))
  long=8*l
  identif = rep(NA,long)
  gr = rep(c('S', 'D', 'F', 'M'),2*l)
  lvl= rep(c(rep('Insecta', 4), rep('Diptera', 4)), l)
  value = rep(NA, long)
  #x=rep(NA,long)
  for (i in 1:nrow(table)){
    S_insecta=8*i-7
    D_insecta=S_insecta+1
    F_insecta=S_insecta+2
    M_insecta=S_insecta+3
    a=S_insecta+4
    b=S_insecta+5
    c=S_insecta+6
    d=8*i
    identif[S_insecta:d] = table[i, 1]
    value[a]=table[i,2]
    value[c]=table[i,4]
    value[d]=table[i,5]
    value[b]=table[i,3]
    value[S_insecta]=table[i,6]
    value[D_insecta]=table[i,7]
    value[F_insecta]=table[i,8]
    value[M_insecta]=table[i,9]
  }
  busco_table = data.frame(cbind(identif, value, gr, lvl)) # new table for graph
  busco_table$val=as.numeric(busco_table$value) 
  return(busco_table)
}

# Get the table to plot
busco_table2=get_busco_table2(info)
# Put correct objects
busco_table2$gr=factor(busco_table2$gr, levels=c('M', 'F', 'D', 'S'), labels=c('M', 'F', 'D', 'S'))
busco_table2$identif=factor(busco_table2$identif, levels=unique(busco_table2$identif))




### Plot tree and busco barplots ----

# Plot tree
t=ggtree(tree)+geom_tiplab(aes(x=x+10), size=3, fontface = "italic")#size=3)

# Plot tree + busco
facet_plot(t, panel='BUSCO Insecta', data=subset(busco_table2, busco_table2$lvl=='Insecta'), geom=geom_bar, mapping=aes(x=val, fill=gr), stat= "identity", position='stack', orientation="y")+
  #labs(fill='BUSCO insecta') +
  scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) + 
  xlim_tree(xlim = c(0,450))+
  theme_tree2()+
  theme(legend.position = 'none', strip.background = element_blank())+
  geom_facet(panel='BUSCO Diptera', data=subset(busco_table2, busco_table2$lvl=='Diptera'), geom=geom_bar, mapping=aes(x=val, fill=gr), stat= "identity", position='stack', orientation="y")+
#  facet_plot(t, panel='BUSCO diptera', data=subset(busco_table2, busco_table2$lvl=='Diptera'), geom=geom_bar, mapping=aes(x=val, fill=gr), stat= "identity", position='stack', orientation="y")+
  #labs(fill='BUSCO insecta') +
  scale_fill_manual(name='BUSCO', values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) + 
  xlim_tree(xlim = c(0,450))+
  theme_tree2()+
  theme(legend.position = 'bottom', strip.background = element_rect(fill = 'black'), strip.text = element_text(colour = 'white', size=12))




### Get table information phylogenomics analysis ----

# Get table 
info_bis=read_excel("~/Desktop/data/Phylogenomics/0_final_list.xlsx", sheet = "for_thesis_R") 
info_bis[is.na(info_bis)] <- ""

# Export LaTeX format table
print(xtable(info_bis), type='latex', include.rownames = FALSE, size=10, file='~/Desktop/table.tex')




### Table BLAST OG phylogenomics ----
blast=data.frame(
         stringsAsFactors = FALSE,
                          OG.ID = c("OG0000128","OG0001113","OG0001405",
                                    "OG0001557","OG0002006"),
                   Aedes.aegyti = c("AAEL004269/AAEL006055/AAEL007002-PD/AAEL009327/AAEL024632-PD/AAEL024632-PA/AAEL019966/AAEL007002-PC","AAEL021100-PA/AAEL021100-PB",
                                    "AAEL024441","AAEL027589","AAEL001513"),
              Anopheles.gambiae = c("AGAP007247/AGAP007248/AGAP004258/AGAP013126/AGAP009626/AGAP002000",
                                    "AGAP000045-PB/AGAP000045-PA","AGAP000189","AGAP000628/AGAP012340",
                                    "AGAP011562"),
        Drosophila.melanogaster = c("FBgn0030897/FBgn0083228/FBgn0036926/FBgn0039380/FBgn0265595/FBgn0265595/FBgn0013303",
                                    "FBgn0024944/FBgn0024944","FBgn0026262",
                                    "FBgn0263974","FBgn0035264"),
                           Name = c("Neuronal calcium-binding protein",
                                    "Octopamine receptor","TBP-associated factor 3","qin",
                                    "WD-repeat protein outer segment 4")
      )
colnames(blast)=c('OG ID', 'Aedes aegypti', 'Anopheles gambiae', 'Drosophila melanogaster', 'Name')

print(xtable(blast), type='latex', include.rownames = FALSE, size=10, file='~/Desktop/table.tex')
















############################################################################
###############             Genomes assemblies              ################
############################################################################

#### Libraries ----
library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(MetBrewer)
library(magick)
library(readxl)
library(ggpubr)

setwd('~/Desktop/')

# Function get busco table graph form 
get_busco <- function(x){
  s=unlist(str_split(x, '[,:%]'))[4]
  d=unlist(str_split(x, '[,:%]'))[7]
  f=unlist(str_split(x, '[,:%]'))[10]
  m=unlist(str_split(x, '[,:%]'))[13]
  busco=c('S', 'D', 'F', 'M')
  values=c(s,d,f,m)
  df=as.data.frame(cbind(values,busco))
  return(df)
}

# Import tables
scaffolds=read_delim("Coding_current/toxo.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
scaffolds_before_hic=read_delim("Coding_current/scaffold_before_hic.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
stats_assemblies=read_excel("data/DNA_work/Toxo/Assemblies_stats.xlsx")

# Plot
s1brev=ggplot(scaffolds[1:100,], aes(x=reorder(scaffold, -length), y=length, fill=met.brewer('Egypt')[3]))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=met.brewer('Egypt')[3])+
  scale_y_continuous(limits = c(0,253000000), breaks = c(0, 1e8, 2.5e8))+
  labs(x='scaffolds', y='length (bases)')+
  theme_classic()+
  theme(legend.position = 'none', axis.text.x = element_blank(),plot.margin = margin(20,0,0,20, 'pt'))

s2brev=ggplot(scaffolds_before_hic[1:100,], aes(x=reorder(scaffold, -length), y=length, fill=met.brewer('Egypt')[1]))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=met.brewer('Egypt')[1])+
  scale_y_continuous(limits = c(0,253000000), breaks = c(0, 1e8, 2.5e8))+
  labs(x='scaffolds', y='length (bases)')+
  theme_classic()+
  theme(legend.position = 'none', axis.text.x = element_blank(),plot.margin = margin(20,0,0,20, 'pt'))

img=image_read('~/Desktop/data/DNA_work/Toxo/HiC/hic_results/2024.08.20.16.24.23.HiCImage.png')
hic_contact_map=image_ggplot(img)+theme(plot.margin = margin(25,0,0,25, 'pt'))

info=ggtexttable(stats_assemblies, rows = NULL, theme = ttheme(base_size = 11))

a=plot_grid(s2brev, s1brev,  nrow=1, labels = c('A)', 'B)'))



busco_t='C:96.9%[S:95.6%,D:1.3%],F:0.9%,M:2.2%,n:3285'
buscoT=get_busco(busco_t)
busco_pT=ggplot(buscoT, aes(x='', y=as.numeric(values), group=busco)) +
  geom_bar(aes(fill = factor(busco, levels=c('M', 'F', 'D', 'S'))), width = 1, stat = "identity") +
  ylab(' ') + 
  xlab(' ')+
  coord_polar("y", start=0)+
  labs(fill='BUSCO Diptera')+
  scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) +
  theme_void()+
  theme(strip.text = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), legend.position = 'right', panel.spacing = unit(10, "lines"), plot.margin = margin(0,0,0,60, 'pt'))

c=plot_grid(hic_contact_map, busco_pT, nrow=1, labels = c('C)', 'D)'))

plot_grid(a,c,info, ncol = 1, rel_heights = c(2,2,1), labels=c('', '', 'E)'))




# Import table
scaffolds_wyem=read_delim("Coding_current/scaffolds_wyeom.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
scaffolded_wyem=read_delim("Coding_current/scaffold_ntlink.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#scaffolded_wyem_nodup=read_delim("Coding_current/nodups.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
scaffolded_on_ref=read_delim("Coding_current/Wyem_final.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

scaffolded_wyem <- scaffolded_wyem[order(scaffolded_wyem$length, decreasing = TRUE),]
#scaffolded_wyem_nodup <- scaffolded_wyem_nodup[order(scaffolded_wyem_nodup$length, decreasing = TRUE),]
scaffolds_wyem <- scaffolds_wyem[order(scaffolds_wyem$length, decreasing = TRUE),]
scaffolded_on_ref <- scaffolded_on_ref[order(scaffolded_on_ref$length, decreasing = TRUE),]

stats_assemblies_wye=read_excel("data/DNA_work/Toxo/Assemblies_stats.xlsx", sheet = 'Feuil2')

# Plot
w1=ggplot(scaffolds_wyem[1:100,], aes(x=reorder(scaffold, -length), y=length, fill=met.brewer('Egypt')[1]))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=met.brewer('Egypt')[1])+
  labs(x='scaffolds', y='length (bases)')+
  scale_y_continuous(limits = c(0,22000000))+
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.position = 'none')

busco_w_init='C:96.7%[S:56.7%,D:40.0%],F:0.7%,M:2.6%,n:3285	'
busco_w_init_p=get_busco(busco_w_init)
busco_pWinit=ggplot(busco_w_init_p, aes(x='', y=as.numeric(values), group=busco)) +
  geom_bar(aes(fill = factor(busco, levels=c('M', 'F', 'D', 'S'))), width = 1, stat = "identity") +
  ylab(' ') + 
  xlab(' ')+
  coord_polar("y", start=0)+
  labs(fill='BUSCO Diptera')+
  scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) +
  theme_void()+
  theme(strip.text = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), legend.position = 'none')

w1.inset <-
  ggdraw() +
  draw_plot(w1) +
  draw_plot(busco_pWinit, x = .6, y = .5, width = .5, height = .5)

w2=ggplot(scaffolded_wyem[1:100,], aes(x=reorder(scaffold, -length), y=length, fill = met.brewer('Egypt')[4]))+
  geom_bar(stat="identity")+
  labs(x='scaffolds', y='length (bases)')+
  scale_y_continuous(limits = c(0,22000000))+
  scale_fill_manual(values=met.brewer('Egypt')[4])+
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.position = 'none')

busco_w_s1='C:96.2%[S:94.4%,D:1.8%],F:1.0%,M:2.8%,n:3285'
busco_w_s1_p=get_busco(busco_w_s1)
busco_pWs1=ggplot(busco_w_s1_p, aes(x='', y=as.numeric(values), group=busco)) +
  geom_bar(aes(fill = factor(busco, levels=c('M', 'F', 'D', 'S'))), width = 1, stat = "identity") +
  ylab(' ') + 
  xlab(' ')+
  coord_polar("y", start=0)+
  labs(fill='BUSCO Diptera')+
  scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) +
  theme_void()+
  theme(strip.text = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), legend.position = 'none')

w2.inset <-
  ggdraw() +
  draw_plot(w2) +
  draw_plot(busco_pWs1, x = .6, y = .5, width = .5, height = .5)


w3=ggplot(scaffolded_wyem[1:100,], aes(x=reorder(scaffold, -length), y=length, fill = met.brewer('Egypt')[4]))+
  geom_bar(stat="identity")+
  labs(x='scaffolds', y='length (bases)')+
  scale_y_continuous(limits = c(0,285000000))+
  scale_fill_manual(values=met.brewer('Egypt')[4])+
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.position = 'none')

w3.inset <-
  ggdraw() +
  draw_plot(w3) +
  draw_plot(busco_pWs1, x = .6, y = .5, width = .5, height = .5)

w4=ggplot(scaffolded_on_ref[1:100,], aes(x=reorder(scaffold, -length), y=length, fill = met.brewer('Egypt')[3]))+
  geom_bar(stat="identity")+
  labs(x='scaffolds', y='length (bases)')+
  scale_y_continuous(limits = c(0,285000000))+
  scale_fill_manual(values=met.brewer('Egypt')[3])+
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.position = 'none')

busco_w_sf='C:96.1%[S:95.1%,D:1.0%],F:1.1%,M:2.8%,n:3285'
busco_w_sf_p=get_busco(busco_w_sf)
busco_pWsf=ggplot(busco_w_sf_p, aes(x='', y=as.numeric(values), group=busco)) +
  geom_bar(aes(fill = factor(busco, levels=c('M', 'F', 'D', 'S'))), width = 1, stat = "identity") +
  ylab(' ') + 
  xlab(' ')+
  coord_polar("y", start=0)+
  labs(fill='BUSCO Diptera')+
  scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) +
  theme_void()+
  theme(strip.text = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), legend.position = 'none')

w4.inset <-
  ggdraw() +
  draw_plot(w4) +
  draw_plot(busco_pWsf, x = .6, y = .5, width = .5, height = .5)


wplots=plot_grid(w1.inset, w2.inset, w3.inset, w4.inset, labels = c('A)', 'B)', 'C)', 'D)'))

info_w=ggtexttable(stats_assemblies_wye, rows = NULL, theme = ttheme(base_size = 11))

plot_grid(wplots,info_w, ncol = 1, labels = c('', 'E)'), rel_heights = c(3,1))




























############################################################################
###############           Transposable elements             ################
############################################################################

#### Libraries ----
library(ggplot2)
library(tidyverse)
library(cowplot)
library(MetBrewer)
library(readxl)
library(ggtree)
library(rotl)
library(stringr)

setwd('~/Desktop/data/DNA_work/')

#### Import data ----
TEs=read_excel("TEs.xlsx")
TEs=TEs[,c(1, 4:8)]
TEs_tib=pivot_longer(TEs, cols = 2:6, names_to = 'TE', values_to = 'val')

#### Tree ----
taxa <- tnrs_match_names(names = TEs$sp)
tree <- tol_induced_subtree(ott_ids = ott_id(taxa))
tips=tree$tip.label
for(i in 1:length(tips)){
  split=unlist(str_split(tips[i], '_'))
  tips[i]=paste(split[1], split[2])
}
tree$tip.label=tips

#### Add table ----

tr=ggtree(tree)+geom_tiplab(aes(x=x+.1), size=4, fontface = "italic")

TE_plot=facet_plot(tr, panel='TEs', data=TEs_tib, geom=geom_bar, mapping=aes(x=as.numeric(val), fill=TE), stat= "identity", position='stack', orientation="y")+
  xlim_tree(xlim = c(0,5))+
  scale_fill_met_d(name = 'Juarez') + 
  theme_tree2()+
  xlim_expand(100, 'TEs')+
  theme(legend.position = 'bottom', strip.background = element_rect(fill = 'black'), strip.text = element_text(colour = 'white', size=12))

facet_widths(TE_plot, c(Tree = .5))
  






















