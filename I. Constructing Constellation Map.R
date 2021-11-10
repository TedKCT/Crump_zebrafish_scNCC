# This is a general script of how to generate the Constellation map
## Library the required packages
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(BSgenome.Drerio.UCSC.danRer11)
library(corrplot)
library(dendextend)
library(umap)
library(parameters)
library(ggrepel)

## Define the time points we're working on, in this case we're only using the first 5 time points
libs = c('1.5dpf', '2dpf', '3dpf', '5dpf', '14dpf')

## Import required rds objects (ATAC-RNA aligned multiome-like Signac objects)
snap.1.5dpf = readRDS('snapATAC.1.5dpf.rds.gz')
snap.2dpf = readRDS('snapATAC.2dpf.rds.gz')
snap.3dpf = readRDS('snapATAC.3dpf.rds.gz')
snap.5dpf = readRDS('snapATAC.5dpf.rds.gz')
snap.14dpf = readRDS('snapATAC.14dpf.rds.gz')

## Cluster/Tissue identification of the 14 dpf object 
DefaultAssay(snap.14dpf) = 'ATAC'
snap.14dpf = FindClusters(snap.14dpf, resolution = 1, algorithm = 3)

## Generate the cluster/tissue specifically enriched peak list
Idents(snap.14dpf) = snap.14dpf$ATAC_snn_res.1
diff.tb = FindAllMarkers(snap.14dpf, only.pos = T, test.use = 'LR', latent.vars = 'peak_region_fragments')
diff.list = lapply(levels(diff.tb$cluster),
                   function(x){
                     diff.tb %>% 
                       filter(cluster == x) %>% 
                       filter(p_val_adj <= 1e-3) %>% 
                       pull(gene)
                   })
names(diff.list) = paste0('lib.14dpf.cluster.', 0:(length(diff.list)-1))

snap.list = list(snap.1.5dpf, snap.2dpf, snap.3dpf, snap.5dpf, snap.14dpf)
names(snap.list) = paste0('lib.', libs)

## Calculate the tissue module score chromatin module score of each tissue based on their enriched peaks
for(obj in seq_along(snap.list)){
  snap.list[[obj]] = AddChromatinModule(snap.list[[obj]], features = diff.list, genome = BSgenome.Drerio.UCSC.danRer11)
}

## Extract the meta data from each Signac object for their module scores
meta.list = lapply(snap.list, function(x){
  x@meta.data %>% 
    dplyr::select(starts_with('lib.14dpf.cluster.'))
}
)
names(meta.list) = paste0('lib.', libs)

## Calculate the hierarchy cluster with every module scores for each library
hc.list = lapply(meta.list, 
                 function(x){
                   d = dist(x %>% as.matrix() %>% t(), method = 'euclidean')
                   hc = hclust(d, method = 'complete')
                   return(hc)
                 })

## Calculate the distances (Robinson-Foulds distance) between each pair of libraries
dend.dist = function(hc.list){
  dend.list = lapply(hc.list, as.dendrogram)
  dend.list = as.dendlist(dend.list)
  names(dend.list) = paste0('lib.', libs)
  
  dist.dend = dist.dendlist(dend.list)
  return(dist.dend)
}

lib.dist = dend.dist(hc.list)

## Calculate the distances between each pair of tissue module scores at every time point
tissue.dist.list = lapply(meta.list, function(x){
  d = dist(x %>% as.matrix %>% t, method = 'euclidean')
  return(as.matrix(d))
})

## Calculate the skewness of every tissue module score
meta.skew = lapply(meta.list,
                   function(x){
                     dfl = Reduce(rbind, apply(x, 2, skewness))
                     df = data.frame(tissue = colnames(x), Skewness = dfl$Skewness, SE = dfl$SE)
                     df
                   })

for(i in seq_along(meta.skew)){
  meta.skew[[i]]$time.point = libs[i]
}

meta.skew = Reduce(rbind, meta.skew)
meta.skew$time.point = factor(meta.skew$time.point, levels = libs)

## Build up the skewed tissue-timepoint table
skew = c(.4, 1, 1, 1, 0)
meta.skewed = data.frame()
passed.clusters = meta.skew$tissue %>% unique()

for(l in 1:length(libs)){
  meta.passed = meta.skew %>% 
    filter(time.point == libs[l] & Skewness >= skew[l] & tissue %in% passed.clusters)
  
  meta.skewed = rbind(meta.skewed, meta.passed)
  passed.clusters = meta.passed$tissue %>% unique()
}

meta.skewed$tissue = regmatches(meta.skewed$tissue, regexpr('[^\\.]*$', meta.skewed$tissue))
meta.skewed$time.tissue = paste0(meta.skewed$time.point, '.', meta.skewed$tissue)
## output meta.skewed table
write.table(meta.skewed, file = 'meta.skewed.tsv', sep = '\t', row.names = T,
            col.names = T, quote = F)

cluster.n = tissue.dist.list[[1]] %>% dim %>% .[1]
lib.n = length(tissue.dist.list)

## Expand the tissue distance
D.tissue = Reduce(cbind, tissue.dist.list)
D.tissue = D.tissue[rep(1:cluster.n, lib.n), ]
D.tissue = (D.tissue + t(D.tissue))/2

## Expand the lib distance
D.lib = lib.dist %>% as.matrix()
D.lib = D.lib[rep(1:lib.n, each = cluster.n), rep(1:lib.n, each = cluster.n)]

## Build up the full distance matrix for each pair of tissue-timepoint
## Estimate the parameter "a" and calculate the final matrix
## a = floor(max(D.tissue)/max(D.lib))/4
a = 12 # optimized after visualization
D = D.tissue + D.lib*a

## Dimensional reduction by UMAP
set.seed(420)
Umap = umap(D, config = umap.defaults, input = 'dist')
Umap.df = Umap$layout %>% as.data.frame()
colnames(Umap.df) = c('x', 'y')
Umap.df = Umap.df %>% 
  mutate(time = rep(libs[1:lib.n], each=cluster.n), 
         tissue = as.character(rep(0:(cluster.n-1), time = lib.n)), 
         time.tissue = paste0(time, '.', tissue)) %>% 
  mutate(time = factor(time, levels = libs))

Umap.df$tissue = factor(Umap.df$tissue, levels = 0:(cluster.n-1))
## output UMAP coordinations
write.table(Umap.df, file = 'Constellation umap.tsv', sep = '\t',
            row.names = T, col.names = T, quote = F)

Umap.df.solid = Umap.df %>% filter(time.tissue %in% meta.skewed$time.tissue)

## ggplot
ggplot(Umap.df, aes(x = x, y = y, color = tissue)) +
  geom_point(aes(shape = time), size = 1.5) +
  scale_shape_manual(values = 1:7) +
  geom_path(aes(group = tissue, color = tissue), linetype = 'dotted') +
  geom_text_repel(data = Umap.df %>% filter(time == '7mpf'), aes(label = tissue), color = 'black') +
  geom_path(data = Umap.df.solid, aes(group = tissue, color = tissue), linetype = 'solid') +
  cowplot::theme_cowplot() +
  guides(color = F)



