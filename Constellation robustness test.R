# Ted: 092721
# Constellation robustness test
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(BSgenome.Drerio.UCSC.danRer11)
library(dendextend)
library(corrplot)
library(dendextend)
library(umap)
library(parameters)
library(ggrepel)
setwd('~/Desktop/scNCC revision')
rds.path = '~/Desktop/Cell competency test 36hpf snATAC/'
libs = c('36hpf', '48hpf', '72hpf', '5dpf', '14dpf')

# First, I need to generate an array of diff list based on different cluster resolutions on 14dpf
atac5 = readRDS(paste0(rds.path, 'mesenchyme2_ATAC_14dpf.rds.gz'))
DefaultAssay(atac5) = 'ATAC'
atac5 = FindClusters(atac5, resolution = .4, algorithm = 3)
atac5 = FindClusters(atac5, resolution = .6, algorithm = 3)
atac5 = FindClusters(atac5, resolution = 0.8, algorithm = 3)
atac5 = FindClusters(atac5, resolution = 1.2, algorithm = 3)
atac5 = FindClusters(atac5, resolution = 1.4, algorithm = 3)

res = c('0.4', '0.6', '0.8', '1', '1.2', '1.4')
DimPlot(atac5, group.by = paste0('ATAC_snn_res.', res), label = T) & NoLegend() & NoAxes()
atac5@meta.data %>% 
  select(starts_with('ATAC_snn_')) %>% 
  summarise_all(.funs = function(x){length(unique(x))})

# Generate enriched peak list
for(r in seq_along(res)){
  Idents(atac5) = atac5@meta.data[, paste0('ATAC_snn_res.', res[r])]
  diff.tb = FindAllMarkers(atac5, only.pos = T, test.use = 'LR', latent.vars = 'peak_region_fragments')
  diff.list = lapply(levels(diff.tb$cluster),
                     function(x){
                       diff.tb %>% 
                         filter(cluster == x) %>% 
                         filter(p_val_adj <= 1e-3) %>% 
                         pull(gene)
                     })
  names(diff.list) = paste0('lib.14dpf.res.', res[r], '.', 0:(length(diff.list)-1))
  saveRDS(diff.list, file = paste0('diff.list.14dpf.res.', res[r], '.rds.gz'))
}

# Constellation 
setwd(rds.path)

atac1 = readRDS('36hpf_forTed.rds.gz')
atac2 = readRDS('mesenchyme_48hpf.rds.gz')
atac3 = readRDS('mesenchyme_ATAC_3dpf.rds.gz')
atac4 = readRDS('mesenchyme_ATAC_5dpf.rds.gz')
atac5 = readRDS('mesenchyme2_ATAC_14dpf.rds.gz')
atac.list = list(atac1, atac2, atac3, atac4, atac5)
names(atac.list) = paste0('lib.', libs)
rm(atac1, atac2, atac3, atac4, atac5);gc()

setwd('~/Desktop/scNCC revision')
# Import enrich peak lists from different resolutions
res = c('0.4', '0.6', '0.8', '1', '1.2', '1.4')
for(r in seq_along(res)){
  assign(paste0('diff.list.', res[r]), readRDS(paste0('diff.list.14dpf.res.', res[r], '.rds.gz')))
}

# Add Chromatin Module score of each enrich peak lists
for(r in seq_along(res)){
  for(atac in seq_along(atac.list)){
    DefaultAssay(atac.list[[atac]]) = 'ATAC'
    atac.list[[atac]] = AddChromatinModule(atac.list[[atac]], 
                                           features = get(paste0('diff.list.', res[r])),
                                           genome = BSgenome.Drerio.UCSC.danRer11)
  }
}

# Extract the metadata with module scores from each resolution
for(r in seq_along(res)){
  cluster.n = length(get(paste0('diff.list.', res[r])))
  assign(paste0('meta.list.', res[r]),
                lapply(atac.list, function(x){
                  x@meta.data %>% 
                    dplyr::select(paste0('lib.14dpf.res.', res[r], '.', 0:(cluster.n-1)))
                  })
         )
}

# Calculate between library distances for each resolution
for(r in seq_along(res)){
  assign(paste0('hc.list.', res[r]),
         lapply(get(paste0('meta.list.', res[r])), 
                function(x){
                  d = dist(x %>% as.matrix() %>% t(), method = 'euclidean')
                  hc = hclust(d, method = 'complete')
                  return(hc)
         }))
}

dend.dist = function(hc.list){
  dend.list = lapply(hc.list, as.dendrogram)
  dend.list = as.dendlist(dend.list)
  names(dend.list) = paste0('lib.', c('36hpf', '48hpf', '72hpf', '5dpf', '14dpf'))
  
  dist.dend = dist.dendlist(dend.list)
  return(dist.dend)
}

for(r in seq_along(res)){
  assign(paste0('lib.dist.', res[r]),
         dend.dist(get(paste0('hc.list.', res[r]))))
}

# Calculate between tissue distance for each resolution
for(r in seq_along(res)){
  assign(paste0('tissue.dist.list.', res[r]),
         lapply(get(paste0('meta.list.', res[r])), function(x){
           d = dist(x %>% as.matrix %>% t, method = 'euclidean')
           return(as.matrix(d))
         }))
}

# Calculate the skewness table
for(r in seq_along(res)){
  assign(paste0('meta.skew.', res[r]),
         lapply(get(paste0('meta.list.', res[r])),
                function(x){
                  dfl = Reduce(rbind, apply(x, 2, skewness))
                  df = data.frame(tissue = colnames(x), Skewness = dfl$Skewness, SE = dfl$SE)
                  df
                })
  )
  
  meta.working = get(paste0('meta.skew.', res[r]))
  for(i in seq_along(meta.working)){
    meta.working[[i]]$time.point = libs[i]
  }
  
  meta.working = Reduce(rbind, meta.working)
  meta.working$time.point = factor(meta.working$time.point, levels = libs)
  
  assign(paste0('meta.skew.', res[r]), meta.working)
  rm(meta.working)
}

# Build up skewness table
skew = c(.4, .8, 1, .9, 0)
for(r in seq_along(res)){
  meta.solid = data.frame()
  meta.working = get(paste0('meta.skew.', res[r]))
  passed.clusters = meta.working$tissue %>% unique()
  
  for(l in length(libs):1){
    meta.passed = meta.working %>% 
      filter(time.point == libs[l] & Skewness >= skew[l] & tissue %in% passed.clusters)
    
    meta.solid = rbind(meta.solid, meta.passed)
    passed.clusters = meta.passed$tissue %>% unique()
  }
  
  meta.solid$tissue = regmatches(meta.solid$tissue, regexpr('[^\\.]*$', meta.solid$tissue))
  meta.solid$time.tissue = paste0(meta.solid$time.point, '.', meta.solid$tissue)
  assign(paste0('meta.solid.', res[r]), meta.solid)
  rm(meta.passed, meta.solid, passed.clusters, meta.working)
}

# Calculate the final distance matrix for umap
# Manual a inputs
As = c(20,20,15.75,12,10,11)
for(r in seq_along(res)){
  cluster.n = get(paste0('tissue.dist.list.', res[r]))[[1]] %>% dim %>% .[1]
  lib.n = length(get(paste0('tissue.dist.list.', res[r])))
  
  # Expand the tissue distance
  D.tissue = Reduce(cbind, get(paste0('tissue.dist.list.', res[r])))
  D.tissue = D.tissue[rep(1:cluster.n, lib.n), ]
  D.tissue = (D.tissue + t(D.tissue))/2
  
  # Expand the lib distance
  D.lib = get(paste0('lib.dist.', res[r])) %>% as.matrix()
  D.lib = D.lib[rep(1:lib.n, each = cluster.n), rep(1:lib.n, each = cluster.n)]
  
  # Estimate the parameter "a" and calculate the final matrix
  # a = floor(max(D.tissue)/max(D.lib))/4
  a = As[r]
  D = D.tissue + D.lib*a
  
  # assign(paste0('Dist.', res[r]), D)
  set.seed(420)
  Umap = umap(D, config = umap.defaults, input = 'dist')
  Umap.df = Umap$layout %>% as.data.frame()
  colnames(Umap.df) = c('x', 'y')
  Umap.df = Umap.df %>% 
    mutate(time = rep(libs[1:lib.n], each=cluster.n), 
           tissue = as.character(rep(0:(cluster.n-1), time = lib.n)), 
           time.tissue = paste0(time, '.', tissue)) %>% 
    mutate(time = factor(time, levels = libs))
  # write.table(Umap.df, 'Constellation umap seed 420.tsv', sep = '\t')
  
  # Umap.df = read.table('Constellation umap seed 420.tsv', sep = '\t', header = T)
  Umap.df$tissue = factor(Umap.df$tissue, levels = 0:(cluster.n-1))
  
  meta.solid = get(paste0('meta.solid.', res[r]))
  Umap.df.solid = Umap.df %>% filter(time.tissue %in% meta.solid$time.tissue)
  
  p = ggplot(Umap.df, aes(x = x, y = y, color = tissue)) +
    geom_point(aes(shape = time), size = 3) +
    geom_path(aes(group = tissue, color = tissue), linetype = 'dashed') +
    geom_text_repel(data = Umap.df %>% filter(time == '14dpf'), aes(label = tissue), color = 'black') +
    geom_path(data = Umap.df.solid, aes(group = tissue, color = tissue), linetype = 'solid') +
    cowplot::theme_cowplot() +
    guides(color = F) +
    ggtitle(paste0('resolution = ', res[r], ', a = ', a))
  
  ggsave(filename = paste0('Constellation with skew adj ', res[r], '.pdf'), plot = p,
         device = 'pdf', units = 'in', width = 7, height = 6)
}

# output meta lists for faster calculation next time
for(r in seq_along(res)){
  saveRDS(get(paste0('meta.list.', res[r])), file = paste0('meta.list.', res[r], '.rds.gz'))
}

