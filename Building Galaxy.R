library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
libs = c('1.5dpf', '2dpf', '3dpf', '5dpf', '14dpf')

# Using the UMAP coordination from the previous pipeline "Constructing Constellation Map"
Umap.df = read.table('Constellation umap.tsv', sep = '\t', header = T)
Umap.df$tissue = factor(Umap.df$tissue, levels = 0:22)

# Using the meta.skew table from the previous pipline
Meta = read.table('meta.skewed.tsv', sep = '\t', header = T, row.names = 1)
Meta$time.point = factor(Meta$time.point, levels = libs)
Meta.filtered = Meta %>% filter(Skewness >= 1 | (time.point == '1.5dpf' & Max.rank <= 10))
## Manually continue the solid line in cluster #13,18, and discontinue the point for #5 at 1.5dpf
Meta.filtered = Meta %>% filter(Skewness >= 1 | 
                                  (time.point == '1.5dpf' & Max.rank <= 10) | 
                                  tissue %in% c('lib.14dpf.13', 'lib.14dpf.18')) %>% 
  filter(!(tissue == 'lib.14dpf.5' & time.point == '1.5dpf'))
Meta.filtered$time.tissue = paste0(Meta.filtered$time.point, '.', gsub('lib.14dpf.', '', Meta.filtered$tissue))

Meta.solid.line = left_join(Meta.filtered, Umap.df, by = c('time.tissue' = 'time.tissue'))

# sigmoid = function(x){
#   y = x - (max(x, na.rm = T) + min(x, na.rm = T))/2
#   sig = 1/(1+exp(-y*10))
#   return(sig)
# }

# Import the correlation test results of motifs generated from "TF correlation test"
DF2 = data.frame()
for(i in seq_along(libs)){
  TF.test.path = paste0('TF.test/TF.test.', libs[i], '/')
  
  for(j in 0:22){
    df = read.csv(paste0(TF.test.path, 'TF.test.lib.14dpf.', j, '.csv'), stringsAsFactors = F)
    colnames(df)[1] = 'Motif.id'
    df = df %>% select(cor, Motif.id) %>% 
      mutate(time.tissue = paste0(libs[i], '.', j))
    
    df = df %>% spread(key = Motif.id, value = cor)
    
    DF2 = rbind(DF2, df)
  }
  
  message(paste0('done import ', libs[i]))
}

Umap.DF2 = inner_join(Umap.df, DF2, by = 'time.tissue')

# Import index for motif.id to motif.name
df.index = read.csv(paste0(TF.test.path, 'TF.test.lib.14dpf.', j, '.csv'), stringsAsFactors = F)
colnames(df.index)[1] = 'Motif.id'





# Import the updated motif
new.lookup = read.csv('Constellation curating TF and Motif homologous look up 053021.csv', na.strings = '', header = T, row.names = 1)
HOM.TF.lookup = new.lookup %>% filter(!is.na(Motif)) %>% select(zebrafish, Motif)
colnames(HOM.TF.lookup) = c('TF', 'Motif')
fish.TFs = HOM.TF.lookup$TF

# Import the correlation test results of gene activities generated from "TF correlation test"
DF3 = data.frame()
for(i in seq_along(libs)){
  gene.test.path = paste0('G.test/G.test.', libs[i], '/')
  
  for(j in 0:22){
    df = read.csv(paste0(gene.test.path, 'G.test.lib.14dpf.', j, '.csv'), stringsAsFactors = F, row.names = 1)
    df = df %>% filter(gene.name %in% fish.TFs) %>% 
      select(cor, gene.name) %>% 
      mutate(time.tissue = paste0(libs[i], '.', j))
    
    df = df %>% spread(key = gene.name, value = cor)
    
    DF3 = rbind(DF3, df)
  }
  
  message(paste0('done import ', libs[i]))
}

Umap.DF3 = inner_join(Umap.DF2, DF3, by = 'time.tissue')

# calculate the coefficient between Motif and TF
# calculate the NA value for both Motif and TF
HOM.TF.lookup = HOM.TF.lookup %>% mutate(coef = NA, Motif.NA = NA, TF.NA = NA)
for(i in 1:nrow(HOM.TF.lookup)){
  motif.id = df.index$Motif.id[which(df.index$Motif.names == HOM.TF.lookup$Motif[i])]
  
  HOM.TF.lookup$Motif.NA[i] = sum(is.na(Umap.DF3[, motif.id]))
  TF.NA = try(sum(is.na(Umap.DF3[, HOM.TF.lookup$TF[i]])), silent = T)
  if(class(TF.NA) != 'try-error'){HOM.TF.lookup$TF.NA[i] = TF.NA}
  
  l = try(lm(Umap.DF3[, motif.id] ~ Umap.DF3[, HOM.TF.lookup$TF[i]]), silent = T)
  if(class(l) == 'try-error'){next}
  coeff = l$coefficients[2] %>% unname()
  HOM.TF.lookup$coef[i] = coeff

}
rownames(HOM.TF.lookup) = HOM.TF.lookup$TF

#
Passed.motif = HOM.TF.lookup %>% filter(coef >= 0) %>% pull(Motif)
Passed.TF = HOM.TF.lookup %>% filter(!is.na(coef) & TF.NA <= 23 & coef >= 0) %>% pull(TF)

# Generate supplementary table 6
Meta.solid.line$TF = Meta.solid.line$motif = NA
for(i in 1:nrow(Meta.solid.line)){
  TF.test.path = paste0('TF.test/TF.test.', Meta.solid.line$time.point[i], '/')
  df = read.csv(paste0(TF.test.path, 'TF.test.lib.14dpf.', Meta.solid.line$tissue.y[i], '.csv'), stringsAsFactors = F)
  colnames(df)[1] = 'Motif.id'
  df.motif = df %>% filter(Motif.names %in% Passed.motif) %>% 
    filter(p.adj <= .05) %>% 
    arrange(desc(cor))
  row.names(df.motif) = df.motif$Motif.names
  
  gene.test.path = paste0('G.test/G.test.', Meta.solid.line$time.point[i], '/')
  df = read.csv(paste0(gene.test.path, 'G.test.lib.14dpf.', Meta.solid.line$tissue.y[i], '.csv'), stringsAsFactors = F, row.names = 1)
  df.TF = df %>% filter(gene.name %in% fish.TFs) %>% 
    select(cor, gene.name) %>% 
    mutate(time.tissue = paste0(Meta.solid.line$time.point[i], '.', Meta.solid.line$tissue.y[i])) %>% 
    arrange(desc(cor))
  rownames(df.TF) = df.TF$gene.name
  
  TF.names = df.TF %>% filter(gene.name %in% Passed.TF) %>% slice_max(order_by = cor, n = 20) %>% pull(gene.name)
  Motif.names = df.motif[unique(HOM.TF.lookup[which(HOM.TF.lookup$TF %in% TF.names), ]$Motif), ] %>% 
    slice_max(order_by = cor, n = 100) %>% 
    filter(!is.na(motif.id) & cor >= 0) %>% 
    pull(Motif.names)
  
  Meta.solid.line$TF[i] = paste(TF.names, collapse = ' | ')
  Meta.solid.line$motif[i] = paste(Motif.names, collapse = ' | ')
  
}

supplementary.tb = Meta.solid.line %>% select(time.point, tissue.x, motif, TF) %>% arrange(tissue.x)
colnames(supplementary.tb)[2] = 'tissue.score'
write.csv(supplementary.tb, file = 'supplementary table 6.csv')

# Build up "Galaxy": supplementary figure 22
DF4 = Umap.DF3[, 6:ncol(Umap.DF3)] %>% as.matrix()
DF4[DF4 < 0] = 0
Umap.DF4 = cbind(Umap.DF3[, 1:5], DF4)

outputed.TF = c()
for(i in 1:nrow(Meta.output)){
  outputed.TF = c(outputed.TF, strsplit(Meta.output$TF[i], split = ' | ', fixed = T) %>% unlist)
}
outputed.TF = unique(outputed.TF)

galaxy.path = 'Constellation UMAP Galaxy/'
for(i in seq_along(outputed.TF)){
  gene = outputed.TF[i]
  motif = HOM.TF.lookup[which(HOM.TF.lookup$TF == gene), 'Motif']
  motif.id = df.index$Motif.id[which(df.index$Motif.names == motif)]
  
  Umap.plot = Umap.DF4 %>% select(x, y, all_of(motif.id), all_of(gene)) %>% arrange((get(gene)))
  p = ggplot(Umap.plot, aes(x = x, y = y)) +
    geom_point(aes(size = get(motif.id), color = get(gene))) +
    scale_color_gradient(low = 'lightgray', high = 'red') +
    geom_path(data = Umap.DF4, aes(group = tissue), linetype = 'dashed', color = 'lightgray') +
    geom_path(data = Meta.solid.line, aes(group = tissue.y), linetype = 'solid', color = 'lightgray') +
    cowplot::theme_cowplot() +
    guides(color = guide_legend(title = 'correlation with gene'), size = guide_legend(title = 'correlation with motif')) +
    ggtitle(paste0(motif, ' and ', gene))
  
  ggsave(filename = paste0(motif, ' and ', gene, '.pdf'), path = galaxy.path, plot = p, device = 'pdf',
         units = 'in', height = 5, width = 8)
}
