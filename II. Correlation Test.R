suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(dplyr))
motif.name = atac.list[[1]][['ATAC']]@motifs@motif.names %>% unlist

for(i in seq_along(atac.list)){
  chromvar.m = atac.list[[i]][['chromvar']]@data
  meta = atac.list[[i]]@meta.data %>% select(starts_with('lib.14dpf.'))
  
  chromvar.mean = Matrix::rowMeans(chromvar.m)
  chromvar.m_avg.m = sweep(chromvar.m, MARGIN = 1, STATS = chromvar.mean, FUN = '-')
  chromvar.cov = Matrix::rowSums(chromvar.m_avg.m^2)
  
  for(arg in 1:23){
    module.score = matrix(meta[, arg])
    
    module.score.mean = mean(module.score)
    module.m_avg.score = module.score - module.score.mean
    
    r.nominator = chromvar.m_avg.m %*% module.m_avg.score %>% as.matrix
    
    module.cov = sum(module.m_avg.score^2)
    r.deliminator = sqrt(chromvar.cov * module.cov) %>% as.matrix
    
    r = r.nominator/r.deliminator
    
    n = length(module.score)
    
    t.r = r*sqrt((n-2)/(1-r^2))
    
    p.r = pmin(pt(t.r, df = n-2, lower.tail = F), pt(t.r, df = n-2, lower.tail = T)) * 2
    
    p.adj = p.adjust(p.r)
    
    TF.test = data.frame(p.value = p.r, cor = r, gene.name = row.names(chromvar.m), p.log = p.log, p.adj = p.adj,
                        row.names = row.names(chromvar.m))
    
    TF.test$Motif.names = motif.name[TF.test$gene.name]
    
    output.path = paste0('TF.test.', libs[i])
    dir.create(path = output.path)
    
    write.csv(TF.test, file = paste0(output.path, '/TF.test.', colnames(meta)[arg], '.csv'), quote = F)
  }
}

for(i in seq_along(atac.list)){
  act.m = atac.list[[i]][['ACTIVITY']]@data
  meta = atac.list[[i]]@meta.data %>% select(starts_with('lib.14dpf.'))

  gene.activity.mean = Matrix::rowMeans(act.m)
  act.m_avg.m = sweep(act.m, MARGIN = 1, STATS = gene.activity.mean, FUN = '-')
  act.cov = Matrix::rowSums(act.m_avg.m^2)
  
  for(arg in 1:23){
    module.score = matrix(meta[, arg])
    
    module.score.mean = mean(module.score)
    module.m_avg.score = module.score - module.score.mean
    
    r.nominator = act.m_avg.m %*% module.m_avg.score %>% as.matrix
    
    module.cov = sum(module.m_avg.score^2)
    r.deliminator = sqrt(act.cov * module.cov) %>% as.matrix
    
    r = r.nominator/r.deliminator
    
    n = length(module.score)
    
    t.r = r*sqrt((n-2)/(1-r^2))
    
    p.r = pmin(pt(t.r, df = n-2, lower.tail = F), pt(t.r, df = n-2, lower.tail = T)) * 2
    
    p.adj = p.adjust(p.r)
    
    G.test = data.frame(p.value = p.r, cor = r, gene.name = row.names(act.m), p.adj = p.adj,
                        row.names = row.names(act.m))
    
    output.path = paste0('G.test.', libs[i])
    dir.create(path = output.path)
    
    write.csv(G.test, file = paste0(output.path, '/G.test.', colnames(meta)[arg], '.csv'), quote = F)
  }
}
