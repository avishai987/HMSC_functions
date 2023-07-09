top_correlated <-function(dataset, genes, threshold,anti_cor = F,n_vargenes =0) {
  require(Seurat)
  require(dplyr)
  markers_expression = FetchData(object = dataset,vars = genes,slot = "data") #get genes expression
  markers_average = rowMeans(markers_expression) %>% as.data.frame() %>% dplyr::rename("average" = 1) #average them
  expression = GetAssayData(object = dataset,assay = "RNA",slot = "data") %>% as.data.frame() #get all genes expression data
  
  if (n_vargenes != 0){ #filter genes to var genes
    if ((VariableFeatures(dataset) %>% length())<n_vargenes) { #if exist vat genes in not enought, compute
      dataset = FindVariableFeatures(object = dataset,nfeatures = n_vargenes)
    }
    var_features = head(VariableFeatures(dataset),n_vargenes)
    expression = expression[rownames(expression) %in% var_features,]
  }
  cor_mat = cor(expression %>% t(), markers_average)%>% as.data.frame() #cor with all genes
  cor_mat = cor_mat[complete.cases(cor_mat),,drop=F]  %>% as.data.frame %>%  dplyr::rename("corr" = 1) #remove rows with NA in at least one col
  if (threshold<1){ #if threshold is based on pearson correlation 
    if(anti_cor == T){top_genes =   cor_mat %>% as.data.frame %>% select(1) %>% dplyr::filter(.< threshold) %>% rownames()}else{
      top_genes =   cor_mat %>% as.data.frame %>% select(1) %>% dplyr::filter(.> threshold) %>% rownames()
    }
  }else{ #if threshold is based on top correlated genes 
    if(anti_cor == T){threshold  = threshold*(-1)}
    top_genes =   cor_mat %>%  top_n(threshold,corr) %>% rownames()
  }
  return(top_genes)
}

top_genes_cor_heatmap <- function(dataset, top_genes) {
  require(pheatmap)
  require(Seurat)
  expression = GetAssayData(object = dataset,assay = "RNA",slot = "data") %>% as.data.frame()
  top_expression = expression %>% dplyr::filter(rownames(expression) %in% top_genes)
  colors <- c(seq(-1,1,by=0.01))
  my_palette <- c("blue",colorRampPalette(colors = c("blue", "white", "red"))
                  (n = length(colors)-3), "red")
  pht = pheatmap(mat = cor(top_expression %>% t(), top_expression %>% t()),color = my_palette, breaks = colors)
  return(pht)
}


enriched_score_umap <- function(dataset,enrich_res, genes,col,distribution = F) {
  rownames(enrich_res) = enrich_res$pathway_name
  enriched_genes = enrich_res[col,"geneID"] %>% strsplit(split = "/") %>% .[[1]] %>% c(.,genes) #add original markers
  enriched_genes_score=apply(dataset@assays$RNA@data[enriched_genes,],2,mean)
  dataset=AddMetaData(dataset,enriched_genes_score,"enriched_genes_score")
  if (distribution == F) {
    print(FeaturePlot(object = dataset,features = "enriched_genes_score"))
  }else{
    data = FetchData(object = dataset,vars = "enriched_genes_score")
    print(
      data %>% 
        ggplot(aes( x=enriched_genes_score)) + 
        geom_density() 
    )
  }
  return(list( score = enriched_genes_score,genes = enriched_genes))
}


calculate_score = function(dataset,myo_genes,lum_genes,lum_threshold =1 , myo_threshold = -1) {
  myoscore=FetchData(object =dataset,vars =  myo_genes,slot = "data") %>% rowMeans()
  lescore=FetchData(object =dataset,vars =  lum_genes,slot = "data") %>% rowMeans()
  correlation = cor(lescore,myoscore) %>% round(digits = 2)
  message("correlation of lum score and myo score:" %>% paste(correlation))
  
  
  original_myo_genes = c("TP63", "TP73", "CAV1", "CDH3", "KRT5", "KRT14", "ACTA2", "TAGLN", "MYLK", "DKK3")
  original_lum_genes = c("KIT", "EHF", "ELF5", "KRT7", "CLDN3", "CLDN4", "CD24", "LGALS3", "LCN2", "SLPI")
  orig_myoscore=FetchData(object =dataset,vars =  original_myo_genes,slot = "data") %>% rowMeans()
  orig_lescore=FetchData(object =dataset,vars =  original_lum_genes,slot = "data") %>% rowMeans()
  correlation_to_original_lum = cor(orig_lescore,lescore) %>% round(digits = 2)
  correlation_to_original_myo = cor(orig_myoscore,myoscore) %>% round(digits = 2)

  message("correlation of lum score and original lum score:" %>% paste(correlation_to_original_lum))
  message("correlation of myo score and original myo score:" %>% paste(correlation_to_original_myo))

  dataset=AddMetaData(dataset,lescore-myoscore,"luminal_over_myo")
  print(
    FeaturePlot(object = dataset,features = "luminal_over_myo")
  )
  data = FetchData(object = dataset,vars = "luminal_over_myo")
  print(
    data %>% 
    ggplot(aes( x=luminal_over_myo)) + 
    geom_density() 
    )
  
lum_cells_num = subset(x = dataset,luminal_over_myo >(lum_threshold)) %>% ncol() /ncol(dataset)
myo_cells_num = subset(x = dataset,luminal_over_myo <(myo_threshold)) %>% ncol()/ncol(dataset)
df = data.frame(cell_type = c("myo_cells","lum_cells"),percentage = c(myo_cells_num,lum_cells_num))
ggplot(data=df, aes(x=cell_type, y=percentage)) +
  geom_bar(stat="identity") + ggtitle("ACC cell types")
}
