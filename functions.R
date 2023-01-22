top_correlated <- function(dataset, genes, threshold,anti_cor = F) {
  require(Seurat)
  markers_expression = FetchData(object = dataset,vars = genes,slot = "data") #get genes expression
  markers_average = rowMeans(markers_expression) %>% as.data.frame() %>% rename("average" = 1) #average them
  cor_mat = cor(expression %>% t(), markers_average)%>% as.data.frame() #cor with all genes
  cor_mat = cor_mat[complete.cases(cor_mat),,drop=F]  %>% as.data.frame %>%  rename("corr" = 1) #remove rows with NA in at least one col
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