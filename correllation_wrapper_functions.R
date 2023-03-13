filter_input_correllation_data <- function(
  ph,
  mta = 0,   # min taxon count
  swpr = 0.2 # sample-wise taxon presence threshold
  
  ) {
  #load("data/ll_deep.rda")
  
  library(phyloseq)
  library(magrittr)
  
  counts = matrix(otu_table(ph))
  
  dim(counts) = otu_table(ph) %>% dim()
  dimnames(counts) = list(
    rownames(otu_table(ph)),
    colnames(otu_table(ph))
  )
  counts %<>% t()
  
  taxa = matrix(tax_table(ph)[,1:6])
  dim(taxa)= tax_table(ph)[,1:6] %>% dim
  dimnames(taxa) = list(
    rownames(tax_table(ph)[,1:6]),
    colnames(tax_table(ph)[,1:6])
  )
  
  if(mta > 0){
    ts = colSums(counts)
    counts = cbind(counts[, ts >= mta ], orphans = rowSums(counts[, ts < mta ]))
    taxa = rbind(taxa[ts >= mta,], orphans = c('Bacteria',rep('orphans', 5)))
  }
  
  counts_el =
    counts[,apply(counts > 0, 2, sum) / dim(counts)[1] > swpr] %>%
    as.matrix
  counts_el[,'orphans'] = 
    counts_el[,'orphans'] + rowSums(counts[,apply(counts > 0, 2, sum) / dim(counts)[1] <= swpr])
  #counts_el = as.matrix(counts_el)
  taxmat_el = taxa[apply(counts > 0, 2, sum) / dim(counts)[1] > swpr,]
  if(grep('Kingdom',colnames(taxmat_el)) == 0){
    colnames(taxmat_el) = c("Kingdom","Phylum","Class","Order","Family","Genus")
  }
  return(list(counts = counts_el, taxa = taxmat_el))
}

# Wrapper function to compute correlations by available methods
compute_corr_graphs <- function(input_data, target_edge_number = 1200, raw_corrs = NULL) {
  
  library(SpiecEasi)
  library(huge)
  library(rMAGMA)
  library(JuliaCall)
  
  library(magrittr)
  #library(igraph)
  
  setwd("/media/repo4T/yuragal_home/LIN/metagenome/Likhoshway_post_anaysis/bashenhaeva/mol_a/data/Zi-LN")
  
  source("utils/magma_fix.R")
  source("inference.R")
  
  if(is.null(raw_corrs)) {
    JULIA_HOME = "/home/yuragal/tools/dev/julia-1.7.0-rc3/bin"
    julia_setup(JULIA_HOME = JULIA_HOME)
  }
  
  edge_number = 0
  
  counts_el = input_data$counts
  taxmat_el = input_data$taxa
  
  ### ZiLN
  if(is.null(raw_corrs)) {
    Z = infer_Z(counts_el)
  } else {
    Z = raw_corrs$zi
  }
  
  compute_graph_with_target_number_of_edges_mb = function(X, target_edge_number, tol = 1) {
    l_min = 0
    l = 0.5
    l_max = 1
    res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
    edge_number = sum(res$path[[1]] != 0) /2
    while(edge_number > target_edge_number + tol || edge_number < target_edge_number - tol){
      if(edge_number > target_edge_number + tol) {
        l_ = (l_max + l) / 2
        l_min = l
        l = l_
      }
      else{
        l_ = (l_min + l) / 2
        l_max = l
        l = l_
      }
      print(l)
      res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
      edge_number = sum(res$path[[1]] != 0) /2
    }
    return(res$path[[1]])
  }
  
  graph.zi = compute_graph_with_target_number_of_edges_mb(Z, target_edge_number)
  
  if(is.null(raw_corrs)) {
    X_se = t(SpiecEasi::clr(counts_el + 1, 1))
  } else {
    X_se = raw_corrs$se
  }
  graph.se = compute_graph_with_target_number_of_edges_mb(X_se, target_edge_number)
  
  ########### Flashweave
  flashweave_wrapper = function(X, alpha = 0.01) {
    fw_input_data = paste(getwd(), 'data', 'data_mat.csv', sep='/')
    
    write.csv(X, file = fw_input_data)
    
    julia_command("using FlashWeave")
    julia_command(paste0('data_path = "', fw_input_data, '"'))
    
    julia_command(paste0('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ', alpha, ')'))
    julia_command(paste0('save_network("data/network.edgelist", netw_results)'))
    
    edge_list = read.table('data/network.edgelist')
    
    d = dim(counts_el)[2]
    node_names = colnames(counts_el)
    
    v1 = factor(edge_list[,1], levels = node_names)
    v2 = factor(edge_list[,2], levels = node_names)
    edges = abs(edge_list[,3])
    S_flashweave = matrix(0, d, d)
    for(i in 1:length(v1)){
      val = edges[i]
      S_flashweave[v1[i], v2[i]] = val
      S_flashweave[v2[i], v1[i]] = val
    }
    return(S_flashweave)
  }
  
  
  compute_graph_with_target_number_of_edges_corr = function(S, target_edge_number, tol = 1) {
    l_min = 0
    l = 0.5
    l_max = 1
    diag(S) = 0
    edge_number = sum(S > l) /2
    while(edge_number > target_edge_number + tol || edge_number < target_edge_number - tol){
      if(edge_number > target_edge_number + tol) {
        l_ = (l_max + l) / 2
        l_min = l
        l = l_
      }
      else{
        l_ = (l_min + l) / 2
        l_max = l
        l = l_
      }
      print(c(l, edge_number))
      edge_number = sum(S > l) /2
      if(l == 0){
        return(S > l)
      }
    }
    print('Final:')
    print(c(l, edge_number))
    return(S > l)
  }
  
  if(is.null(raw_corrs)) {
    S_flashweave = flashweave_wrapper(counts_el, 0.01)
  } else {
    S_flashweave = raw_corrs$fw
  }
  graph.flashweave = compute_graph_with_target_number_of_edges_corr(S_flashweave, target_edge_number)
  
  ########## Sparcc 
  
  if(is.null(raw_corrs)) {
    scc = sparcc(counts_el)
  } else {
    scc = raw_corrs$sc
  }
  graph.sparcc = compute_graph_with_target_number_of_edges_corr(abs(scc$Cor), target_edge_number)
  
  ######## MAGMA
  
  if(is.null(raw_corrs)) {
    res_norm <- magma.norm(data=counts_el, distrib="ZINB", seq_depth="GMPR")
  } else {
    res_norm = raw_corrs$mg
  }
  graph.magma = compute_graph_with_target_number_of_edges_mb(res_norm$z, target_edge_number)
  
  return(list(
    raw_graphs = list(
      se = graph.se,
      mg = graph.magma,
      zi = graph.zi,
      sc = graph.sparcc,
      fw = graph.flashweave
    ),
    raw_corrs = list(
      se = X_se,
      mg = res_norm,
      zi = Z,
      sc = scc,
      fw = S_flashweave
    )
  ))
}

##### Draw Venn
draw_venn = function(graphs, fname = NULL){
  library(VennDiagram)
  
  z1 = graphs$se != 0
  z2 = graphs$mg != 0
  z3 = graphs$zi != 0
  z4 = graphs$sc != 0
  z5 = graphs$fw != 0
  
  graph = list(z1, z2, z3, z4, z5)
  
  intersections = double(10 + 10 + 5)
  cex = double(5 + 10 + 10 + 5 + 1)
  
  n = 5
  l = 0
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      l = l + 1
      intersections[l] = sum(graph[[i]] & graph[[j]]) / 2
    }
  }
  
  for(i in 1:(n - 2)){
    for(j in (i + 1):(n - 1)){
      for(k in (j + 1):n) {
        l = l + 1
        intersections[l] = sum(graph[[i]] & graph[[j]] & graph[[k]]) / 2
      }
    }
  }
  
  for(i in 1:(n - 3)){
    for(j in (i + 1):(n - 2)){
      for(k in (j + 1):(n - 1)) {
        for(m in (k + 1):n){
          l = l + 1
          intersections[l] = sum(graph[[i]] & graph[[j]] & graph[[k]] & graph[[m]]) / 2
        }
      }
    }
  }
  a1 = sum(graph[[1]]) / 2
  a2 = sum(graph[[2]]) / 2
  a3 = sum(graph[[3]]) / 2
  a4 = sum(graph[[4]]) / 2
  a5 = sum(graph[[5]]) / 2
  a12345 = sum(graph[[1]] & graph[[2]] & graph[[3]] & graph[[4]] & graph[[5]]) / 2
  
  if(!is.null(fname)){
    pdf(file=fname, height=11, width=11)
  }
  plot.new()
  venn.plot <- draw.quintuple.venn(
    area1 = a1, 
    area2 = a2, 
    area3 = a3, 
    area4 = a4, 
    area5 = a5, 
    n12 = intersections[1], 
    n13 = intersections[2], 
    n14 = intersections[3], 
    n15 = intersections[4], 
    n23 = intersections[5],
    n24 = intersections[6],
    n25 = intersections[7],
    n34 = intersections[8],
    n35 = intersections[9],
    n45 = intersections[10],
    n123 = intersections[11],
    n124 = intersections[12],
    n125 = intersections[13],
    n134 = intersections[14],
    n135 = intersections[15],
    n145 = intersections[16], 
    n234 = intersections[17],
    n235 = intersections[18],
    n245 = intersections[19],
    n345 = intersections[20],
    n1234 = intersections[21],
    n1235 = intersections[22],
    n1245 = intersections[23],
    n1345 = intersections[24],
    n2345 = intersections[25],
    n12345 = a12345,
    category = c("Spiec-Easi", "Magma", "ZiLN", "Sparcc", "Flashweave"),
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), 
    cat.cex = 1.5,
    margin = 0.07,
    ind = TRUE
  );
  if(!is.null(fname)){
    dev.off()
  }
  return(c(a1, a2, a3, a4, a5, intersections, a12345))
}

count_shared_edges <- function(
  CGS,
  r = 3 #number of methods to calculate shared edges
){
  library(igraph)
  library(magrittr)
  
  if(! r %in% 2:3){
    warning('Only can compute shared edges for r %in% c(2:3).')
    return()
  }
  
  combs = gtools::combinations(n = names(CGS) %>% length, r = r, v = names(CGS))
  if(r == 3) {
    she = apply(combs, 1, function(x) {
      l = E(CGS[[ x[1] ]]$gr %s% CGS[[ x[2] ]]$gr %s% CGS[[ x[3] ]]$gr) %>% length
      return(c(x, l))
    })
  }
  
  if(r == 2) {
    she = apply(combs, 1, function(x) {
      l = E(CGS[[ x[1] ]]$gr %s% CGS[[ x[2] ]]$gr) %>% length
      return(c(x, l))
    })
  }
  
  return(t(she))
}

# Merge co-occurence graphs generated by specific methods
# Draw the final graph colored by:
# -- taxons of a rank
# -- by mean abundance and LFC for specific contrast

merge_graphs_and_draw_assort_contrast <- function(
  corr_graphs,
  corr_input,
  taxon_rank = 'Phylum',
  method_names = c('mg','se','zi'), # graphs to merge
  ds2 = NULL,
  raw_counts = NULL,
  contrast = 'lifestyle',
  lfc_range = c(-5,5),
  alpha = 0.05,
  lfcThreshold = log2(1.5),
  show.vertex.names = FALSE,
  lo = NA,
  return.results = FALSE,
  fn = NULL,
  plot.phylotypes = FALSE,
  plot.de = FALSE,
  new.plot = TRUE,
  plt = 'RdBy',
  vertex_size = 3,
  ...
  ){
  
  library(igraph)
  library(magrittr)
  library(DESeq2)
  
  # construct graph with edges shared between target graphs/methods
  # remove vetrices with no edges
  gr = corr_graphs[[ method_names[1] ]]$gr
  for(m in 2:length(method_names)){
    gr = gr %s% corr_graphs[[ method_names[m] ]]$gr
  }
  todel = which(degree(gr) == 0 | V(gr)$name == 'orphans')
  if(length(todel) > 0){
    gr = delete.vertices(gr, todel)
  }
  
  if(is.null(dim(lo))){
    lo = layout_with_fr(gr)
  }
  ceb <- cluster_edge_betweenness(gr)
  
  lbls = corr_input$taxa[V(gr)$name, taxon_rank] %>% unique() %>% sort
  cols = pals::kelly(lbls %>% length)
  names(cols) = lbls
  cols
  V(gr)$color = cols[corr_input$taxa[V(gr)$name, taxon_rank]]
  
  # get abundance data from deseq2 results
  if(is.null(ds2)){
    ds2 = estimate_deseq2(raw_counts, contrast[1])
  }
  
  if(length(contrast) == 1){
    ds2_res = results(
      ds2,
      lfcThreshold = lfcThreshold,
      alpha = alpha,
      #contrast = c(contrast, ds2_cc@colData[ , contrast ] %>% levels())
    )
  } else {
    ds2_res = results(
      ds2,
      lfcThreshold = lfcThreshold,
      alpha = alpha,
      contrast = c(contrast[1], paste0(contrast[2], contrast[-(1:2)]))
    )
  }
  lfcs = ds2_res[V(gr)$name, c('baseMean', 'log2FoldChange','padj')]
  names(lfcs) = c('abund','lfc','padj')
  lfcs$lfc %<>% ifelse(
    abs(.) <= max(lfc_range),
    .,
    ifelse(. > lfc_range[2], lfc_range[2], lfc_range[1]))
  lfcs$lfc[lfcs$padj >= alpha] = NA
  
  colfunc = scales::col_numeric(plt, lfc_range, reverse = ifelse(plt == 'RdBy', TRUE, FALSE))
  
  if(show.vertex.names){
    vertex.names = V(gr)$name
  } else {
    #vertex.names = rep(NA, length(V(gr)) )
    vertex.names = NA
  }
  
  if(!is.null(fn)){
    width = 24
    height = 12
    if(!plot.phylotypes){
      width = width/2
    }
    
    svg(filename = fn, width = width, height = height)
  }
  
  if(new.plot){
    par(mfrow=c(1, sum(plot.phylotypes, plot.de)))
  }
  
  if(plot.phylotypes){
    plot(
      gr,
      layout = lo,
      #vertex.label = ifelse(V(gr)$name == 'OTU155','OTU155',NA),
      vertex.label = vertex.names,
      vertex.label.cex = 0.2,
      vertex.label.color = 'green4',
      vertex.size = vertex_size,
      main = sprintf(
        '%s, deg > 0, color by %s',
        paste(c('mg','se','zi'), collapse = ' & '),
        taxon_rank
        ),
      edge.color = 'gray'
    )
    
    legend(
      x = -1.3,
      y = ifelse(is.null(fn), -0.8, -1),
      legend = names(cols),
      pt.bg = cols,
      pch = 21,
      ncol = 3,
      x.intersp = 0.3,
      y.intersp = ifelse(is.null(list(...)$y.intersp), ifelse(is.null(fn), 0.15, 1), list(...)$y.intersp),
      text.width = 0.35,
      bty = 'n'
    )
  }
  
  if(plot.de){
    if(length(contrast) == 1){
      colby = ds2@colData[ , contrast ] %>% levels()
    } else {
      colby = paste0(contrast[2], contrast[-(1:2)])
    }
    colby = paste(colby, collapse = '/')
    
    plot(
      gr,
      layout = lo,
      vertex.label = vertex.names,
      vertex.label.cex = 0.2,
      vertex.label.color = 'green4',
      vertex.color = colfunc( lfcs[V(gr)$name, 'lfc'] ),
      vertex.size = 2 + log1p(lfcs[V(gr)$name, 'abund']),
      main = sprintf('color by %s, size by mean abundance', colby),
      edge.color = 'gray'
    )
    
    lfc_seq = lfc_range[1]:lfc_range[2]
    heatlbls = seq(-5, 5, l = length(lfc_seq))
    heatcol = as.raster(matrix(colfunc(lfc_seq), nrow = 1))
    
    rasterImage(heatcol, -0.5, -1.2, 0.5, -1.1)
    text(
      x = seq(-0.5, 0.5, l = length(lfc_seq)),
      y=-1.25,
      labels = heatlbls
    )
    
    heatlbls = rep(NA, length(heatlbls))
    tcont = sub(
      paste0('log2 fold change \\(MLE\\): ', contrast[1], ' '),
      '',
      ds2_res@elementMetadata$description[2], perl = TRUE
    ) %>% 
      strsplit(' ') %>% unlist %>% rev
    
    heatlbls[1] = tcont[1]
    heatlbls[length(heatlbls)] = tcont[3]
    heatlbls[which(lfc_seq == 0)] = 'Log2FC'
    
    text(
      x = seq(-0.5, 0.5, l = length(lfc_seq)),
      y=ifelse(is.null(fn), -1.35, -1.3),
      labels = heatlbls,
      cex = 1.2
    )
    
    szs = lfcs[V(gr)$name, 'abund']
    szsl = log1p(szs)
    szcat = sprintf('%.0f', seq(min(szsl), max(szsl), by = (max(szsl) - min(szsl))/4) %>% exp(.) - 1)
    legend(
      x=0.6,
      y=ifelse(is.null(fn), -0.8, -0.9),
      legend = szcat,
      pt.cex = (3+log1p(as.integer(szcat)))/3,
      col='black',
      pch=21,
      x.intersp = ifelse(is.null(fn), 0.3, 1.5),
      y.intersp = ifelse(is.null(fn), 0.3, 1.5),
      text.width = 0.1,
      bty = 'n',
      xpd = TRUE
    )
  }
  
  if(!is.null(fn)){
    dev.off()
  }
  
  if(return.results){
    return(list(graph = gr, layout = lo, deseq = ds2, deseq_res = ds2_res))
  }
}

estimate_deseq2 <- function(ph, contrast){
  library(DESeq2)
  library(phyloseq)
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  ds2 = phyloseq_to_deseq2(ph, formula(paste0('~', contrast)))
  geoMeans = apply(counts(ds2), 1, gm_mean)
  ds2 = estimateSizeFactors(ds2, geoMeans = geoMeans)
  return(DESeq(ds2, fitType="local"))
}
