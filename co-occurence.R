# Prepare datasets to infer correlations by fastspar
#
# NB: Drop samples 43-46 first as possible outliers
#
# 1. Generate OTU-count tables for different groups of samples and amplicon types
# a) 16S FL/PA vs Apr/Mar/Jun :
#   - FL-Apr, FL-Mar, FL-Jun
#   - PA-Apr, PA-Mar, PA-Jun
#   - FL-Apr-Mar and PA-Apr-Mar combined
# b) 18S Apr vs Mar vs Jun
#   - Apr, May, Jun
#   - Apr-Mar combined
#   

source('/media/repo4T/yuragal_home/LIN/metagenome/Likhoshway_post_anaysis/bashenhaeva/mol_a/data/Zi-LN/correllation_wrapper_functions.R')
source('~/Desktop/LIN/metagenome/scripts_functions/analysis_functions.R')
ipak()

tph = ph$bact$raw %>% subset_samples(!id %in% 43:46)
sample_data(tph)$month %<>% relevel(., ref = 'Mar')
corr_input = filter_input_correllation_data(tph, mta = 50, swpr = 0.25)

# filter_stats = cbind(raw = sample_sums(tph), filtered = rowSums(corr_input$counts)) %>% as.data.frame()
# filter_stats$missed_pct = (1-filter_stats$filtered/filter_stats$raw)*100
# 
# View(filter_stats)

coocc = list()

coocc$bact = compute_corr_graphs(corr_input, target_edge_number = 1200)
#corr_graphs.1000 = compute_corr_graphs(corr_input, target_edge_number = 1000, raw_corrs = coocc$bactcorrs)
par(mfrow=c(1,1))
draw_venn(coocc$bact$raw_graphs)

library(igraph)

###########################################
# Compare co-occurrence graphs generated by all methods
# list of co-occurrence graphs -- CGS:
###########################################
# se -- SpiecEasi
# mg -- MAGMA
# zi -- Zi-LN
# sc -- SparcCC
# fw -- FlashWeave

CGS = lapply(names(coocc$bact$raw_graphs), function(x) {
  z = coocc$bact$raw_graphs[[x]] != 0
  g = graph.adjacency(z, mode = 'undirected')
  V(g)$name = colnames(corr_input$counts)
  return(
    list(
      mx = z,
      gr = g,
      coords = layout_with_fr(g)
    )
  )
})

names(CGS) = names(coocc$bact$raw_graphs)
coocc$bact$binary_graphs = CGS

#CGS = coocc$bact$binary_graphs

apply(count_shared_edges(CGS, r = 3), 1, function(x) cat(x,'\n'))
# fw mg sc 59 
# fw mg se 205 
# fw mg zi 249 
# fw sc se 66 
# fw sc zi 73 
# fw se zi 242 
# mg sc se 175 
# mg sc zi 146 
# mg se zi 494 <-- max
# sc se zi 157 

# plot(
#   CGS$zi$gr - which(degree(CGS$zi$gr) == 0),
#   #CGS$zi$gr - edges_by_weight[[2]] - edges_by_weight[[3]] - edges_by_weight[[4]] - edges_by_weight[[5]],
#   #CGS$zi$gr - edges_by_weight[[3]] - edges_by_weight[[4]] - edges_by_weight[[5]],
#   layout=CGS$zi$coords[which(degree(CGS$zi$gr) > 0),],
#   # ceb,
#   # gr,
#   # layout = lo,
#   vertex.label=NA,
#   vertex.color = labels$colors$V[corr_input$taxa[V(CGS$zi$gr - which(degree(CGS$zi$gr) == 0)),'Phylum']],
#   #vertex.color = membership(ceb),
#   vertex.size=3,
#   #edge.color = ifelse(E(CGS[[x]]) %in% grays, 'red', 'gray'),
#   #edge.color = labels$colors$E[labels$corr],
#   main = 'Zi-LN CEB',
#   edge.color = 'gray'
# )

# Calculate graph, layout, and deseq
graph_cc_merged = merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input,
  raw_counts = tph,
  #ds2 = ds2$all,
  contrast = 'lifestyle',
  lo = graph_cc_merged$layout,
  show.vertex.names = FALSE,
  return.results = TRUE,
  plot.phylotypes = TRUE
)

# save merged graph
# layout is important as one need to use constant layout
# when comparing differential abundance patterns
# coocc$bact$merged_mg_se_zi = graph_cc_merged

graph_cc_merged = coocc$bact$merged_mg_se_zi

# Use the calculated results to redraw a graph
merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input,
  ds2 = graph_cc_merged$deseq,
  contrast = 'lifestyle',
  show.vertex.names = FALSE,
  lo = graph_cc_merged$layout,
  plot.phylotypes = TRUE,
  #fn = '16S_co-occurence_and_diff_abundance_with_OTU_ids.svg'
  fn = '16S_co-occurence_and_diff_abundance.svg'
)


# plot FL/PA contrast by month
sample_data(tph)$group = factor(paste0(sample_data(tph)$month, sample_data(tph)$lifestyle))
graph_cc_merged_by_month = merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input,
  raw_counts = tph,
  contrast = c('group','Apr','PA','FL'),
  show.vertex.names = FALSE,
  lo = graph_cc_merged$layout,
  return.results = TRUE,
  plot.phylotypes = FALSE
)

svg('16S_co-occurence_and_diff_abundance_by_month.svg', width = 18, height = 6)
par(mfrow=c(1,3))

sapply(graph_cc_merged_by_month$deseq@colData$month %>% levels, function(x){
  merge_graphs_and_draw_assort_contrast(
    corr_graphs = CGS,
    corr_input = corr_input,
    ds2 = graph_cc_merged_by_month$deseq,
    contrast = c('group',x,'PA','FL'),
    show.vertex.names = FALSE,
    lo = graph_cc_merged$layout,
    return.results = FALSE,
    plot.phylotypes = FALSE
  )
})
dev.off()

# plot FL/PA contrast by month
sample_data(tph)$group = factor(paste0(sample_data(tph)$lifestyle, sample_data(tph)$month))
graph_cc_merged_by_month_FL = merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input,
  raw_counts = tph,
  contrast = c('group','FL','Apr','Jul'),
  show.vertex.names = FALSE,
  lo = graph_cc_merged$layout,
  return.results = TRUE,
  plot.phylotypes = FALSE,
  new.plot = FALSE
)

svg('16S_co-occurence_and_diff_abundance_by_month_and_lifestyle1.svg', width = 18, height = 12)
par(mfrow=c(2,3))

apply(gtools::combinations(1:3, n=3, r=2), 1, function(x){
  mlvl = sample_data(tph)$month %>% levels
  merge_graphs_and_draw_assort_contrast(
    corr_graphs = CGS,
    corr_input = corr_input,
    ds2 = graph_cc_merged_by_month_FL$deseq,
    contrast = c('group','FL', mlvl[ x[2] ], mlvl[ x[1] ]),
    show.vertex.names = FALSE,
    lo = graph_cc_merged$layout,
    return.results = FALSE,
    plot.phylotypes = FALSE,
    plt = 'PiYG',
    new.plot = FALSE,
    plot.de = TRUE
  )
})


apply(gtools::combinations(1:3, n=3, r=2), 1, function(x){
  mlvl = sample_data(tph)$month %>% levels
  merge_graphs_and_draw_assort_contrast(
    corr_graphs = CGS,
    corr_input = corr_input,
    ds2 = graph_cc_merged_by_month_FL$deseq,
    contrast = c('group','PA', mlvl[ x[2] ], mlvl[ x[1] ]),
    show.vertex.names = FALSE,
    lo = graph_cc_merged$layout,
    return.results = FALSE,
    plot.phylotypes = FALSE,
    plt = 'PiYG',
    new.plot = FALSE,
    plot.de = TRUE
  )
})
dev.off()

# sapply(3:5, function(y){
#   #cl = c('magenta','green','blue', 'red')
#   cl = c('green','blue', 'red')
#   plot(
#     #CGS$zi$gr %s% edges_by_weight[[y]],
#     #layout=CGS$zi$coords,
#     gr %s% edges_by_weight[[y]],
#     vertex.label=NA,
#     #vertex.color = labels$colors$V[labels$taxa],
#     vertex.color = labels$colors$V[corr_input$taxa[V(gr),'Phylum']],
#     vertex.size=3,
#     #edge.color = ifelse(E(CGS[[x]]) %in% grays, 'red', 'gray'),
#     #edge.color = labels$colors$E[labels$corr],
#     #main = 'Zi-LN',
#     #edge.color = add.alpha(cl[y-1], 1),
#     edge.color = add.alpha(cl[y-2], 1),
#     add = TRUE
#     #main = x
#   )
# })
# 
# legend(x=-1.2, y=-1, legend=names(cols), pt.bg = cols, pch=21, ncol = 3, x.intersp = 0.3, y.intersp = 0.2, text.width = 0.45)
# 
# # build graph with edge weight being a number of co-occurrence methods supporting the edge
# g2u = graph.adjacency(CGS$se$mx + CGS$mg$mx + CGS$zi$mx + CGS$sc$mx + CGS$fw$mx, mode = 'undirected',  weighted = TRUE)
# V(g2u)$name = colnames(corr_input$counts)
# #gall = g2u - E(g2u)[[weight < 5]]
# 
# edges_by_weight = lapply(1:5, function(x) {
#   g2u - E(g2u)[[weight != x]]
# })
# 
# labels = list()
# labels$taxa = corr_input$taxa[,'Phylum']
# #labels$colors$V = viridis::inferno(labels$taxa %>% unique %>% length)
# labels$colors$V = pals::kelly(labels$taxa %>% unique %>% length)
# names(labels$colors$V) = labels$taxa %>% unique
# #labels$colors$E = c('gray','red')
# #labels$corr = as.numeric(E(g1) %in% E(g2a))+1
# 
# # Draw co-occurrence graphs generated by all methods 
# # with edges colored by number of methods supporting the edge:
# # red   -- supported by all methods
# # blue  -- supported by any four
# # green -- supported by any three
# # gray  -- any other edges inferred by particular method
# #
# # Graph layout is constant and optimal for ZI-LN
# 
# par(mfrow=c(2,3))
# 
# sapply(names(CGS), function(x) {
#   #bgr = CGS[[x]]$gr
#   bgr = CGS[[x]]$gr - edges_by_weight[[1]] - edges_by_weight[[2]] - edges_by_weight[[3]]
#   # ivs = which(degree(bgr) == 0)
#   # bgr = delete.vertices(bgr, ivs)
#   # lo = layout_with_fr(bgr)
#   
#   #bgr = CGS[[x]]$gr
#   plot(
#     #CGS$zi,
#     bgr - edges_by_weight[[3]] - edges_by_weight[[4]] - edges_by_weight[[5]],
#     #layout=CGS$zi$coords,
#     layout=CGS[[x]]$coords,
#     #layout=lo,
#     vertex.label=NA,
#     vertex.color = labels$colors$V[corr_input$taxa[V(bgr),'Phylum']],
#     #vertex.color = labels$colors$V[labels$taxa],
#     vertex.size=3,
#     #edge.color = ifelse(E(CGS[[x]]) %in% grays, 'red', 'gray'),
#     #edge.color = labels$colors$E[labels$corr],
#     #main = 'Zi-LN',
#     edge.color = 'gray',
#     main = x
#   )
#   
#   sapply(3:5, function(y){
#     cl = c('green','blue', 'red')
#     plot(
#       bgr %s% edges_by_weight[[y]],
#       layout=CGS[[x]]$coords,
#       #layout=CGS$zi$coords,
#       #layout=lo,
#       vertex.label=NA,
#       #vertex.color = labels$colors$V[labels$taxa],
#       vertex.color = labels$colors$V[corr_input$taxa[V(bgr),'Phylum']],
#       vertex.size=3,
#       #edge.color = ifelse(E(CGS[[x]]) %in% grays, 'red', 'gray'),
#       #edge.color = labels$colors$E[labels$corr],
#       #main = 'Zi-LN',
#       edge.color = add.alpha(cl[y-2], 1),
#       add = TRUE
#       #main = x
#     )
#   })
# })


gr = coocc$bact$merged_mg_se_zi$graph
lo = coocc$bact$merged_mg_se_zi$layout
ceb = cluster_edge_betweenness(gr)



write_xls(
  cbind(group = membership(ceb), corr_input$taxa[V(graph_cc_merged$graph)$name,]),
  fn = 'co-occurrence.xls',
  sheet = 'Bacteria_tax',
  append = TRUE
)

svg('16S_group20.svg', width = 12, height = 8)
par(mfrow=c(1,2))
sapply(ceb[[20]], function(x) plotCounts(graph_cc_merged_by_month$deseq, x, intgroup = 'group'))
dev.off()

sapply(c('OTU24'), function(x) plotCounts(graph_cc_merged_by_month$deseq, x, intgroup = 'group'))

#############################################################
# Build and analyze co-occurrence of Eukaryotes
#############################################################
tph = ph$euk$raw %>% subset_samples(!id %in% 43:46)
sample_data(tph)$month %<>% relevel(., ref = 'Mar')
corr_input1 = filter_input_correllation_data(tph, mta = 100, swpr = 0.25) #%>% str

clph = tax_glom(tph, 'Class')
taxa_names(clph) = tax_table(clph)[,'Class']

df = data.frame(pch = as.numeric(taxa_sums(clph)))
df[,2:3] = tax_table(clph)[,2:3]
View(df)

cln = c(tax_table(clph)[,'Class'])

coocc$euk1 = compute_corr_graphs(corr_input1, target_edge_number = 500)

par(mfrow=c(1,1))
draw_venn(coocc$bact$raw_graphs)
draw_venn(coocc$euk$raw_graphs)

CGS = lapply(names(coocc$euk1$raw_graphs), function(x) {
  z = coocc$euk1$raw_graphs[[x]] != 0
  g = graph.adjacency(z, mode = 'undirected')
  V(g)$name = colnames(corr_input1$counts)
  return(
    list(
      mx = z,
      gr = g,
      coords =layout_with_fr(g)
    )
  )
})

names(CGS) = names(coocc$euk1$raw_graphs)
# coocc$euk1$binary_graphs = CGS

plot(eigen(as_adjacency_matrix(coocc$euk$merged_mg_se_zi$graph))$values, type="b")
abline(h=1,col="red", lty = 3)

coocc$euk$merged_mg_se_zi$graph

apply(count_shared_edges(CGS, r = 3), 1, function(x) cat(x,'\n'))
# fw mg sc 48 
# fw mg se 104 
# fw mg zi 115 
# fw sc se 56 
# fw sc zi 56 
# fw se zi 111 
# mg sc se 103 
# mg sc zi 109 
# mg se zi 291 <-- max
# sc se zi 104 

gr = coocc$euk$binary_graphs$zi$gr
lo = coocc$euk$binary_graphs$zi$coords
ceb = cluster_edge_betweenness(gr)

gr = coocc$euk1$merged_mg_se_zi$graph
lo = coocc$euk1$merged_mg_se_zi$layout
lo = layout_with_fr(gr)
lo = layout_with_kk(gr)
lo = layout_with_mds(gr)
ceb = cluster_edge_betweenness(gr)
sg = spinglass.community(gr)

merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input,
  ds2 = graph_cc_merged$deseq,
  contrast = 'month',
  show.vertex.names = FALSE,
  lo = lo,
  #lo = graph_cc_merged$layout,
  plot.phylotypes = TRUE,
  y.intersp = 0.25,
  new.plot = FALSE
  #fn = '16S_co-occurence_and_diff_abundance_with_OTU_ids.svg'
  #fn = '16S_co-occurence_and_diff_abundance.svg'
)



par(mfrow=c(1,1))
plot(
  #CGS$zi$gr - which(degree(CGS$zi$gr) == 0),
  #CGS$zi$gr - edges_by_weight[[2]] - edges_by_weight[[3]] - edges_by_weight[[4]] - edges_by_weight[[5]],
  #CGS$zi$gr - edges_by_weight[[3]] - edges_by_weight[[4]] - edges_by_weight[[5]],
  #layout=CGS$zi$coords[which(degree(CGS$zi$gr) > 0),],
  ceb,
  gr,
  layout = lo,
  vertex.label= membership(ceb),
  vertex.label.cex = 0.4,
  vertex.label.color = 'black',
  #vertex.color = labels$colors$V[corr_input$taxa[V(CGS$zi$gr - which(degree(CGS$zi$gr) == 0)),'Phylum']],
  #vertex.color = membership(ceb),
  vertex.size=3,
  vertex.color= membership(ceb),
  #edge.color = ifelse(E(CGS[[x]]) %in% grays, 'red', 'gray'),
  #edge.color = labels$colors$E[labels$corr],
  main = 'mg & se & se CEB',
  edge.color = 'gray'
)

# Calculate graph, layout, and deseq
graph_cc_merged = merge_graphs_and_draw_assort_contrast(
  corr_graphs = CGS,
  corr_input = corr_input1,
  raw_counts = tph,
  #ds2 = ds2$all,
  contrast = 'month',
  #lo = graph_cc_merged$layout,
  show.vertex.names = FALSE,
  return.results = TRUE,
  plot.phylotypes = TRUE,
  plot.de = FALSE
)

# save merged graph
# layout is important as one need to use constant layout
# when comparing differential abundance patterns
# coocc$euk1$merged_mg_se_zi = graph_cc_merged

graph_cc_merged = coocc$euk1$merged_mg_se_zi
# CGS = coocc$euk1$binary_graphs

svg('18S_co-occurence_and_diff_abundance_by_month.svg', width = 12, height = 12)
par(mfrow=c(2,2))

merge_graphs_and_draw_assort_contrast(
  # corr_graphs = CGS,
  corr_graphs = coocc$euk1$binary_graphs,
  corr_input = corr_input1,
  # ds2 = graph_cc_merged$deseq,
  ds2 = coocc$euk1$merged_mg_se_zi$deseq,
  contrast = 'month',
  show.vertex.names = FALSE,
  #lo = lo,
  lo = graph_cc_merged$layout,
  plot.phylotypes = TRUE,
  y.intersp = 0.25,
  new.plot = FALSE,
  vertex_size = 5
  #fn = '16S_co-occurence_and_diff_abundance_with_OTU_ids.svg'
  #fn = '16S_co-occurence_and_diff_abundance.svg'
)

apply(gtools::combinations(1:3, n=3, r=2), 1, function(x){
  mlvl = sample_data(tph)$month %>% levels
  merge_graphs_and_draw_assort_contrast(
    corr_graphs = CGS,
    corr_input = corr_input1,
    ds2 = graph_cc_merged$deseq,
    contrast = c('month','', mlvl[ x[2] ], mlvl[ x[1] ]),
    show.vertex.names = FALSE,
    lo = graph_cc_merged$layout,
    return.results = FALSE,
    plot.phylotypes = FALSE,
    plot.de = TRUE,
    new.plot = FALSE,
    plt = 'PiYG'
  )
})
dev.off()

write_xls(
  cbind(group = membership(ceb), corr_input1$taxa[V(graph_cc_merged$graph)$name,]),
  fn = 'co-occurrence.xls',
  sheet = 'Eukarya',
  #append = TRUE
)

svg('18S_group1.svg', width = 12, height = 4)
par(mfrow=c(2,6))
sapply(ceb[[1]], function(x) plotCounts(graph_cc_merged$deseq, x, intgroup = 'month'))
dev.off()

svg('18S_group4.svg', width = 10, height = 4)
par(mfrow=c(2,5))
sapply(ceb[[4]], function(x) plotCounts(graph_cc_merged$deseq, x, intgroup = 'month'))
dev.off()

svg('18S_group5.svg', width = 12, height = 6)
par(mfrow=c(3,6))
sapply(ceb[[5]], function(x) plotCounts(graph_cc_merged$deseq, x, intgroup = 'month'))
dev.off()

