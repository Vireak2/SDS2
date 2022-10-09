library(linkcomm)
library(modMax)

# Test with undirected network
sample_network = RandomSBM(number_nodes=25, mode = 'Disassortative', number_blocks = 3) # Set mode = 'Assortative', 'Disassortative', 'Mixed'
heatmap(sample_network$adjacency, Colv = NA, Rowv = NA,)
igraph = sample_network$graph



### Link community detection with link density ###
lm <- getLinkCommunities(sample_network$edges)  # a tendance à bloquer a D = 1
# Couper le dendrogram là ou la modularité E est maximisée? newLinkCommsAt(lm, cutat = 0.5)
lm
plot(lm, type = "graph")
# Node community centrality scores:
node_centrality <- getCommunityCentrality(lm)
node_centrality
# Community connectedness scores:
corectness <- getCommunityConnectedness(lm)
corectness




### Node clustering
c = cluster_louvain(igraph)
membership(c)
sizes(c)
plot(igraph, vertex.label=NA, vertex.color=membership(c), vertex.size=10)




### Link community detection with line graph modularity (C, D and E)###
lc = graph_from_adjacency_matrix(sample_network$line_graph, mode = 'undirected')
ld = graph_from_adjacency_matrix(sample_network$D_line_graph, 'undirected', weighted = TRUE)
le = graph_from_adjacency_matrix(sample_network$E1_line_graph, 'undirected', weighted = TRUE)
cc = cluster_louvain(lc)
cd = cluster_louvain(ld, weights = E(ld)$weight)
ce = cluster_louvain(le, weights = E(le)$weight)
membership(cc)
membership(cd)
membership(ce)
sizes(cc)
sizes(cd)
sizes(ce)
plot(igraph, vertex.label=NA, edge.color=membership(cc), vertex.size=3, edge.width=2,
     main = 'link communities with C method', vertex.color = 'gray')
plot(igraph, vertex.label=NA, edge.color=membership(cd), vertex.size=3, edge.width=2,
     main = 'link communities with D method', vertex.color = 'gray')
plot(igraph, vertex.label=NA, edge.color=membership(ce), vertex.size=3, edge.width=2,
     main = 'link communities with E method', vertex.color = 'gray')

# E method tends to find less communities

