library(linkcomm)
library(modMax)
library(network)
library(plot.matrix)
library(KneeArrower)


# Test with undirected network
sample_network = RandomSBM(number_nodes=100,
                           mode = 'Assortative', #Try mode = 'Disassortative', 'Mixed'
                           number_blocks = 3)
plot(sample_network$adjacency)
igraph = sample_network$graph



### Link community detection with link density ###
lm <- getLinkCommunities(sample_network$edges)
lm
plot(lm, type = "graph")

# Node community centrality scores:
node_centrality <- getCommunityCentrality(lm)
node_centrality
# Community connectedness scores:
corectness <- getCommunityConnectedness(lm)
corectness


# Often lm cut a 1 and gives no information, new_lm cut a bit lower and gives a better result
new_lm = newLinkCommsAt(lm, cutat = lm$pdmax*.99)
new_lm
plot(new_lm, type = "graph")




### Node clustering
c = cluster_louvain(igraph)
membership(c)
sizes(c)
plot(igraph, vertex.label=NA, vertex.color=membership(c), vertex.size=10)




### Link community detection with line graph modularity (C, D and E)###

# igraphs objects with adjacency=C, D and E1
lc = graph_from_adjacency_matrix(sample_network$line_graph, mode = 'undirected')
ld = graph_from_adjacency_matrix(sample_network$D_line_graph, 'undirected', weighted = TRUE)
le = graph_from_adjacency_matrix(sample_network$E1_line_graph, 'undirected', weighted = TRUE)

# modularity maximisation on the new line graphs
cc = cluster_louvain(lc)
cd = cluster_louvain(ld, weights = E(ld)$weight)
ce = cluster_louvain(le, weights = E(le)$weight)

# memberships
membership(cc)
membership(cd)
membership(ce)

# plots of colored edges
plot(igraph, vertex.label=NA, edge.color=membership(cc), vertex.size=3, edge.width=2,
     main = 'link communities with C method', vertex.color = 'gray')
plot(igraph, vertex.label=NA, edge.color=membership(cd), vertex.size=3, edge.width=2,
     main = 'link communities with D method', vertex.color = 'gray')
plot(igraph, vertex.label=NA, edge.color=membership(ce), vertex.size=3, edge.width=2,
     main = 'link communities with E method', vertex.color = 'gray')


# Further plots for E
v_membership_e = labelize(sample_network$node_edge, membership(ce), mode='most')
plot(igraph,
     vertex.label=NA,
     vertex.color=v_membership_e,
     vertex.size=10,
     edge.width=1,
     edge.color=membership(ce),
     main = 'vertex clustering E method')
colored_e = Color_Adjacency(sample_network$adjacency, sample_network$edges,
                            membership(ce))
plot(colored_e, col = topo.colors)


# mixed membership plot
mixed_membership_e = labelize(sample_network$node_edge, membership(ce), mode='multi')
plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_e,
     vertex.size=10,
     edge.width=1,
     edge.color=membership(ce),
     col = topo.colors,
     main = 'mixed membership E method')




# compare modularity
v_membership_linkcomm = as.integer(new_lm$nodeclusters[,2]) # did not implement the modularity on lincomm yet
v_membership_c = labelize(sample_network$node_edge, membership(cc), mode='most')
v_membership_d = labelize(sample_network$node_edge, membership(cd), mode='most')

modularity(igraph, membership(c)) # Benchmark
modularity(igraph, v_membership_c)
modularity(igraph, v_membership_d)
modularity(igraph, v_membership_e)











# plot for overlaping
bp = matrix(c(1, .2, 0, .2,0.4,.2, 0,.2,1), 3,3)
overlap = SBM(c(5,5,5), bp)
igraph = overlap$graph
c = cluster_louvain(igraph)
membership(c)
sizes(c)
plot(igraph, vertex.label=NA, vertex.color=membership(c), vertex.size=10)
# plot for nested
np =  matrix(c(.5, .5, 0.1, .5,1,.1, .1,.1,.3), 3,3)
nested = SBM(c(5,5,10), np)
igraph = nested$graph



