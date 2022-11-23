library(linkcomm)
library(modMax)
library(network)
library(plot.matrix)
library(KneeArrower)
library(xtable)
library(combinat)
# Test with undirected network
#n= 256
#K = 2,5,10
sample_network = RandomSBM(number_nodes=256,
                           mode = 'Disassortative', #Try mode = 'Disassortative', 'Mixed', 'Overlap'
                           number_blocks = 2)
#sample_network = readRDS("256n2K_assortative.rds")
#saveRDS(sample_network, file = "256n2K_disassortative.rds")
colored_adj = sample_network$colored_adjacency
plot(colored_adj,
     col = topo.colors,
     main = 'Adjacency matrix')
igraph = sample_network$graph



### Link community detection with link density ###
lm <- getLinkCommunities(sample_network$edges,  hcmethod = 'single')
lm
plot(lm, type = "graph")

# Node community centrality scores:
node_centrality <- getCommunityCentrality(lm)
node_centrality
# Community connectedness scores:
corectness <- getCommunityConnectedness(lm)
corectness
lm$pdmax

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

# Further plots for E
colored_c = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(cc))
colored_d = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(cd))
colored_e = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(ce))
plot(colored_c, col = topo.colors)
plot(colored_d, col = topo.colors)
plot(colored_e, col = topo.colors)


# mixed membership
mixed_membership_c = labelize(sample_network$node_edge, membership(cc), mode='multi')
mixed_membership_d = labelize(sample_network$node_edge, membership(cd), mode='multi')
mixed_membership_e = labelize(sample_network$node_edge, membership(ce), mode='multi')

plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_c,
     vertex.size=10,
     edge.width=1,
     edge.color=membership(cc),
     col = topo.colors,
     main = 'mixed membership C method')

plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_d,
     vertex.size=10,
     edge.width=1,
     edge.color=membership(cd),
     col = topo.colors,
     main = 'mixed membership D method')

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






expected = expected_communities(sample_network, TRUE)
norm(expected - number_communities(mixed_membership_c), type="2")
norm(expected - number_communities(mixed_membership_d), type="2")
norm(expected - number_communities(mixed_membership_e), type="2")



#misclassification
misclassification_edges(sample_network, colored_e)





#Tiebow
tiebow=matrix(c(0,1,1,0,0, 1,0,1,0,0, 1,1,0,1,1, 0,0,1,0,1, 0,0,1,1,0), 5, 5)
tiebow_edges = edges_from_adjacency(tiebow)
tiebow_igraph = graph_from_adjacency_matrix(tiebow, 'undirected')
plot(tiebow_igraph,
     vertex.size = 30,
     edge.width=7,
     vertex.color = 'white')
tiebow_linkcomm = getLinkCommunities(tiebow_edges,  hcmethod = 'single')
plot(tiebow_linkcomm, type='graph')


#dendrogram
tiebow_jaccard = Link_similarities(tiebow_igraph,selfloops = FALSE)$LS_weights
edgesnames= c('(1,2)', '(1,3)', '(2,3)', '(3,4)', '(3,5)', '(4,5)')
tiebow_jaccard = data.frame(tiebow_jaccard, row.names = edgesnames)
colnames(tiebow_jaccard) = edgesnames
as.dist(tiebow_jaccard)
tiebow_dist = as.dist(-tiebow_jaccard)
plot(hclust(tiebow_dist),
     main = '',
     axes = FALSE)

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



