library(linkcomm)
library(modMax)
library(network)
library(plot.matrix)
library(KneeArrower)
library(xtable)
library(combinat)
library(LaplacesDemon)
library(ggplot2)
library(RColorBrewer)
library(Matrix)

#n= 256
#K = 2,5,10
#sample_network = RandomSBM(number_nodes=100,
#                           mode = 'Overlap', #Try mode = 'Disassortative', 'Mixed', 'Overlap'
#                           number_blocks = 3,
#                           linegraphs = FALSE)

sample_network = readRDS("SDS2/256n2K_assortative.rds")
#saveRDS(sample_network, file = "256n4K_overlap.rds")

plot_adjacency(sample_network$colored_adjacency)

igraph = sample_network$graph
line_graph = line.graph(sample_network$graph)
N_misclass = 8

### LINKCOMM ###
lm <- getLinkCommunities(sample_network$edges,  hcmethod = 'single')
lm
plot(lm, type = "graph",)

# Node community centrality scores:
node_centrality <- getCommunityCentrality(lm)
node_centrality
# Community connectedness scores:
corectness <- getCommunityConnectedness(lm)
corectness
lm$pdmax

membership_linkcomm = memberize_linkcomm(lm$clusters, sample_network$number_edges)
mixed_membership_linkcomm = labelize(sample_network$node_edge, membership_linkcomm, mode = 'multi')
colored_linkcomm = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcomm)
plot_adjacency(colored_linkcomm)

### Linkcomm mod ###
bestlm = bestLinkcomm(lm, 25, line_graph, sample_network)
plot(bestlm, type = "graph")
bestlm$pdmax
membership_linkcommmod = memberize_linkcomm(bestlm$clusters, sample_network$number_edges)
mixed_membership_linkcommmod = labelize(sample_network$node_edge, membership_linkcommmod, mode = 'multi')
colored_linkcommmod = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcommmod)
plot_adjacency(colored_linkcommmod)

### IMPERIAL ###
# igraphs objects with adjacency = C, D and E1
lc = graph_from_adjacency_matrix(get_linegraph(sample_network, type = 'C'), mode = 'undirected')
cc = cluster_louvain(lc)
rm(lc)


ld = graph_from_adjacency_matrix(get_linegraph(sample_network, type = 'D'), 'undirected', weighted = TRUE)
cd = cluster_louvain(ld, weights = E(ld)$weight)
rm(ld)


le = graph_from_adjacency_matrix(get_linegraph(sample_network, type = 'E'), 'undirected', weighted = TRUE)
ce = cluster_louvain(le, weights = E(le)$weight)
rm(le)


n_cc = max(membership(cc))
n_cd = max(membership(cd))
n_ce = max(membership(ce))

# Adjacency plots
colored_c = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(cc))
colored_d = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(cd))
colored_e = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership(ce))

plot_adjacency(colored_c)
plot_adjacency(colored_d)
plot_adjacency(colored_e)

misclass_c = misclassification_edges(sample_network, colored_c, N_misclass)[[2]]
misclass_d = misclassification_edges(sample_network, colored_d, N_misclass)[[2]]
misclass_e = misclassification_edges(sample_network, colored_e, N_misclass)[[2]]

rm(colored_c); rm(colored_d); rm(colored_e)

# mixed membership
mixed_membership_c = labelize(sample_network$node_edge, membership(cc), mode='multi')
mixed_membership_d = labelize(sample_network$node_edge, membership(cd), mode='multi')
mixed_membership_e = labelize(sample_network$node_edge, membership(ce), mode='multi')

plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_c,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(cc),
     main = 'mixed membership C method')

plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_d,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(cd),
     main = 'mixed membership D method')

plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_e,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(ce),
     main = 'mixed membership E method')


### LINKSCAN* ###
LS_clustering = LinkScan(sample_network,
                         mu = 0.7,
                         n_eps = 20,
                         n_critic = 5,
                         poly_approx = FALSE,
                         star= TRUE,
                         alpha = 0,
                         beta = 1)

epsilon = LS_clustering$epsilon

# Plot Linkscan adjacency
colored_linkscan = Color_Adjacency(sample_network$adjacency, sample_network$edges, LS_clustering$membership)
plot_adjacency(colored_linkscan)

# Plot graph
linkscan_membership = labelize(sample_network$node_edge, LS_clustering$membership, mode='multi')
plot(sample_network$graph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=linkscan_membership,
     vertex.size=1,
     edge.width=1,
     edge.color=LS_clustering$membership,
     main = 'mixed membership linkscan')






### Comparisons ###


# Modularity
mod_linkcomm = modularity(line_graph, membership_linkcomm+1)
mod_linkcommmod = modularity(line_graph, membership_linkcommmod+1)
mod_c = modularity(line_graph, membership(cc))
mod_d = modularity(line_graph, membership(cd))
mod_e = modularity(line_graph, membership(ce))
mod_linkscan = modularity(line_graph, LS_clustering$membership)
modularities = c(mod_linkcomm, mod_c, mod_d, mod_e, mod_linkscan) 

# Expected
expected = expected_communities(sample_network, TRUE)
exp_linkcomm = norm(expected - number_communities(mixed_membership_linkcomm), type="2")
exp_linkcommmod = norm(expected - number_communities(mixed_membership_linkcommmod), type="2")
exp_c = norm(expected - number_communities(mixed_membership_c), type="2")
exp_d = norm(expected - number_communities(mixed_membership_d), type="2")
exp_e = norm(expected - number_communities(mixed_membership_e), type="2")
exp_linkscan = norm(expected - number_communities(linkscan_membership), type="2")
expectations = c(exp_linkcomm, exp_c, exp_d, exp_e, exp_linkscan) 

# misclassification approx
misclass_linkcomm = misclassification_edges(sample_network, colored_linkcomm, N_misclass)[[2]]
misclass_linkcommmod = misclassification_edges(sample_network, colored_linkcommmod, N_misclass)[[2]]
misclass_c
misclass_d
misclass_e
misclass_linkscan = misclassification_edges(sample_network, colored_linkscan, N_misclass)[[2]]
missclassifications = c(misclass_linkcomm, misclass_c, misclass_d, misclass_e, misclass_linkscan)

measurement = data.frame(Methods = c('Linkcomm', 'Imperial C','Imperial D', 'Imperial E', 'LinkScan'),
                         Modularity = modularities,
                         Expectation_L2error = expectations,
                         Approx_misclassification = missclassifications)
measurement
#saveRDS(measurement, file = "SDS2/10K_mixed_measures.rds")
# In case
#saveRDS(measurement, file = "10K_mixed_measures.rds")

c(mod_linkcommmod, exp_linkcommmod, misclass_linkcommmod)

y = readRDS("SDS2/10K_mixed_measures.rds")
