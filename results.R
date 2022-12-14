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
library(CliquePercolation)

#min 20 edges
#misclass on vertices

#K = 2,5,10
number_nodes = 100
sample_network = RandomSBM(number_nodes,
                           mode = 'Assortative', #Try mode = 'Disassortative', 'Mixed', 'Overlap'
                           number_blocks = 3,
                           linegraphs = FALSE)

sample_network = readRDS("networks/2K_assortative.rds")
#saveRDS(sample_network, file = "networks/10K_mixed.rds")
plot_adjacency(sample_network$colored_adjacency)

igraph = sample_network$graph
line_graph = line.graph(sample_network$graph)
N_misclass = 8


### LINKCOMM ###
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

membership_linkcomm = memberize_linkcomm(lm$clusters, sample_network$number_edges)
mixed_membership_linkcomm = labelize(sample_network$node_edge, membership_linkcomm, mode = 'multi')
colored_linkcomm = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcomm)
plot_adjacency(colored_linkcomm)

### Linkcomm mod ###
bestlm = bestLinkcomm(lm, 50, line_graph, sample_network)
plot(bestlm, type = "graph")
bestlm$pdmax
membership_linkcommmod = memberize_linkcomm(bestlm$clusters, sample_network$number_edges)
mixed_membership_linkcommmod = labelize(sample_network$node_edge, membership_linkcommmod, mode = 'multi')
colored_linkcommmod = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcommmod)
plot_adjacency(colored_linkcommmod)



### Clique Percolation ###
cp = cpAlgorithm(sample_network$adjacency, k=6, method = 'unweighted')
membership_cp = memberize_linkcomm(cp$list.of.communities.labels, sample_network$number_nodes)




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
membership_c = as.double(membership(cc))
membership_d = as.double(membership(cd))
membership_e = as.double(membership(ce))


# Adjacency plots
colored_c = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership_c)
colored_d = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership_d)
colored_e = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            membership_e)

plot_adjacency(colored_c)
plot_adjacency(colored_d)
plot_adjacency(colored_e)

rm(colored_c); rm(colored_d); rm(colored_e)

# mixed membership
mixed_membership_c = labelize(sample_network$node_edge, membership_c, mode='multi')
mixed_membership_d = labelize(sample_network$node_edge, membership_d, mode='multi')
mixed_membership_e = labelize(sample_network$node_edge, membership_e, mode='multi')

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

LS_clustering$epsilon
membershipLS = LS_clustering$membership
# Plot Linkscan adjacency
colored_linkscan = Color_Adjacency(sample_network$adjacency, sample_network$edges, membershipLS)
plot_adjacency(colored_linkscan)
linkscan_membership = labelize(sample_network$node_edge, membershipLS, mode='multi')

# Plot graph
plot(sample_network$graph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=linkscan_membership,
     vertex.size=4,
     edge.width=1,
     edge.color=membershipLS,
     main = 'mixed membership linkscan')


# Save the clusterings
memberships = data.frame(membership_linkcomm,
                         membership_linkcommmod,
                         membership_c,
                         membership_d,
                         membership_e,
                         membershipLS)
write.csv(memberships, paste("networks/",".rds", sep = name), row.names=FALSE)
rm(memberships)
### Comparisons ###


# Modularity
mod_linkcomm = modularity(line_graph, membership_linkcomm+1)
mod_linkcommmod = modularity(line_graph, membership_linkcommmod+1)
mod_c = modularity(line_graph, membership_c)
mod_d = modularity(line_graph, membership_d)
mod_e = modularity(line_graph, membership_e)
mod_linkscan = modularity(line_graph, membershipLS)
modularities = c(mod_linkcomm,mod_linkcommmod, mod_c, mod_d, mod_e, mod_linkscan) 

# Expected
expected = expected_communities(sample_network, TRUE)
exp_linkcomm = norm(expected - number_communities(mixed_membership_linkcomm), type="2")
exp_linkcommmod = norm(expected - number_communities(mixed_membership_linkcommmod), type="2")
exp_c = norm(expected - number_communities(mixed_membership_c), type="2")
exp_d = norm(expected - number_communities(mixed_membership_d), type="2")
exp_e = norm(expected - number_communities(mixed_membership_e), type="2")
exp_linkscan = norm(expected - number_communities(linkscan_membership), type="2")
expectations = c(exp_linkcomm,exp_linkcommmod, exp_c, exp_d, exp_e, exp_linkscan) 

# misclassification approx
misclass_linkcomm = misclassification_edges(sample_network, colored_linkcomm, N_misclass)[[2]]
misclass_linkcommmod = misclassification_edges(sample_network, colored_linkcommmod, N_misclass)[[2]]
misclass_c = misclassification_edges(sample_network, colored_c, N_misclass)[[2]]
misclass_d = misclassification_edges(sample_network, colored_d, N_misclass)[[2]]
misclass_e = misclassification_edges(sample_network, colored_e, N_misclass)[[2]]
misclass_linkscan = misclassification_edges(sample_network, colored_linkscan, N_misclass)[[2]]
misclassifications = c(misclass_linkcomm,misclass_linkcommmod, misclass_c, misclass_d, misclass_e, misclass_linkscan)

# misclassifications on vertex using diagonal prediction
v_mis_LC = misclassification_vertex(sample_network, colored_linkcomm)
v_mis_LCmod = misclassification_vertex(sample_network, colored_linkcommmod)
v_mis_c = misclassification_vertex(sample_network, colored_c)
v_mis_d = misclassification_vertex(sample_network, colored_d)
v_mis_e = misclassification_vertex(sample_network, colored_e)
v_mis_LS = misclassification_vertex(sample_network, colored_linkscan)
v_missclassifications=c(v_mis_LC, v_mis_LCmod,v_mis_c, v_mis_d, v_mis_e, v_mis_LS) 

measurement = data.frame(Methods = c('Linkcomm','Linkcomm + mod', 'Imperial C','Imperial D', 'Imperial E', 'LinkScan*'),
                         Modularity = modularities,
                         Expectation_L2error = expectations,
                         Approx_misclassification = misclassifications,
                         vertex_misclassifications = v_missclassifications)
measurement

saveRDS(measurement, file = "measures/10K_mixed.rds")
# In case
#saveRDS(measurement, file = "10K_mixed_measures.rds")

y = read.csv("measures_csv/10K_mixed.csv")
xtable(y)
