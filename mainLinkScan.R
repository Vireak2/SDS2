library(linkcomm)
library(modMax)
library(network)
library(plot.matrix)
library(KneeArrower)
library(akmedoids)


### Test Linkscan ###

sample_network = RandomSBM(number_nodes=100,
                           mode = 'Assortative',
                           number_blocks = 2)

plot(sample_network$adjacency)
igraph = sample_network$graph


LS_clustering = LinkScan(sample_network,
                         mu = 0.7,
                         n_eps = 20,
                         n_critic = 4,
                         poly_approx = FALSE,
                         star= TRUE,
                         alpha = 0,
                         beta = 1)

LS_clustering$epsilon
# Plot of the edges communities
linkscan_membership = labelize(sample_network$node_edge, LS_clustering$membership, mode='multi')
plot(sample_network$graph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=linkscan_membership,
     vertex.size=10,
     edge.width=1,
     edge.color=LS_clustering$membership,
     main = 'mixed membership linkscan')


# plot the adjacency matrix
colored = Color_Adjacency(sample_network$adjacency, sample_network$edges, LS_clustering$membership)
plot(colored, main = 'LinkScan', col = topo.colors)


expected = expected_communities(sample_network, TRUE)
norm(expected - number_communities(linkscan_membership), type="2")
