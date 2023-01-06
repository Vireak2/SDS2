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
library(pdfCluster)
library(NMI)


# This script create a SBM sample_network and proceed to test Linkcomm, the Imperial methods and LinkSCAN
# All plots of graphs and adjacency matrices are saved in the folder name
# All networks are saved in the "Networks" folder
# All edge membership vectors are saved in a dataframe in the "memberships" folder
# All edge measurement are saved in a dataframe in the "measurements" folder

# Network creation
number_nodes = 256
sample_network = RandomSBM(number_nodes,
                           inter = c(0.01, log(number_nodes)*3/number_nodes),
                           intra = c(log(number_nodes)*10/number_nodes, 0.5),
                           mode = 'Overlap', #Try mode = 'Disassortative', 'Mixed', 'Overlap'
                           number_blocks = 3,
                          linegraphs = FALSE)

# Vector of names
cname= c('2K_Assortative', '5K_Assortative', '10K_Assortative',
         '2K_Disa', '5K_Disa',
         '5K_mixed', '10K_mixed',
         '3K_Overlap')

name = cname[] #Choose the name you want

# Save the network in the "Networks" folder
saveRDS(sample_network, file = paste("Networks/",".rds", sep = name))

# Load the network
sample_network = readRDS(paste("Networks/",".rds", sep = name))

igraph = sample_network$graph
line_graph = line.graph(sample_network$graph)

# Save the plot of the adjacency 
png(file=paste(name, '/adj.png', sep=''), width=600, height=600)
plot_adjacency(sample_network$colored_adjacency)
dev.off()


### LINKCOMM ###
# Clustering
png(file=paste(name, '/LC.png', sep=''), width=600, height=600)
lm <- getLinkCommunities(sample_network$edges,  hcmethod = 'single')
dev.off()
lm

#Plot
png(file=paste(name, '/LCplot.png', sep=''), width=600, height=600)
plot(lm, type = "graph")
dev.off()

# Membership cleaning
membership_linkcomm = clean_membership(memberize_linkcomm(lm$clusters, sample_network$number_edges))
mixed_membership_linkcomm = labelize(sample_network$node_edge, membership_linkcomm, mode = 'multi')
colored_linkcomm = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcomm)

# Plot of adjacency
png(file=paste(name, '/LCadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_linkcomm)
dev.off()


### Linkcomm mod ###
# New cutoff point
bestlm = bestLinkcomm(lm, 50, line_graph, sample_network)

png(file=paste(name, '/LCmodplot.png', sep=''), width=600, height=600)
plot(bestlm, type = "graph")
dev.off()

bestlm$pdmax
membership_linkcommmod = clean_membership(memberize_linkcomm(bestlm$clusters, sample_network$number_edges))
mixed_membership_linkcommmod = labelize(sample_network$node_edge, membership_linkcommmod, mode = 'multi')
colored_linkcommmod = Color_Adjacency(sample_network$adjacency, sample_network$edges, membership_linkcommmod)

png(file=paste(name, '/LCmodadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_linkcommmod)
dev.off()


### IMPERIAL ###
# lc, ld, le are igraphs objects with adjacency = C, D and E1
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

# Membership cleaning
membership_c = clean_membership(as.double(membership(cc)))
membership_d = clean_membership(as.double(membership(cd)))
membership_e = clean_membership(as.double(membership(ce)))


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

png(file=paste(name, '/Cadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_c)
dev.off()
png(file=paste(name, '/Dadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_d)
dev.off()
png(file=paste(name, '/Eadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_e)
dev.off()

# mixed membership
mixed_membership_c = labelize(sample_network$node_edge, membership_c, mode='multi')
mixed_membership_d = labelize(sample_network$node_edge, membership_d, mode='multi')
mixed_membership_e = labelize(sample_network$node_edge, membership_e, mode='multi')

png(file=paste(name, '/Dplot.png', sep=''), width=600, height=600)
plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_c,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(cc),
     main = 'mixed membership C method')
dev.off()

png(file=paste(name, '/Dplot.png', sep=''), width=600, height=600)
plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_d,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(cd),
     main = 'mixed membership D method')
dev.off()

png(file=paste(name, '/Eplot.png', sep=''), width=600, height=600)
plot(igraph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=mixed_membership_e,
     vertex.size=5,
     edge.width=1,
     edge.color=membership(ce),
     main = 'mixed membership E method')
dev.off()

### LINKSCAN* ###
# LinkSCAN algorithm
png(file=paste(name, '/pls.png', sep=''), width=600, height=600)
LS_clustering = LinkScan(sample_network,
                         mu = 0.5,
                         n_eps = 20,
                         n_critic = 5,
                         poly_approx = FALSE,
                         star= TRUE,
                         alpha = 0,
                         beta = 1,
                         measure='NMI')
dev.off()

LS_clustering$epsilon
membershipLS = clean_membership(LS_clustering$membership)
write.csv(membershipLS, paste(name, '/lsmembership.csv', sep=''), row.names = FALSE)

# Plot Linkscan adjacency
colored_linkscan = Color_Adjacency(sample_network$adjacency, sample_network$edges, membershipLS)
png(file=paste(name, '/LSadj.png', sep=''), width=600, height=600)
plot_adjacency(colored_linkscan)
dev.off()

linkscan_membership = labelize(sample_network$node_edge, membershipLS, mode='multi')
# Plot graph
png(file=paste(name, '/LSplot.png', sep=''), width=600, height=600)
plot(sample_network$graph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=linkscan_membership,
     vertex.size=4,
     edge.width=1,
     edge.color=membershipLS,
     main = 'mixed membership linkscan')
dev.off()



# Save the clusterings
memberships = data.frame(membership_linkcomm,
                         membership_linkcommmod,
                         membership_c,
                         membership_d,
                         membership_e,
                         membershipLS)

write.csv(memberships, paste("Memberships/",".csv", sep = name), row.names=FALSE)


### Comparisons ###
# Modularity
mod_linkcomm = modularity(line_graph, membership_linkcomm+1)
mod_linkcommmod = modularity(line_graph, membership_linkcommmod+1)
mod_c = modularity(line_graph, membership_c+1)
mod_d = modularity(line_graph, membership_d+1)
mod_e = modularity(line_graph, membership_e+1)
mod_linkscan = modularity(line_graph, membershipLS+1)
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
#Set the misclassification bound
N_misclass = 8
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

# Adjusted Rand Index
sbm_membership = colored_to_membership(sample_network)
ARI_LC = adj.rand.index(sbm_membership, membership_linkcomm)
ARI_LCmod = adj.rand.index(sbm_membership, membership_linkcommmod)
ARI_c = adj.rand.index(sbm_membership, membership_c)
ARI_d = adj.rand.index(sbm_membership, membership_d)
ARI_e = adj.rand.index(sbm_membership, membership_e)
ARI_LS = adj.rand.index(sbm_membership, membershipLS)
ARIs = c(ARI_LC, ARI_LCmod, ARI_c, ARI_d, ARI_e, ARI_LS)


# NMI
n_edges = sample_network$number_edges
df_mem_sbm = data.frame(1:n_edges, sbm_membership)
df = cbind(1:n_edges, memberships)
NMI_LC = as.double(NMI(df_mem_sbm, df[,c(1,2)]))
NMI_LCmod = as.double(NMI(df_mem_sbm, df[,c(1,3)]))
NMI_c = as.double(NMI(df_mem_sbm, df[,c(1,4)]))
NMI_d = as.double(NMI(df_mem_sbm, df[,c(1,5)]))
NMI_e = as.double(NMI(df_mem_sbm, df[,c(1,6)]))
NMI_LS = as.double(NMI(df_mem_sbm, df[,c(1,7)]))
NMIs = c(NMI_LC, NMI_LCmod, NMI_c, NMI_d, NMI_e, NMI_LS)

measurement = data.frame(Methods = c('Linkcomm','Linkcomm + mod', 'Imperial C','Imperial D', 'Imperial E', 'LinkScan*'),
                         Modularity = modularities,
                         Expectation_L2error = expectations,
                         Approx_misclassification = misclassifications,
                         ARI = ARIs,
                         NMI = NMIs)
measurement
write.csv(measurement, file = paste("Measurements/",".csv", sep = name), row.names = FALSE)

