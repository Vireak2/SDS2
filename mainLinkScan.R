library(linkcomm)
library(modMax)
library(network)
library(plot.matrix)
library(KneeArrower)
library(akmedoids)


### test Linkscan ###

sample_network = RandomSBM(number_nodes=100,
                           mode = 'Assortative',
                           number_blocks = 3)

plot(sample_network$adjacency)

mu = 0.7

# Get the LS graph and adjacency
LinkSpace = Link_similarities(sample_network$graph)
LS_graph = LinkSpace$LS_graph
LS_adjacency = LinkSpace$LS_weights


### Try to find a epsilon ###
# 100mu-percentile lowest similarity
PLS = similarity_chart(LS_adjacency, mu = mu)
index = 1:length(PLS)
plot(index, PLS)

# find cutoff on the plot with the discrete x and y

eps_1 = findCutoff(PLS, index, method = 'first', .9) #Try eps_1 =elbow_point(PLS, index)
eps_1
#Plot the point
points(eps_1$x, eps_1$y, col="red", pch=20, cex = 2)


# Use cubic approximation to get a smoother function
approxPLS = lm(PLS~index + I(index^2) + I(index^3))
approxPLS_line = predict(approxPLS)[index]
# plot the curve
lines(approxPLS_line, col="green")
eps_cubic = findCutoff(index, approxPLS_line, method = 'first', .1)
eps_cubic
# plot the the new found epsilon
points(eps_cubic$x, eps_cubic$y, col="blue", pch=20, cex = 2)




### Structural Clustering ###
epsilon = eps_cubic$y # set your epsilon as you want
t = StructuralClustering(sample_network$graph, epsilon, mu = mu) #Main algorithm, try with different epsilons

length(t$NeutralCluster)  # See how big is the neutral cluster if too big => epsilon too big

# Plot of the edges communities
plot(sample_network$graph,
     vertex.label=NA,
     edge.color=t$membership,
     vertex.size=3,
     edge.width=2,
     main = 'linkscan',
     vertex.color = 'gray')


# further plots
linkscan_membership = labelize(sample_network$node_edge, t$membership, mode='multi')


plot(sample_network$graph,
     vertex.label=NA,
     vertex.shape="pie",
     vertex.pie=linkscan_membership,
     vertex.size=10,
     edge.width=1,
     main = 'mixed membership linkscan')












# Use own function to find the knee points of PLS
plot(PLS)
lines(approxPLS_line, col="green")
knees = find_knee(index, approxPLS_line, 8)
knees3 = find_knee(index, PLS, 1)

eps2_x = knees$x
eps2_y = knees$y

eps3_x = knees3$x
eps3_y = knees3$y

points(eps2_x, eps2_y, col="gold", pch=20, cex = 1)
points(eps3_x, eps3_y, col="orange", pch=20, cex = 2)




#Try to normalize the graph
normalizedPLS=(PLS-mean(PLS))/sd(PLS)
normalizedIndex = (1:length(PLS) - mean(1:length(PLS))) / sd(1:length(PLS) )
plot( normalizedIndex, normalizedPLS)
