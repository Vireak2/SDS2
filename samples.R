library(linkcomm)

# Test with undirected, assortative network and 5 communities
# Edges probabilities with P~U(0,0.2) for inter community and P~U(0.6,1) for intra community. And number of nodes ~ Pois(30)
blockprob = matrix(runif(25,0,0.2), 5, 5)
blockprob[lower.tri(blockprob)] <- t(blockprob)[lower.tri(blockprob)]
diag(blockprob) = runif(5,0.6,1)
nodesdistr = rpois(5,30)

# Build the network
sample_network = SBM(nodesdistr, blockprob)

# Link community detection
lm <- getLinkCommunities(sample_network$edges)
lm
plot(lm, type = "graph")


# Test with undirected, disassortative network and 5 communities
blockprob = matrix(runif(25,0.6,1), 5, 5)
blockprob[lower.tri(blockprob)] <- t(blockprob)[lower.tri(blockprob)]
diag(blockprob) = runif(5,0,0.2)
nodesdistr = rpois(5,30)

# Link community detection
lm <- getLinkCommunities(sample_network$edges)
lm
plot(lm, type = "graph")
