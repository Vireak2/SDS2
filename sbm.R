SBM <- function(nodes_partition, block_prob, directed=FALSE, loops=FALSE){
  # Given a vector that set the number of nodes on each block and the matrix of edges probabilities
  # build a stochastic block model including the number of nodes, nodes partition,matrix of edges probabilities, adjacency matrix, degree vector and edge set
  # It also return the line graph adjacency matrix and the 'D' and 'E' versions of it
  # The function works for directed graphs, for undirected graph the lower triangular matrix of edges probabilities is considered
  
  # Define the number of community and total number of nodes 
  number_community = length(nodes_partition)
  number_nodes = sum(nodes_partition)
  
  # Declare the adjacency matrix
  adjacency = matrix(0, number_nodes, number_nodes)
  
  # Indices that correspond to the beginning and ending of nodes communities in the adjacency matrix
  communities_index_b = head( cumsum(c(1,nodes_partition)), -1)
  communities_index_e = cumsum(nodes_partition)
  
  # Computation of the adjacency matrix block by block when the graph is undirected
  if (directed == FALSE) {
    for (c1 in 1:number_community) {
      for (c2 in c1:number_community) {
        nc1=nodes_partition[c1]; nc2=nodes_partition[c2]
        adj_c1c2 = matrix(rbinom(nc1*nc2, 1, block_prob[[c1,c2]]), nc1, nc2)
        adjacency[communities_index_b[c1]:communities_index_e[c1], communities_index_b[c2]:communities_index_e[c2]] = adj_c1c2
      }
    }
    adjacency[lower.tri(adjacency)] <- t(adjacency)[lower.tri(adjacency)]
  }
  
  # Computation of the adjacency matrix block by block when the graph is directed
  if (directed == TRUE) {
    for (c1 in 1:number_community) {
      for (c2 in 1:number_community) {
        nc1=nodes_partition[c1]; nc2=nodes_partition[c2]
        adj_c1c2 = matrix(rbinom(nc1*nc2, 1, block_prob[[c1,c2]]), nc1, nc2)
        adjacency[communities_index_b[c1]:communities_index_e[c1], communities_index_b[c2]:communities_index_e[c2]] = adj_c1c2
      }
    }    
  }
  
  if (!loops){
    diag(adjacency)=0
  }
  
  
  graph = graph_from_adjacency_matrix(adjacency, 'undirected')
  
  # Compute degrees vector
  degrees = colSums(adjacency)
  
  # Declare the edge set
  edges = matrix(0, 0, 2)
  
  # Computation of the edge set undirected
  if (directed == FALSE) {
    for (i in 1:number_nodes) {
      for (j in i:number_nodes) {
        if (adjacency[[i,j]] != 0){
          edges = rbind(edges, c(i,j))
        }
      }
    }
  }
  
  # Computation of the edge set directed
  if (directed == TRUE) {
    for (i in 1:number_nodes) {
      for (j in 1:number_nodes) {
        if (adjacency[[i,j]] != 0){
          edges = rbind(edges, c(i,j))
        }
      }
    }
  }
  
  
  number_edges = dim(edges)[1]
  
  # Compute the Matrix C, the line graph matrix
  # Compute the weight penalties
  node_edge = matrix(0, number_nodes, number_edges)
  for (k in 1:number_edges) {
    e = edges[k,]
    node_edge[e[1], k] = 1
    node_edge[e[2], k] = 1
  }
  
  l_graph = matrix(0, number_edges, number_edges)
  for (alpha in 1:(number_edges-1)) {
    for (beta in (alpha+1):number_edges) {
      l_graph[alpha, beta] = sum(node_edge[,alpha]*node_edge[,beta])
    }
  }
  l_graph[lower.tri(l_graph)] <- t(l_graph)[lower.tri(l_graph)]
  
  # Compute the Matrix D, the weighted line graph matrix
  weight_penalty = degrees
  for (k in 1:number_nodes) {
    if (degrees[k] != 1){
      weight_penalty[k] = 1/(degrees[k]-1)
    }
  }
  # Computation of D
  D_l_graph = matrix(0, number_edges, number_edges)
  for (alpha in 1:(number_edges-1)) {
    for (beta in (alpha+1):number_edges) {
      D_l_graph[alpha, beta] = sum(weight_penalty*node_edge[,alpha]*node_edge[,beta])
    }
  }
  D_l_graph[lower.tri(D_l_graph)] <- t(D_l_graph)[lower.tri(D_l_graph)]
  
  
  # Compute the Matrix E, the corrected weighted node-edge matrix with self loops
  # Compute the E weight penalties
  E_weight_penalty = degrees
  for (k in 1:number_nodes) {
    if (degrees[k] != 0){
      E_weight_penalty[k] = 1/degrees[k]
    }
  }
  
  # Computation of E
  E_l_graph = matrix(0, number_edges, number_edges)
  for (alpha in 1:number_edges) {
    for (beta in alpha:number_edges) {
      E_l_graph[alpha, beta] = sum(E_weight_penalty*node_edge[,alpha]*node_edge[,beta])
    }
  }
  E_l_graph[lower.tri(E_l_graph)] <- t(E_l_graph)[lower.tri(E_l_graph)]
  
  E1_l_graph = E_l_graph %*% E_l_graph - E_l_graph
  
  # Return the properties of the network
  return_list = list(number_nodes, nodes_partition, block_prob, adjacency, graph, degrees, edges, number_edges, node_edge, l_graph, D_l_graph, E_l_graph, E1_l_graph)
  names(return_list) = c('number_nodes', 'nodes_partition', 'block_prob', 'adjacency','graph', 'degrees', 'edges','number_edges', 'node_edge', 'line_graph', 'D_line_graph', 'E_line_graph', 'E1_line_graph')
  return(return_list)
}



















RandomSBM <- function(number_nodes = 100, inter=c(0.01, 0.2), intra = c(0.4,0.6), mode = 'Assortative',
                      number_blocks=floor(sqrt(number_nodes)), prob_mix = 0.5, mode_node = 'random'){

  if(mode == 'Assortative'){
  blockprob = matrix(runif(number_blocks*number_blocks,inter[1], inter[2]), number_blocks, number_blocks)
  blockprob[lower.tri(blockprob)] <- t(blockprob)[lower.tri(blockprob)]
  diag(blockprob) = runif(number_blocks, intra[1], intra[2])
  }
  
  
  if(mode == 'Disassortative'){
    blockprob = matrix(runif(number_blocks*number_blocks,intra[1], intra[2]), number_blocks, number_blocks)
    blockprob[lower.tri(blockprob)] <- t(blockprob)[lower.tri(blockprob)]
    diag(blockprob) = runif(number_blocks, inter[1], inter[2])
  }
  
  
  if(mode == 'Mixed'){
    blockprob = matrix(0, number_blocks, number_blocks)
    temp = matrix(rbinom(number_blocks*number_blocks, 1, prob_mix), number_blocks, number_blocks)
    for (i in 1:number_blocks) {
      for (j in 1:i) {
        if(temp[i,j] == 1){
          blockprob[i,j]= runif(1, intra[1], intra[2])
        }
        else{
          blockprob[i,j]= runif(1, inter[1], inter[2])
        }
      }
    }
    blockprob[upper.tri(blockprob)] <- t(blockprob)[upper.tri(blockprob)]
  }
  
  
  if(mode_node == 'random'){
  nodesdistr = rpois(number_blocks, number_nodes/number_blocks)
  }
  else{
    nodesdistr = replicate(n_blocks, floor(number_nodes/number_blocks))
  }
  sample_network = SBM(nodesdistr, blockprob)
  
  return(sample_network)
  
  }

