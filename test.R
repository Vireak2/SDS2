overlap_coefficient<- function(x, y){
  # Compute the overlap quality between two 0/1 vectors
  # Formula: O(x,y) = |x^y|/min(|x|,|y|)
  
  similarity = t(x)%*%y / (min(sum(x), sum(y)))
  return(similarity)
  
}


overlap_quality<- function(adjacency){
  adjacency = ceiling(adjacency)
  n_nodes = dim(adjacency)[1]
  edges = matrix(0, 0, 2)
  edges = edges_from_adjacency(adjacency)
  
  # Matrix of the gamma vector that is (the set of vertices adjacent to v) + v
  gamma = adjacency + diag(1, dim(adjacency))
  n_edges = dim(edges)[1]
  
  ### Computation of the weights ###
  #Declare the weighted adjacency
  OC = matrix(0, n_edges, n_edges)
  
  
  # Compute the weights of
  # the adjacent edges of the type e1 = (i, k) and e2 = (k, j) or (j, k)
  for(e1 in 1:(n_edges-1)){
    w1 = edges[e1,1]
    k = edges[e1,2]
    
    for(e2 in (e1+1):n_edges){ # Look for edges stored after e1 in the edge set
      if(edges[e2,2] == k){ #if e2 = (k, j)
        w2 = edges[e2,1]
        similarity = overlap_coefficient(gamma[w1,], gamma[w2,])
        OC[e1, e2] = similarity
        OC[e2, e1] = similarity
      }
      
      if(edges[e2,1] == k){ #if e2 = (j, k)
        w2 = edges[e2,2]
        similarity = overlap_coefficient(gamma[w1,], gamma[w2,])
        OC[e1, e2] = similarity
        OC[e2, e1] = similarity
      }
    }
  }
  
  for(e1 in 1:(n_edges-1)){
    w1 = edges[e1,2]
    k = edges[e1,1]
    
    for(e2 in (e1+1):n_edges){
      if(edges[e2,2] == k){
        w2 = edges[e2,1]
        similarity = overlap_coefficient(gamma[w1,], gamma[w2,])
        OC[e1, e2] = similarity
        OC[e2, e1] = similarity
      }
      
      if(edges[e2,1] == k){
        w2 = edges[e2,2]
        similarity = overlap_coefficient(gamma[w1,], gamma[w2,])
        OC[e1, e2] = similarity
        OC[e2, e1] = similarity
      }
    }
  }
  return(OC)
}


expected_communities <- function(sample_network, as_vertex=TRUE){
  blockprob = sample_network$block_prob
  nodes_partition = sample_network$nodes_partition
  n_blocks = length(nodes_partition)
  expected = rep(0, n_blocks)
  number_nodes = sum(nodes_partition)
  blocks_index = 1:n_blocks
  pmatrix = matrix(0, n_blocks, n_blocks)
  # Indices that correspond to the beginning and ending of nodes communities in the adjacency matrix
  
  
  for (c1 in 1:n_blocks) {
    for (c2 in 1:n_blocks) {
      nc1=nodes_partition[c1]; nc2=nodes_partition[c2]
      p_c1c2 = (1-blockprob[c1,c2])^nc2 #Probability that a vertex of c1 has no edges to c2
      pmatrix[c1,c2] = p_c1c2
     }
   }
  
  expectations = c()
  for (i in blocks_index) {
    Ei=0
    for(k in 1:(n_blocks-1)){
    c_ik = 0 # Probability that vertex in community i has k communities
    kcomb = combn(blocks_index,k)
    possibilities = (dim(kcomb))[2]
      for(comb_index in 1:possibilities){
        comi = kcomb[,comb_index]
        c_ik=c_ik+prod(pmatrix[i,-comi])*prod(1-pmatrix[i,comi])
      }
    Ei = Ei + c_ik*k
    }
    comi = 1:n_blocks
    c_ik=c_ik+prod(pmatrix[i,-comi])*prod(1-pmatrix[i,comi])
    Ei = Ei + c_ik*k
    expectations = append(expectations, Ei)
  }
  
  if(as_vertex){
    E = c()
    for(c in blocks_index){
      Ec = rep(expectations[c], nodes_partition[c])
      E = append(E, Ec)
    }
    expectations = E
  }
  
  return(expectations)
}


number_communities <- function(mixed_membership){
  n_vertex = length(mixed_membership)
  n_communities = rep(0, n_vertex)
  
  for(i in 1:n_vertex){
    n_communities[i] = sum(mixed_membership[[i]] != 0)
  }
  return(n_communities)
}


misclassification_edges <- function(sample_network,
                                    colored_clustering,
                                    n_com = 0){
  
  n_edges = sample_network$number_edges
  colored_adj = sample_network$colored_adjacency
  n_vertex = dim(colored_adj)[1]
  
  
  communities_obs = 1:max(colored_adj)
  communities_predi = 1:max(colored_clustering)
  
  obs_count = c()
  predi_count = c()
  
  for(c in communities_obs){
    obs_count = append(obs_count, sum(colored_adj == c))
  }
  
  for(c in communities_predi){
    predi_count = append(predi_count, sum(colored_clustering == c))
  }
  
  obs_sort = sort(obs_count, index.return = TRUE, decreasing =TRUE)[[2]]
  predi_sort = sort(predi_count, index.return = TRUE, decreasing = TRUE)[[2]]
  
  
  
  
  if(length(communities_predi) > length(communities_obs)){
    matrix1 = colored_clustering
    com_sort = predi_sort
    com = communities_predi
    matrix2 = colored_adj
    com_sort2 = obs_sort
    n_max = length(communities_predi)
    if(n_com == 0){
    n_com = length(communities_obs)
    }
    
  }else{
    
    matrix1 = colored_adj
    com_sort = obs_sort
    com = communities_obs
    matrix2 = colored_clustering
    com_sort2 = predi_sort
    n_max = length(communities_obs)
    if(n_com == 0){
      n_com = length(communities_predi)
    }
  }
  
  colored_matrix = matrix(match(matrix1,
                                c(0,com_sort[1:n_com]),
                                nomatch = n_max+1),
                          n_vertex,n_vertex)-1
  
  colored_matrix2 = matrix(match(matrix2,
                                 c(0,com_sort2[1:n_com]),
                                nomatch = n_max+2),
                          n_vertex,n_vertex)-1
  
  colored_matrix2 = colored_matrix2[upper.tri(colored_matrix2, diag = TRUE)]
  

  permutations = permn(1:n_com)
  misclass_rate = c()
  for (iter in 1:length(permutations)) {
    permut_i = permutations[[iter]]
    temp_colored = colored_matrix
      
    for(i in 1:n_com){
      temp_colored[colored_matrix==i] = permut_i[i]
    }
    
    temp_colored = temp_colored[upper.tri(temp_colored, diag = TRUE)]
    misclass_i = sum((temp_colored - colored_matrix2)!=0) / n_edges
    misclass_rate = append(misclass_rate, misclass_i)
  }
  
  return_list = list(misclass_rate,min(misclass_rate))
  return(return_list)
}