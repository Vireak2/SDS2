overlap_coefficient<- function(x, y){
  # Compute the overlap quality between two 0/1 vectors
  # Formula: O(x,y) = |x^y|/min(|x|,|y|)
  
  similarity = t(x)%*%y / (min(sum(x), sum(y)))
  return(similarity)
  
}


overlap_quality<- function(colored_adjacency){
  n_block = max(colored_adjacency)
  OQ = matrix(1,n_block, n_block)
  
  for (i in 1:n_block) {
    comi = rowSums(colored_adjacency == i)!=0
    for (j in i+1:n_block) {
      comj = rowSums(colored_adjacency == j)!=0
      OQ[i,j] = overlap_coefficient(comi, comj)
    }
  }
  OQ[lower.tri(OQ)] = t(OQ)[lower.tri(OQ)]
  return(OQ)
}


clean_membership <- function(membership, threshold = 8){
  mem_table = table(membership)
  new_membership_table = mem_table[mem_table >= threshold]
  communities = sort(as.double(names(new_membership_table)))
  communities=communities[communities!=0]
  new_n_com = length(communities)
  temp = sapply(membership, function(x){if(x %in% communities){return(x)}else{return(0)}})
  new_membership=temp
  
  for(i in 1:new_n_com){
    new_membership[temp == communities[i]] = i
  }
  
  return(new_membership)
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
                                    n_com = 8){
  
  n_edges = sample_network$number_edges
  colored_adj = sample_network$colored_adjacency
  n_vertex = dim(colored_adj)[1]
  
  sorted_predi_coms = sort(table(colored_clustering[colored_clustering!=0]), decreasing = TRUE)
  sorted_actual_coms = sort(table(colored_adj[colored_adj!=0]), decreasing = TRUE)
  
  n_predi_coms = min(length(sorted_predi_coms),n_com)
  n_actual_coms = min(length(sorted_actual_coms),n_com)
  
  most_predi_coms = as.double(names(sorted_predi_coms[1:n_predi_coms]))
  most_actual_coms = as.double(names(sorted_actual_coms[1:n_actual_coms]))
  
  matrix1 = matrix(0,n_vertex,n_vertex)
  matrix2 = matrix(0,n_vertex,n_vertex)
  
  for (i in 1:n_predi_coms) {
    matrix1[colored_clustering == most_predi_coms[i]] = i
  }
  
  for (i in 1:n_actual_coms) {
    matrix2[colored_adj == most_actual_coms[i]] = i
  }
  
  n_edges_considered = max(sum(matrix1 != 0), sum(matrix2 != 0))/2
  
  permutations = permn(1:n_predi_coms)
  misclass_rate = c()
  for (iter in 1:length(permutations)) {
    permut_iter = permutations[[iter]]
    temp_colored = matrix1
    
    for(i in 1:n_predi_coms){
      temp_colored[which(matrix1==i)] = permut_iter[i]
    }
    
    differences = (temp_colored - matrix2)!=0
    differences = differences[upper.tri(differences, diag = TRUE)]
    misclass_i = sum(differences) / n_edges_considered
    misclass_rate = append(misclass_rate, misclass_i)
  }
  
  return_list = list(misclass_rate,min(misclass_rate))
  return(return_list)
  
}


misclassification_vertex <- function(sample_network,
                                    colored_clustering){
  
  n_vertex = sample_network$number_nodes
  nodes_partition = sample_network$nodes_partition
  colored_adj = sample_network$colored_adjacency
  communities_index_b = head( cumsum(c(1,nodes_partition)), -1)
  communities_index_e = cumsum(nodes_partition)
  actual_v_class = rep(0,n_vertex)
  predi_v_class_raw = rep(0,n_vertex)
  
  for(index in 1:length(communities_index_b)){
    interval = communities_index_b[index]:communities_index_e[index]
    actual_v_class[interval] = index
    
    predicted_values = c(colored_clustering[interval,interval][colored_clustering[interval,interval]!=0])
    most_common_predi = as.integer(names(sort(table(predicted_values), decreasing =TRUE)[1]))
    predi_v_class_raw[interval] = most_common_predi
  }
  
  n_unique_values = length(unique(predi_v_class_raw))
  predi_v_class_raw = predi_v_class_raw + n_unique_values
  predi_v_class = rep(0,n_vertex)
  for(j in 1:n_unique_values){
    predi_v_class[predi_v_class_raw==unique(predi_v_class_raw)[j]] = j
  }
  
  n_com = max(unique(predi_v_class))
  permutations = permn(1:n_com)
  misclassifications = rep(1,length(permutations))
  
  for(iter in 1:length(permutations)){
    permut_i = permutations[[iter]]
    temp_predict = predi_v_class
    
    for(i in 1:n_com){
      temp_predict[predi_v_class==i] = permut_i[i]
    }
    
    misclass_i = sum((temp_predict - actual_v_class)!=0)/n_vertex
    misclassifications[iter] = misclass_i
  }
  
  return(min(misclassifications)) 
}



plot_adjacency <- function(colored_adj){
  n_val = max(colored_adj)
  to_plot = t(apply(colored_adj, 2, rev))
  colors = c( 'white', colorRampPalette((brewer.pal(min(c(9,n_val)), name = 'Set1')))(n_val) )
  
  image(t(to_plot),
        col = colors,
        useRaster = TRUE,
        axes = FALSE)
  box()
}

bestLinkcomm = function(lm, n=50, linegraph, sample_network){
  cutoffpoints = seq(from = .9*lm$pdmax, to= 1.1*lm$pdmax, length.out = 50)
  modularities = c(0)
  for(point in cutoffpoints){
    newlm = newLinkCommsAt(lm, cutat=point)
    membership_linkcomm = memberize_linkcomm(newlm$clusters, sample_network$number_edges)
    mod = modularity(linegraph, membership_linkcomm+1)
    if(mod > max(modularities)){
      bestlm = newlm
    }
    modularities = append(modularities,mod)
  }
  return(bestlm)
}
