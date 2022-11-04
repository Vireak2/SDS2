Jaccard <- function(x, y){
  # Compute the Jaccard similarity between two 0/1 vectors
  # Formula: sigma(x,y) = |x^y|/|xUy|
  
  similarity = t(x)%*%y / (sum((x+y)!=0))
  return(similarity)
  
}



Link_similarities <- function(graph, selfloops = TRUE){
  # Given a graph return the weighted line graph matrix 
  # (adjacency of the LinkSpace graph)
  # where the weights are given by the Jaccard similarity
  # of the two endpoints of the adjacent edges
  # Self loops if we consider the similarity are define when both edges are the same
  
  adjacency = as.matrix(as_adjacency_matrix(graph))
  n_nodes = dim(adjacency)[1]
  
  # Declare the edge set 
  edges = matrix(0, 0, 2)
  
  # Construct the edge set
  for (i in 1:n_nodes) {
    for (j in i:n_nodes) {
      if (adjacency[i,j] != 0){
        edges = rbind(edges, c(i,j))
      }
    }
  }
  
  
  # Matrix of the gamma vector that is (the set of vertices adjacent to v) + v
  gamma = adjacency + diag(1, dim(adjacency))
  n_edges = dim(edges)[1]
  
  
  
  ### Computation of the weights ###
  #Declare the weighted adjacency
  LS_weights = matrix(0, n_edges, n_edges)
  
  
  # Compute the weights of
  # the adjacent edges of the type e1 = (i, k) and e2 = (k, j) or (j, k)
  for(e1 in 1:(n_edges-1)){
    w1 = edges[e1,1]
    k = edges[e1,2]
    
    for(e2 in (e1+1):n_edges){ # Look for edges stored after e1 in the edge set
      if(edges[e2,2] == k){ #if e2 = (k, j)
        w2 = edges[e2,1]
        similarity = Jaccard(gamma[w1,], gamma[w2,])
        LS_weights[e1, e2] = similarity
        LS_weights[e2, e1] = similarity
      }
      
      if(edges[e2,1] == k){ #if e2 = (j, k)
        w2 = edges[e2,2]
        similarity = Jaccard(gamma[w1,], gamma[w2,])
        LS_weights[e1, e2] = similarity
        LS_weights[e2, e1] = similarity
      }
    }
  }
  
  
  
  # Compute the weights of 
  # adjacent edges of the type (k, i) + (k, j) or (j, k)
  for(e1 in 1:(n_edges-1)){
    w1 = edges[e1,2]
    k = edges[e1,1]
    
    for(e2 in (e1+1):n_edges){
      if(edges[e2,2] == k){
        w2 = edges[e2,1]
        similarity = Jaccard(gamma[w1,], gamma[w2,])
        LS_weights[e1, e2] = similarity
        LS_weights[e2, e1] = similarity
      }
      
      if(edges[e2,1] == k){
        w2 = edges[e2,2]
        similarity = Jaccard(gamma[w1,], gamma[w2,])
        LS_weights[e1, e2] = similarity
        LS_weights[e2, e1] = similarity
      }
    }
  }
  
  
  if(selfloops){LS_weights= LS_weights + diag(1, dim(LS_weights)[1], dim(LS_weights)[2])}
  
  # Return the LinkScan weights matrix
  # and igraph object with the corresponding adjacency
  LS_graph = graph_from_adjacency_matrix(LS_weights, mode = 'undirected', weighted = TRUE)
  return_list = list(LS_graph, LS_weights)
  names(return_list) = c('LS_graph', 'LS_weights')
  return(return_list)
}





StructuralClustering <- function(graph, epsilon, mu=0.7){
  # Compute the selections function and the epsilon neighbors of the LinkSpace graph
  # Then apply the Algorithm 3 'StructuralClustering' with parameters epsilon and mu
  
  # Get the LinkSpace weights matrix and LinkScan Igraph
  LinkSpace = Link_similarities(graph)
  adjacency_weigthed = LinkSpace$LS_weights
  n_vertex = dim(adjacency_weigthed)[1]
  
  LS_graph = LinkSpace$LS_graph
  
  #Declare the vector of degree (needed for the selection functions)
  degree = degree(LS_graph)
  
  # Declare the vector of selection functions
  selection = rep(0, n_vertex)
  
  # Declare the matrix that store the epsilon-Neighbors of the LS-vertex in the rows
  N_eps = data.frame(row.names = 'vertex')
  
  # Calculate the selection function and epsilon neighbors for each LS-vertices
  for (v in 1:n_vertex) {
    Nv = adjacency_weigthed[v,]
    Nv_eps = which(Nv>epsilon) # epsilon neighbors of v
    
    if (length(Nv_eps) == 0){ # if v is isolated
      selection[v] = 0
      Nv_eps = 0
    }else{
      fv = sum((Nv>epsilon))/degree[v] # formula for the selection functions
      selection[v] = fv
    }
    
    N_eps = rbind(N_eps, Nv_eps)  #add the row of epsilon neighbors corresponding to v
  }
  
  
  
  ### StructuralClustering ###
  #Initialisation
  
  # Declare the data frame that contains the clusters on the rows
  clusters=data.frame(row.names = 'cluster')
  
  # Initial set of unclassified LS-vertices and neutral set
  unclassified=1:n_vertex
  neutral = c()
  
  # Cluster ID start at 1+1=2 (1 is for the neutral community)
  clusterID = 1
  
  # vector of membership of the LS-vertices (edges of original graph)
  membership = rep(0, n_vertex)
  
  while(length(unclassified) !=0 ){
    v = unclassified[1]
    if((selection[v] > mu)){ # if v is a core node
      print(v)
      clusterID=clusterID+1 # generate a new clusterID
      cluster_v = c(v) # generate a new cluster
      unclassified = unclassified[unclassified != v ] # remove v from unclassified
      Q = unique(as.integer(N_eps[v,])) #add N_eps(v) to a queue Q; uniqe to avoid redondancy in Q
      
        while(length(Q) != 0 ) { #While Q!=0
          y = Q[1] # y = first node in Q
          R = unique(as.integer(N_eps[y,])) #R = DirReach(v)
                
            for (x in R) {
              
              if ( ((x %in% unclassified) | (x %in% neutral)) ){ #If x is unclassified or neutral
                cluster_v = append(cluster_v,x) #add x to the current cluster
                membership[x] = clusterID
                neutral = neutral[neutral != x] #remove x from neutral cluster
              }
              
              if(x %in% unclassified){  #If x is unclassified
                Q = append(Q,x) # add x to the end of Q
              }
              
              unclassified = unclassified[unclassified != x ] # remove x from unclassied
            }
            Q = Q[-1] # remove y from Q
        }
      cluser_v = unique(cluster_v)  # avoid redondancy in the cluster of v (although should not happen)
      clusters = rbind(clusters, cluster_v) # add the cluster of v in the dataframe
      
    }else{ # if v is not a core node
      neutral = append(neutral, v) # add v in neutral
      unclassified = unclassified[unclassified != v ] # remove v from unclassified
    }
  }
  
  neutral = append(neutral, unclassified) # add the forgotten node in neutral; never happens normally
  membership[neutral] = 1 # neutral community takes clusterID = 1
  clusters = rbind(clusters, neutral) # add neutral cluster to the dataframe
  
  
  
  StructClust=list(LS_graph,adjacency_weigthed, N_eps, selection, clusters, neutral, membership)
  names(StructClust)=c('LS_graph', 'Weights','EpsilonNeighbors', 'SelectionFunction', 'Clusters', 'NeutralCluster', 'membership')
  return(StructClust)
}

