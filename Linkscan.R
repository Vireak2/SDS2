Jaccard <- function(x, y){
  # Compute the Jaccard similarity between two 0/1 vectors
  # Formula: sigma(x,y) = |x^y|/|xUy|
  
  similarity = t(x)%*%y / (sum((x+y)!=0))
  return(similarity)
  
}



Link_similarities <- function(graph, selfloops = FALSE, star= FALSE, alpha = 0, beta = 1){
  # Given a graph return the weighted line graph matrix 
  # (adjacency of the LinkSpace graph)
  # where the weights are given by the Jaccard similarity
  # of the two endpoints of the adjacent edges
  # Self loops if we consider the similarity are define when both edges are the same (do not use)
  # alpha and beta are used in the sampling step of linkscan* alpha = 0 sets alpha to 2<d>
  
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
  
  
  
  
  
  
  ### Sampling ###
  if(star){
    n_LS_vertices = dim(LS_weights)[1]
    if(alpha == 0){
      alpha = mean(degree(graph, mode = 'out', loops =TRUE))*2
    }
    
    LS_degrees = degree(LS_graph, mode = 'out', loops =TRUE)
    LS_star_weights = LS_weights
    
    for(v in 1:n_LS_vertices){
      nv = alpha + beta*log(LS_degrees[v])
      if(nv < LS_degrees[v]){
        neighbors = which(LS_weights[v,] != 0)
        sampled = sample(neighbors, size = floor(nv))
        LS_star_weights[v, -sampled] = 0
        LS_star_weights[-sampled, v] = 0
      }
    }
    LS_weights = LS_star_weights
    LS_graph = graph_from_adjacency_matrix(LS_weights, mode = 'undirected', weighted = TRUE)
  }
  
  
  
  return_list = list(LS_graph, LS_weights)
  names(return_list) = c('LS_graph', 'LS_weights')
  return(return_list)
}




StructuralClustering <- function(adjacency_weigthed, LS_graph, epsilon, mu=0.7){
  # Compute the selections function and the epsilon neighbors of the LinkSpace graph
  # Then apply the Algorithm 3 'StructuralClustering' with parameters epsilon and mu
  
  # Get the LinkSpace weights matrix and LinkScan Igraph
  n_vertex = dim(adjacency_weigthed)[1]
  
  #Declare the vector of degree (needed for the selection functions)
  degree = degree(LS_graph, mode = 'out', loops = TRUE)
  
  # Declare the vector of selection functions
  selection = rep(0, n_vertex)
  
  # Declare the matrix that store the epsilon-Neighbors of the LS-vertex in the rows
  N_eps = list()
  
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
    
    N_eps = append(N_eps, list(Nv_eps))  #add the row of epsilon neighbors corresponding to v
  }
  
  
  # Initial set of unclassified LS-vertices and neutral set
  unclassified=1:n_vertex
  neutral = c()
  
  # Set of core nodes
  core_nodes = unclassified[selection > mu]
  
  
  # Create the DIRREACH set
  dir.reach = list()
  for (v in unclassified) {
    if(v %in% core_nodes){
      dir.reach_v = as.integer(N_eps[[v]])
    }else{
      dir.reach_v = c()
    }
    dir.reach = append(dir.reach, list(dir.reach_v))
  }
  
  
  ### StructuralClustering ###
  #Initialisation
  
  # Declare the data frame that contains the clusters on the rows
  clusters=data.frame(row.names = 'cluster')
  
  
  # Cluster ID start at 2 (1 is for the neutral community)
  clusterID = 1
  
  # vector of membership of the LS-vertices (edges of original graph)
  membership = rep(0, n_vertex)
  
  while(length(unclassified) !=0 ){
    v = unclassified[1]
    if((selection[v] > mu)){ # if v is a core node
      clusterID=clusterID+1 # generate a new clusterID
      cluster_v = c(v) # generate a new cluster
      membership[v] = clusterID
      unclassified = unclassified[unclassified != v ] # remove v from unclassified
      Q = dir.reach[[v]] #add N_eps(v) to a queue Q
      
      while(length(Q) != 0 ) { #While Q!=0
        y = Q[1] # y = first node in Q
        R = dir.reach[[y]]#R = DirReach(v)
        
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
  
  
  
  StructClust=list(LS_graph, N_eps, selection, clusters, neutral, membership, epsilon)
  names(StructClust)=c('LS_graph','EpsilonNeighbors', 'SelectionFunction', 'Clusters', 'NeutralCluster', 'membership', 'epsilon')
  return(StructClust)
}



similarity_chart <- function(LS_weigths, fraction_nodes = 1, mu=0.7){
  # Return the sorted vector of 100mu-th percentile lowest similarity of each nodes
  # fraction node denotes the fraction of nodes to sample 
  
  adjacency = LS_weigths
  n_nodes = dim(adjacency)[1]
  
  
  # Sampling
  if (fraction_nodes != 1){
    nodes = sample(1:n_nodes, floor(fraction_nodes * n_nodes))
    adjacency = adjacency[nodes, nodes]
    n_nodes = dim(adjacency)[1]
  }
  
  # Vector of 100mu% lowest similarities
  mu_PLS = c()
  
  for (i in 1:n_nodes){
    S_i = adjacency[,i]
    N_i = S_i[which(S_i != 0)]
    mu_percentile = floor(mu * length(N_i))
    v_i =sort(N_i, decreasing = FALSE)[mu_percentile]  #Dont know if decreasing or not
    mu_PLS = append(mu_PLS, v_i)
  }
  
  mu_PLS = sort(mu_PLS)
  
  return(mu_PLS)
}




find_knee <- function(x, y, n_max, percent_inter = 0.05){
  # Returns the points of the function y=f(x) that maximize the amplitude of the 1st derivative 
  # n_max : number of maximums returned
  # percent_inter : radius of the interval where the others maximums can be selected
  # as a fraction of the total size of x
  
  deriv1 = diff(y)/diff(x) # first derivative
  n_indices = length(x)
  radius_inter = floor(percent_inter * n_indices)
  
  xmax=c()
  ymax=c()
  
  deriv_domain = deriv1
  for (k in 1:n_max) {
    # look for the max of |f'(x)| in the domain
    max_deriv1 = max( abs( deriv_domain ) )
    
    i_max = which.max(abs(deriv_domain))
    x_i = x[i_max]
    y_i = y[i_max]
    
    l_bound = max(i_max - radius_inter, 1)
    u_bound = min(i_max + radius_inter, length(deriv_domain))
    deriv_domain[(l_bound:u_bound)] = 0
    
    xmax = append(xmax, x_i)
    ymax = append(ymax, y_i)
    
  }
  return_list = list(xmax, ymax)
  names(return_list) = c('x', 'y')
  return(return_list)
}


LinkScan <- function(sample_network,
                     mu =0.7,
                     n_eps = 20,
                     n_critic = 5,
                     poly_approx = FALSE,
                     star= TRUE,
                     alpha = 0,
                     beta = 1,
                     measure='NMI'){
  # Apply the LinkScan algorithm on the SBM sample_network
  # set mu, alpha and beta as the parameters of LinkScan (if alpha = 0, alpha = 2<d>)
  # n_eps is the number of sampled epsilon to try
  # n_critic is the number of points that maximize the curvature of the curve of similarity to try
  # star = TRUE to use the LinkSCAN*
  # measure = "ARI", "NMI", "modularity" is the quantity we wish to maximize after each try of epsilon
  
  
  # Get the LS graph and adjacency
  print('Computation of the LinkSpace weights and sampling...')
  linksim = Link_similarities(sample_network$graph, selfloops = FALSE, star, alpha, beta)
  
  ### Try to find a epsilon ###
  # 100mu-percentile lowest similarity
  PLS = similarity_chart(linksim$LS_weights, mu = mu)
  index = 1:length(PLS)
  plot(index, PLS, xlab = 'LinkScan vertex index', ylab = 'mu-percentile lowest similarity')
  
  
  if(poly_approx){
    # Use cubic approximation to get a smoother function
    approxPLS = spline(index, PLS)
    # plot the curve
    lines(approxPLS_line, col="green")
    # plot the the new found epsilon
    knees = find_knee(index, approxPLS_line, n_max=n_critic)
    eps_poly_x = knees$x
    eps_poly_y = knees$y
    points(eps_poly_x, eps_poly_y, col="blue", pch=20, cex = 2)
  }
  
  
  # Use own function to find the knee points of PLS
  knees3 = find_knee(index, PLS, n_max=n_critic)
  eps_x = knees3$x
  eps_y = knees3$y
  points(eps_x, eps_y, col="orange", pch=20, cex = 2)
  M = length(eps_x)
  
  
  # Choose epsilons in the interval using gaussian mixture
  possible_eps = rnormm(n_eps, rep(1,M)/M, eps_y, sigma = rep(1,M)/20)
  possible_eps = c(abs(possible_eps),eps_y)
  possible_eps = sort(unique(possible_eps))
  
  # Set the progress bar
  print('Structural clustering...')
  pb = txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(possible_eps)+1, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")
  progress=0  
  setTxtProgressBar(pb, progress)

  # Objects to compare the measures
  if(measure=='modularity'){
    line_graph = line.graph(sample_network$graph)
  }
  if(measure=='ARI'){
    sbm_membership = colored_to_membership(sample_network)
  }
  if(measure=='NMI'){
    number_edges = sample_network$number_edges
    df_mem_sbm = data.frame(1:number_edges,colored_to_membership(sample_network))
  }
  Measures=c(0)
  
  
  
  # Try all the epsilons
  for (epsilon in possible_eps) {
    LS_clustering = StructuralClustering(linksim$LS_weights, linksim$LS_graph, epsilon, mu = mu)
    
    if(measure=='modularity'){
    meas = modularity(line_graph, LS_clustering$membership)
    }
    if(measure=='ARI'){
      meas = adj.rand.index(sbm_membership, LS_clustering$membership)
    }
    if(measure=='NMI'){
      df_ls = data.frame(1:number_edges,LS_clustering$membership)
      meas = as.double(NMI(df_mem_sbm, df_ls))
    }
    
    if(meas >= max(Measures)){
      best_clustering = LS_clustering
    }
    Measures = append(Measures, meas)
    
    progress = progress+1
    setTxtProgressBar(pb, progress)
  }
  
  # Get the best one
  epsilon = best_clustering$epsilon
  setTxtProgressBar(pb, progress+1)
  close(pb)
  return(best_clustering)
}

