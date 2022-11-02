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




find_knee <- function(x, y, n_intervals, threshold = 0){
  # Returns the points of the function y=f(x)
  # that maximize the amplitude of the 1st derivative 
  # n_intervals: cut x into n_intervals equal size pieces and return the maximum on each intervals
  # threshold: keep the points that correspond to a value of |f'(x)| >threshold
  # (Will change thow threshold works later, might implement the n-maximums)
  
  deriv1 = diff(y)/diff(x) # first derivative
  deriv2 = diff(deriv1)/diff(x[-1]) # second derivative (useless here)
  indices = length(x)
  chunk_size = floor(length(x)/n_intervals)
  
  xmax=c()
  ymax=c()
  
  for (i in 1:n_intervals) {
    # look for the max of |f'(x)| in the i-th interval
    max_deriv1 = max( abs( deriv1[((i-1)*chunk_size+1) : min(i*chunk_size, length(deriv2))] ) )
    
    x_i = which.max(abs(deriv1[((i-1)*chunk_size+1):min(i*chunk_size, length(deriv2))])) + chunk_size*(i-1)
    y_i = y[x_i]
    
    if(max_deriv1 >threshold ){
    xmax = append(xmax, x_i)
    ymax = append(ymax, y_i)
    }
    
  }
  return_list = list(xmax, ymax)
  names(return_list) = c('x', 'y')
  return(return_list)
}