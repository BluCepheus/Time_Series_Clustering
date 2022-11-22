library(tidyverse)
library(dtwclust)

# Hierarchical clustering function.
# Input is normalized dataset (matrix).
# Output is table of hierarchical clusterings (tibble).
# Function creates 48 clusterings using various distance measures,
# centroid calculation methods and linkage methods.
h_clust_fun <- function(dataset) {
  
  # Initializing tibble of clusterings
  clust_tib <- tibble(
    distance = rep(
      c('dtw_basic', 'lbk', 'lbi', 'sbd', rep('gak', 2)), each = 8
    ),
    centroid = rep(
      c(
        rep('dba', 3), 'shape_extraction', 'dba', 'shape_extraction'
      ), each = 8
    ),
    linkage = rep(
      c(
        'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty',
        'centroid', 'median'
      ), 6
    )
  )
  
  clustering_list = vector('list', length(48))
  
  for (i in 1:48) {
    if (clust_tib$centroid[i] == 'dba') {
      tsclust(dataset,
              type = 'hierarchical',
              k = 2L,
              distance = clust_tib$distance[i],
              centroid = dba,
              control = hierarchical_control(method = clust_tib$linkage[i]),
              window.size = 5L) ->
        clustering_list[[i]]
    } else {
      tsclust(dataset,
              type = 'hierarchical',
              k = 2L,
              distance = clust_tib$distance[i],
              centroid = shape_extraction,
              control = hierarchical_control(method = clust_tib$linkage[i]),
              window.size = 5L) ->
        clustering_list[[i]]
    }
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 48,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  clust_tib$clustering <- clustering_list
  
  clust_tib %>%
    add_column(
      code = paste(
        clust_tib$distance,
        clust_tib$centroid,
        clust_tib$linkage,
        sep = '_'
      ),
      .before = 'distance'
    ) %>%
    return()
}

# Function to compute internal CVIs for hierarchical clusterings
# Input is table of hierarchical clusterings (tibble) - output of h_clust_fun.
# Output is table with internal CVIs for hierarchical clusterings (tibble).
int_cvi_h_clust_fun <- function(h_clust_tib) {
  
  # Initializing tibble of CVIs
  cvi_tib <- tibble(
    code = h_clust_tib$code,
    distance = h_clust_tib$distance,
    centroid = h_clust_tib$centroid,
    linkage = h_clust_tib$linkage,
    Silhouette = rep(NA, 48),
    Score_Function = rep(NA, 48),
    Calinski_Harabasz = rep(NA, 48),
    Davies_Bouldin = rep(NA, 48),
    Mod_Davies_Bouldin = rep(NA, 48),
    Dunn = rep(NA, 48),
    COP = rep(NA, 48)
  )
  
  for (i in 1:48) {
    cvi(
      h_clust_tib$clustering[[i]],
      type = c('Sil', 'SF', 'CH', 'DB', 'DBstar', 'D', 'COP')
    ) %>% t() ->
      cvi_tib[i, 4:10]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 48,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }

  cvi_tib %>% return()
}

# Function to transform internal CVIs table to rank table
# Input is table with internal CVIs for hierarchical clusterings (tibble) - output of int_cvi_h_clust_fun.
# Output is table with ranking of internal CVIs for hierarchical clusterings (tibble).
int_rank_h_clust_fun <- function(int_cvi_h_clust_tib) {
  # Initializing tibble of ranks
  rank_tib <- tibble(
    code = int_cvi_h_clust_tib$code,
    distance = int_cvi_h_clust_tib$distance,
    centroid = int_cvi_h_clust_tib$centroid,
    linkage = int_cvi_h_clust_tib$linkage,
    Silhouette = rank(int_cvi_h_clust_tib$Silhouette),
    Score_Function = rank(int_cvi_h_clust_tib$Score_Function),
    Calinski_Harabasz = rank(int_cvi_h_clust_tib$Calinski_Harabasz),
    Davies_Bouldin = rank(-int_cvi_h_clust_tib$Davies_Bouldin),
    Mod_Davies_Bouldin = rank(-int_cvi_h_clust_tib$Mod_Davies_Bouldin),
    Dunn = rank(int_cvi_h_clust_tib$Dunn),
    COP = rank(-int_cvi_h_clust_tib$COP)
  )
  
  rank_tib %>%
    add_column(Total = rowSums(.[-c(1:4)])) %>%
    arrange(desc(Total)) %>%
    return()
}

# Function to calculate inertion of clusterings
# Input is normalized dataset (matrix) and table of hierarchical clusterings (tibble) - output of h_clust_fun.
# Output is table with inertion of each clustering for number of groups between 2 and 20.
inert_h_clust_fun <- function(dataset, h_clust_tib) {
  # Subfunction to calculate distance matrices
  # Input is normalized dataset (matrix) and method from c('dtw_basic', 'lbk', 'lbi', 'sbd', 'gak')
  dist_mat_fun <- function(dataset, method) {
    if (!(method %in% c('lbk', 'lbi'))) {
      dataset %>%
        scale() %>%
        proxy::dist(method = method) %>%
        return()
    } else {
      dataset %>%
        scale() %>%
        proxy::dist(method = method, window.size = 5L) %>%
        return()
    }
  }
  
  dist_mat_fun(dataset, 'dtw_basic') -> dtw_basic_dist
  dist_mat_fun(dataset, 'lbk') -> lbk_dist
  dist_mat_fun(dataset, 'lbi') -> lbi_dist
  dist_mat_fun(dataset, 'sbd') -> sbd_dist
  dist_mat_fun(dataset, 'gak') -> gak_dist
  
  
  clust_tib <- tibble(
    distance = rep(NA, 912),
    centroid = rep(NA, 912),
    linkage = rep(NA, 912),
    clusters = rep(NA, 912),
    inertion = rep(NA, 912)
  )
  
  n = 1
  for (i in 1:48) {
    if (h_clust_tib[i, 1] == 'dtw_basic') {
      current_dist <- dtw_basic_dist
    } else if (h_clust_tib[i, 1] == 'lbk') {
      current_dist <- lbk_dist
    } else if (h_clust_tib[i, 1] == 'lbi') {
      current_dist <- lbi_dist
    } else if (h_clust_tib[i, 1] == 'sbd') {
      current_dist <- sbd_dist
    } else if (h_clust_tib[i, 1] == 'gak') {
      current_dist <- gak_dist
    }
    
    for (k in 2:20) {
      clust_tib[n, 1] <- h_clust_tib[i, 1]
      clust_tib[n, 2] <- h_clust_tib[i, 2]
      clust_tib[n, 3] <- h_clust_tib[i, 3]
      clust_tib[n, 4] <- k
      
      1 - (
        ClustGeo::withindiss(
          stats::as.dist(current_dist),
          part = cutree(hierarchical_clustering$clustering[[i]], k = k)
        ) /
          ClustGeo::inertdiss(
            stats::as.dist(current_dist)
          )
      ) -> clust_tib[n, 5]
      
      n = n + 1
    }
  }
  
  clust_tib %>%
    add_column(
      code = paste(
        clust_tib$distance,
        clust_tib$centroid,
        clust_tib$linkage,
        sep = '_'
      ),
      .before = 'distance'
    ) %>%
    return()
}

# Function to compute external CVIs for hierarchical clusterings
# Input is row from h_clust_fun output tibble with chosen clustering and part of h_clust_fun output tibble with other clusterings to compare.
# Output is table with external CVIs for hierarchical clusterings (tibble).
ext_cvi_h_clust_fun <- function(chosen_clustering, compare_clusterings) {
  cvi_tib <- tibble(
    code = rep(compare_clusterings$code, each = 19),
    distance = rep(compare_clusterings$distance, each = 19),
    centroid = rep(compare_clusterings$centroid, each = 19),
    linkage = rep(compare_clusterings$linkage, each = 19),
    clusters = rep(2:20, nrow(compare_clusterings)),
    Rand = rep(NA, 19 * nrow(compare_clusterings)),
    Adjusted_Rand = rep(NA, 19 * nrow(compare_clusterings)),
    Jaccard = rep(NA, 19 * nrow(compare_clusterings)),
    Fowlkes_Mallows = rep(NA, 19 * nrow(compare_clusterings)),
    Variation_of_Information = rep(NA, 19 * nrow(compare_clusterings))
  )
  
  for (i in 1:(19 * nrow(compare_clusterings))) {
    compare_clusterings %>%
      filter(code == cvi_tib$code[i]) ->
      current_clustering
    
    cvi(
      cutree(chosen_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      cutree(current_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      type = 'RI'
    ) -> cvi_tib$Rand[i]
    cvi(
      cutree(chosen_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      cutree(current_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      type = 'ARI'
    ) -> cvi_tib$Adjusted_Rand[i]
    cvi(
      cutree(chosen_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      cutree(current_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      type = 'J'
    ) -> cvi_tib$Jaccard[i]
    cvi(
      cutree(chosen_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      cutree(current_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      type = 'FM'
    ) -> cvi_tib$Fowlkes_Mallows[i]
    cvi(
      cutree(chosen_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      cutree(current_clustering$clustering[[1]], k = cvi_tib$clusters[i]),
      type = 'VI'
    ) -> cvi_tib$Variation_of_Information[i]
  }
  
  cvi_tib %>% return()
}

#############################################################################

# Partitional clustering function.
# Input is normalized dataset (matrix).
# Output is table of partitional clusterings (tibble).
# Function creates 72 clusterings using various distance measures,
# centroid calculation methods and number of clusters from 2 to 10.
p_clust_fun <- function(dataset) {
  
  # Initializing tibble of clusterings
  clust_tib <- tibble(
    distance = rep(
      c('dtw_basic', 'lbk', 'lbi', 'sbd', rep('gak', 4)), each = 9
    ),
    centroid = rep(
      c(
        rep('dba', 3), 'shape', 'dba', 'pam', 'mean', 'median'
      ), each = 9
    ),
    clusters = rep(2:10, 8)
  )
  
  clustering_list = vector('list', length(72))
  
  for (i in 1:72) {
    
    tsclust(
      dataset,
      type = 'partitional',
      k = clust_tib[i, 3],
      distance = clust_tib$distance[i],
      centroid = clust_tib$centroid[i],
      control = partitional_control(
        pam.precompute = T,
        iter.max = 1000L,
        pam.sparse = T
      ),
      window.size = 5L
    ) -> clustering_list[[i]]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 72,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  clust_tib$clustering <- clustering_list
  
  clust_tib %>%
    add_column(
      code = paste(
        clust_tib$distance,
        clust_tib$centroid,
        clust_tib$clusters,
        sep = '_'
      ),
      .before = 'distance'
    ) %>%
    return()
}

# Function to compute internal CVIs for partitional clusterings
# Input is table of partitional clusterings (tibble) - output of p_clust_fun
# Output is table with internal CVIs for partitional clusterings (tibble).
int_cvi_p_clust_fun <- function(p_clust_tib) {
  
  # Initializing tibble of CVIs
  cvi_tib <- tibble(
    code = p_clust_tib$code,
    distance = p_clust_tib$distance,
    centroid = p_clust_tib$centroid,
    clusters = p_clust_tib$clusters,
    Silhouette = rep(NA, 72),
    Score_Function = rep(NA, 72),
    Calinski_Harabasz = rep(NA, 72),
    Davies_Bouldin = rep(NA, 72),
    Mod_Davies_Bouldin = rep(NA, 72),
    Dunn = rep(NA, 72),
    COP = rep(NA, 72)
  )
  
  for (i in 1:72) {
    cvi(
      p_clust_tib$clustering[[i]],
      type = c('Sil', 'SF', 'CH', 'DB', 'DBstar', 'D', 'COP')
    ) %>% t() ->
      cvi_tib[i, 5:11]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 72,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  cvi_tib %>% return()
}

# Function to transform internal CVIs table to rank table
# Input is table with internal CVIs for partitional clusterings (tibble) - output of int_cvi_p_clust_fun.
# Output is table with ranking of internal CVIs for partitional clusterings (tibble).
int_rank_p_clust_fun <- function(int_cvi_p_clust_tib) {
  # Initializing tibble of ranks
  rank_tib <- tibble(
    code = int_cvi_p_clust_tib$code,
    distance = int_cvi_p_clust_tib$distance,
    centroid = int_cvi_p_clust_tib$centroid,
    clusters = int_cvi_p_clust_tib$clusters,
    Silhouette = rank(int_cvi_p_clust_tib$Silhouette),
    Score_Function = rank(int_cvi_p_clust_tib$Score_Function),
    Calinski_Harabasz = rank(int_cvi_p_clust_tib$Calinski_Harabasz),
    Davies_Bouldin = rank(-int_cvi_p_clust_tib$Davies_Bouldin),
    Mod_Davies_Bouldin = rank(-int_cvi_p_clust_tib$Mod_Davies_Bouldin),
    Dunn = rank(int_cvi_p_clust_tib$Dunn),
    COP = rank(-int_cvi_p_clust_tib$COP)
  )
  
  rank_tib %>%
    add_column(Total = rowSums(.[-c(1:4)])) %>%
    arrange(desc(Total)) %>%
    return()
}

#############################################################################

# TADPole clustering function.
# Input is normalized dataset (matrix).
# Output is table of partitional clusterings (tibble).
# Function creates 54 clusterings using dtw_basic distance, dba centroids,
# lbk / lbi lower bounds, number of clusters from 2 to 10 and cutoff distance in c(20, 50, 100).
t_clust_fun <- function(dataset) {
  
  # Initializing tibble of clusterings
  clust_tib <- tibble(
    distance = rep('dtw_basic', 54),
    centroid = rep('dba', 54),
    lower_bound = rep(c('lbk', 'lbi'), each = 27),
    clusters = rep(2:10, each = 3, times = 2),
    cutoff = rep(c(100, 250, 500), 18)
  )
  
  clustering_list = vector('list', length(18))
  
  for (i in 1:54) {
    tsclust(
      dataset,
      type = 'tadpole',
      k = clust_tib$clusters[i],
      centroid = dba,
      control = tadpole_control(
        dc = clust_tib$cutoff[i],
        window.size = 5L,
        lb = clust_tib$lower_bound[i]
      ),
      window.size = 5L
    ) -> clustering_list[[i]]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 54,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  clust_tib$clustering <- clustering_list
  
  clust_tib %>%
    add_column(
      code = paste(
        clust_tib$distance,
        clust_tib$centroid,
        clust_tib$lower_bound,
        clust_tib$clusters,
        clust_tib$cutoff,
        sep = '_'
      ),
      .before = 'distance'
    ) %>%
    return()
}

# Function to compute internal CVIs for TADPole clusterings
# Input is table of TADPole clusterings (tibble) - output of t_clust_fun
# Output is table with internal CVIs for TADPole clusterings (tibble).
int_cvi_t_clust_fun <- function(t_clust_tib) {
  
  # Initializing tibble of CVIs
  cvi_tib <- tibble(
    code = t_clust_tib$code,
    distance = t_clust_tib$distance,
    centroid = t_clust_tib$centroid,
    lower_bound = t_clust_tib$lower_bound,
    clusters = t_clust_tib$clusters,
    cutoff = t_clust_tib$cutoff,
    Silhouette = rep(NA, 54),
    Score_Function = rep(NA, 54),
    Calinski_Harabasz = rep(NA, 54),
    Davies_Bouldin = rep(NA, 54),
    Mod_Davies_Bouldin = rep(NA, 54),
    Dunn = rep(NA, 54),
    COP = rep(NA, 54)
  )
  
  for (i in 1:54) {
    cvi(
      t_clust_tib$clustering[[i]],
      type = c('Sil', 'SF', 'CH', 'DB', 'DBstar', 'D', 'COP')
    ) %>% t() ->
      cvi_tib[i, 7:13]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 54,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  cvi_tib %>% return()
}

# Function to transform internal CVIs table to rank table
# Input is table with internal CVIs for TADPole clusterings (tibble) - output of int_cvi_p_clust_fun.
# Output is table with ranking of internal CVIs for TADPole clusterings (tibble).
int_rank_t_clust_fun <- function(int_cvi_t_clust_tib) {
  # Initializing tibble of ranks
  rank_tib <- tibble(
    code = int_cvi_t_clust_tib$code,
    distance = int_cvi_t_clust_tib$distance,
    centroid = int_cvi_t_clust_tib$centroid,
    lower_bound = int_cvi_t_clust_tib$lower_bound,
    clusters = int_cvi_t_clust_tib$clusters,
    cutoff = int_cvi_t_clust_tib$cutoff,
    Silhouette = rank(int_cvi_t_clust_tib$Silhouette),
    Score_Function = rank(int_cvi_t_clust_tib$Score_Function),
    Calinski_Harabasz = rank(int_cvi_t_clust_tib$Calinski_Harabasz),
    Davies_Bouldin = rank(-int_cvi_t_clust_tib$Davies_Bouldin),
    Mod_Davies_Bouldin = rank(-int_cvi_t_clust_tib$Mod_Davies_Bouldin),
    Dunn = rank(int_cvi_t_clust_tib$Dunn),
    COP = rank(-int_cvi_t_clust_tib$COP)
  )
  
  rank_tib %>%
    add_column(Total = rowSums(.[-c(1:6)])) %>%
    arrange(desc(Total)) %>%
    return()
}

#############################################################################

# Fuzzy clustering function.
# Input is normalized dataset (matrix).
# Output is table of partitional clusterings (tibble).
# Function creates 270 clusterings using various distance measures,
# two centroid calculation methods, number of clusters from 2 to 10 and
# convergence criterion in c(.001, .01, 1)
f_clust_fun <- function(dataset) {
  
  # Initializing tibble of clusterings
  clust_tib <- tibble(
    distance = rep(c('dtw_basic', 'lbk', 'lbi', 'sbd', 'gak'), each = 54),
    centroid = rep( c('fcm', 'fcmdd'), each = 27, times = 5),
    clusters = rep(2:10, times = 30),
    delta = rep(c(.001, .01, 1), each = 9, times = 10)
  )
  
  clustering_list = vector('list', length(18))
  
  for (i in 1:270) {
    tsclust(
      dataset,
      type = 'fuzzy',
      k = clust_tib$clusters[i],
      distance = clust_tib$distance[i],
      centroid = clust_tib$centroid[i],
      control = fuzzy_control(
        fuzziness = 2,
        iter.max = 1000,
        delta = clust_tib$delta[i]
      ),
      window.size = 5L
    ) -> clustering_list[[i]]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 270,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  clust_tib$clustering <- clustering_list
  
  clust_tib %>%
    add_column(
      code = paste(
        clust_tib$distance,
        clust_tib$centroid,
        clust_tib$clusters,
        sep = '_'
      ),
      .before = 'distance'
    ) %>%
    return()
}

# Function to compute internal CVIs for fuzzy clusterings
# Input is table of fuzzy clusterings (tibble) - output of f_clust_fun
# Output is table with internal CVIs for fuzzy clusterings (tibble).
int_cvi_f_clust_fun <- function(f_clust_tib) {
  
  # Initializing tibble of CVIs
  cvi_tib <- tibble(
    code = f_clust_tib$code,
    distance = f_clust_tib$distance,
    centroid = f_clust_tib$centroid,
    clusters = f_clust_tib$clusters,
    delta = f_clust_tib$delta,
    MPC = rep(NA, 270),
    Kwon = rep(NA, 270),
    Tang = rep(NA, 270),
    SC = rep(NA, 270),
    PBMF = rep(NA, 270)
  )
  
  for (i in 1:270) {
    cvi(
      f_clust_tib$clustering[[i]],
      type = c('MPC', 'K', 'T', 'SC', 'PBMF')
    ) %>% t() ->
      cvi_tib[i, 6:10]
    
    # Progress bar
    txtProgressBar(min = 0,
                   max = 270,
                   style = 3,
                   width = 50,
                   char = ':') %>%
      setTxtProgressBar(., i)
  }
  
  cvi_tib %>% return()
}

# Function to transform internal CVIs table to rank table
# Input is table with internal CVIs for fuzzy clusterings (tibble) - output of int_cvi_f_clust_fun.
# Output is table with ranking of internal CVIs for fuzzy clusterings (tibble).
int_rank_f_clust_fun <- function(int_cvi_f_clust_tib) {
  # Initializing tibble of ranks
  rank_tib <- tibble(
    code = int_cvi_f_clust_tib$code,
    distance = int_cvi_f_clust_tib$distance,
    centroid = int_cvi_f_clust_tib$centroid,
    clusters = int_cvi_f_clust_tib$clusters,
    delta = int_cvi_f_clust_tib$delta,
    MPC = rank(int_cvi_f_clust_tib$MPC),
    Kwon = rank(-int_cvi_f_clust_tib$Kwon),
    Tang = rank(-int_cvi_f_clust_tib$Tang),
    SC = rank(int_cvi_f_clust_tib$SC),
    PBMF = rank(int_cvi_f_clust_tib$PBMF)
  )
  
  rank_tib %>%
    add_column(Total = rowSums(.[-c(1:5)])) %>%
    arrange(desc(Total)) %>%
    return()
}













