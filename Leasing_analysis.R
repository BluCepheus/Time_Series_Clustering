setwd('C:\\Users\\U157766\\Desktop\\GitHub\\Leasing')

library(tidyverse)
library(dtwclust)

source('Leasing_functions.R')

# Importing and normalizing leasing dataset

read.csv('leasing_dataset.csv',
         sep = ';',
         dec = ',',
         header = T,
         row.names = 1
        ) %>%
  filter(
    !grepl(
      'BAWAG|BOS|BPS|CLIF|GBG|GRENKE|KBC|KOPEX|LE|LEASEPLAN|LHI|PSA|TOYOTA|VBRES',
      rownames(.)
    ) # Removing companies which reported in less than 20% of analysed quarters
  ) %>%
  zscore() ->
  leasing_dataset

#############################################################################

# Computional complexity of calculating the distance matrix

time_tib <- tibble(
  Measure_of_distance = c(
    'dtw', 'dtw2', 'dtw_basic', 'sdtw', 'lbk', 'lbi', 'sbd', 'gak'
  ),
  Time_elapsed = rep(NA, 8)
)

i = 1
for (
  method in c('dtw', 'dtw2', 'dtw_basic', 'sdtw', 'lbk', 'lbi', 'sbd', 'gak')
) {
  system.time(
    proxy::dist(
      leasing_dataset, method = method, window.size = 5L
    )
  ) ->
  time_tib$Time_elapsed[i]
  i = i + 1
}

time_tib %>%
  ggplot(., aes(x = Measure_of_distance, y = Time_elapsed)) +
  geom_bar(stat = 'identity') +
  labs(
    title = 'Time elapsed during calculation of the distance matrix',
    subtitle = 'breakdown by measure of distance',
    x = 'Measure of distance',
    y = 'Time elapsed in seconds',
    caption = 'Source: own study'
  ) +
  scale_y_continuous(
    labels = scales::number
  ) +
  theme_bw() ->
  time_bar_plot

#############################################################################

# Hierarchical clustering

# Calculating 48 hierarchical clusterings
hierarchical_clustering <- h_clust_fun(leasing_dataset)

# Internal CVIs
h_clust_i_cvi <- int_cvi_h_clust_fun(hierarchical_clustering)

# Internal CVIs ranking
h_clust_i_cvi_rank <- int_rank_h_clust_fun(h_clust_i_cvi)

# Inertion of hierarchical clusterings with different number of clusters
inert_h_clust_tib <- inert_h_clust_fun(leasing_dataset, hierarchical_clustering)

inert_h_clust_tib %>%
  filter(code %in% head(h_clust_i_cvi_rank$code, 4)) %>%
  ggplot(., aes(
    x = clusters,
    y = inertion,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Inertia of best hierarchical clusterings',
    subtitle = 'breakdown by clustering',
    y = 'Inertia',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  inert_h_clust_plot_1

inert_h_clust_tib %>%
  filter(
    code %in% filter(
      inert_h_clust_tib,
      inertion > .8,
      clusters <= 11
    )$code
  ) %>%
  ggplot(., aes(
    x = clusters,
    y = inertion,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Inertia of best lbk_dba clusterings',
    subtitle = 'breakdown by clustering',
    y = 'Inertia',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  inert_h_clust_plot_2

# External CVIs
h_clust_e_cvi <- ext_cvi_h_clust_fun(
  hierarchical_clustering %>% filter(code == 'lbk_dba_average'),
  hierarchical_clustering %>% filter(
    code %in% c(
      'dtw_basic_dba_mcquitty',
      'dtw_basic_dba_average',
      'lbk_dba_median',
      'lbk_dba_ward.D',
      'lbk_dba_ward.D2',
      'lbk_dba_mcquitty',
      'lbk_dba_complete'
    )
  )
)

h_clust_e_cvi %>%
  ggplot(., aes(
    x = clusters,
    y = Rand,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Comparison of chosen clusterings vs lbk_dba_average',
    subtitle = 'according to Rand index',
    y = 'Rand index',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  rand_h_clust_plot

h_clust_e_cvi %>%
  ggplot(., aes(
    x = clusters,
    y = Adjusted_Rand,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Comparison of chosen clusterings vs lbk_dba_average',
    subtitle = 'according to adjusted Rand index',
    y = 'Adjusted Rand index',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  adj_rand_h_clust_plot

h_clust_e_cvi %>%
  ggplot(., aes(
    x = clusters,
    y = Jaccard,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Comparison of chosen clusterings vs lbk_dba_average',
    subtitle = 'according to Jaccard index',
    y = 'Jaccard index',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  jaccard_h_clust_plot

h_clust_e_cvi %>%
  ggplot(., aes(
    x = clusters,
    y = Fowlkes_Mallows,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Comparison of chosen clusterings vs lbk_dba_average',
    subtitle = 'according to Fowlkes-Mallows index',
    y = 'Fowlkes-Mallows index',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  fowlkes_mallows_h_clust_plot

h_clust_e_cvi %>%
  ggplot(., aes(
    x = clusters,
    y = Variation_of_Information,
    color = code,
    fill = code
  )) +
  geom_line() +
  geom_point() +
  labs(
    title = 'Comparison of chosen clusterings vs lbk_dba_average',
    subtitle = 'according to Variation of Information index',
    y = 'Variation of Information index',
    x = 'Number of clusters',
    caption = 'Source: own study'
  ) +
  scale_color_discrete(name = 'Clusterings') +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() ->
  variation_information_h_clust_plot

# Dendrograms

hierarchical_clustering$clustering[[13]] %>%
  as.dendrogram() %>%
  dendextend::set('labels_col', k = 11) %>%
  dendextend::set('labels_cex', .7) %>%
  dendextend::set('branches_k_color', k = 11) %>%
  plot(type = 'rectangle', main = 'lbk_dba_average clustering with 11 groups')

hierarchical_clustering$clustering[[13]] %>%
  as.dendrogram() %>%
  dendextend::rect.dendrogram(k = 11, border = 1, lty = 2, lwd = 1)

hierarchical_clustering$clustering[[13]] %>% plot(hang = -1)
rect.hclust(hierarchical_clustering$clustering[[13]], k = 11)

# Resulting clustering

hierarchical_res_tib <- tibble(
  cluster = cutree(hierarchical_clustering$clustering[[13]], k = 11),
  company = leasing_dataset %>% rownames()
) %>%
  arrange(cluster)

#############################################################################

# Partitional clustering

# Calculating 72 partitional clusterings
partitional_clustering <- p_clust_fun(leasing_dataset)

# Internal CVIs
p_clust_i_cvi <- int_cvi_p_clust_fun(partitional_clustering)

# Internal CVIs ranking
p_clust_i_cvi_rank <- int_rank_p_clust_fun(p_clust_i_cvi)

p_clust_i_cvi %>%
  filter(code %in% head(p_clust_i_cvi_rank$code, 5))

# Comparison of lbk_dba_average_11 and gak_mean_8
cvi(
  cutree(hierarchical_clustering$clustering[[13]], k = 11),
  partitional_clustering$clustering[[61]]@cluster,
  type = 'external'
)

# Resulting clustering
partitional_res_tib <- tibble(
  cluster = partitional_clustering$clustering[[61]]@cluster,
  company = leasing_dataset %>% rownames()
) %>%
  arrange(cluster)

#############################################################################

# TADPole clustering

# Calculating 90 TADPole clusterings
tadpole_clustering <- t_clust_fun(leasing_dataset)

# Internal CVIs
t_clust_i_cvi <- int_cvi_t_clust_fun(tadpole_clustering)

# Internal CVIs ranking
t_clust_i_cvi_rank <- int_rank_t_clust_fun(t_clust_i_cvi)

# dtw_basic_dba_lbk_10 vs dtw_basic_dba_lbi_8 
cvi(
  tadpole_clustering$clustering[[25]]@cluster,
  tadpole_clustering$clustering[[46]]@cluster,
  type = 'external'
)

# Comparison of chosen clusterings
crisp_clust_e_cvi <- tibble(
  Clustering = c(
    rep('lbk_dba_average_11', 3), rep('gak_mean_8', 2), 'dtw_basic_dba_lbk_10'
  ),
  Comparison = c(
    'gak_mean_8', 'dtw_basic_dba_lbk_10', 'dtw_basic_dba_lbi_8',
    'dtw_basic_dba_lbk_10', rep('dtw_basic_dba_lbi_8', 2)
  ),
  ARI = rep(NA, 6),
  RI = rep(NA, 6),
  J = rep(NA, 6),
  FM = rep(NA, 6),
  VI = rep(NA, 6)
)

cvi(
  cutree(hierarchical_clustering$clustering[[13]], k = 11),
  partitional_clustering$clustering[[61]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[1, 3:7]

cvi(
  cutree(hierarchical_clustering$clustering[[13]], k = 11),
  tadpole_clustering$clustering[[25]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[2, 3:7]

cvi(
  cutree(hierarchical_clustering$clustering[[13]], k = 11),
  tadpole_clustering$clustering[[46]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[3, 3:7]

cvi(
  partitional_clustering$clustering[[61]]@cluster,
  tadpole_clustering$clustering[[25]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[4, 3:7]

cvi(
  partitional_clustering$clustering[[61]]@cluster,
  tadpole_clustering$clustering[[46]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[5, 3:7]

cvi(
  tadpole_clustering$clustering[[25]]@cluster,
  tadpole_clustering$clustering[[46]]@cluster,
  type = c('RI', 'ARI', 'J', 'FM', 'VI')
) %>% t() -> crisp_clust_e_cvi[6, 3:7]

# Resulting clustering
tadpole_res_tib <- tibble(
  cluster = tadpole_clustering$clustering[[46]]@cluster,
  company = leasing_dataset %>% rownames()
) %>%
  arrange(cluster)

#############################################################################

# Fuzzy clustering

# Calculating 270 fuzzy clusterings
fuzzy_clustering <- f_clust_fun(leasing_dataset)

# Internal CVIs
f_clust_i_cvi <- int_cvi_f_clust_fun(fuzzy_clustering)

# Internal CVIs ranking
f_clust_i_cvi_rank <- int_rank_f_clust_fun(f_clust_i_cvi)

f_clust_i_cvi %>%
  filter(code %in% head(f_clust_i_cvi_rank$code, 4))

# External CVIs
f_clust_e_cvi <- tibble(
  Clustering = c(
    rep('lbi_fcmdd_4', 3), rep('lbi_fcmdd_3', 2), 'sbd_fcmdd_3'
  ),
  Comparison = c(
    'lbi_fcmdd_3', 'sbd_fcmdd_3', 'gak_fcmdd_3', 'sbd_fcmdd_3', rep('gak_fcmdd_3', 2)
  ),
  ARI = rep(NA, 6),
  RI = rep(NA, 6),
  VI = rep(NA, 6),
  NMIM = rep(NA, 6)
)

cvi(
  fuzzy_clustering$clustering[[147]]@fcluster,
  fuzzy_clustering$clustering[[155]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[1, 3:6]

cvi(
  fuzzy_clustering$clustering[[147]]@fcluster,
  fuzzy_clustering$clustering[[200]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[2, 3:6]

cvi(
  fuzzy_clustering$clustering[[147]]@fcluster,
  fuzzy_clustering$clustering[[254]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[3, 3:6]

cvi(
  fuzzy_clustering$clustering[[155]]@fcluster,
  fuzzy_clustering$clustering[[200]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[4, 3:6]

cvi(
  fuzzy_clustering$clustering[[155]]@fcluster,
  fuzzy_clustering$clustering[[254]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[5, 3:6]

cvi(
  fuzzy_clustering$clustering[[200]]@fcluster,
  fuzzy_clustering$clustering[[254]]@fcluster,
  type = c('external')
) %>% t() -> f_clust_e_cvi[6, 3:6]

# Resulting clustering
fuzzy_res_tib_1 <- tibble(
  cluster = fuzzy_clustering$clustering[[147]]@fcluster,
  company = leasing_dataset %>% rownames()
) %>%
  arrange(cluster)

fuzzy_res_tib_2 <- tibble(
  cluster = fuzzy_clustering$clustering[[254]]@fcluster,
  company = leasing_dataset %>% rownames()
) %>%
  arrange(cluster)
