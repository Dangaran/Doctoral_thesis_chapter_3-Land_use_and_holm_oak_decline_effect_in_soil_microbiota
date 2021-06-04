#### Summary network results ####
results_summary_network <- function(network) {
  round(c(length(V(network)), # Number of OTUs
          clique_num(network), # Larquest_clique
          length(cliques(network, min = 3, max = 3)), # Count number of 3 cliques
          length(cliques(network, min = 4, max = 4)), # Count number of 4 cliques
          length(cliques(network, min = 5, max = 5)), # Count number of 5 cliques
          diameter(network, directed = FALSE, unconnected = TRUE), # Diameter
          mean_distance(network, unconnected = TRUE), # Avg_path_lenght
          length(E(network))/length(V(network)), # Connectance
          transitivity(network, type = "global", vids = NULL, weights = NULL, isolates = c("NaN", "zero"))), 3) # Clustering_coefficient
}