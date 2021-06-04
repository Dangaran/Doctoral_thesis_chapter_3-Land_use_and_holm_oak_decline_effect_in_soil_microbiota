#### Good's coverage function ####
remove.bad.coverage.samples <- function(physeq, cov.threshold = NA){
  require('phyloseq') # phyloseq
  require('ggplot2') # ggplot2
  require('dplyr') # dplyr
  
  # Good function (https://rdrr.io/github/jfq3/QsRutils/src/R/goods.R):
  # Calculates Good's coverage from a community data matrix with samples as rows and OTUs as columns.
  # A table with the headings number of singletons, number of sequences, and Good's coverage for each sample in rows.
  # references Good, I. J. 1953. The Population Frequencies of Species and the Estimation of Population Parameters. Biometrika 40:237-264.
  
  goods <-
    function(com){
      no.seqs <- rowSums(com)
      sing <- com==1
      no.sing <- apply(sing, 1, sum)
      goods <- 100*(1-no.sing/no.seqs)
      goods.sum <- cbind(no.sing, no.seqs, goods)
      goods.sum <- as.data.frame(goods.sum)
      return(goods.sum)
    }
  
  otu.table <- otu_table(physeq) # Extract otu table
  otu.table.transpose <- t(otu.table) # Transpose otu table
  good.coverage <- goods(otu.table.transpose) # Use the formula
  good.coverage.sort <- good.coverage[order(good.coverage$goods),] # Order the samples from lower to highes coverage
  
  samples_85 <- sum(good.coverage$goods < 85) # Number of samples lost with 85% of good's coverage
  samples_80 <- sum(good.coverage$goods < 80) # Number of samples lost with 80% of good's coverage
  samples_75 <- sum(good.coverage$goods < 75) # Number of samples lost with 75% of good's coverage
  
  library(ggplot2) # Plot results
  plot_coverage <- ggplot(good.coverage,aes(row.names(good.coverage), goods)) + 
    geom_point() + labs(x = "Samples", y = "% Good's coverage") +
    geom_hline(yintercept = 85, color = "green") +
    geom_hline(yintercept = 80, color = "red")
  
  if(!missing(cov.threshold) && is.numeric(cov.threshold) == TRUE) {
    good.coverage$SampleID <- rownames(good.coverage)
    
    high_coverage_samples_data <- filter(good.coverage, goods > cov.threshold)
    
    samples_to_keep <- high_coverage_samples_data$SampleID
    
    phy_remove <- prune_samples(samples_to_keep, physeq)
    
    return(phy_remove)
    
  } else {
    print(paste("Samples lost with 85% of good's coverage:", samples_85,sep = " "))
    print(paste("Samples lost with 80% of good's coverage:", samples_80,sep = " "))
    print(paste("Samples lost with 75% of good's coverage:", samples_75,sep = " "))
    output <- list("coverage_data" = good.coverage.sort, "coverage_plot" = plot_coverage)
    
    return(output)
  }
}