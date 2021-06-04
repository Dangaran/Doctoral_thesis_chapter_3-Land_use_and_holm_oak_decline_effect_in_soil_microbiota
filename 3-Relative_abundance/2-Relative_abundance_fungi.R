##%######################################################%##
#                                                          #
####           Calculate relative abundance             ####
#                                                          #
##%######################################################%##
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(reshape2)

# Forests ----
# First agglomerate by phylum to define the phyla you want to keep.
fungi_damaged_for_merged <- merge_samples(fungi_damaged_for, "Type") # Agglomerate phyla by management
physeqPhylum_for = tax_glom(fungi_damaged_for_merged, "Phylum")
physeqPhylumRA_for = transform_sample_counts(physeqPhylum_for, function(x) x/sum(x))
physeqPhylumRAF_for = filter_taxa(physeqPhylumRA_for, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abuncance

# Phyla with more than 1% of relative abundande
keepPhyla_for = get_taxa_unique(physeqPhylumRAF_for, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_for <- data.frame(tax_table(physeqPhylum_for)[,1:2])
tax_table_for$Phylum <- ifelse(as.vector(tax_table_for[,2]) %in% keepPhyla_for == TRUE, as.vector(tax_table_for[,2]), "Other")
tax_table(physeqPhylum_for)[,2] <- tax_table_for$Phylum 

fungi_phyla_for <- physeqPhylum_for %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


fungi_phyla_for$Sample <- factor(fungi_phyla_for$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
fungi_phyla_for$Phylum <- factor(fungi_phyla_for$Phylum, levels = c("p__Ascomycota", "p__Basidiomycota", "p__Mortierellomycota", "p__Mucoromycota", "Other"))


ggplot(fungi_phyla_for, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of fungil communities in Forest according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))

# Calculate significance difference between tree health
physeqPhylum_for

relab_fungi_for <-  transform_sample_counts(fungi_damaged_for, function(x) x / sum(x)*100)
relab_fungi_for = tax_glom(relab_fungi_for, "Phylum") # Unir por filo
relab_fungi_for

# ANOVA
relab_fungi_for_otu_tab <- data.frame(otu_table(relab_fungi_for))
rownames(relab_fungi_for_otu_tab) == rownames(tax_table(relab_fungi_for)) # Check order of rows

relab_fungi_for_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_for)[,2]) # Create Phyla column

library(reshape2)
relab_fungi_for_sig <- melt(relab_fungi_for_otu_tab)
relab_fungi_for_sig$Type <- ifelse(grepl('S$', relab_fungi_for_sig$variable), "Healthy", ifelse(grepl('D$', relab_fungi_for_sig$variable), "Affected", ifelse(grepl('M$', relab_fungi_for_sig$variable), "Dead", ifelse(grepl('C$', relab_fungi_for_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_fungi_for_sig, mean)

# Function to calculate significance per phylum
significant_phylum <- function(physeq, Phylum){
  phylum <- subset(physeq, Phyla %in% c(Phylum))
  anova <- aov(value ~ Type, phylum)
  
  anova_results <- summary(anova)
  tukey_results <- TukeyHSD(anova)
  output <- list(anova_results=anova_results, tukey_results=tukey_results)
  output
  
}

significant_phylum(relab_fungi_for_sig, "p__Ascomycota") # NO SIG
significant_phylum(relab_fungi_for_sig, "p__Basidiomycota") # NO SIG
significant_phylum(relab_fungi_for_sig, "p__Mortierellomycota") # NO SIG
significant_phylum(relab_fungi_for_sig, "p__Mucoromycota") # NO SIG


# Dehesa -----
# First agglomerate by phylum to define the phyla you want to keep.
fungi_damaged_deh_merged <- merge_samples(fungi_damaged_deh, "Type") # Agglomerate phyla by management

physeqPhylum_deh = tax_glom(fungi_damaged_deh_merged, "Phylum")
physeqPhylumRA_deh = transform_sample_counts(physeqPhylum_deh, function(x) x/sum(x))
physeqPhylumRAF_deh = filter_taxa(physeqPhylumRA_deh, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abuncance

# Phyla with more than 1% of relative abundande
keepPhyla_deh = get_taxa_unique(physeqPhylumRAF_deh, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_deh <- data.frame(tax_table(physeqPhylum_deh)[,1:2])
tax_table_deh$Phylum <- ifelse(as.vector(tax_table_deh[,2]) %in% keepPhyla_deh == TRUE, as.vector(tax_table_deh[,2]), "Other")
tax_table(physeqPhylum_deh)[,2] <- tax_table_deh$Phylum 

fungi_phyla_deh <- physeqPhylum_deh %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

fungi_phyla_deh$Sample <- factor(fungi_phyla_deh$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
fungi_phyla_deh$Phylum <- factor(fungi_phyla_deh$Phylum, levels = c("p__Ascomycota", "p__Basidiomycota", "p__Mucoromycota", "Other"))

ggplot(fungi_phyla_deh, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of fungil communities in Dehesa according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))


# Calculate significance difference between tree health
physeqPhylum_deh

relab_fungi_deh <-  transform_sample_counts(fungi_damaged_deh, function(x) x / sum(x)*100)
relab_fungi_deh = tax_glom(relab_fungi_deh, "Phylum") # Unir por filo
relab_fungi_deh

# ANOVA
relab_fungi_deh_otu_tab <- data.frame(otu_table(relab_fungi_deh))
rownames(relab_fungi_deh_otu_tab) == rownames(tax_table(relab_fungi_deh)) # Check order of rows

relab_fungi_deh_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_deh)[,2]) # Create Phyla column

library(reshape2)
relab_fungi_deh_sig <- melt(relab_fungi_deh_otu_tab)
relab_fungi_deh_sig$Type <- ifelse(grepl('S$', relab_fungi_deh_sig$variable), "Healthy", ifelse(grepl('D$', relab_fungi_deh_sig$variable), "Affected", ifelse(grepl('M$', relab_fungi_deh_sig$variable), "Dead", ifelse(grepl('C$', relab_fungi_deh_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_fungi_deh_sig, mean)


significant_phylum(relab_fungi_deh_sig, "p__Ascomycota") # Sig: H-BS 
significant_phylum(relab_fungi_deh_sig, "p__Basidiomycota") # Sig: H-BS 
significant_phylum(relab_fungi_deh_sig, "p__Mucoromycota") # Sig: A-BS 


# Open woodlands -------
# First agglomerate by phylum to define the phyla you want to keep.
fungi_damaged_opw_merged <- merge_samples(fungi_damaged_opw, "Type") # Agglomerate phyla by management

physeqPhylum_opw = tax_glom(fungi_damaged_opw_merged, "Phylum")
physeqPhylumRA_opw = transform_sample_counts(physeqPhylum_opw, function(x) x/sum(x))
physeqPhylumRAF_opw = filter_taxa(physeqPhylumRA_opw, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abuncance

# Phyla with more than 1% of relative abundande
keepPhyla_opw = get_taxa_unique(physeqPhylumRAF_opw, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_opw <- data.frame(tax_table(physeqPhylum_opw)[,1:2])
tax_table_opw$Phylum <- ifelse(as.vector(tax_table_opw[,2]) %in% keepPhyla_opw == TRUE, as.vector(tax_table_opw[,2]), "Other")
tax_table(physeqPhylum_opw)[,2] <- tax_table_opw$Phylum 

fungi_phyla_opw <- physeqPhylum_opw %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

fungi_phyla_opw$Sample <- factor(fungi_phyla_opw$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
fungi_phyla_opw$Phylum <- factor(fungi_phyla_opw$Phylum, levels = c("p__Ascomycota", "p__Basidiomycota", "p__Mortierellomycota", "p__Mucoromycota", "Other"))

ggplot(fungi_phyla_opw, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of fungil communities in Open woodland according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))

# Calculate significance difference between tree health
physeqPhylum_opw

relab_fungi_opw <-  transform_sample_counts(fungi_damaged_opw, function(x) x / sum(x)*100)
relab_fungi_opw = tax_glom(relab_fungi_opw, "Phylum") # Unir por filo
relab_fungi_opw

# ANOVA
relab_fungi_opw_otu_tab <- data.frame(otu_table(relab_fungi_opw))
rownames(relab_fungi_opw_otu_tab) == rownames(tax_table(relab_fungi_opw)) # Check order of rows

relab_fungi_opw_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_opw)[,2]) # Create Phyla column

library(reshape2)
relab_fungi_opw_sig <- melt(relab_fungi_opw_otu_tab)
relab_fungi_opw_sig$Type <- ifelse(grepl('S$', relab_fungi_opw_sig$variable), "Healthy", ifelse(grepl('D$', relab_fungi_opw_sig$variable), "Affected", ifelse(grepl('M$', relab_fungi_opw_sig$variable), "Dead", ifelse(grepl('C$', relab_fungi_opw_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_fungi_opw_sig, mean)


significant_phylum(relab_fungi_opw_sig, "p__Ascomycota") # NO SIG
significant_phylum(relab_fungi_opw_sig, "p__Basidiomycota") # NO SIG
significant_phylum(relab_fungi_opw_sig, "p__Mortierellomycota") # NO SIG
significant_phylum(relab_fungi_opw_sig, "p__Mucoromycota") # NO SIG

