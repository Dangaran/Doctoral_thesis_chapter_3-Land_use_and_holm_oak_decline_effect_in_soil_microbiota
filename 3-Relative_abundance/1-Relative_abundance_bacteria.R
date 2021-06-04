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
illumina_damaged_for_merged <- merge_samples(illumina_damaged_for, "Type") # Agglomerate phyla by management
physeqPhylum_for = tax_glom(illumina_damaged_for_merged, "Phylum")
physeqPhylumRA_for = transform_sample_counts(physeqPhylum_for, function(x) x/sum(x))
physeqPhylumRAF_for = filter_taxa(physeqPhylumRA_for, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abundance

# Phyla with more than 1% of relative abundance
keepPhyla_for = get_taxa_unique(physeqPhylumRAF_for, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_for <- data.frame(tax_table(physeqPhylum_for)[,1:2])
tax_table_for$Phylum <- ifelse(as.vector(tax_table_for[,2]) %in% keepPhyla_for == TRUE, as.vector(tax_table_for[,2]), "Other")
tax_table_for$Phylum <- ifelse(as.vector(tax_table_for[,2]) == "Crenarchaeota", "Other",as.vector(tax_table_for[,2]))

tax_table(physeqPhylum_for)[,2] <- tax_table_for$Phylum 
tax_table(physeqPhylum_for)[,1] <- "Bacteria"


illumina_phyla_for <- physeqPhylum_for %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


illumina_phyla_for$Sample <- factor(illumina_phyla_for$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
illumina_phyla_for$Phylum <- factor(illumina_phyla_for$Phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Verrucomicrobia", "Other"))


ggplot(illumina_phyla_for, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of bacterial communities in Forest according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))

# Calculate significance difference between tree health
physeqPhylum_for

relab_bact_for <-  transform_sample_counts(illumina_damaged_for, function(x) x / sum(x)*100)
relab_bact_for = tax_glom(relab_bact_for, "Phylum") 
relab_bact_for

# ANOVA
relab_bact_for_otu_tab <- data.frame(otu_table(relab_bact_for))
rownames(relab_bact_for_otu_tab) == rownames(tax_table(relab_bact_for)) # Check order of rows

relab_bact_for_otu_tab$Phyla <- as.vector(tax_table(relab_bact_for)[,2]) # Create Phyla column


relab_bact_for_sig <- melt(relab_bact_for_otu_tab)
relab_bact_for_sig$Type <- ifelse(grepl('S$', relab_bact_for_sig$variable), "Healthy", ifelse(grepl('D$', relab_bact_for_sig$variable), "Affected", ifelse(grepl('M$', relab_bact_for_sig$variable), "Dead", ifelse(grepl('C$', relab_bact_for_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_bact_for_sig, mean)

# Function to calculate significance per phylum
significant_phylum <- function(physeq, Phylum){
  phylum <- subset(physeq, Phyla %in% c(Phylum))
  anova <- aov(value ~ Type, phylum)
  
  anova_results <- summary(anova)
  tukey_results <- TukeyHSD(anova)
  output <- list(anova_results=anova_results, tukey_results=tukey_results)
  output
  
}

significant_phylum(relab_bact_for_sig, "Acidobacteria") # Sig: H-BS
significant_phylum(relab_bact_for_sig, "Actinobacteria") # Sig: H-BS // H-D
significant_phylum(relab_bact_for_sig, "Bacteroidetes") # NO SIG
significant_phylum(relab_bact_for_sig, "Chloroflexi") # NO SIG
significant_phylum(relab_bact_for_sig, "Firmicutes") # NO SIG
significant_phylum(relab_bact_for_sig, "Gemmatimonadetes") # NO SIG
significant_phylum(relab_bact_for_sig, "Planctomycetes") # Sig: H-BS
significant_phylum(relab_bact_for_sig, "Proteobacteria") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_for_sig, "Verrucomicrobia") # Sig: H-BS


# Dehesa ----
# First agglomerate by phylum to define the phyla you want to keep.
illumina_damaged_deh_merged <- merge_samples(illumina_damaged_deh, "Type") # Agglomerate phyla by management

physeqPhylum_deh = tax_glom(illumina_damaged_deh_merged, "Phylum")
physeqPhylumRA_deh = transform_sample_counts(physeqPhylum_deh, function(x) x/sum(x))
physeqPhylumRAF_deh = filter_taxa(physeqPhylumRA_deh, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abundance

# Phyla with more than 1% of relative abundance
keepPhyla_deh = get_taxa_unique(physeqPhylumRAF_deh, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_deh <- data.frame(tax_table(physeqPhylum_deh)[,1:2])
tax_table_deh$Phylum <- ifelse(as.vector(tax_table_deh[,2]) %in% keepPhyla_deh == TRUE, as.vector(tax_table_deh[,2]), "Other")
tax_table(physeqPhylum_deh)[,2] <- tax_table_deh$Phylum 
tax_table(physeqPhylum_deh)[,1] <- "Bacteria"


illumina_phyla_deh <- physeqPhylum_deh %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


illumina_phyla_deh$Sample <- factor(illumina_phyla_deh$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
illumina_phyla_deh$Phylum <- factor(illumina_phyla_deh$Phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Verrucomicrobia", "Other"))


ggplot(illumina_phyla_deh, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of bacterial communities in Dehesa according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))

# Calculate significance difference between tree health
physeqPhylum_deh

relab_bact_deh <-  transform_sample_counts(illumina_damaged_deh, function(x) x / sum(x)*100)
relab_bact_deh = tax_glom(relab_bact_deh, "Phylum")
relab_bact_deh

# ANOVA
relab_bact_deh_otu_tab <- data.frame(otu_table(relab_bact_deh))
rownames(relab_bact_deh_otu_tab) == rownames(tax_table(relab_bact_deh)) # Check order of rows

relab_bact_deh_otu_tab$Phyla <- as.vector(tax_table(relab_bact_deh)[,2]) # Create Phyla column


relab_bact_deh_sig <- melt(relab_bact_deh_otu_tab)
relab_bact_deh_sig$Type <- ifelse(grepl('S$', relab_bact_deh_sig$variable), "Healthy", ifelse(grepl('D$', relab_bact_deh_sig$variable), "Affected", ifelse(grepl('M$', relab_bact_deh_sig$variable), "Dead", ifelse(grepl('C$', relab_bact_deh_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_bact_deh_sig, mean)

significant_phylum(relab_bact_deh_sig, "Acidobacteria") # Sig: H-BS // A-BS
significant_phylum(relab_bact_deh_sig, "Actinobacteria") # NO SIG
significant_phylum(relab_bact_deh_sig, "Bacteroidetes") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_deh_sig, "Chloroflexi") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_deh_sig, "Firmicutes") # NO SIG
significant_phylum(relab_bact_deh_sig, "Gemmatimonadetes") # NO SIG
significant_phylum(relab_bact_deh_sig, "Planctomycetes") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_deh_sig, "Proteobacteria") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_deh_sig, "Verrucomicrobia") # NO SIG


# Open woodland ----
# First agglomerate by phylum to define the phyla you want to keep.
illumina_damaged_opw_merged <- merge_samples(illumina_damaged_opw, "Type") # Agglomerate phyla by management

physeqPhylum_opw = tax_glom(illumina_damaged_opw_merged, "Phylum")
physeqPhylumRA_opw = transform_sample_counts(physeqPhylum_opw, function(x) x/sum(x))
physeqPhylumRAF_opw = filter_taxa(physeqPhylumRA_opw, function(x) mean(x) > 0.01, TRUE) # Phyla with more than 1% of abuncance

# Phyla with more than 1% of relative abundande
keepPhyla_opw = get_taxa_unique(physeqPhylumRAF_opw, "Phylum")

# subset to just the phyla that you want, using the original phyloseq object
tax_table_opw <- data.frame(tax_table(physeqPhylum_opw)[,1:2])
tax_table_opw$Phylum <- ifelse(as.vector(tax_table_opw[,2]) %in% keepPhyla_opw == TRUE, as.vector(tax_table_opw[,2]), "Other")
tax_table(physeqPhylum_opw)[,2] <- tax_table_opw$Phylum 
tax_table(physeqPhylum_opw)[,1] <- "Bacteria"


illumina_phyla_opw <- physeqPhylum_opw %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {100*x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


illumina_phyla_opw$Sample <- factor(illumina_phyla_opw$Sample, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
illumina_phyla_opw$Phylum <- factor(illumina_phyla_opw$Phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Verrucomicrobia", "Other"))


ggplot(illumina_phyla_opw, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + labs(x = "Tree health", y = "Relative Abundance (%)", title = "Relative abundance of bacterial communities in Open woodland according to tree health")+
  geom_text(aes(label = round(Abundance, 2)), size = 3, position = position_stack(vjust = 0.5)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA))

# Calculate significance difference between tree health
physeqPhylum_opw

relab_bact_opw <-  transform_sample_counts(illumina_damaged_opw, function(x) x / sum(x)*100)
relab_bact_opw = tax_glom(relab_bact_opw, "Phylum") # Unir por filo
relab_bact_opw

# ANOVA
relab_bact_opw_otu_tab <- data.frame(otu_table(relab_bact_opw))
rownames(relab_bact_opw_otu_tab) == rownames(tax_table(relab_bact_opw)) # Check order of rows

relab_bact_opw_otu_tab$Phyla <- as.vector(tax_table(relab_bact_opw)[,2]) # Create Phyla column


relab_bact_opw_sig <- melt(relab_bact_opw_otu_tab)
relab_bact_opw_sig$Type <- ifelse(grepl('S$', relab_bact_opw_sig$variable), "Healthy", ifelse(grepl('D$', relab_bact_opw_sig$variable), "Affected", ifelse(grepl('M$', relab_bact_opw_sig$variable), "Dead", ifelse(grepl('C$', relab_bact_opw_sig$variable), "Bare soil", "Other"))))
aggregate(value ~ Phyla + Type, relab_bact_opw_sig, mean)


significant_phylum(relab_bact_opw_sig, "Acidobacteria") # NO SIG
significant_phylum(relab_bact_opw_sig, "Actinobacteria") # NO SIG
significant_phylum(relab_bact_opw_sig, "Bacteroidetes") # NO SIG
significant_phylum(relab_bact_opw_sig, "Chloroflexi") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_opw_sig, "Firmicutes") # Sig: A-BS 
significant_phylum(relab_bact_opw_sig, "Gemmatimonadetes") # NO SIG
significant_phylum(relab_bact_opw_sig, "Planctomycetes") # Sig: A-BS // D-BS
significant_phylum(relab_bact_opw_sig, "Proteobacteria") # Sig: H-BS // A-BS // D-BS
significant_phylum(relab_bact_opw_sig, "Verrucomicrobia") # NO SIG

