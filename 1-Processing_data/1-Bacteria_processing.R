library(phyloseq)
library(ggplot2)
source('../accessory_functions/remove_bad_coverage_samples.R')

##%######################################################%##
#                                                          #
####                Recommended websites                ####
#                                                          #
##%######################################################%##

# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# https://github.com/microbiome/microbiome/tree/master/R
# http://genoweb.toulouse.inra.fr/~formation/15_FROGS/5-June2016/FROGS_phyloseq_23062016.pdf
# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html

##%######################################################%##
#                                                          #
####                 Prepare input file                 ####
#                                                          #
##%######################################################%##

# Import biom file
illumina_16s <- import_biom("Data/Bacteria/otu_table.json", # https://twbattaglia.gitbooks.io/introduction-to-qiime/content/phyloseq.html
                                 treefilename = "Data/Bacteria/rep_set.tre",
                                 parseFunction = parse_taxonomy_greengenes) # Reason to use that parseFunction --> https://github.com/joey711/phyloseq/issues/302 (michberr 3/3/14)

# Add sample variables
Sample.data <- read.csv(file = "Data/Bacteria/database.csv", row.names = 1, header = TRUE, sep = ";")
Sample.data$Manejo_Conteo <- factor(Sample.data$Manejo_Conteo, levels = c("Dehesa", "Open woodland", "Forest"))
Sample.data$Type <- factor(Sample.data$Type, levels = c("Healthy", "Affected", "Dead", "Bare soil"))

sample_data(illumina_16s) <- Sample.data


##%######################################################%##
#                                                          #
####                Sequencing filters                  ####
#                                                          #
##%######################################################%##

# Good's coverage: Singletons/total reads (First filter)
remove_seq <- remove.bad.coverage.samples(illumina_16s)
remove_seq$coverage_plot

illumina_16s_good <- remove.bad.coverage.samples(illumina_16s, 80)

# Bokulich et al filter
sum(sample_sums(illumina_16s_good))*0.00005

# Remove taxa with less than 10 reads (Second filter)
illumina_16s_10seq <-  prune_taxa(taxa_sums(illumina_16s_good) > 206, illumina_16s_good) #Cambiar a 10

##%######################################################%##
#                                                          #
####      Info about physeq data after filtering        ####
#                                                          #
##%######################################################%##

illumina_16s_10seq


# Number of sequences per sample
sort(sample_sums(illumina_16s_10seq))
hist(sample_sums(illumina_16s_10seq))

# Mean, median, minimum and maximum reads per sample
summary(sample_sums(illumina_16s_10seq))
