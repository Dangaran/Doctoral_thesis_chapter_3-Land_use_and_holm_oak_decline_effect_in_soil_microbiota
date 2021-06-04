library(phyloseq)
library(ggplot2)
source('../accessory_functions/remove_bad_coverage_samples.R')

##%######################################################%##
#                                                          #
####                Prepare input file                  ####
#                                                          #
##%######################################################%##

illumina_fungi <- import_biom("Data/Fungi/otu_table.json", parseFunction = parse_taxonomy_default) 
colnames(tax_table(illumina_fungi)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Specie")

# Add sample variables
Sample.data <- read.csv(file = "Data/Bacteria/database.csv", row.names = 1, header = TRUE, sep = ";")
Sample.data$Manejo_Conteo <- factor(Sample.data$Manejo_Conteo, levels = c("Dehesa", "Open woodland", "Forest"))
Sample.data$Type <- factor(Sample.data$Type, levels = c("Healthy", "Affected", "Dead", "Bare soil"))

sample_data(illumina_fungi) <- Sample.data


##%######################################################%##
#                                                          #
####                 Sequencing filters                 ####
#                                                          #
##%######################################################%##

# Good's coverage: Singletons/total reads (First filter)
remove_seq <- remove.bad.coverage.samples(illumina_fungi)
remove_seq$coverage_plot

illumina_fungi_good <- remove.bad.coverage.samples(illumina_fungi,85)



#### Remove non-fungi ####
mycota <- subset_taxa(illumina_fungi_good, Phylum=="p__Basidiomycota" | Phylum=="p__Ascomycota" | Phylum=="p__Chytridiomycota" | Phylum=="p__Mucoromycota"
                      | Phylum=="p__Mortierellomycota" | Phylum=="p__GS19" | Phylum=="p__Rozellomycota" | Phylum=="p__Entomophthoromycota"
                      | Phylum=="p__Monoblepharomycota" | Phylum=="p__Calcarisporiellomycota" | Phylum=="p__Zoopagomycota" | Phylum=="p__Glomeromycota"
                      | Phylum=="p__Oomycota" | Phylum=="p__Kickxellomycota" | Phylum=="p__Olpidiomycota" | Phylum=="p__Neocallimastigomycota"
                      | Phylum=="p__Entorrhizomycota" | Phylum=="p__Blastocladiomycota" | Phylum=="p__Aphelidiomycota" | Phylum=="p__Bacillariophyta")

hist(sample_sums(mycota))
sum(sample_sums(mycota))
sort(sample_sums(mycota))

# Rename species s__unidentified 
tax_table(mycota)@.Data[,7] <- ifelse(tax_table(mycota)@.Data[,7] == "s__unidentified", sprintf("s__unidentified%d", 1:6625), tax_table(mycota)@.Data[,7])

mycota_glom <- tax_glom(mycota, taxrank = "Specie")


# Bokulich et al filter
sum(sample_sums(mycota_glom))*0.00005


# Remove taxa with less than 10 reads (Second filter) http://fungal-sequencing-methods-discussion.blogspot.dk/p/funguild.html (last sentence)
mycota_glom_filt <-  prune_taxa(taxa_sums(mycota_glom) > 253, mycota_glom)


# Remove samples with less than 1000 sequcens
mycota_glom_filt <- prune_samples(sample_sums(mycota_glom_filt)>2000, mycota_glom_filt)


test <- data.frame(tax_table(mycota_glom))
write.csv(test, "tax_table_glom.csv")
mycota_damaged <- subset_samples(mycota_glom, Treatment == "Damaged")

##%######################################################%##
#                                                          #
####      Info about physeq data after filtering        ####
#                                                          #
##%######################################################%##

sort(sample_sums(mycota_glom))
hist(sample_sums(mycota_glom))

# Mean, median, minimum and maximum reads per sample
summary(sample_sums(mycota_glom))
