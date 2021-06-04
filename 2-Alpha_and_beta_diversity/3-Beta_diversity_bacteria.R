##%######################################################%##
#                                                          #
####                  Beta diversity                    ####
#                                                          #
##%######################################################%##
# For each land-use
library(vegan)
illumina_16s_good_damaged

illumina_damaged_raref <- rarefy_even_depth(illumina_16s_good_damaged)

ordu <-  ordinate(illumina_damaged_raref, "NMDS", "wunifrac")
plot_ordination(illumina_damaged_raref, ordu, color="Manejo_Conteo")  + 
  stat_ellipse(type = "t", linetype = 2)

# Dehesa vs Forest
illumina_damaged_DehvsFor <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo != "Open woodland")

sort(sample_sums(illumina_damaged_DehvsFor))

illumina_damaged_DehvsFor_raref <- rarefy_even_depth(illumina_damaged_DehvsFor)

illumina_damaged_DehvsFor_raref_OTU <- phyloseq::distance(illumina_damaged_DehvsFor_raref, method = "wunifrac")
illumina_damaged_DehvsFor_raref_data <- data.frame(sample_data(illumina_damaged_DehvsFor_raref))


# Andonis 
set.seed(123)
adonis(illumina_damaged_DehvsFor_raref_OTU ~ Manejo_Conteo, illumina_damaged_DehvsFor_raref_data,permutations = 999)

# Dehesa vs Open woodland
illumina_damaged_Dehvsopen <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo != "Forest")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(illumina_damaged_Dehvsopen))

illumina_damaged_Dehvsopen_raref <- rarefy_even_depth(illumina_damaged_Dehvsopen)

illumina_damaged_Dehvsopen_raref_OTU <-  phyloseq::distance(illumina_damaged_Dehvsopen_raref, method = "wunifrac")
illumina_damaged_Dehvsopen_raref_data <- data.frame(sample_data(illumina_damaged_Dehvsopen_raref))


# Andonis 
set.seed(123)
adonis(illumina_damaged_Dehvsopen_raref_OTU ~ Manejo_Conteo, illumina_damaged_Dehvsopen_raref_data,permutations = 999)


# Forest vs Open woodland
illumina_damaged_forvsopen <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo != "Dehesa")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(illumina_damaged_forvsopen))

illumina_damaged_forvsopen_raref <- rarefy_even_depth(illumina_damaged_forvsopen)

illumina_damaged_forvsopen_raref_OTU <- phyloseq::distance(illumina_damaged_forvsopen_raref, method = "wunifrac")
illumina_damaged_forvsopen_raref_data <- data.frame(sample_data(illumina_damaged_forvsopen_raref))


# Andonis 
set.seed(123)
adonis(illumina_damaged_forvsopen_raref_OTU ~ Manejo_Conteo, illumina_damaged_forvsopen_raref_data,permutations = 999)


# Beta diversidad para cada manejo por salud ----

# Forest ----
illumina_damaged_for <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Forest")
illumina_damaged_for_bs <- subset_samples(illumina_damaged_for, Type != "Bare soil")


# Healthy vs Affected
bacteria_forest_HvsA <- subset_samples(illumina_damaged_for_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_forest_HvsA))

set.seed(123)
bacteria_forest_HvsA_raref <- rarefy_even_depth(bacteria_forest_HvsA)

set.seed(123)
bacteria_forest_HvsA_raref_OTU <- phyloseq::distance(bacteria_forest_HvsA_raref, method = "bray")
mean(bacteria_forest_HvsA_raref_OTU)

bacteria_forest_HvsA_raref_data <- data.frame(sample_data(bacteria_forest_HvsA_raref))

# Andonis with strata = Site

set.seed(123)
adonis(bacteria_forest_HvsA_raref_OTU ~ Type, strata = bacteria_forest_HvsA_raref_data$Site, data = bacteria_forest_HvsA_raref_data,  permutations = 999)


# Healthy vs Dead
bacteria_forest_HvsD <- subset_samples(illumina_damaged_for_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_forest_HvsD))

set.seed(123)
bacteria_forest_HvsD_raref <- rarefy_even_depth(bacteria_forest_HvsD)

set.seed(123)
bacteria_forest_HvsD_raref_OTU <- phyloseq::distance(bacteria_forest_HvsD_raref, method = "bray")
mean(bacteria_forest_HvsD_raref_OTU)
bacteria_forest_HvsD_raref_data <- data.frame(sample_data(bacteria_forest_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_forest_HvsD_raref_OTU ~ Type, strata = bacteria_forest_HvsD_raref_data$Site, data = bacteria_forest_HvsD_raref_data,  permutations = 999)


# Affected vs Dead
bacteria_forest_AvsD <- subset_samples(illumina_damaged_for_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_forest_AvsD))

set.seed(123)
bacteria_forest_AvsD_raref <- rarefy_even_depth(bacteria_forest_AvsD)

set.seed(123)
bacteria_forest_AvsD_raref_OTU <- phyloseq::distance(bacteria_forest_AvsD_raref, method = "bray")
mean(bacteria_forest_AvsD_raref_OTU)
bacteria_forest_AvsD_raref_data <- data.frame(sample_data(bacteria_forest_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_forest_AvsD_raref_OTU ~ Type, strata = bacteria_forest_AvsD_raref_data$Site, data = bacteria_forest_AvsD_raref_data,  permutations = 999)


#### Healthy vs Bare soil 
bacteria_forest_HvsBS <- subset_samples(illumina_damaged_for, Type != "Affected")
bacteria_forest_HvsBS <- subset_samples(bacteria_forest_HvsBS, Type != "Dead")

sort(sample_sums(bacteria_forest_HvsBS))

set.seed(123)
bacteria_forest_HvsBS_raref <- rarefy_even_depth(bacteria_forest_HvsBS)

set.seed(123)
bacteria_forest_HvsBS_raref_OTU <- phyloseq::distance(bacteria_forest_HvsBS_raref, method = "bray")
mean(bacteria_forest_HvsBS_raref_OTU)
bacteria_forest_HvsBS_raref_data <- data.frame(sample_data(bacteria_forest_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_forest_HvsBS_raref_OTU ~ Type, strata = bacteria_forest_HvsBS_raref_data$Site, data = bacteria_forest_HvsBS_raref_data,  permutations = 999)


#### Affected vs Bare soil 
bacteria_forest_AvsBS <- subset_samples(illumina_damaged_for, Type != "Healthy")
bacteria_forest_AvsBS <- subset_samples(bacteria_forest_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_forest_AvsBS))

set.seed(123)
bacteria_forest_AvsBS_raref <- rarefy_even_depth(bacteria_forest_AvsBS)

set.seed(123)
bacteria_forest_AvsBS_raref_OTU <- phyloseq::distance(bacteria_forest_AvsBS_raref, method = "bray")
mean(bacteria_forest_AvsBS_raref_OTU)
bacteria_forest_AvsBS_raref_data <- data.frame(sample_data(bacteria_forest_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_forest_AvsBS_raref_OTU ~ Type, strata = bacteria_forest_AvsBS_raref_data$Site, data = bacteria_forest_AvsBS_raref_data,  permutations = 999)


#### Dead vs Bare soil 
bacteria_forest_DvsBS <- subset_samples(illumina_damaged_for, Type != "Healthy")
bacteria_forest_DvsBS <- subset_samples(bacteria_forest_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_forest_DvsBS))

set.seed(123)
bacteria_forest_DvsBS_raref <- rarefy_even_depth(bacteria_forest_DvsBS)

set.seed(123)
bacteria_forest_DvsBS_raref_OTU <- phyloseq::distance(bacteria_forest_DvsBS_raref, method = "bray")
mean(bacteria_forest_DvsBS_raref_OTU)
bacteria_forest_DvsBS_raref_data <- data.frame(sample_data(bacteria_forest_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_forest_DvsBS_raref_OTU ~ Type, strata = bacteria_forest_DvsBS_raref_data$Site, data = bacteria_forest_DvsBS_raref_data,  permutations = 999)


# Dehesa ----
illumina_damaged_deh <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Dehesa")
illumina_damaged_deh_bs <- subset_samples(illumina_damaged_deh, Type != "Bare soil")


# Healthy vs Affected
bacteria_Dehesa_HvsA <- subset_samples(illumina_damaged_deh_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_Dehesa_HvsA))

set.seed(123)
bacteria_Dehesa_HvsA_raref <- rarefy_even_depth(bacteria_Dehesa_HvsA)

set.seed(123)
bacteria_Dehesa_HvsA_raref_OTU <- phyloseq::distance(bacteria_Dehesa_HvsA_raref, method = "bray")
mean(bacteria_Dehesa_HvsA_raref_OTU)
bacteria_Dehesa_HvsA_raref_data <- data.frame(sample_data(bacteria_Dehesa_HvsA_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_HvsA_raref_OTU ~ Type, strata = bacteria_Dehesa_HvsA_raref_data$Site, data = bacteria_Dehesa_HvsA_raref_data,  permutations = 999)


# Healthy vs Dead
bacteria_Dehesa_HvsD <- subset_samples(illumina_damaged_deh_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_Dehesa_HvsD))

set.seed(123)
bacteria_Dehesa_HvsD_raref <- rarefy_even_depth(bacteria_Dehesa_HvsD)

set.seed(123)
bacteria_Dehesa_HvsD_raref_OTU <- phyloseq::distance(bacteria_Dehesa_HvsD_raref, method = "bray")
mean(bacteria_Dehesa_HvsD_raref_OTU)
bacteria_Dehesa_HvsD_raref_data <- data.frame(sample_data(bacteria_Dehesa_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_HvsD_raref_OTU ~ Type, strata = bacteria_Dehesa_HvsD_raref_data$Site, data = bacteria_Dehesa_HvsD_raref_data,  permutations = 999)


# Affected vs Dead
bacteria_Dehesa_AvsD <- subset_samples(illumina_damaged_deh_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_Dehesa_AvsD))

set.seed(123)
bacteria_Dehesa_AvsD_raref <- rarefy_even_depth(bacteria_Dehesa_AvsD)

set.seed(123)
bacteria_Dehesa_AvsD_raref_OTU <- phyloseq::distance(bacteria_Dehesa_AvsD_raref, method = "bray")
mean(bacteria_Dehesa_AvsD_raref_OTU)
bacteria_Dehesa_AvsD_raref_data <- data.frame(sample_data(bacteria_Dehesa_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_AvsD_raref_OTU ~ Type, strata = bacteria_Dehesa_AvsD_raref_data$Site, data = bacteria_Dehesa_AvsD_raref_data,  permutations = 999)


#### Healthy vs Bare soil 
bacteria_Dehesa_HvsBS <- subset_samples(illumina_damaged_deh, Type != "Affected")
bacteria_Dehesa_HvsBS <- subset_samples(bacteria_Dehesa_HvsBS, Type != "Dead")

sort(sample_sums(bacteria_Dehesa_HvsBS))
set.seed(123)
bacteria_Dehesa_HvsBS_raref <- rarefy_even_depth(bacteria_Dehesa_HvsBS)

set.seed(123)
bacteria_Dehesa_HvsBS_raref_OTU <- phyloseq::distance(bacteria_Dehesa_HvsBS_raref, method = "bray")
mean(bacteria_Dehesa_HvsBS_raref_OTU)
bacteria_Dehesa_HvsBS_raref_data <- data.frame(sample_data(bacteria_Dehesa_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_HvsBS_raref_OTU ~ Type, strata = bacteria_Dehesa_HvsBS_raref_data$Site, data = bacteria_Dehesa_HvsBS_raref_data,  permutations = 999)


#### Affected vs Bare soil 
bacteria_Dehesa_AvsBS <- subset_samples(illumina_damaged_deh, Type != "Healthy")
bacteria_Dehesa_AvsBS <- subset_samples(bacteria_Dehesa_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_Dehesa_AvsBS))

set.seed(123)
bacteria_Dehesa_AvsBS_raref <- rarefy_even_depth(bacteria_Dehesa_AvsBS)

set.seed(123)
bacteria_Dehesa_AvsBS_raref_OTU <- phyloseq::distance(bacteria_Dehesa_AvsBS_raref, method = "bray")
mean(bacteria_Dehesa_AvsBS_raref_OTU)
bacteria_Dehesa_AvsBS_raref_data <- data.frame(sample_data(bacteria_Dehesa_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_AvsBS_raref_OTU ~ Type, strata = bacteria_Dehesa_AvsBS_raref_data$Site, data = bacteria_Dehesa_AvsBS_raref_data,  permutations = 999)


#### Dead vs Bare soil 
bacteria_Dehesa_DvsBS <- subset_samples(illumina_damaged_deh, Type != "Healthy")
bacteria_Dehesa_DvsBS <- subset_samples(bacteria_Dehesa_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_Dehesa_DvsBS))

set.seed(123)
bacteria_Dehesa_DvsBS_raref <- rarefy_even_depth(bacteria_Dehesa_DvsBS)

set.seed(123)
bacteria_Dehesa_DvsBS_raref_OTU <- phyloseq::distance(bacteria_Dehesa_DvsBS_raref, method = "bray")
mean(bacteria_Dehesa_DvsBS_raref_OTU)
bacteria_Dehesa_DvsBS_raref_data <- data.frame(sample_data(bacteria_Dehesa_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_Dehesa_DvsBS_raref_OTU ~ Type, strata = bacteria_Dehesa_DvsBS_raref_data$Site, data = bacteria_Dehesa_DvsBS_raref_data,  permutations = 999)


# Open woodland ----
illumina_damaged_opw <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Open woodland")
illumina_damaged_opw_bs <- subset_samples(illumina_damaged_opw, Type != "Bare soil")


# Healthy vs Affected
bacteria_openwood_HvsA <- subset_samples(illumina_damaged_opw_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_openwood_HvsA))

set.seed(123)
bacteria_openwood_HvsA_raref <- rarefy_even_depth(bacteria_openwood_HvsA)

set.seed(123)
bacteria_openwood_HvsA_raref_OTU <- phyloseq::distance(bacteria_openwood_HvsA_raref, method = "bray")
mean(bacteria_openwood_HvsA_raref_OTU)
bacteria_openwood_HvsA_raref_data <- data.frame(sample_data(bacteria_openwood_HvsA_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_HvsA_raref_OTU ~ Type, strata = bacteria_openwood_HvsA_raref_data$Site, data = bacteria_openwood_HvsA_raref_data,  permutations = 999)


# Healthy vs Dead
bacteria_openwood_HvsD <- subset_samples(illumina_damaged_opw_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_openwood_HvsD))

set.seed(123)
bacteria_openwood_HvsD_raref <- rarefy_even_depth(bacteria_openwood_HvsD)

set.seed(123)
bacteria_openwood_HvsD_raref_OTU <- phyloseq::distance(bacteria_openwood_HvsD_raref, method = "bray")
mean(bacteria_openwood_HvsD_raref_OTU)
bacteria_openwood_HvsD_raref_data <- data.frame(sample_data(bacteria_openwood_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_HvsD_raref_OTU ~ Type, strata = bacteria_openwood_HvsD_raref_data$Site, data = bacteria_openwood_HvsD_raref_data,  permutations = 999)


# Affected vs Dead
bacteria_openwood_AvsD <- subset_samples(illumina_damaged_opw_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_openwood_AvsD))

set.seed(123)
bacteria_openwood_AvsD_raref <- rarefy_even_depth(bacteria_openwood_AvsD)

set.seed(123)
bacteria_openwood_AvsD_raref_OTU <- phyloseq::distance(bacteria_openwood_AvsD_raref, method = "bray")
mean(bacteria_openwood_AvsD_raref_OTU)
bacteria_openwood_AvsD_raref_data <- data.frame(sample_data(bacteria_openwood_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_AvsD_raref_OTU ~ Type, strata = bacteria_openwood_AvsD_raref_data$Site, data = bacteria_openwood_AvsD_raref_data,  permutations = 999)


#### Healthy vs Bare soil 
bacteria_openwood_HvsBS <- subset_samples(illumina_damaged_opw, Type != "Affected")
bacteria_openwood_HvsBS <- subset_samples(bacteria_openwood_HvsBS, Type != "Dead")

sort(sample_sums(bacteria_openwood_HvsBS))

set.seed(123)
bacteria_openwood_HvsBS_raref <- rarefy_even_depth(bacteria_openwood_HvsBS)

set.seed(123)
bacteria_openwood_HvsBS_raref_OTU <- phyloseq::distance(bacteria_openwood_HvsBS_raref, method = "bray")
mean(bacteria_openwood_HvsBS_raref_OTU)
bacteria_openwood_HvsBS_raref_data <- data.frame(sample_data(bacteria_openwood_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_HvsBS_raref_OTU ~ Type, strata = bacteria_openwood_HvsBS_raref_data$Site, data = bacteria_openwood_HvsBS_raref_data,  permutations = 999)


#### Affected vs Bare soil 
bacteria_openwood_AvsBS <- subset_samples(illumina_damaged_opw, Type != "Healthy")
bacteria_openwood_AvsBS <- subset_samples(bacteria_openwood_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_openwood_AvsBS))

set.seed(123)
bacteria_openwood_AvsBS_raref <- rarefy_even_depth(bacteria_openwood_AvsBS)

set.seed(123)
bacteria_openwood_AvsBS_raref_OTU <- phyloseq::distance(bacteria_openwood_AvsBS_raref, method = "bray")
mean(bacteria_openwood_AvsBS_raref_OTU)
bacteria_openwood_AvsBS_raref_data <- data.frame(sample_data(bacteria_openwood_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_AvsBS_raref_OTU ~ Type, strata = bacteria_openwood_AvsBS_raref_data$Site, data = bacteria_openwood_AvsBS_raref_data,  permutations = 999)


#### Dead vs Bare soil 
bacteria_openwood_DvsBS <- subset_samples(illumina_damaged_opw, Type != "Healthy")
bacteria_openwood_DvsBS <- subset_samples(bacteria_openwood_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(bacteria_openwood_DvsBS))

set.seed(123)
bacteria_openwood_DvsBS_raref <- rarefy_even_depth(bacteria_openwood_DvsBS)

set.seed(123)
bacteria_openwood_DvsBS_raref_OTU <- phyloseq::distance(bacteria_openwood_DvsBS_raref, method = "bray")
mean(bacteria_openwood_DvsBS_raref_OTU)
bacteria_openwood_DvsBS_raref_data <- data.frame(sample_data(bacteria_openwood_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(bacteria_openwood_DvsBS_raref_OTU ~ Type, strata = bacteria_openwood_DvsBS_raref_data$Site, data = bacteria_openwood_DvsBS_raref_data,  permutations = 999)

