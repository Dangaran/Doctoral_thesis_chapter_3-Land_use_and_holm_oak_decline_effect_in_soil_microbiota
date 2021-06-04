##%######################################################%##
#                                                          #
####                  Beta diversity                    ####
#                                                          #
##%######################################################%##
# For each land-use
library(vegan)
fungi_good_damaged

fungi_damaged_raref <- rarefy_even_depth(fungi_good_damaged)

ordu <-  ordinate(fungi_damaged_raref, "PCoA", "bray", weighted=TRUE)
plot_ordination(fungi_damaged_raref, ordu, color="Manejo_Conteo")  + 
  stat_ellipse(type = "t", linetype = 2)

# Dehesa vs Forest
fungi_damaged_DehvsFor <- subset_samples(fungi_good_damaged, Manejo_Conteo != "Open woodland")

sort(sample_sums(fungi_damaged_DehvsFor))

fungi_damaged_DehvsFor_raref <- rarefy_even_depth(fungi_damaged_DehvsFor)

fungi_damaged_DehvsFor_raref_OTU <- veganotu(fungi_damaged_DehvsFor_raref)
fungi_damaged_DehvsFor_raref_data <- data.frame(sample_data(fungi_damaged_DehvsFor_raref))


# Andonis 
set.seed(123)
adonis(fungi_damaged_DehvsFor_raref_OTU ~ Manejo_Conteo, fungi_damaged_DehvsFor_raref_data,permutations = 999)

# Dehesa vs Open woodland
fungi_damaged_Dehvsopen <- subset_samples(fungi_good_damaged, Manejo_Conteo != "Forest")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_damaged_Dehvsopen))

fungi_damaged_Dehvsopen_raref <- rarefy_even_depth(fungi_damaged_Dehvsopen)

fungi_damaged_Dehvsopen_raref_OTU <- veganotu(fungi_damaged_Dehvsopen_raref)
fungi_damaged_Dehvsopen_raref_data <- data.frame(sample_data(fungi_damaged_Dehvsopen_raref))


# Andonis 
set.seed(123)
adonis(fungi_damaged_Dehvsopen_raref_OTU ~ Manejo_Conteo, fungi_damaged_Dehvsopen_raref_data,permutations = 999)


# Forest vs Open woodland
fungi_damaged_forvsopen <- subset_samples(fungi_good_damaged, Manejo_Conteo != "Dehesa")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_damaged_forvsopen))

fungi_damaged_forvsopen_raref <- rarefy_even_depth(fungi_damaged_forvsopen)

fungi_damaged_forvsopen_raref_OTU <- veganotu(fungi_damaged_forvsopen_raref)
fungi_damaged_forvsopen_raref_data <- data.frame(sample_data(fungi_damaged_forvsopen_raref))


# Andonis 
set.seed(123)
adonis(fungi_damaged_forvsopen_raref_OTU ~ Manejo_Conteo, fungi_damaged_forvsopen_raref_data,permutations = 999)


# Beta diversity by land-use and tree health ----

# Forest ----
fungi_damaged_for <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Forest")
fungi_damaged_for_bs <- subset_samples(fungi_damaged_for, Type != "Bare soil")


# Healthy vs Affected
fungi_forest_HvsA <- subset_samples(fungi_damaged_for_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_forest_HvsA))

set.seed(123)
fungi_forest_HvsA_raref <- rarefy_even_depth(fungi_forest_HvsA)

set.seed(123)
fungi_forest_HvsA_raref_OTU <- phyloseq::distance(fungi_forest_HvsA_raref, method = "bray")
mean(fungi_forest_HvsA_raref_OTU)
fungi_forest_HvsA_raref_data <- data.frame(sample_data(fungi_forest_HvsA_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_HvsA_raref_OTU ~ Type, strata = fungi_forest_HvsA_raref_data$Site, data = fungi_forest_HvsA_raref_data,  permutations = 999, method = "bray")


# Healthy vs Dead
fungi_forest_HvsD <- subset_samples(fungi_damaged_for_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_forest_HvsD))

set.seed(123)
fungi_forest_HvsD_raref <- rarefy_even_depth(fungi_forest_HvsD)

set.seed(123)
fungi_forest_HvsD_raref_OTU <- phyloseq::distance(fungi_forest_HvsD_raref, method = "bray")
mean(fungi_forest_HvsD_raref_OTU)
fungi_forest_HvsD_raref_data <- data.frame(sample_data(fungi_forest_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_HvsD_raref_OTU ~ Type, strata = fungi_forest_HvsD_raref_data$Site, data = fungi_forest_HvsD_raref_data,  permutations = 999, method = "bray")


# Affected vs Dead
fungi_forest_AvsD <- subset_samples(fungi_damaged_for_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_forest_AvsD))

set.seed(123)
fungi_forest_AvsD_raref <- rarefy_even_depth(fungi_forest_AvsD)

set.seed(123)
fungi_forest_AvsD_raref_OTU <- phyloseq::distance(fungi_forest_AvsD_raref, method = "bray")
mean(fungi_forest_AvsD_raref_OTU)
fungi_forest_AvsD_raref_data <- data.frame(sample_data(fungi_forest_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_AvsD_raref_OTU ~ Type, strata = fungi_forest_AvsD_raref_data$Site, data = fungi_forest_AvsD_raref_data,  permutations = 999, method = "bray")


#### Healthy vs Bare soil 
fungi_forest_HvsBS <- subset_samples(fungi_damaged_for, Type != "Affected")
fungi_forest_HvsBS <- subset_samples(fungi_forest_HvsBS, Type != "Dead")

sort(sample_sums(fungi_forest_HvsBS))

set.seed(123)
fungi_forest_HvsBS_raref <- rarefy_even_depth(fungi_forest_HvsBS)

set.seed(123)
fungi_forest_HvsBS_raref_OTU <- phyloseq::distance(fungi_forest_HvsBS_raref, method = "bray")
mean(fungi_forest_HvsBS_raref_OTU)
fungi_forest_HvsBS_raref_data <- data.frame(sample_data(fungi_forest_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_HvsBS_raref_OTU ~ Type, strata = fungi_forest_HvsBS_raref_data$Site, data = fungi_forest_HvsBS_raref_data,  permutations = 999, method = "bray")


#### Affected vs Bare soil 
fungi_forest_AvsBS <- subset_samples(fungi_damaged_for, Type != "Healthy")
fungi_forest_AvsBS <- subset_samples(fungi_forest_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_forest_AvsBS))

set.seed(123)
fungi_forest_AvsBS_raref <- rarefy_even_depth(fungi_forest_AvsBS)

set.seed(123)
fungi_forest_AvsBS_raref_OTU <- phyloseq::distance(fungi_forest_AvsBS_raref, method = "bray")
mean(fungi_forest_AvsBS_raref_OTU)
fungi_forest_AvsBS_raref_data <- data.frame(sample_data(fungi_forest_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_AvsBS_raref_OTU ~ Type, strata = fungi_forest_AvsBS_raref_data$Site, data = fungi_forest_AvsBS_raref_data,  permutations = 999, method = "bray")


#### Dead vs Bare soil 
fungi_forest_DvsBS <- subset_samples(fungi_damaged_for, Type != "Healthy")
fungi_forest_DvsBS <- subset_samples(fungi_forest_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_forest_DvsBS))

set.seed(123)
fungi_forest_DvsBS_raref <- rarefy_even_depth(fungi_forest_DvsBS)

set.seed(123)
fungi_forest_DvsBS_raref_OTU <- phyloseq::distance(fungi_forest_DvsBS_raref, method = "bray")
mean(fungi_forest_DvsBS_raref_OTU)
fungi_forest_DvsBS_raref_data <- data.frame(sample_data(fungi_forest_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_forest_DvsBS_raref_OTU ~ Type, strata = fungi_forest_DvsBS_raref_data$Site, data = fungi_forest_DvsBS_raref_data,  permutations = 999, method = "bray")


# Dehesa ----
fungi_damaged_deh <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Dehesa")
fungi_damaged_deh_bs <- subset_samples(fungi_damaged_deh, Type != "Bare soil")


# Healthy vs Affected
fungi_Dehesa_HvsA <- subset_samples(fungi_damaged_deh_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_Dehesa_HvsA))

set.seed(123)
fungi_Dehesa_HvsA_raref <- rarefy_even_depth(fungi_Dehesa_HvsA)

set.seed(123)
fungi_Dehesa_HvsA_raref_OTU <- phyloseq::distance(fungi_Dehesa_HvsA_raref, method = "bray")
mean(fungi_Dehesa_HvsA_raref_OTU)
fungi_Dehesa_HvsA_raref_data <- data.frame(sample_data(fungi_Dehesa_HvsA_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_HvsA_raref_OTU ~ Type, strata = fungi_Dehesa_HvsA_raref_data$Site, data = fungi_Dehesa_HvsA_raref_data,  permutations = 999, method = "bray")


# Healthy vs Dead
fungi_Dehesa_HvsD <- subset_samples(fungi_damaged_deh_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_Dehesa_HvsD))

set.seed(123)
fungi_Dehesa_HvsD_raref <- rarefy_even_depth(fungi_Dehesa_HvsD)

set.seed(123)
fungi_Dehesa_HvsD_raref_OTU <- phyloseq::distance(fungi_Dehesa_HvsD_raref, method = "bray")
mean(fungi_Dehesa_HvsD_raref_OTU)
fungi_Dehesa_HvsD_raref_data <- data.frame(sample_data(fungi_Dehesa_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_HvsD_raref_OTU ~ Type, strata = fungi_Dehesa_HvsD_raref_data$Site, data = fungi_Dehesa_HvsD_raref_data,  permutations = 999, method = "bray")


# Affected vs Dead
fungi_Dehesa_AvsD <- subset_samples(fungi_damaged_deh_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_Dehesa_AvsD))

fungi_Dehesa_AvsD_raref <- rarefy_even_depth(fungi_Dehesa_AvsD)

set.seed(123)
fungi_Dehesa_AvsD_raref_OTU <- phyloseq::distance(fungi_Dehesa_AvsD_raref, method = "bray")
mean(fungi_Dehesa_AvsD_raref_OTU)
fungi_Dehesa_AvsD_raref_data <- data.frame(sample_data(fungi_Dehesa_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_AvsD_raref_OTU ~ Type, strata = fungi_Dehesa_AvsD_raref_data$Site, data = fungi_Dehesa_AvsD_raref_data,  permutations = 999, method = "bray")


#### Healthy vs Bare soil 
fungi_Dehesa_HvsBS <- subset_samples(fungi_damaged_deh, Type != "Affected")
fungi_Dehesa_HvsBS <- subset_samples(fungi_Dehesa_HvsBS, Type != "Dead")

sort(sample_sums(fungi_Dehesa_HvsBS))

set.seed(123)
fungi_Dehesa_HvsBS_raref <- rarefy_even_depth(fungi_Dehesa_HvsBS)

set.seed(123)
fungi_Dehesa_HvsBS_raref_OTU <- phyloseq::distance(fungi_Dehesa_HvsBS_raref, method = "bray")
mean(fungi_Dehesa_HvsBS_raref_OTU)
fungi_Dehesa_HvsBS_raref_data <- data.frame(sample_data(fungi_Dehesa_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_HvsBS_raref_OTU ~ Type, strata = fungi_Dehesa_HvsBS_raref_data$Site, data = fungi_Dehesa_HvsBS_raref_data,  permutations = 999, method = "bray")


#### Affected vs Bare soil 
fungi_Dehesa_AvsBS <- subset_samples(fungi_damaged_deh, Type != "Healthy")
fungi_Dehesa_AvsBS <- subset_samples(fungi_Dehesa_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_Dehesa_AvsBS))

set.seed(123)
fungi_Dehesa_AvsBS_raref <- rarefy_even_depth(fungi_Dehesa_AvsBS)

set.seed(123)
fungi_Dehesa_AvsBS_raref_OTU <- phyloseq::distance(fungi_Dehesa_AvsBS_raref, method = "bray")
mean(fungi_Dehesa_AvsBS_raref_OTU)
fungi_Dehesa_AvsBS_raref_data <- data.frame(sample_data(fungi_Dehesa_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_AvsBS_raref_OTU ~ Type, strata = fungi_Dehesa_AvsBS_raref_data$Site, data = fungi_Dehesa_AvsBS_raref_data,  permutations = 999, method = "bray")


#### Dead vs Bare soil 
fungi_Dehesa_DvsBS <- subset_samples(fungi_damaged_deh, Type != "Healthy")
fungi_Dehesa_DvsBS <- subset_samples(fungi_Dehesa_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_Dehesa_DvsBS))

set.seed(123)
fungi_Dehesa_DvsBS_raref <- rarefy_even_depth(fungi_Dehesa_DvsBS)

set.seed(123)
fungi_Dehesa_DvsBS_raref_OTU <- phyloseq::distance(fungi_Dehesa_DvsBS_raref, method = "bray")
mean(fungi_Dehesa_DvsBS_raref_OTU)
fungi_Dehesa_DvsBS_raref_data <- data.frame(sample_data(fungi_Dehesa_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_Dehesa_DvsBS_raref_OTU ~ Type, strata = fungi_Dehesa_DvsBS_raref_data$Site, data = fungi_Dehesa_DvsBS_raref_data,  permutations = 999, method = "bray")


# Open woodland ----
fungi_damaged_opw <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Open woodland")
fungi_damaged_opw_bs <- subset_samples(fungi_damaged_opw, Type != "Bare soil")


# Healthy vs Affected
fungi_openwood_HvsA <- subset_samples(fungi_damaged_opw_bs, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_openwood_HvsA))

set.seed(123)
fungi_openwood_HvsA_raref <- rarefy_even_depth(fungi_openwood_HvsA)

set.seed(123)
fungi_openwood_HvsA_raref_OTU <- phyloseq::distance(fungi_openwood_HvsA_raref, method = "bray")
mean(fungi_openwood_HvsA_raref_OTU)
fungi_openwood_HvsA_raref_data <- data.frame(sample_data(fungi_openwood_HvsA_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_HvsA_raref_OTU ~ Type, strata = fungi_openwood_HvsA_raref_data$Site, data = fungi_openwood_HvsA_raref_data,  permutations = 999, method = "bray")


# Healthy vs Dead
fungi_openwood_HvsD <- subset_samples(fungi_damaged_opw_bs, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_openwood_HvsD))

set.seed(123)
fungi_openwood_HvsD_raref <- rarefy_even_depth(fungi_openwood_HvsD)

set.seed(123)
fungi_openwood_HvsD_raref_OTU <- phyloseq::distance(fungi_openwood_HvsD_raref, method = "bray")
mean(fungi_openwood_HvsD_raref_OTU)
fungi_openwood_HvsD_raref_data <- data.frame(sample_data(fungi_openwood_HvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_HvsD_raref_OTU ~ Type, strata = fungi_openwood_HvsD_raref_data$Site, data = fungi_openwood_HvsD_raref_data,  permutations = 999, method = "bray")


# Affected vs Dead
fungi_openwood_AvsD <- subset_samples(fungi_damaged_opw_bs, Type != "Healthy")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_openwood_AvsD))

set.seed(123)
fungi_openwood_AvsD_raref <- rarefy_even_depth(fungi_openwood_AvsD)

set.seed(123)
fungi_openwood_AvsD_raref_OTU <- phyloseq::distance(fungi_openwood_AvsD_raref, method = "bray")
mean(fungi_openwood_AvsD_raref_OTU)
fungi_openwood_AvsD_raref_data <- data.frame(sample_data(fungi_openwood_AvsD_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_AvsD_raref_OTU ~ Type, strata = fungi_openwood_AvsD_raref_data$Site, data = fungi_openwood_AvsD_raref_data,  permutations = 999, method = "bray")


#### Healthy vs Bare soil 
fungi_openwood_HvsBS <- subset_samples(fungi_damaged_opw, Type != "Affected")
fungi_openwood_HvsBS <- subset_samples(fungi_openwood_HvsBS, Type != "Dead")

sort(sample_sums(fungi_openwood_HvsBS))

set.seed(123)
fungi_openwood_HvsBS_raref <- rarefy_even_depth(fungi_openwood_HvsBS)

set.seed(123)
fungi_openwood_HvsBS_raref_OTU <- phyloseq::distance(fungi_openwood_HvsBS_raref, method = "bray")
mean(fungi_openwood_HvsBS_raref_OTU)
fungi_openwood_HvsBS_raref_data <- data.frame(sample_data(fungi_openwood_HvsBS_raref))


# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_HvsBS_raref_OTU ~ Type, strata = fungi_openwood_HvsBS_raref_data$Site, data = fungi_openwood_HvsBS_raref_data,  permutations = 999, method = "bray")


#### Affected vs Bare soil 
fungi_openwood_AvsBS <- subset_samples(fungi_damaged_opw, Type != "Healthy")
fungi_openwood_AvsBS <- subset_samples(fungi_openwood_AvsBS, Type != "Dead")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_openwood_AvsBS))

set.seed(123)
fungi_openwood_AvsBS_raref <- rarefy_even_depth(fungi_openwood_AvsBS)

set.seed(123)
fungi_openwood_AvsBS_raref_OTU <- phyloseq::distance(fungi_openwood_AvsBS_raref, method = "bray")
mean(fungi_openwood_AvsBS_raref_OTU)
fungi_openwood_AvsBS_raref_data <- data.frame(sample_data(fungi_openwood_AvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_AvsBS_raref_OTU ~ Type, strata = fungi_openwood_AvsBS_raref_data$Site, data = fungi_openwood_AvsBS_raref_data,  permutations = 999, method = "bray")


#### Dead vs Bare soil 
fungi_openwood_DvsBS <- subset_samples(fungi_damaged_opw, Type != "Healthy")
fungi_openwood_DvsBS <- subset_samples(fungi_openwood_DvsBS, Type != "Affected")

# Remove taxa with less than 10 reads (Second filter)
sort(sample_sums(fungi_openwood_DvsBS))

set.seed(123)
fungi_openwood_DvsBS_raref <- rarefy_even_depth(fungi_openwood_DvsBS)

set.seed(123)
fungi_openwood_DvsBS_raref_OTU <- phyloseq::distance(fungi_openwood_DvsBS_raref, method = "bray")
mean(fungi_openwood_DvsBS_raref_OTU)
fungi_openwood_DvsBS_raref_data <- data.frame(sample_data(fungi_openwood_DvsBS_raref))

# Andonis with strata = Site
set.seed(123)
adonis(fungi_openwood_DvsBS_raref_OTU ~ Type, strata = fungi_openwood_DvsBS_raref_data$Site, data = fungi_openwood_DvsBS_raref_data,  permutations = 999, method = "bray")


