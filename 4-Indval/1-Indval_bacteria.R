##%######################################################%##
#                                                          #
####                Specialist species                  ####
#                                                          #
##%######################################################%##
library(labdsv)
library(phyloseq)

# Select the Management to study
illumina_16s_good_damaged

# Make rarefaction to compare between samples
set.seed(123)
illumina_damaged_man_raref <- rarefy_even_depth(illumina_16s_good_damaged)

# Change taxa names
taxa_names(illumina_damaged_man_raref) <- paste("OTU", seq(1,length(rownames(tax_table(illumina_damaged_man_raref)))), sep="_")
illumina_damaged_man_otu_tab <- veganotu(illumina_damaged_man_raref)
illumina_damaged_man_data <- data.frame(sample_data(illumina_damaged_man_raref))
illumina_damaged_man_tax <- data.frame(tax_table(illumina_damaged_man_raref)[,1:6])

# Calculate indval (specialist species) according to tree health (https://www.rdocumentation.org/packages/labdsv/versions/1.8-0/topics/indval)
set.seed(123)
indval_man <- indval(illumina_damaged_man_otu_tab,illumina_damaged_man_data$Manejo_Conteo)

summary(indval_man, p=0.05)
nrow(illumina_damaged_man_tax)
# Select significant values
sum(table(which(indval_man$pval <= 0.05)))
indval.sig_man <- names(which(indval_man$pval <= 0.05))

# Select indval higher than 0.3 for each health tree so we can know how many specialist we have (https://www.nature.com/articles/ismej2012168#abstract)
indval.df_man <- indval_man$indval

indval.df_man$OTU_ID <- ifelse(match(rownames(indval.df_man), indval.sig_man), indval.sig_man, "Undefined")
indval.df_man <- na.omit(indval.df_man)
nrow(indval.df_man)

dehesa.inval_man <- subset(indval.df_man, indval.df_man[ , 1] >= 0.3)
openwood.inval_man <- subset(indval.df_man, indval.df_man[ , 2] >= 0.3)
forest.inval_man <- subset(indval.df_man, indval.df_man[ , 3] >= 0.3)


# Recount
(nrow(dehesa.inval_man)/nrow(illumina_damaged_man_tax))*100
(nrow(openwood.inval_man)/nrow(illumina_damaged_man_tax))*100
(nrow(forest.inval_man)/nrow(illumina_damaged_man_tax))*100

# Extract proportion of indval by land-use
dehesa.inval_man$Phylum <- as.vector(illumina_damaged_man_tax[match(rownames(dehesa.inval_man), rownames(illumina_damaged_man_tax)), 2])
openwood.inval_man$Phylum <- as.vector(illumina_damaged_man_tax[match(rownames(openwood.inval_man), rownames(illumina_damaged_man_tax)), 2])
forest.inval_man$Phylum <- as.vector(illumina_damaged_man_tax[match(rownames(forest.inval_man), rownames(illumina_damaged_man_tax)), 2])

dehesa.inval_man
openwood.inval_man
forest.inval_man

dehesa_inval_count = aggregate(OTU_ID ~ Phylum, dehesa.inval_man, length)
dehesa_inval_count$OTU_ID = (dehesa_inval_count$OTU_ID/nrow(dehesa.inval_man))*100

openwood_inval_count = aggregate(OTU_ID ~ Phylum, openwood.inval_man, length)
openwood_inval_count$OTU_ID = (openwood_inval_count$OTU_ID/nrow(openwood.inval_man))*100

forest_inval_count = aggregate(OTU_ID ~ Phylum, forest.inval_man, length)
forest_inval_count$OTU_ID = (forest_inval_count$OTU_ID/nrow(forest.inval_man))*100

# Calculate indval to forests ----
# Select the Management to study
illumina_damaged_for <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Forest")

# Make rarefaction to compare between samples
set.seed(123)
illumina_damaged_for_raref <- rarefy_even_depth(illumina_damaged_for)

# Change taxa names
taxa_names(illumina_damaged_for_raref) <- paste("OTU", seq(1,length(rownames(tax_table(illumina_damaged_for_raref)))), sep="_")
illumina_damaged_for_otu_tab <- veganotu(illumina_damaged_for_raref)
illumina_damaged_for_data <- data.frame(sample_data(illumina_damaged_for_raref))
illumina_damaged_for_tax <- data.frame(tax_table(illumina_damaged_for_raref)[,1:6])

# Calculate indval (specialist species) according to tree health (https://www.rdocumentation.org/packages/labdsv/versions/1.8-0/topics/indval)
set.seed(123)
indval_for <- indval(illumina_damaged_for_otu_tab,illumina_damaged_for_data$Type)

summary(indval_for, p=0.05)
nrow(illumina_damaged_for_tax)
# Select significant values
sum(table(which(indval_for$pval <= 0.05)))
indval.sig_for <- names(which(indval_for$pval <= 0.05))

# Select indval higher than 0.3 for each health tree so we can know how many specialist we have (https://www.nature.com/articles/ismej2012168#abstract)
indval.df_for <- indval_for$indval

indval.df_for$OTU_ID <- ifelse(match(rownames(indval.df_for), indval.sig_for), indval.sig_for, "Undefined")
indval.df_for <- na.omit(indval.df_for)
nrow(indval.df_for)

healthy.inval_for <- subset(indval.df_for, indval.df_for[ , 1] >= 0.3)
affected.inval_for <- subset(indval.df_for, indval.df_for[ , 2] >= 0.3)
dead.inval_for <- subset(indval.df_for, indval.df_for[ , 3] >= 0.3)
baresoil.inval_for <- subset(indval.df_for, indval.df_for[ , 4] >= 0.3)

# Recount
(nrow(healthy.inval_for)/nrow(illumina_damaged_for_tax))*100
(nrow(affected.inval_for)/nrow(illumina_damaged_for_tax))*100
(nrow(dead.inval_for)/nrow(illumina_damaged_for_tax))*100
(nrow(baresoil.inval_for)/nrow(illumina_damaged_for_tax))*100

# Extract proportion of indval by tree health
healthy.inval_for$Phylum <- as.vector(illumina_damaged_for_tax[match(rownames(healthy.inval_for), rownames(illumina_damaged_for_tax)), 2])
affected.inval_for$Phylum <- as.vector(illumina_damaged_for_tax[match(rownames(affected.inval_for), rownames(illumina_damaged_for_tax)), 2])
dead.inval_for$Phylum <- as.vector(illumina_damaged_for_tax[match(rownames(dead.inval_for), rownames(illumina_damaged_for_tax)), 2])
baresoil.inval_for$Phylum <- as.vector(illumina_damaged_for_tax[match(rownames(baresoil.inval_for), rownames(illumina_damaged_for_tax)), 2])

healthy.inval_for
affected.inval_for
dead.inval_for
baresoil.inval_for


# Calculate indval to dehesas ----
# Select the Management to study
illumina_damaged_deh <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Dehesa")

# Make rarefaction to compare between samples
set.seed(123)
illumina_damaged_deh_raref <- rarefy_even_depth(illumina_damaged_deh)

# Change taxa names
taxa_names(illumina_damaged_deh_raref) <- paste("OTU", seq(1,length(rownames(tax_table(illumina_damaged_deh_raref)))), sep="_")
illumina_damaged_deh_otu_tab <- veganotu(illumina_damaged_deh_raref)
illumina_damaged_deh_data <- data.frame(sample_data(illumina_damaged_deh_raref))
illumina_damaged_deh_tax <- data.frame(tax_table(illumina_damaged_deh_raref)[,1:6])

# Calculate indval (specialist species) according to tree health (https://www.rdocumentation.org/packages/labdsv/versions/1.8-0/topics/indval)
set.seed(123)
indval_deh <- indval(illumina_damaged_deh_otu_tab,illumina_damaged_deh_data$Type)

summary(indval_deh, p=0.05)
nrow(illumina_damaged_deh_tax)
# Select significant values
sum(table(which(indval_deh$pval <= 0.05)))
indval.sig_deh <- names(which(indval_deh$pval <= 0.05))

# Select indval higher than 0.3 for each health tree so we can know how many specialist we have (https://www.nature.com/articles/ismej2012168#abstract)
indval.df_deh <- indval_deh$indval

indval.df_deh$OTU_ID <- ifelse(match(rownames(indval.df_deh), indval.sig_deh), indval.sig_deh, "Undefined")
indval.df_deh <- na.omit(indval.df_deh)
nrow(indval.df_deh)

healthy.inval_deh <- subset(indval.df_deh, indval.df_deh[ , 1] >= 0.3)
affected.inval_deh <- subset(indval.df_deh, indval.df_deh[ , 2] >= 0.3)
dead.inval_deh <- subset(indval.df_deh, indval.df_deh[ , 3] >= 0.3)
baresoil.inval_deh <- subset(indval.df_deh, indval.df_deh[ , 4] >= 0.3)

# Recount
(nrow(healthy.inval_deh)/nrow(illumina_damaged_deh_tax))*100
(nrow(affected.inval_deh)/nrow(illumina_damaged_deh_tax))*100
(nrow(dead.inval_deh)/nrow(illumina_damaged_deh_tax))*100
(nrow(baresoil.inval_deh)/nrow(illumina_damaged_deh_tax))*100

# Extract proportion of indval by tree health
healthy.inval_deh$Phylum <- as.vector(illumina_damaged_deh_tax[match(rownames(healthy.inval_deh), rownames(illumina_damaged_deh_tax)), 2])
affected.inval_deh$Phylum <- as.vector(illumina_damaged_deh_tax[match(rownames(affected.inval_deh), rownames(illumina_damaged_deh_tax)), 2])
dead.inval_deh$Phylum <- as.vector(illumina_damaged_deh_tax[match(rownames(dead.inval_deh), rownames(illumina_damaged_deh_tax)), 2])
baresoil.inval_deh$Phylum <- as.vector(illumina_damaged_deh_tax[match(rownames(baresoil.inval_deh), rownames(illumina_damaged_deh_tax)), 2])

healthy.inval_deh
affected.inval_deh
dead.inval_deh
baresoil.inval_deh


# Calculate indval to Open woodlands----
# Select the Management to study
illumina_damaged_opw <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Open woodland")

# Make rarefaction to compare between samples
set.seed(123)
illumina_damaged_opw_raref <- rarefy_even_depth(illumina_damaged_opw)

# Change taxa names
taxa_names(illumina_damaged_opw_raref) <- paste("OTU", seq(1,length(rownames(tax_table(illumina_damaged_opw_raref)))), sep="_")
illumina_damaged_opw_otu_tab <- veganotu(illumina_damaged_opw_raref)
illumina_damaged_opw_data <- data.frame(sample_data(illumina_damaged_opw_raref))
illumina_damaged_opw_tax <- data.frame(tax_table(illumina_damaged_opw_raref)[,1:6])

# Calculate indval (specialist species) according to tree health (https://www.rdocumentation.org/packages/labdsv/versions/1.8-0/topics/indval)
set.seed(123)
indval_opw <- indval(illumina_damaged_opw_otu_tab,illumina_damaged_opw_data$Type)

summary(indval_opw, p=0.05)
nrow(illumina_damaged_opw_tax)
# Select significant values
sum(table(which(indval_opw$pval <= 0.05)))
indval.sig_opw <- names(which(indval_opw$pval <= 0.05))

# Select indval higher than 0.3 for each health tree so we can know how many specialist we have (https://www.nature.com/articles/ismej2012168#abstract)
indval.df_opw <- indval_opw$indval

indval.df_opw$OTU_ID <- ifelse(match(rownames(indval.df_opw), indval.sig_opw), indval.sig_opw, "Undefined")
indval.df_opw <- na.omit(indval.df_opw)
nrow(indval.df_opw)

healthy.inval_opw <- subset(indval.df_opw, indval.df_opw[ , 1] >= 0.3)
affected.inval_opw <- subset(indval.df_opw, indval.df_opw[ , 2] >= 0.3)
dead.inval_opw <- subset(indval.df_opw, indval.df_opw[ , 3] >= 0.3)
baresoil.inval_opw <- subset(indval.df_opw, indval.df_opw[ , 4] >= 0.3)

# Recount
(nrow(healthy.inval_opw)/nrow(illumina_damaged_opw_tax))*100
(nrow(affected.inval_opw)/nrow(illumina_damaged_opw_tax))*100
(nrow(dead.inval_opw)/nrow(illumina_damaged_opw_tax))*100
(nrow(baresoil.inval_opw)/nrow(illumina_damaged_opw_tax))*100

# Extract proportion of indval by tree health
healthy.inval_opw$Phylum <- as.vector(illumina_damaged_opw_tax[match(rownames(healthy.inval_opw), rownames(illumina_damaged_opw_tax)), 2])
affected.inval_opw$Phylum <- as.vector(illumina_damaged_opw_tax[match(rownames(affected.inval_opw), rownames(illumina_damaged_opw_tax)), 2])
dead.inval_opw$Phylum <- as.vector(illumina_damaged_opw_tax[match(rownames(dead.inval_opw), rownames(illumina_damaged_opw_tax)), 2])
baresoil.inval_opw$Phylum <- as.vector(illumina_damaged_opw_tax[match(rownames(baresoil.inval_opw), rownames(illumina_damaged_opw_tax)), 2])

healthy.inval_opw
affected.inval_opw
dead.inval_opw
baresoil.inval_opw

