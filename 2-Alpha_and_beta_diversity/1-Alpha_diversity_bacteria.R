##%######################################################%##
#                                                          #
####                  Alpha diversity                   ####
#                                                          #
##%######################################################%##
# By land-use 
illumina_16s_good_damaged <- subset_samples(illumina_16s_10seq, Treatment == "Damaged")

# Select between shannon, simpson and evenness
bact_alphadiv <- plot_anova_diversity(illumina_16s_good_damaged, method = c("shannon"), 
                                      grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
aggregate(value ~ Manejo_Conteo,bact_alphadiv$data,mean)

# By tree status for each land-use
illumina_damaged_for <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Forest")
illumina_damaged_deh <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Dehesa")
illumina_damaged_opw <- subset_samples(illumina_16s_good_damaged, Manejo_Conteo == "Open woodland")

# Forest
bact_dam_alphadiv_for <- plot_anova_diversity(illumina_damaged_for, method = c("shannon"), 
                                              grouping_column = "Type", pValueCutoff = 0.05)
aggregate(value ~ Type,bact_dam_alphadiv_for$data,mean)

# Dehesa
bact_dam_alphadiv_deh <- plot_anova_diversity(illumina_damaged_deh, method = c("shannon"), 
                                              grouping_column = "Type", pValueCutoff = 0.05)
aggregate(value ~ Type,bact_dam_alphadiv_deh$data,mean)

# Open woodland
bact_dam_alphadiv_opw <- plot_anova_diversity(illumina_damaged_opw, method = c("shannon"), 
                                              grouping_column = "Type", pValueCutoff = 0.05)
aggregate(value ~ Type,bact_dam_alphadiv_opw$data,mean)
