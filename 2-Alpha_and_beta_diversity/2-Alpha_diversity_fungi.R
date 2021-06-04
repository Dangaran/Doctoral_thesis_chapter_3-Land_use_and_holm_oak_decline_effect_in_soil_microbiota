##%######################################################%##
#                                                          #
####                  Alpha diversity                   ####
#                                                          #
##%######################################################%##
# By land-use 
fungi_good_damaged <- subset_samples(mycota_glom_filt, Treatment == "Damaged")

# Select between shannon, simpson and evenness
bact_alphadiv <- plot_anova_diversity(fungi_good_damaged, method = c("evenness"), 
                                      grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
aggregate(value ~ Manejo_Conteo,bact_alphadiv$data,mean)

# By tree status for each land-use
fungi_damaged_for <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Forest")
fungi_damaged_deh <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Dehesa")
fungi_damaged_opw <- subset_samples(fungi_good_damaged, Manejo_Conteo == "Open woodland")

# Forest
bact_dam_alphadiv_for <- plot_anova_diversity(fungi_damaged_for, method = c("evenness"), 
                                              grouping_column = "Type", pValueCutoff = 0.05)
aggregate(value ~ Type,bact_dam_alphadiv_for$data,mean)

# Dehesa
bact_dam_alphadiv_deh <- plot_anova_diversity(fungi_damaged_deh, method = c("evenness"), 
                                              grouping_column = "Type", pValueCutoff = 0.12)
aggregate(value ~ Type,bact_dam_alphadiv_deh$data,mean)

# Open woodland
bact_dam_alphadiv_opw <- plot_anova_diversity(fungi_damaged_opw, method = c("evenness"), 
                                              grouping_column = "Type", pValueCutoff = 0.05)
aggregate(value ~ Type,bact_dam_alphadiv_opw$data,mean)
