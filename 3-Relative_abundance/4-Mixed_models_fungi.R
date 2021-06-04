##%######################################################%##
#                                                          #
####                Mixed model function                ####
#                                                          #
##%######################################################%##

mixed_model_phylum <- function(data_physeq, otu_table_melt_phyla, Phylum){
  phylum <- subset(otu_table_melt_phyla, Phyla %in% c(Phylum))
  data_physeq$Phylum <- phylum$value
  
  library(nlme)
  library(piecewiseSEM)
  mixed_model <- lme(Phylum ~ Porcentaje_defol, random = ~1|Site, data = data_physeq, na.action=na.exclude)
  summary_mixed <- summary(mixed_model)
  p1 <- plot(mixed_model) ## how the residuals bahave - heterocedasticidad
  norm <- shapiro.test(residuals(mixed_model))
  
  library(piecewiseSEM)
  r_sq <- rsquared(mixed_model)
  
  
  output <- list(data= data_physeq,
                 model = mixed_model,
                 summary = summary_mixed, 
                 cheked_res = p1,
                 checked_norm = norm,
                 r_square = r_sq)
  output
}


##%######################################################%##
#                                                          #
####                       Forest                       ####
#                                                          #
##%######################################################%##

# Prepare data
relab_fungi_for <-  transform_sample_counts(fungi_damaged_for, function(x) x / sum(x)*100)
relab_fungi_for = tax_glom(relab_fungi_for, "Phylum")
relab_fungi_for

relab_fungi_for_otu_tab <- data.frame(otu_table(relab_fungi_for))
rownames(relab_fungi_for_otu_tab) == rownames(tax_table(relab_fungi_for)) # Check order of rows
relab_fungi_for_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_for)[,2]) # Create Phyla column
relab_fungi_for_data <- data.frame(sample_data(relab_fungi_for))

library(reshape2)
relab_fungi_for_melt <- melt(relab_fungi_for_otu_tab)

# Input for mixex_model function
relab_fungi_for_data
relab_fungi_for_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_fungi_for_data, relab_fungi_for_melt, "p__Mortierellomycota")


##%######################################################%##
#                                                          #
####                       Dehesa                       ####
#                                                          #
##%######################################################%##

# Prepare data
relab_fungi_deh <-  transform_sample_counts(fungi_damaged_deh, function(x) x / sum(x)*100)
relab_fungi_deh = tax_glom(relab_fungi_deh, "Phylum")
relab_fungi_deh

relab_fungi_deh_otu_tab <- data.frame(otu_table(relab_fungi_deh))
rownames(relab_fungi_deh_otu_tab) == rownames(tax_table(relab_fungi_deh)) # Check order of rows
relab_fungi_deh_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_deh)[,2]) # Create Phyla column
relab_fungi_deh_data <- data.frame(sample_data(relab_fungi_deh))

library(reshape2)
relab_fungi_deh_melt <- melt(relab_fungi_deh_otu_tab)

# Input for mixex_model function
relab_fungi_deh_data
relab_fungi_deh_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_fungi_deh_data, relab_fungi_deh_melt, "p__Mucoromycota")


##%######################################################%##
#                                                          #
####                   Open woodlands                   ####
#                                                          #
##%######################################################%##

# Prepare data
relab_fungi_opw <-  transform_sample_counts(fungi_damaged_opw, function(x) x / sum(x)*100)
relab_fungi_opw = tax_glom(relab_fungi_opw, "Phylum") # Unir por filo
relab_fungi_opw

relab_fungi_opw_otu_tab <- data.frame(otu_table(relab_fungi_opw))
rownames(relab_fungi_opw_otu_tab) == rownames(tax_table(relab_fungi_opw)) # Check order of rows
relab_fungi_opw_otu_tab$Phyla <- as.vector(tax_table(relab_fungi_opw)[,2]) # Create Phyla column
relab_fungi_opw_data <- data.frame(sample_data(relab_fungi_opw))

library(reshape2)
relab_fungi_opw_melt <- melt(relab_fungi_opw_otu_tab)

# Input for mixex_model function
relab_fungi_opw_data
relab_fungi_opw_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_fungi_opw_data, relab_fungi_opw_melt, "p__Mucoromycota")


##%######################################################%##
#                                                          #
####                  Plots by phylum                   ####
#                                                          #
##%######################################################%##

# Ascomycota ----
p1_fung <- ggplot(test3$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 48.93582, slope = -0.00666, linetype = "Forest")) + # p__Ascomycota (NO SIG)
  geom_abline(aes(intercept = 36.99640, slope = 0.17332, linetype = "Dehesa (*)")) + # p__Ascomycota 
  geom_abline(aes(intercept = 31.54421, slope = 0.13342, linetype = "Open woodland (*)")) + # p__Ascomycota 
  ylim(0,60) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa (*)" = "dotted",
                                                      "Open woodland (*)" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Ascomycota")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))

# Basidiomycota ----
p2_fung <- ggplot(test3$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 45.58893, slope = -0.00840, linetype = "Forest")) + # p__Basidiomycota (NO SIG)
  geom_abline(aes(intercept = 59.81499, slope = -0.18976, linetype = "Dehesa (*)")) + # p__Basidiomycota 
  geom_abline(aes(intercept = 63.68352, slope = -0.13382, linetype = "Open woodland (0.077)")) + # p__Basidiomycota (NO SIG)
  ylim(0,65) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa (*)" = "dotted",
                                                      "Open woodland (0.077)" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Basidiomycota")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))


# Mortierellomycota ----
p3_fung <- ggplot(test3$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 2.5407048, slope = 0.0046405, linetype = "Forest")) + # p__Mortierellomycota (NO SIG)
  geom_abline(aes(intercept = 0, slope = 0, linetype = "Dehesa")) + # p__Basidiomycota 
  geom_abline(aes(intercept = 2.9696128, slope = -0.0040003, linetype = "Open woodland")) + # p__Mortierellomycota (NO SIG)
  ylim(0,10) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa" = "dotted",
                                                      "Open woodland" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Mortierellomycota")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))


# Mucoromycota ----
p4_fung <- ggplot(test3$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 2.8772777, slope = 0.0088222, linetype = "Forest")) + # p__Mucoromycota (NO SIG)
  geom_abline(aes(intercept = 1.178014, slope = 0.0148556, linetype = "Dehesa (*)")) + # p__Mucoromycota (NO SIG) TRANSFORMAR LOG
  geom_abline(aes(intercept = 1.7093859, slope = 0.0029219, linetype = "Open woodland")) + # p__Mucoromycota (NO SIG) TRANSFORMAR
  ylim(0,10) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa (*)" = "dotted",
                                                      "Open woodland" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Mucoromycota")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))

library(gridExtra)
grid.arrange(p1_fung, p2_fung, p3_fung, p4_fung, nrow = 2)






