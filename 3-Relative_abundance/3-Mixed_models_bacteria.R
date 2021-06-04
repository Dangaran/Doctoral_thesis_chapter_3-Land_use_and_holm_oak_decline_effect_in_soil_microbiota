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
relab_bact_for <-  transform_sample_counts(illumina_damaged_for, function(x) x / sum(x)*100)
relab_bact_for = tax_glom(relab_bact_for, "Phylum")
relab_bact_for

relab_bact_for_otu_tab <- data.frame(otu_table(relab_bact_for))
rownames(relab_bact_for_otu_tab) == rownames(tax_table(relab_bact_for)) # Check order of rows
relab_bact_for_otu_tab$Phyla <- as.vector(tax_table(relab_bact_for)[,2]) # Create Phyla column
relab_bact_for_data <- data.frame(sample_data(relab_bact_for))

library(reshape2)
relab_bact_for_melt <- melt(relab_bact_for_otu_tab)

# Input for mixex_model function
relab_bact_for_data
relab_bact_for_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_bact_for_data, relab_bact_for_melt, "Bacteroidetes")


##%######################################################%##
#                                                          #
####                       Dehesa                       ####
#                                                          #
##%######################################################%##

# Prepare data
relab_bact_deh <-  transform_sample_counts(illumina_damaged_deh, function(x) x / sum(x)*100)
relab_bact_deh = tax_glom(relab_bact_deh, "Phylum")
relab_bact_deh

relab_bact_deh_otu_tab <- data.frame(otu_table(relab_bact_deh))
rownames(relab_bact_deh_otu_tab) == rownames(tax_table(relab_bact_deh)) # Check order of rows
relab_bact_deh_otu_tab$Phyla <- as.vector(tax_table(relab_bact_deh)[,2]) # Create Phyla column
relab_bact_deh_data <- data.frame(sample_data(relab_bact_deh))

library(reshape2)
relab_bact_deh_melt <- melt(relab_bact_deh_otu_tab)

# Input for mixex_model function
relab_bact_deh_data
relab_bact_deh_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_bact_deh_data, relab_bact_deh_melt, "Acidobacteria")

##%######################################################%##
#                                                          #
####                   Open woodlands                   ####
#                                                          #
##%######################################################%##

# Prepare data
relab_bact_opw <-  transform_sample_counts(illumina_damaged_opw, function(x) x / sum(x)*100)
relab_bact_opw = tax_glom(relab_bact_opw, "Phylum")
relab_bact_opw

relab_bact_opw_otu_tab <- data.frame(otu_table(relab_bact_opw))
rownames(relab_bact_opw_otu_tab) == rownames(tax_table(relab_bact_opw)) # Check order of rows
relab_bact_opw_otu_tab$Phyla <- as.vector(tax_table(relab_bact_opw)[,2]) # Create Phyla column
relab_bact_opw_data <- data.frame(sample_data(relab_bact_opw))

library(reshape2)
relab_bact_opw_melt <- melt(relab_bact_opw_otu_tab)

# Input for mixex_model function
relab_bact_opw_data
relab_bact_opw_melt

# Check each phylum to obtain slope and intercept
data_with_model <- mixed_model_phylum(relab_bact_opw_data, relab_bact_opw_melt, "Verrucomicrobia")


##%######################################################%##
#                                                          #
####                  Plots by phylum                   ####
#                                                          #
##%######################################################%##

# Acidobacteria ----
library(ggplot2)
p1_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 20.534036, slope = 0.028934, linetype = "Forest (*)")) + # Forest
  geom_abline(aes(intercept = 21.461076, slope = 0.027055, linetype = "Dehesa")) + # Dehesa (NO SIG)
  geom_abline(aes(intercept = 22.653012, slope = 0.018122, linetype = "Open woodland")) + # Open woodland (NO SIG) TRANSFORMAR
  ylim(0,30) +
  scale_linetype_manual(name = "Land-use", values = c("Forest (*)" = "solid",
                                                  "Dehesa" = "dotted",
                                                  "Open woodland" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Acidobacteria")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(size = 0.4, 
        linetype = "solid"), panel.background = element_rect(fill = NA))

# Actinobacteria ----
p2_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 14.877278, slope = -0.042337, linetype = "Forest (*)")) + # Actinobacteria
  geom_abline(aes(intercept = 8.067926, slope = -0.013530, linetype = "Dehesa (*)")) + # Actinobacteria TRANSFORMAR LOG
  geom_abline(aes(intercept = 11.073021, slope = -0.026385, linetype = "Open woodland (*)")) + # Actinobacteria TRANSFORMAR LOG
  ylim(0,20) +
  scale_linetype_manual(name = "Land-use", values = c("Forest (*)" = "solid",
                                                      "Dehesa (*)" = "dotted",
                                                      "Open woodland (*)" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Actinobacteria")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))

# Bacteroidetes ----
p3_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 11.546233, slope = 0.008283, linetype = "Forest")) + # Bacteroidetes (NO SIG) TRANSFORMAR
  geom_abline(aes(intercept = 12.307456, slope = 0.007866, linetype = "Dehesa")) + # Bacteroidetes (NO SIG)
  geom_abline(aes(intercept = 10.983815, slope = 0.017509, linetype = "Open woodland")) + # Bacteroidetes (NO SIG)
  ylim(0,20) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa" = "dotted",
                                                      "Open woodland" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Bacteroidetes")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))


# Planctomycetes ----
p4_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 5.346657, slope = 0.006237, linetype = "Forest")) + # Planctomycetes (NO SIG) TRANSFORMAR LOG
  geom_abline(aes(intercept = 5.463937, slope = 0.006035, linetype = "Dehesa")) + # Planctomycetes (NO SIG)
  geom_abline(aes(intercept = 5.702437, slope = 0.007052, linetype = "Open woodland (*)")) + # Planctomycetes
  ylim(0,10) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa" = "dotted",
                                                      "Open woodland (*)" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Planctomycetes")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))

# Proteobacteria ----
p5_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 35.54835, slope = -0.02085, linetype = "Forest")) + # Proteobacteria (NO SIG)
  geom_abline(aes(intercept = 36.2064, slope = -0.0464, linetype = "Dehesa (*)")) + # Proteobacteria
  geom_abline(aes(intercept = 38.46473, slope = -0.03578, linetype = "Open woodland (*)")) + # Proteobacteria
  ylim(0,50) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa (*)" = "dotted",
                                                      "Open woodland (*)" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Proteobacteria")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))


# Verrucomicrobia ----
p6_bact <- ggplot(data_with_model$data, aes(Porcentaje_defol, Phylum)) +
  geom_abline(aes(intercept = 7.243340, slope = 0.016906, linetype = "Forest")) + # Verrucomicrobia (NO SIG)
  geom_abline(aes(intercept = 9.121607, slope = 0.002570, linetype = "Dehesa")) + # Verrucomicrobia (NO SIG)
  geom_abline(aes(intercept = 6.481797, slope = 0.014149, linetype = "Open woodland")) + # Verrucomicrobia (NO SIG)
  ylim(0,20) +
  scale_linetype_manual(name = "Land-use", values = c("Forest" = "solid",
                                                      "Dehesa" = "dotted",
                                                      "Open woodland" = "dashed")) +
  labs(x = "Defoliation Degree (%)", y = "Relative abundance (%)", title = "Verrucomicrobia")+
  geom_blank() + theme(plot.subtitle = element_text(vjust = 1), 
                       plot.caption = element_text(vjust = 1), 
                       axis.line = element_line(size = 0.4, 
                                                linetype = "solid"), panel.background = element_rect(fill = NA))


# Join all plots
library(gridExtra)
grid.arrange(p1_bact, p2_bact, p3_bact, p4_bact, p5_bact, p6_bact, nrow = 2)

