library(data.table)
library(phyloseq)
library(vegan)
library(ggplot2)

##%######################################################%##
#                                                          #
####                     Load data                      ####
#                                                          #
##%######################################################%##

data_site <- fread("Data/Bacteria/database.csv")

##%######################################################%##
#                                                          #
####                   Preparing data                   ####
#                                                          #
##%######################################################%##

data_site[,.N, by = .(Treatment)]
data_damag <- data_site[Treatment == "Damaged"]

## Reordering data_damag$Manejo_Conteo
data_damag$Manejo_Conteo <- factor(data_damag$Manejo_Conteo,
  levels = c("Forest", "Open woodland", "Dehesa")
)

# pH, SOC, N per land-use
data_damag[, .(pH_avg = round(mean(pH, na.rm = T),2), 
               SOC_avg =  round(mean(C_Org_kg_m, na.rm = T),2), 
               N_mineral_avg =  round(mean(N_mineral_inicial_mg_m, na.rm = T),2)), by = .(Manejo_Conteo)]

data_damag[, .(pH_se = round(st_error_NA(pH),2), 
               SOC_se =  round(st_error_NA(C_Org_kg_m),2), 
               N_mineral_se =  round(st_error_NA(N_mineral_inicial_mg_m),2)), by = .(Manejo_Conteo)]

##%######################################################%##
#                                                          #
####                  Statistics ANOVA                  ####
#                                                          #
##%######################################################%##

# By Land-use
pH_aov <- aov(pH ~ Manejo_Conteo, data_damag)
SOC_aov <- aov(C_Org_kg_m ~ Manejo_Conteo, data_damag)
Min_N_aov <- aov(N_mineral_inicial_mg_m ~ Manejo_Conteo, data_damag)

TukeyHSD(pH_aov, conf.level=0.95)
TukeyHSD(SOC_aov, conf.level=0.95)
TukeyHSD(Min_N_aov, conf.level=0.95)


# pH, SOC, N per land-use and defol
data_damag[, .(pH_avg = round(mean(pH, na.rm = T),2), 
               SOC_avg =  round(mean(C_Org_kg_m, na.rm = T),2), 
               N_mineral_avg =  round(mean(N_mineral_inicial_mg_m, na.rm = T),2)), by = .(Manejo_Conteo, Type)]

data_damag[, .(pH_se = round(st_error_NA(pH),2), 
               SOC_se =  round(st_error_NA(C_Org_kg_m),2), 
               N_mineral_se =  round(st_error_NA(N_mineral_inicial_mg_m),2)), by = .(Manejo_Conteo, Type)]

# By Land-use and tree health
# Forests
pH_for_aov <- aov(pH ~ Type, data_damag[Manejo_Conteo == "Forest"])
SOC_for_aov <- aov(C_Org_kg_m ~ Type, data_damag[Manejo_Conteo == "Forest"])
Min_N_for_aov <- aov(N_mineral_inicial_mg_m ~ Type, data_damag[Manejo_Conteo == "Forest"])

TukeyHSD(pH_for_aov, conf.level=0.95)
TukeyHSD(SOC_for_aov, conf.level=0.95)
TukeyHSD(Min_N_for_aov, conf.level=0.95)

# Open woodlands
pH_opw_aov <- aov(pH ~ Type, data_damag[Manejo_Conteo == "Open woodland"])
SOC_opw_aov <- aov(C_Org_kg_m ~ Type, data_damag[Manejo_Conteo == "Open woodland"])
Min_N_opw_aov <- aov(N_mineral_inicial_mg_m ~ Type, data_damag[Manejo_Conteo == "Open woodland"])

TukeyHSD(pH_opw_aov, conf.level=0.95)
TukeyHSD(SOC_opw_aov, conf.level=0.95)
TukeyHSD(Min_N_opw_aov, conf.level=0.95)

# Dehesa
pH_deh_aov <- aov(pH ~ Type, data_damag[Manejo_Conteo == "Dehesa"])
SOC_deh_aov <- aov(C_Org_kg_m ~ Type, data_damag[Manejo_Conteo == "Dehesa"])
Min_N_deh_aov <- aov(N_mineral_inicial_mg_m ~ Type, data_damag[Manejo_Conteo == "Dehesa"])

TukeyHSD(pH_deh_aov, conf.level=0.95)
TukeyHSD(SOC_deh_aov, conf.level=0.95)
TukeyHSD(Min_N_deh_aov, conf.level=0.95)

