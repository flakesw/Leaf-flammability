#*******************************************************************************
# Data Preparation
# Sam Flake, 18 January 2023
# 
# This script brings in raw data from experimental burns of leaf litter,
# measurements of leaf dimensions from the lab, measurements of leaf dimensions
# obtained from 3D scans, and leaf nutrient data. It outputs clean datasets
# at the leaf level, individual level, burn level, and species level
#*******************************************************************************
library("tidyverse")
library("units")

# Process burn-level data -----------------------------------------------------
burn_sum <- read.csv("./raw data/BurnSummary19July2019_NEW.csv")
burn_sum$IndCode <- paste0(burn_sum$Species, burn_sum$sample)

#remove leaves which didn't burn at all
#burn_sum <- burn_sum[burn_sum$Fireline.intensity.kWm.1 > 0 & burn_sum$PercentCombustion > 0.5, ]

burn_sum$success <- ifelse(burn_sum$Fireline.intensity.kWm.1 > 0, TRUE, FALSE)

#remove small leaves
burn_sum <- subset(burn_sum, !(Species %in% c("ANFA", "STRO", "JACA")))


burn_sum$litter_spec_vol <- (burn_sum$Average.Litter.Depth/100) #m3 kg-1
burn_sum$litter_spec_vol <- burn_sum$litter_spec_vol * 1000 #convert to cm3 g-1
burn_sum$litter_bulk_density <- 1/burn_sum$litter_spec_vol 
burn_sum$log_litter_spec_vol <- log(burn_sum$litter_spec_vol) 

# Process individual leaf data -------------------------------------------------
leaf_data <- read.csv("./raw data/Leaf Scans.csv", skip = 1, header = TRUE, stringsAsFactors = FALSE)
dry_leaf_data <- read.csv("./raw data/Dry Leaf Measurments.csv")
cloud_data <- read.csv("./clean data/cloud_data_per_leaf.csv", header = TRUE)

leaf_data <- left_join(leaf_data[, -c(8, 9, 10, 15, 17, 21:24)], cloud_data[, -c(1:3)], by = c("IndLeafCode"))

# leaf_data <- read.csv("./clean data/all_leaf_data_with_sla_and_3d_sf_2020-07-01.csv")
leaf_data$Ind <- strsplit(leaf_data$IndLeafCode, "-") %>% 
  map(1) %>% #pluck() the first element of each element of the list of names & leaf numbers
  unlist() %>%
  as.factor()
leaf_data$dryHull <- leaf_data$dryHull/1000 #convert from mm3 to cm3

#calculate leaf volumes (for Hypothesis 1)

leaf_data$particle_volume <- leaf_data$Areacm2 * leaf_data$FreshTH /10 #in cm3
leaf_data$volume_ratio <- leaf_data$freshHull / leaf_data$particle_volume
leaf_data$hull_specific_vol <- leaf_data$dryHull / leaf_data$DryWeight #units are cm3 g-1
leaf_data$prismatic_volume <- leaf_data$Areacm2 * leaf_data$FreshH #cm3
leaf_data$rectangular_volume <- leaf_data$FreshH * leaf_data$FreshL * leaf_data$FreshW #cm3

#calculate leaf sizes and curls
leaf_data$S1 <- leaf_data$Areacm2^0.5
leaf_data$C1 <- leaf_data$DryH / (leaf_data$S1)
leaf_data$C1H <- ((leaf_data$dryHull/1000)^(1/3))/(leaf_data$S1)
leaf_data$C1Hfresh <- ((leaf_data$freshHull/1000)^(1/3))/(leaf_data$S1)
leaf_data$C1H_new <- (leaf_data$dryHull/1000)/(leaf_data$Areacm2^(3/2))

leaf_data$S2 <- (leaf_data$DryL * leaf_data$DryW)^0.5
leaf_data$C2 <- leaf_data$DryH / leaf_data$S2
leaf_data$C2H <- ((leaf_data$dryHull/1000)^(1/3))/leaf_data$S2
leaf_data$C2Hfresh <- ((leaf_data$freshHull/1000)^(1/3))/(leaf_data$FreshL * leaf_data$FreshW)^0.5
leaf_data$C2H_new <- (leaf_data$dryHull/1000)/(leaf_data$Areacm2^(3/2))

leaf_data$S3 <- leaf_data$DryL
leaf_data$C3 <- leaf_data$DryH / leaf_data$S3
leaf_data$C3H <- ((leaf_data$dryHull/1000)^(1/3))/leaf_data$S3
leaf_data$C3Hfresh <- ((leaf_data$freshHull/1000)^(1/3))/leaf_data$FreshL

# process leaf nutrient data ---------------------------------------------------

nutrients <- read.csv("./raw data/W.Hoffmann-Flake7322_022820.csv", stringsAsFactors = FALSE)
nutrients <- nutrients[-1, ]
nutrients[159, 12] <- NA
nutrients[, 5:12] <- apply(nutrients[, 5:12], 2, as.numeric)

unique(burn_sum$Species)[!(unique(burn_sum$Species) %in% nutrients$Code)]
# 71 of 91 species present in the nutrient analysis
#Missing: CASYf, ZETU, ASTO, TAAU, ERTO, AESE, ERDE, HEAP, MAGU, CRFL, MISE, XYAR, FAMO,
# BYLA, PORA, SYRE, DAFA, BYIN, GUVI, CHLA

nutrients_agg <- nutrients %>%
  group_by(Code) %>%
  dplyr::summarise(Carbon = mean(Carbon, na.rm = TRUE),
                   Nitrogen = mean(Nitrogen, na.rm = TRUE),
                   Al.Conc. = mean(Al.Conc., na.rm = TRUE),
                   All_ash = sum(P.Conc., K.Conc., Ca.Conc., Mg.Conc., Fe.Conc., Al.Conc.))

# Aggregate and join data ---------------------------------------------------------------
leaf_data$IndCode <- leaf_data$IndLeafCode %>% 
  strsplit(., "-") %>%
  map(1) %>%
  unlist()

#aggregate leaf data to species level
leaf_sp_level <- leaf_data %>%
  group_by(Species) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = function(x) mean(x, na.rm = TRUE)))

#aggregate data to burn level, using species-level leaf data
burn_level <- left_join(burn_sum, nutrients_agg, c("Species" = "Code")) %>%
  left_join(., leaf_sp_level, by = c("Species"))

#calculate packing ratio -- number of leaves * leaf particle volume / total fuel volume
burn_level$packing_ratio <- (150/burn_level$DryWeight * burn_level$particle_volume)/(60*25*burn_level$Average.Litter.Depth)

#bring in trait data from Flake et al. 2021 Functional ecology
#TODO: needed?
traits <- read.csv("./raw data/clean_species2020-06-09.csv")
burn_level <- left_join(burn_level, traits[, c("FG", "Code")], by = c("Species" = "Code"))

# Create species-level data ----------------------------------------------------
burn_sp <- burn_level %>%
  group_by(Species) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean),
            FG = FG[1],
            n_success = sum(success),
            prop_success = n_success / 3,
            accumulator = ifelse(Al.Conc. > 1000, TRUE, FALSE))

accum_sp <- burn_sp[which(burn_sp$accumulator), "Species"]
burn_level$accumulator <- ifelse(burn_level$Al.Conc. > 1000, TRUE, FALSE)



write.csv(leaf_data, "./clean data/leaf_level_data.csv")
write.csv(leaf_sp_level, "./clean data/leaf_sp_level_data.csv")
write.csv(burn_level, "./clean data/burn_level_data.csv")
write.csv(burn_sp, "./clean data/burn_sp_level_data.csv")

