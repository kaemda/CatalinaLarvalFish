#LOAD PACKAGES -----
#to see all functions ::

library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(sf)
library(vegan)
library(ggrepel)
library(xlsx)
library(factoextra)
library(FactoMineR)
library(indicspecies)
library(taxize)

# set working directory

setwd("C://KDale/Projects/CatalinaLarvalFish/")
theme_set(theme_classic(base_size = 14))
familyColors = c("Atherinidae" = "#B64A51",
                 "Blennioidei" = "#F86F49",
                 "Chaenopsidae" = "#E6C895",
                 "Engraulidae" = "#5E8564",
                 "Gobiidae" = "#83BFAC",
                 "Labrisomidae" = "#465b77",
                 "Pomacentridae" = "#361968",
                 "Other fishes" = "gray70")

# IMPORT DATA -----------------------------------------------------------------
# Raw datasets
fish.data=read_xlsx("Data/Catalina_Ichthyoplankton_Raw.xlsx")

# Light trapping metadata
trap.data = read_xlsx("Data/Light_trap_data.xlsx", sheet = "AllYears_Summary")

# Merge data and light trapping metadata
fish.data = merge(x = fish.data, y = trap.data, by.x = c("Sample"), by.y= c("Sample_ID"), all.x = TRUE) %>% 
  pivot_longer(., cols = c("Preflexion", "Flexion", "Postflexion", "Mid-transformation","UnknownAge"), names_to = "AgeClass", values_to = "N") %>% 
  subset(., AgeClass != "Mid-transformation") %>% 
  mutate(., N_mL = N / `Zooplankton volume`) %>% 
  mutate(., AgeClass = factor(AgeClass, levels = c("Preflexion", "Flexion", "Postflexion", "UnknownAge"))) %>% 
  mutate(AgeClassNum = as.numeric(as.factor(AgeClass)), Grouping = factor(Grouping, levels = c("Other fishes", "Atherinidae", "Blennioidei", "Chaenopsidae", "Engraulidae", "Gobiidae", "Labrisomidae", "Pomacentridae")))

write.xlsx(x = fish.data, "Data/Catalina_Icthyoplankton_Long.xlsx")

# Octopus data
octoData <- read_xlsx("Data/Octopus_Data.xlsx")
octo.data = merge(x = octoData, y = trap.data, by.x = c("Sample_ID"), by.y= c("Sample_ID"), all.x = TRUE)

# Calculate community matrix
communityMatrix <- group_by_at(fish.data, c("Sample", "MPA", "Year", "Season", "Substrate", "Season", "Family")) %>%
  subset(., Family != "Unknown") %>% 
  summarize(., N_total = max(N)) %>%
  pivot_wider(., names_from = Family, values_from = N_total) %>% 
  replace(is.na(.), 0)

totalPerSample <- group_by_at(fish.data, c("Sample", "MPA", "Site", "Year", "Season", "Substrate", "Grouping")) %>%
  summarize(., N = sum(N)) %>% 
  group_by_at(., c("Sample")) %>%
  summarize(., total = sum(N))

permanovaMatrix <- group_by_at(fish.data, c("Sample", "MPA", "Site", "Year", "Season", "Substrate", "Grouping")) %>%
  summarize(., N = sum(N)) %>% 
  merge(., totalPerSample) %>% 
  mutate(., percent = N/total) %>% 
  pivot_wider(., id_cols = -c(N, total), names_from = Grouping, values_from = percent) %>% 
  replace(is.na(.), 0)

# DENSITIES ------------------------------------------------------------------
densities <- function() {
  
  variables = c("MPA", "Season", "Substrate", "Year")
  responses = c("Num_Fish", "Fish_n_ml", "Num_Ceph", "Fish_n_ml")
  labels = c("N larvae","Larvae/mL","N paralarvae","Paralarvae/mL")
  plotlist = vector(mode = "list", length = length(variables) * length(responses))
  list.index = 1
  trap.data$Year <- as.factor(trap.data$Year)
  
  for (i in 1:length(variables)) {
    if (variables[i] == "Year") {
      data = subset(trap.data, Season == "Summer")
    } else {
      data = trap.data
    }
    for (j in 1:length(responses)) {
      t.test(data = data, as.formula(paste0(responses[j] ," ~ ", variables[i]))) %>% print()
      
      plotlist[[list.index]] <- ggplot(data = trap.data, aes_string(x = variables[i], y = responses[i])) +
        geom_boxplot(notch = FALSE, outlier.shape = 16, outlier.size = 3, lwd = 1, fatten = 1, fill = "gray50") 
      #labs(y = labels[j], x= "")
      list.index = list.index + 1
    }
  }
  
  #pdf("Figures/Fig#_SiteDensities/Fig#_Densities.pdf", width = 8, height = 8)
  print(cowplot::plot_grid(plotlist = plotlist, nrow = 4))
  #dev.off()
}

# DIVERSITY -------------------------------------------------------------------------
diversityIndices <- function () {
  
  variables = c("MPA", "Season", "Substrate")
  
  # Calculate effective number of species (inverse Simpson)
  communityMatrix$effective_N_species = vegan::diversity(x = communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)], index = "invsimpson") 
  communityMatrix <- mutate(communityMatrix, effective_N_species = replace(effective_N_species, is.infinite(effective_N_species), 0))
  
  # Paired t-test for richness: Are there more species?
  boxplot(specnumber(communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)]) ~ communityMatrix$MPA)
  t.test(specnumber(communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)]) ~ communityMatrix$MPA) %>% print()
  t.test(specnumber(communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)]) ~ communityMatrix$Substrate) %>% print()
  t.test(specnumber(communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)]) ~ communityMatrix$Season) %>% print()
  
  # Paired t-test for effective number of species
  t.test(communityMatrix$effective_N_species ~ communityMatrix$MPA) %>% print() # NO DIFFERENCE WITH PROTECTION
  t.test(communityMatrix$effective_N_species ~ communityMatrix$Season) %>% print() # NO DIFFERENCE IN SEASON
  t.test(communityMatrix$effective_N_species ~ communityMatrix$Substrate) %>% print() # NO DIFFERENCE IN SUBSTRATE (p = 0.05)
  
  summersOnly <- subset(communityMatrix, Season == "Summer")
  t.test(summersOnly$effective_N_species ~ summersOnly$Year) %>% print() # NO DIFFERENCE IN SUBSTRATE (p = 0.05)
  
}

# PERMANOVAS ----------------------------------------------------------------
permanovas <- function(variable) {
  
  formula = as.formula(paste('permanovaMatrix[,which(colnames(permanovaMatrix) == "Gobiidae"):ncol(permanovaMatrix)] ~ permanovaMatrix$', variable))
  
  adonis2(data = permanovaMatrix, formula = formula) %>% print()
  
}

# FISH AGE AVERAGES ----------------------------------------------------------
ageAverages <- function(variable, data) {
  
  # Construct initial dataframe to hold results
  familyAgeAverages <-
    data.frame(
      family = unique(data$Grouping),
      average_groupOne = 0,
      average_groupTwo = 0,
      NObsgroupOne = 0,
      NObsgroupTwo = 0,
      stdDev_groupOne = 0,
      stdDev_groupTwo = 0,
      p.value = 0,
      df = 0,
      t = 0
    )
  
  # For each family, get average age for the variable of interest
  for (i in 1:nrow(familyAgeAverages)) {
    
    # Remove fish of unknown age
    data = subset(data, Grouping == familyAgeAverages$family[i] & AgeClassNum != 4) 
    
    # Count number of each age class value
    data = group_by_at(data, c(variable, "AgeClassNum")) %>% 
      summarize(., N = sum(N))
    
    # Create a dataframe consisting of repeated age values (one line for each individual in that age category)
    average_groupOne <- bind_rows(data.frame(seq = seq.int(data$AgeClassNum[1], data$AgeClassNum[1], length.out = data$N[1])),
                                  data.frame(seq = seq.int(data$AgeClassNum[2], data$AgeClassNum[2], length.out = data$N[2])),
                                  data.frame(seq = seq.int(data$AgeClassNum[3], data$AgeClassNum[3], length.out = data$N[3])))
    
    # Calculate average age, number of observations, and the standard deviation for group one
    familyAgeAverages$average_groupOne[i] = mean(average_groupOne$seq)
    familyAgeAverages$NObsgroupOne[i] = length(average_groupOne$seq)
    familyAgeAverages$stdDev_groupOne[i] = sd(average_groupOne$seq)
    
    # Do the same thing for group two
    average_groupTwo <- bind_rows(data.frame(seq = seq.int(data$AgeClassNum[4], data$AgeClassNum[4], length.out = data$N[4])),
                                  data.frame(seq = seq.int(data$AgeClassNum[5], data$AgeClassNum[5], length.out = data$N[5])),
                                  data.frame(seq = seq.int(data$AgeClassNum[6], data$AgeClassNum[6], length.out = data$N[6])))
    familyAgeAverages$average_groupTwo[i] = mean(average_groupTwo$seq)
    familyAgeAverages$NObsgroupTwo[i] = length(average_groupTwo$seq)
    familyAgeAverages$stdDev_groupTwo[i] = sd(average_groupTwo$seq)
    
    # Only if both groups have >1 observation and the means are not the same (which would return an error),
    # calculate a two-tailed t-test and add results to the dataframe
    if(length(average_groupOne$seq) > 1 & length(average_groupTwo$seq) > 1 & mean(average_groupOne$seq) != mean(average_groupTwo$seq)) {
      print(familyAgeAverages$family[i])
      ttest <- t.test(average_groupOne$seq, average_groupTwo$seq)
      print(ttest)
      familyAgeAverages$p.value[i] = ttest$p.value
      familyAgeAverages$df[i] = ttest$parameter
      familyAgeAverages$t[i] = ttest$statistic
    } else {
      familyAgeAverages$p.value[i] = NA
      familyAgeAverages$df[i] = NA
      familyAgeAverages$t[i] = NA
    }
  }
  
  familyAgeAverages <- mutate(familyAgeAverages, family = factor(family, levels = c("Other fishes", "Atherinidae", "Blennioidei", "Chaenopsidae", "Engraulidae", "Gobiidae", "Labrisomidae", "Pomacentridae"))) %>% 
    pivot_longer(., cols = c("average_groupOne", "average_groupTwo", "stdDev_groupOne", "stdDev_groupTwo"), names_to = c(".value", "group"), names_sep = "_")
  
  ageAverageFigures(data = familyAgeAverages, variable = variable, varNames = ungroup(unique(data[,1])))
  
  write.xlsx(x = familyAgeAverages, file = paste0("Results/AgeAverage_t-tests_", variable, ".xlsx"))
  
  return(familyAgeAverages)
}


# PARALARVAE A:M T-TESTS ----------------------------------------------------------------------

paralarvae_age_tests <- function() {
  
  t.test(octo.data$Arm_mantle_ratio ~ octo.data$MPA) %>% print() # Older in Non MPA
  t.test(octo.data$Arm_mantle_ratio ~ octo.data$Substrate)%>% print() # Older in sandy
  
}

# FIG 1 MAPS ---------------------------------------------------------------------------
makeMaps <- function() {
  Catalina = read_sf("C://KDale/GIS/Catalina.shp")
  
  locations=read_xlsx("Data/Light_trap_locations.xlsx", sheet = 1)
  
  locations_sf = locations %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326)
  locations_sf_t = st_transform(locations_sf, crs=2163)
  
  ## Catalina map
  #jpeg(filename = "Catalina_map.jpg",width = 5,height = 5, "in",res= 300)
  pdf("Catalina_map.pdf", width = 5, height = 5)
  print(ggplot() +
          geom_sf(data = Catalina[1,]) +
          geom_sf_text(data = Catalina, aes(label = "Catalina")) +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          geom_sf(data = locations_sf_t) +
          ggsn::scalebar(data = Catalina, border.size = 0.5, location = "bottomleft", model = "WGS84", transform = T, dist_unit = "km", dist = 5, st.size=3, height=0.025,  box.color = "gray20", box.fill = c("gray20", "white"), st.color = "gray20") +
          labs(x = "", y = "") +
          geom_text_repel(data = locations, aes(x = longitude, y = latitude, label = `Site_Name`),
                          fontface = "bold", nudge_x = c(-0.025, 0.02, 0.05), nudge_y = c(-0.03, 0.03, 0)))
  dev.off()
  
  ## Catalina cropped
  catalina_cropped= st_crop(Catalina, xmin = -118.48, xmax = -118.52,
                            ymin = 33.425, ymax = 33.45,
                            expand=FALSE)
  
  pdf("Catalina_map_cropped.pdf", width = 5, height = 5)
  print(ggplot() +
          geom_sf(data = catalina_cropped) +
          geom_sf_text(data = catalina_cropped, aes(label = "Catalina"), fontface = "bold", 
                       nudge_x = c(0.01), nudge_y = c(-0.004)) +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          geom_sf(data = locations_sf) +
          ggsn::scalebar(data = California_crop, model = "WGS84", transform = T, dist_unit = "km", dist = 50, st.size=3, height=0.02, location = "bottomleft", box.color = "gray40", box.fill = c("gray40", "white"), st.color = "gray40") +
          geom_text_repel(data = locations, aes(x = longitude, y = latitude, label = `Site_Name`),
                          fontface = "bold", nudge_x = c(-0.008, 0.002, 0.008), nudge_y = c(-0.002, 0.005, 0.005)) +
          ggtitle("Light Trap Locations"))
  dev.off()
  
  ## Southern California map
  
  California= read_sf("C://KDale/GIS/North_South_America.shp")
  locations_sf = locations %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326)
  California_crop=  st_crop(California, xmin = -117, xmax = -121.28,
                            ymin = 32.56, ymax = 35.0,
                            expand=FALSE)
  pdf(file = "Figures/Fig1_Map/Southern_California.pdf",width = 5,height = 5)
  print(ggplot() +
          geom_sf(data = California_crop) +
          geom_sf_text(data = California_crop, aes(label ="  ")) +
          xlab(expression(paste("longitude"))) +
          ylab(expression(paste("latitude "))) +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          ggtitle("Southern California"))
  dev.off()
  
  # Locations of other surveys
  northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp") %>% st_union() %>% st_transform(., crs = "EPSG:4326")
  
  locations_otherSurveys = read_xlsx("Data/Other Surveys/OtherSurveys_Location.xlsx", sheet = 1) %>%  st_as_sf(coords = c("longitude", "latitude"), crs=4326)
  locations_otherSurveys$survey = factor(locations_otherSurveys$survey, levels = c(c("Light traps", "Avendaño-Ibarra et al. 2004", "Stephens et al. 1984", "SMURFs", "CalCOFI")))
  California_crop=  st_crop(California, xmin = -117.5, xmax = -121.28,
                            ymin = 33, ymax = 34.5,
                            expand=FALSE)
  pdf("Figures/Fig#_Othersurveys_comparison/Fig#_Locations.pdf", width = 5, height = 7)
  print(ggplot() +
          geom_sf(data = northAmerica, fill = "gray70") +
          ggsn::scalebar(data = locations_otherSurveys, border.size = 0.5, st.dist = 0.03, model = "WGS84", transform = T, dist_unit = "km", dist = 100, st.size=3, height=0.02, location = "topright", box.color = "gray20", box.fill = c("gray20", "white"), st.color = "gray20") +
          geom_sf(data = locations_otherSurveys, mapping = aes(fill = survey), pch = 21, size = 3, color = "gray50", alpha = 0.8) +
          scale_fill_manual("Survey", values = c("Light traps" = "black", "SMURFs" = "goldenrod","CalCOFI" = "deepskyblue3", "Stephens et al. 1984" = "indianred4", "Avendaño-Ibarra et al. 2004" = "tomato")) +
          theme_classic() +
          xlim(min = -121.28, max= -111) +
          ylim(min = 24, max = 34.5) +
          labs(x = "", y = "") +
          theme(axis.text.x = element_text(angle = 60,hjust=1), legend.position = c(0.3, 0.17), legend.text = (element_text(size = 12))))
  dev.off()  
}

# FIG 3: FAMILIES -----------------------------------------------------------
familySummary <- function(data) {
  
  # Order families by decreasing total abundance
  data$Family <- as.factor(data$Family)
  familySummary <- data %>% group_by(., Family) %>% summarise(N = sum(N))
  familySummary$Family=factor(x = familySummary$Family,levels = familySummary$Family[order(familySummary$N,decreasing = TRUE)])
  
  pdf("Figures/Fig3_Families/Fig3_Family_summary_wide.pdf", width = 9, height = 5)
  print(ggplot(data = familySummary)+
          geom_col(mapping = aes(x= Family,y= log(N)))+
          labs(x="Taxa: ",y="Log (N fish)")+
          theme_classic(base_size = 16)+
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          scale_y_continuous(expand = expansion(mult = c(0,.1))))
  dev.off()
}

# FIG 4: AGE AVERAGES ---------------------------------------------------
ageAverageFigures <- function(data, variable, varNames) {
  
  file = paste0("Figures/Fig4_AgeAverages/Fig4_AgeAverages_", variable, ".pdf")
  
  pdf(file = file, width = 6, height = 5)
  print(ggplot(subset(data, !is.na(p.value) & family != "Unknown" & family != "Other fishes"), mapping = aes_string(x = "group", y = "average", group = "family", col = "family")) +
          geom_point(size = 3, position = position_dodge(width = .3)) +
          geom_line(lwd = 1, lty = 2, position = position_dodge(width = .3)) +
          geom_linerange(aes(x = group, ymin = average-stdDev, ymax = average+stdDev), lwd = 1, position = position_dodge(width = .3)) +
          labs(x = "", y = "Average age") +
          scale_x_discrete(labels= varNames[[variable]]) +
          scale_color_manual("Taxonomic grouping", values = familyColors) +
          ggtitle(variable))
  dev.off()
}

# FIG 5 A:M RATIO ----------------------------------------------------------------------
paralarvae_age_boxplots <- function() {
  
  a <- ggplot(octo.data, aes(y= Arm_mantle_ratio, x = MPA)) +
    geom_boxplot(color = "black", fill = "gray50", lwd = 1, fatten = 1, outlier.shape = 16, outlier.size = 3) +
    theme_classic(base_size = 14) +
    ggtitle("A. Protection status") +
    labs(x = "", y = "Arm:mantle length")
  
  b <- ggplot(octo.data, aes(y= Arm_mantle_ratio, x = Substrate)) +
    geom_boxplot(color = "black", fill = "gray50", lwd = 1, fatten = 1, outlier.shape = 16, outlier.size = 3) +
    theme_classic(base_size = 14) +
    ggtitle("B. Substrate") +
    labs(x = "", y = "")
  
  pdf("Figures/Fig5_A-M Ratio/Fig5_A-M Ratio.pdf", width = 7, height = 5)
  print(cowplot::plot_grid(a,b,nrow = 1))
  dev.off()
}

# FIG # COMPARE CALCOFI --------------------------------------------------------
compareDatasets <- function() {
  
  # SMURF data - raw
  # Available at https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_UCSB_Fish_Recruitment.1.3
  smurf.raw <- read_xlsx("Data/Other Surveys/PISCO_UCSB_subtidal_recruitment_fish_data.1.2.xlsx", sheet =  1)
  smurf.spp <- read_xlsx("Data/Other Surveys/PISCO_UCSB_subtidal_recruitment_fish_spp_list.1.2.xlsx", sheet =  1)
  smurf.sites <- read_xlsx("Data/Other Surveys/PISCO_UCSB_subtidal_recruitment_site_list.1.3.xlsx", sheet =  1)
  
  smurf.dat <- merge(smurf.raw, smurf.spp) %>% merge(., smurf.sites) %>% subset(total_fish_collected > 0) %>% 
    rename(total_larvae = total_fish_collected, latitude = LAT_WGS84, longitude = LONG_WGS84) %>% 
    subset(latitude < 34.1) %>% # Select only Channel Islands stations
    mutate(survey = "SMURFs") %>% 
    group_by_at(c("family", "survey")) %>% 
    summarize(total_larvae = sum(total_larvae)) %>% 
    ungroup() %>% 
    mutate(percentage_abs = total_larvae/sum(total_larvae)) %>% 
    mutate(percentage = percentage_abs * -1)
  
  # # CalCOFI data - raw
  # calcofi.coast <- read.xlsx(file = "Data/Other Surveys/CalCOFI_Line90_Stations27.7-28.xlsx", sheetName = "Summary") %>% 
  #   mutate(survey = "CalCOFI Sta 27.7, 28 (nearshore)")
  # calcofi.cat <- read.xlsx(file = "Data/Other Surveys/CalCOFI_Line90_Stations35-37.xlsx", sheetName = "Summary") %>% 
  #   mutate(survey = "CalCOFI Sta 35, 37 (Catalina)")
  # 
  # # Get family information for each scientific name
  # calcofi.taxa <- unique(calcofi.coast$scientific_name, calcofi.cat$scientific_name) %>% subset(!is.na(.))
  # families <- tax_name(calcofi.taxa, get = "family")
  # 
  # # Summarize by family and add percentage of total
  # calcofi.coast <-  calcofi.coast %>% 
  #   merge(., families[c("query", "family")], by.x = "scientific_name", by.y = "query") %>% 
  #   group_by_at(c("family", "survey")) %>% summarize(total_larvae = sum(total_larvae)) %>% 
  #   ungroup() %>% 
  #   mutate(percentage = total_larvae/sum(total_larvae))
  # 
  # # Summarize by family and add percentage of total
  # calcofi.cat <- calcofi.cat %>% 
  #   merge(., families[c("query", "family")], by.x = "scientific_name", by.y = "query") %>% 
  #   group_by_at(c("family", "survey")) %>% summarize(total_larvae = sum(total_larvae)) %>% 
  #   ungroup() %>% 
  #   mutate(percentage = total_larvae/sum(total_larvae))
  # 
  # # Combine CalCOFI datasets
  # calcofi.family.summary <- bind_rows(calcofi.coast, calcofi.cat) %>% 
  #   mutate(percentage_abs = percentage, total_larvae_log_abs = log(total_larvae)) %>% 
  #   mutate(percentage = percentage_abs*-1, total_larvae_log = total_larvae_log_abs*-1)
  # 
  # # Read/write CalCOFI data
  # write.xlsx(x = comparison, file = "Data/Other Surveys/CalCOFI_family_summary.xlsx")
  calcofi.family.summary <- read_xlsx("Data/Other Surveys/CalCOFI_family_summary.xlsx")
  
  # Stephens et al 1984 data
  stephens <- read_xlsx("Data/Other Surveys/Stephens_1984.xlsx", sheet = 1) %>% 
    group_by_at("family") %>%  
    summarize(total_larvae = sum(total_larvae)) %>% 
    mutate(survey = "Stephens et al. 1984",  percentage_abs = total_larvae/sum(total_larvae)) %>% 
    mutate(total_larvae_log = log(total_larvae), percentage = percentage_abs*-1, total_larvae_log_abs = log(total_larvae))
  
  # Avendano-Ibarra et al 2004
  avendano <- read_xlsx("Data/Other Surveys/Avendano-Ibarra_2004.xlsx", sheet = 1) %>% 
    group_by_at("family") %>% 
    summarize(percentage_abs = sum(percentage)) %>% 
    mutate(survey = "Avenda\u00f1o-Ibarra et al. 2004") %>% 
    mutate(percentage = percentage_abs*-1)
  
  # Summarize light trap data
  family.summary <- fish.data %>% group_by_at("Family") %>% summarize(total_larvae = sum(N)) %>% 
    rename(family = Family) %>% 
    mutate(survey = "Light traps", percentage = total_larvae/sum(total_larvae), percentage_abs = total_larvae/sum(total_larvae)) %>% 
    mutate(total_larvae_log = log(total_larvae), total_larvae_log_abs = log(total_larvae))
  
  # Combine datasets
  comparison <- bind_rows(calcofi.family.summary, smurf.dat, stephens, avendano) %>% 
    subset(., !is.na(family) & family != "Unknown" & percentage_abs > 0.005) # Subset to families that make up > 0.1% of catch
  comparison$survey = factor(comparison$survey, levels = c("Light traps", "Avenda\u00f1o-Ibarra et al. 2004", "Stephens et al. 1984",   "SMURFs", "CalCOFI Sta 35, 37 (Catalina)", "CalCOFI Sta 27.7, 28 (nearshore)"))
  comparison$family = as.factor(comparison$family)
  comparison$family = reorder(comparison$family, comparison$percentage_abs, FUN = sum)
  
  # Bar chart showing percentages for each survey type
  pdf(file = "Figures/Fig#_Othersurveys_comparison/Fig#_Othersurveys_comparison.pdf", width = 8.5, height = 11)
  print(ggplot() +
          geom_col(comparison, mapping = aes(y = family, x = percentage, fill = survey), width = .80) +
          scale_fill_manual("Survey", values = c("black", "tomato", "indianred4", "goldenrod","deepskyblue1", "deepskyblue4")) +
          #scale_color_manual("Survey", values = c("black","tomato", "indianred4",  "goldenrod","deepskyblue1", "deepskyblue4")) +
          labs(x = "Percentage of total catch", y = "") +
          theme_classic(base_size = 16) +
          scale_x_continuous(breaks = c(-1.5, -1.0, -0.5, 0, 0.5), labels = c(1.5, 1, 0.5, 0, 0.5)) +
          theme(axis.text = element_text(size = 14), legend.text=element_text(size=14), axis.text.y = element_text(size = 14), legend.position = c(0.27,0.1)))
  dev.off()
  
}

# SUPP FIG: COMMUNITY COMPOSITIONS ---------------------------------------------
communityComposition <- function(variable) {
  
  # Summarize by general taxonomic categories
  pivot_by_grouping <- group_by_at(fish.data, c(variable, "Grouping")) %>% 
    summarize(N = sum(N))
  
  # Get total abundance
  totals <- group_by_at(pivot_by_grouping, variable) %>% summarize(N_total = sum(N))
  
  # Calculate percentage of total for each grouping
  pivot_by_grouping <- merge(pivot_by_grouping, totals, by = variable) %>% 
    mutate(., freq = N/N_total)
  
  plot <- ggplot(pivot_by_grouping) +
    geom_col(mapping = aes_string(x = variable, y = "freq", fill = "Grouping")) +
    theme_classic(base_size = 12) +
    ggtitle("All fish larvae", subtitle =   variable) +
    labs(y = "Proportion", x = "") +
    scale_fill_manual("Taxonomic grouping", values = familyColors)
  
  pdf(file = paste0("Figures/SupplementalFig_Proportions/All_larvae_", variable, ".pdf"), width = 4.15, height = 4.25)
  print(plot)
  dev.off()
}

# RUN FUNCTIONS ----------------------------------------------------------------------
makeMaps()
densities()
familySummary(data = fish.data)
diversityIndices()
paralarvae_age_tests()
paralarvae_age_boxplots()

variables = c("MPA","Substrate", "Season", "Year")
for (i in 1:length(variables)) {
  permanovas(variable = variables[i])
  communityComposition(variable = variables[i])
}

# Average ages across categories
familyAgeAvgMPA <- ageAverages(variable = "MPA", data = fish.data)
familyAgeAvgSeasons <- ageAverages(variable = "Season", data = fish.data)
familyAgeAvgSubstrate <- ageAverages(variable = "Substrate", data = fish.data)

# SANDBOX ---------------------------------------------------------------------


write.xlsx(x= comparison, file = "Data/Other Surveys/CalCOFI_comparison.xlsx")
# PCAs
pca <- prcomp(permanovaMatrix[,8:ncol(permanovaMatrix)], scale. = TRUE, center = TRUE)
#pca <- prcomp(communityMatrix[,which(colnames(communityMatrix) == "Gobiidae"):ncol(communityMatrix)], scale. = TRUE, center = TRUE)

# MCAs
mca.dat <- communityMatrix[1:nrow(communityMatrix), 2:6] %>% ungroup()
mca <- MCA(X = mca.dat)

fviz_mca_biplot(res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal())

summary(pca)
fviz_pca_biplot(pca, geom = "point", habillage = permanovaMatrix$MPA, addEllipses = TRUE)
#fviz_pca_biplot(pca, geom = "point", habillage = communityMatrix$N_S, addEllipses = TRUE)

permanovaMatrix <- bind_cols(permanovaMatrix, pca$scores[,1:2])

sites <- read.xlsx("Data/Light_trap_locations.xlsx", sheetIndex =  1)

light <- merge(light, sites)

source("C://KDale/Projects/Phenology/Code/linkRoms.R")

test <- linkroms(tows = trap.data)

