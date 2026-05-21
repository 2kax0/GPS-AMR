library(tidyverse)
library(scatterpie)
library(ggpubr)
library(meta)
library(maps)
library(ggpattern)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggnewscale)
library(ggpattern)

filter <- dplyr::filter

# set working directory
setwd("C:/Users/Max Choong/OneDrive/Desktop/GPS-AMR/")

#-------------------------------------------- read in files -------------------------------------------

# read in colour schemes for serotypes
colour <- read.csv(file = "serotype_colours.csv")

# read in metadata
metadata <- read.csv(file = "monocle-metadata - Copy.csv", na.strings = c("","NA"))

# tidy antibiotics for plotting
antibiotic_tidying <- read.csv(file = "antibiotic.csv")
#----------------------------------------- metadata cleanup -----------------------------------------

# fill in missing continent data
metadata$Continent <- ifelse(metadata$Country %in% "NIGERIA" & is.na(metadata$Continent), "AFRICA",
                             ifelse(metadata$Country %in% "THAILAND" & is.na(metadata$Continent), "ASIA", metadata$Continent))

# filter out pre-1990 sequences
### majority of pre-1990 sequences were from Norway
metadata <- metadata |>
  filter(Year > 1989)

# country name tidying
metadata$Country <- str_to_title(metadata$Country)

metadata$Country_mapping <- ifelse(metadata$Country %in% "Russian Federation", "Russia",
                           ifelse(metadata$Country %in% "Trinidad And Tobago", "Trinidad and Tobago",
                                  ifelse(metadata$Country %in% "United States", "United States of America",
                                         ifelse(metadata$Country %in% "The Gambia", "Gambia", metadata$Country))))

metadata$GPSC <- ifelse(is.na(metadata$GPSC), "unknown", metadata$GPSC)

metadata$Vaccine_period_clean <- gsub("\\-.*", "", metadata$Vaccine_period)

metadata$Vaccine_period_clean <- factor(metadata$Vaccine_period_clean, 
                                            levels = c("PREPCV", "POSTPCV7", "POSTPCV10", "POSTPCV13"))

library(usethis)
create_github_token()
gitcreds::gitcreds_set()
#----------------------------------------- maps -----------------------------------------
################# need to set consistent scale


# initialise map and removing antarctica
world_map <- map_data("world")
world_map <- subset(world_map, region != "Antarctica")

# tally table of number of isolates from each country, then rename for joining
metadata_tally <- metadata |>
  group_by(Country_mapping) |>
  tally() |>
  rename(
    sovereignt = Country_mapping
  )

# heatmap
### potential issue: trinidad and tobago joined together in dataset but split in map data

# initialise worldmap
world <- ne_countries(scale = "medium", returnclass = "sf")

# join tally with world data
world_annotated <- left_join(world, metadata_tally, by = "sovereignt")

# remove antarctica from world data
world_annotated <- world_annotated |>
  filter(!sovereignt %in% "Antarctica")

# world map
world_plot <- ggplot(data = world_annotated) +
  geom_sf(fill = "white") +
  new_scale_fill() +
  geom_sf(aes(fill = log10(n))) +
  scale_fill_gradientn(colours = c("blue", "steelblue", "pink", "red"), na.value = "white", limits = c(0, 4.2)) +
  theme_minimal() +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

world_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "world plot.svg", fsep = "/")

ggsave(world_plot, filename = world_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")

# zoomed in maps for europe
europe <- world_annotated |>
  filter(continent %in% "Europe" & !sovereignt %in% "Russia")

eu_world_plot <- ggplot(data = europe) +
  geom_sf(fill = "white") +
  new_scale_fill() +
  geom_sf(aes(fill = log10(n))) +
  scale_fill_gradientn(colours = c("blue", "steelblue", "pink", "red"), na.value = "white", limits = c(0, 4.2)) +  
  theme_minimal() +
  xlim(-40, 42) +
  ylim(35, 73) +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

eu_world_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "eu world plot.svg", fsep = "/")

ggsave(eu_world_plot, filename = eu_world_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")

# zoomed in map of sub-saharan africa
west_africa <- world_annotated |>
  filter(region_wb %in% "Sub-Saharan Africa")

ssa_world_plot <- ggplot(data = west_africa) +
  geom_sf(fill = "white") +
  new_scale_fill() +
  geom_sf(aes(fill = log10(n))) +
  scale_fill_gradientn(colours = c("blue", "steelblue", "pink", "red"), na.value = "white", limits = c(0, 4.2)) +  
  theme_minimal()+
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

ssa_world_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "ssa world plot.svg", fsep = "/")

ggsave(ssa_world_plot, filename = ssa_world_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")

# zoomed in map of south and southeast asia
south_asia <- world_annotated |>
  filter(subregion %in% c("Southern Asia", "South-Eastern Asia"))

sea_world_plot <- ggplot(data = south_asia) +
  geom_sf(fill = "white") +
  new_scale_fill() +
  geom_sf(aes(fill = log10(n))) +
  scale_fill_gradientn(colours = c("blue", "steelblue", "pink", "red"), na.value = "white", limits = c(0, 4.2)) +  
  theme_minimal()+
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

sea_world_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "sea world plot.svg", fsep = "/")

ggsave(sea_world_plot, filename = sea_world_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")

#----------------------------------------- vaccine reference table ----------------------------------------

# read in table of serotypes covered by each vaccine
vaccine_table <- read_csv(file = "vaccine_chart.csv")

# create a column with values as 1
vaccine_table$value <- 1

# reorder factors
vaccine_table$vaccine <- factor(vaccine_table$vaccine, levels = c("PCV20", 
                                                                  "PCV15",
                                                                  "PCV13",
                                                                  "PCV10_PNEUMOSIL",
                                                                  "PCV10_GSK",
                                                                  "PCV7"))

# labels
vaccine_labels <- c("PCV20", 
                    "PCV15",
                    "PCV13",
                    "PCV10 (PNEUMOSIL)",
                    "PCV10 (GSK)",
                    "PCV7")

# vaccine reference table
vaccine_table |>
  complete(serotype, vaccine) |>
  ggplot(aes(x = serotype, y = vaccine, fill = value)) +
  geom_tile(colour = "black") +
  theme_minimal() +
  scale_fill_gradient(na.value = "white", guide = "none") +
  scale_y_discrete(labels = vaccine_labels) +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 12)
  )

pcv7 <- vaccine_table |>
  filter(vaccine %in% "PCV7")

pcv7_serotypes <- pcv7$serotype

pcv10 <- vaccine_table |>
  filter(vaccine %in% "PCV10_GSK")

pcv10_serotypes <- pcv10$serotype

pcv13 <- vaccine_table |>
  filter(vaccine %in% "PCV13")

pcv13_serotypes <- pcv13$serotype

pcv20 <- vaccine_table |>
  filter(vaccine %in% "PCV20")

pcv20_serotypes <- pcv20$serotype


#----------------------------------------- gpsc/serotype plot -----------------------------------------
# tally table of gpscs in the dataset
GPSC_tally <- metadata |>
  group_by(GPSC) |>
  tally()

# 40 most common gpsc in the dataset
GPSC_40 <- GPSC_tally[order(GPSC_tally$n, decreasing = TRUE),] |>
  filter(!is.na(GPSC) & !GPSC %in% "unknown") |>
  slice_head(n = 40) 

# save as vector
n_40 <- unique(GPSC_40$GPSC)

# classify as others
tidy <- metadata |>
  mutate(
    GPSC_tidy = ifelse(GPSC %in% n_40, GPSC, "NA")
  ) |>
  filter(!In_silico_serotype %in% "_" & !is.na(GPSC_tidy)) |>
  relocate(GPSC_tidy, .after = GPSC)


# read in file containing tidied serotypes
clean_serotype <- read_csv(file = "tidied_serotypes.csv")

# left join to match serotype with tidied serotype
tidy <- left_join(tidy, clean_serotype, by = "In_silico_serotype") |>
  relocate(serotype_clean, .after = In_silico_serotype)

# replace NAs in clean serotype column with in silico serotype
tidy$serotype_clean <- ifelse(is.na(tidy$serotype_clean), tidy$In_silico_serotype, tidy$serotype_clean)


# tally table of serotypes in the 40 most abundant gpsc 
graph_dataset <- tidy |>
  group_by(GPSC_tidy, serotype_clean) |>
  tally()

# read in serotype colours
serotype_colours <- read_csv(file = "serotype_colours.csv")

# join in colours
graph_dataset <- left_join(graph_dataset, serotype_colours, by = "serotype_clean")

serotype_colors_vector <- setNames(graph_dataset$In_silico_serotype__colour, graph_dataset$serotype_clean)

# annotate vaccine types
graph_dataset$vaccine <- ifelse(graph_dataset$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")

# barchart
graph_dataset |>
  filter(!GPSC_tidy %in% "NA") |>
  ggplot() +
  (aes(x = as.factor(GPSC_tidy), y = n, fill = serotype_clean, pattern = vaccine)) +
  geom_col_pattern(
    pattern_spacing = 0.01, 
    pattern_fill = "black", 
    pattern_size = 0.5, 
    pattern_linetype = 0.5,
    pattern_key_scale_factor = 0.5,
    pattern_orientation = "vertical", show.legend = TRUE) +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(override.aes = list(pattern = "none"), ncol = 3)) +
  scale_pattern_manual(
    name = "Vaccine",
    values = c("NVT" = "stripe", "VT" = "none")
  )  +
  labs(
    title = "Serotype and antibiotic resistance within the 40 most abundant GPSCs",
    x = "GPSC",
    y = "Number of isolates",
    fill = "Serotype"
  )

continent_year_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "continent year plot.svg", fsep = "/")

ggsave(continent_year_plot, filename = continent_year_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")


#----------------------------------------- year/manifest plot ----------------------------------------

# isolates with >=15 or >5 will be labelled as NA 
# re-classify and tidying up age and clinical manifestation
metadata_age_manifest <- tidy |>
  mutate(
    Age_year_num = as.integer(Age_years),
    Age_category = case_when(
      Age_years %in% c("<=5", "<1", ">=2 & <5", ">2<=5", "<=2", "<5", ">2<=5", "0-2") ~ "< 5",
      Age_years %in% c(">15<=24", ">24<=44", ">44<=65", ">5<=15", "15-44", "23-64", "34-48", "45-64", "49-64", "5-14", "3-5", "6-12", "4-12") ~ ">=5 - <=65",
      Age_years %in% c(">65", ">=65", "65-78", "79-98") ~ "> 65",
      Age_years %in% "_" & Age_months %in% "_" & Age_days %in% "_" ~ "NA",
      !is.na(Age_year_num) & Age_year_num < 5 ~ "< 5",
      !is.na(Age_year_num) & Age_year_num >= 5 & Age_year_num <= 65 ~ ">=5 - <=65",
      !is.na(Age_year_num) & Age_year_num > 65 ~ "> 65",
      is.na(Age_year_num) & !is.na(Age_months) ~ "< 5",
      is.na(Age_year_num) & is.na(Age_months) & !is.na(Age_days) ~ "< 5",
      TRUE ~ NA
    ),
    Manifest_category = ifelse(Manifestation %in% "UNKNOWN" & Clinical_manifestation %in% c("ARTHRITIS", "SEPTIC ARTHRITIS", "MENINGITIS", "SEPSIS"), "IPD", Manifestation)
  ) |>
  relocate(Manifest_category, .after = Manifestation) 

# annotate vaccine types
metadata_age_manifest$vaccine <- ifelse(metadata_age_manifest$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")


# table with continent manifestation for plotting
continent_year_tally <- metadata_age_manifest |>
  filter(!is.na(Year) & !is.na(Continent) & !Year %in% "_" & !Manifest_category %in% "UNKNOWN") |>
  group_by(Continent, Year, Manifest_category) |>
  tally() |>
  mutate(
    join = paste0(Continent, " ", Manifest_category)
  )

# relevel factors
continent_year_tally$Manifest_category <- factor(continent_year_tally$Manifest_category, levels = c("CARRIAGE", "DIS_CAR", "IPD", "NONINVASIVE_DISEASE", "UNKNOWN"))

# manifest, time of collection and continent plot
continent_year_plot <- continent_year_tally |>
  filter(as.integer(Year) > 1989 & !Continent %in% "CENTRAL AMERICA") |>
  ggplot(aes(x = as.factor(Year), y = n, fill = Manifest_category)) +
  geom_bar(stat = "identity") +
  facet_grid(Continent ~ ., scales = "free_y") +
  scale_fill_manual(values = c("#BA55D3", "#2ECC71", "#00BFFF", "#FF8C00"),
                    labels = c("Carriage", "Diseased carriage", "IPD", "Non IPD")) +
  labs(
    title = "Clinical manifestation of isolates collected in each region",
    x = "Year",
    y = "Number of isolates",
    fill = "Clinical \nManifestation"
  ) +
  geom_vline(xintercept = "2000", linetype = "dashed") +
  geom_vline(xintercept = "2010", linetype = "dashed") +
  geom_vline(xintercept = "2020", linetype = "dashed") +
  theme(
    legend.title = element_text(size = 16),
    title = element_text(size = 18),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.key.size = unit(0.8, "cm")
  )

# join in colours 
metadata_age_manifest <- left_join(metadata_age_manifest, serotype_colours, by = "serotype_clean") |>
  rename(
    serotype_colour = In_silico_serotype__colour.y 
  )

continent_year_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "continent year plot.svg", fsep = "/")

ggsave(continent_year_plot, filename = continent_year_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")

#----------------------------------------- panels ----------------------------------------

# vancomycin not included as R isolates were not found
# no CHL resistance found for GPSC1
# join serotype/gpsc table for CLI, COT, ERY, PEN, RIF, TET
antibiotic <- c("WGS_PEN_SIR_Meningitis", "WGS_AMO_SIR", "WGS_TAX_SIR_Meningitis", "WGS_CFT_SIR_Meningitis",
                "WGS_MER_SIR", "WGS_CFX_SIR", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR",
                "WGS_CHL_SIR", "WGS_RIF_SIR")


Continent <- c("AFRICA", "NORTH AMERICA", "LATIN AMERICA", "EUROPE", "ASIA", "OCEANIA")

world <- map_data("world")
world <- subset(world, region != "Antarctica")

setwd("C:/Users/Max Choong/OneDrive/Desktop/GPS-AMR/figures")


metadata_age_manifest[33:70] <- lapply(metadata_age_manifest[33:70], gsub, pattern = "FLAG", replacement = "S", fixed = TRUE)

metadata_age_manifest <- left_join(metadata_age_manifest, serotype_colours, by = "serotype_clean")


serotype_colors_vector <- setNames(metadata_age_manifest$In_silico_serotype__colour, metadata_age_manifest$serotype_clean)

metadata_age_manifest <- metadata_age_manifest |>
  filter(Vaccine_period_clean %in% c("POSTPCV10", "POSTPCV13") & Manifest_category %in% "IPD")

### request from steph
### change for full report


### panel and forest plot for loops
for(i in antibiotic){
  
  #filter out isolates that failed resistance prediction
  tidy_dataset <- metadata_age_manifest |> filter(!!((sym(i))) %in% c("S", "I", "R"))
  
  ### create new directory for each looped variable
  newdir <- paste0(i, " plots")

  dir.create(newdir)
  
  ### data tidying, relevel factors
  tidy_dataset$Age_category <- factor(tidy_dataset$Age_category, levels = c("< 5", ">=5 - <=65", "> 65"))
  tidy_dataset$Manifest_category <- factor(tidy_dataset$Manifest_category, levels = c("CARRIAGE", "DIS_CAR", "IPD", "NONINVASIVE_DISEASE", "UNKNOWN"))
  tidy_dataset[[i]] <- factor(tidy_dataset[[i]], levels = c("S", "I", "R"))
  
  tidy_dataset$Manifest_category_tidy <- ifelse(tidy_dataset$Manifest_category %in% "CARRIAGE", "Carriage",
                                                ifelse(tidy_dataset$Manifest_category %in% "DIS_CAR", "Diseased carriage",
                                                       ifelse(tidy_dataset$Manifest_category %in% "IPD", "IPD",
                                                              ifelse(tidy_dataset$Manifest_category %in% "NONINVASIVE_DISEASE", "Non-invasive disease",
                                                                     ifelse(tidy_dataset$Manifest_category %in% "UNKNOWN", "Unknown", "Unknown")))))
  
  ### prep for age group bar plot
  age_tally <- tidy_dataset |>
    filter(!is.na(Age_category)) |>
    filter(!is.na(Manifest_category_tidy)) |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(Age_category, !!((sym(i))), Manifest_category_tidy) |>
    tally()
  
  total_tally <- tidy_dataset |>
    filter(!is.na(Age_category)) |>
    filter(!is.na(Manifest_category_tidy)) |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(Age_category) |>
    tally()
  
  join <- left_join(age_tally, total_tally, by = "Age_category")
  
  ### define path to save csv
  age_manifest_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/age manifest table.csv", sep = ""), fsep = "/")
  
  ### csv file for number of resistant isolates in each age group for each disease manifestation
  write_csv(join, file = age_manifest_path)
  
  ### calculate percentage
  age_tally_manifest <- join |>
    summarise(
      res_manifest = ifelse(!.data[[i]] == "R", "S", as.character(Manifest_category_tidy)),
      percentage = n.x/n.y,
      .groups = "keep"
    ) 
  
  ### manual relevelling of factors
  age_tally_manifest$res_manifest <- factor(age_tally_manifest$res_manifest, levels = c("S", "Carriage", "Diseased carriage", "IPD", "Non-invasive disease", "Unknown"))
  age_tally_manifest$Age_category <- factor(age_tally_manifest$Age_category, levels = c("< 5", ">=5 - <=65", "> 65"))
  
  ### manifestation in resistant isolates within each age group
  age_manifest <- age_tally_manifest |>
    filter(!res_manifest %in% "S") |>
    ggplot(aes(x = Age_category, y = percentage, fill = res_manifest)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#BA55D3", "#2ECC71", "#00BFFF", "#FF8C00", "gray")) +
    scale_y_continuous(labels = scales::percent) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      plot.margin = margin(0, 0, 0, 5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    ) +
    labs(
      fill = "Clinical \nmanifestation",
      y = "Percentage",
      x = "Age category"
    ) 
  
  ### scatterpie prep
  
  grouping <- tidy_dataset |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(!!((sym(i)))) |>
    tally() |>
    pivot_wider(names_from = !!((sym(i))), values_from = n)
  
  if("I" %in% names(grouping)){
    
    
    pivot_SIR <- metadata_age_manifest |>
      filter(!is.na(Continent)) |>
      filter(!!((sym(i))) %in% c("S", "I", "R")) |>
      group_by(Continent, !!((sym(i)))) |>
      tally() |>
      pivot_wider(names_from = !!((sym(i))), values_from = n)
    
    pivot_SIR[is.na(pivot_SIR)] = 0
    
    pivot_plot <- pivot_SIR |>
      mutate(
        Longitude = as.numeric(case_when(
          Continent == "AFRICA" ~ 21.890909,
          Continent == "NORTH AMERICA" ~ -100.538118,
          Continent == "LATIN AMERICA" ~ -53.604523,
          Continent == "EUROPE" ~ 15.738565,
          Continent == "ASIA" ~ 110.308879,
          Continent == "OCEANIA" ~ 150.562787
        )),
        Latitude = as.numeric(case_when(
          Continent == "AFRICA" ~ 6.208920,
          Continent == "NORTH AMERICA" ~ 49.249542,
          Continent == "LATIN AMERICA" ~ -12.581874,
          Continent == "EUROPE" ~ 50.890969,
          Continent == "ASIA" ~ 41.691447,
          Continent == "OCEANIA" ~ -28.553076
        )))
    
    pivot_plot$Total <- rowSums(pivot_plot[, !colnames(pivot_plot) %in% c("Continent", "Longitude", "Latitude")], na.rm = TRUE)
  } else
  {
    
    
    pivot_SIR <- metadata_age_manifest |>
      filter(!is.na(Continent)) |>
      filter(!!((sym(i))) %in% c("S", "R")) |>
      group_by(Continent, !!((sym(i)))) |>
      tally() |>
      pivot_wider(names_from = !!((sym(i))), values_from = n)
    
    pivot_SIR[is.na(pivot_SIR)] = 0
    
    pivot_plot <- pivot_SIR |>
      mutate(
        Longitude = as.numeric(case_when(
          Continent == "AFRICA" ~ 21.890909,
          Continent == "NORTH AMERICA" ~ -100.538118,
          Continent == "LATIN AMERICA" ~ -53.604523,
          Continent == "EUROPE" ~ 15.738565,
          Continent == "ASIA" ~ 110.308879,
          Continent == "OCEANIA" ~ 150.562787
        )),
        Latitude = as.numeric(case_when(
          Continent == "AFRICA" ~ 6.208920,
          Continent == "NORTH AMERICA" ~ 49.249542,
          Continent == "LATIN AMERICA" ~ -12.581874,
          Continent == "EUROPE" ~ 50.890969,
          Continent == "ASIA" ~ 41.691447,
          Continent == "OCEANIA" ~ -28.553076
        )))
    
    pivot_plot$Total <- rowSums(pivot_plot[, !colnames(pivot_plot) %in% c("Continent", "Longitude", "Latitude")], na.rm = TRUE)
  }
  
  if(!"I" %in% names(pivot_plot)){
    plot <- ggplot() +
      geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
      coord_fixed(1.3) +
      geom_scatterpie(data = pivot_plot, aes(x = Longitude, y = Latitude, r = 13), cols = c("S", "R"), color = NA) + 
      coord_equal() +
      scale_fill_manual(values = c("#0069EC", "#FF2722")) +
      labs(
        x = "Longitude",
        y = "Latitude",
        fill = "Resistance"
      ) +
      theme(
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position.inside = c(.1,.35),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(-20, 0, 0, 0)
      ) 
    
  } else
    plot <- ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
    coord_fixed(1.3) +
    geom_scatterpie(data = pivot_plot, aes(x = Longitude, y = Latitude, r = 13), cols = c("S", "I", "R"), color = NA) + 
    coord_equal() +
    scale_fill_manual(values = c("#0069EC", "#F797B1", "#FF2722")) +
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Resistance"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.4, "cm"),
      legend.position.inside = c(.1,.35),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(-20, 0, 0, 0)
    ) 
  
  SIR_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/SIR table.csv", sep = ""), fsep = "/")
  write_csv(pivot_plot, file = SIR_path)
  
  map_manifest <- ggarrange(plot, age_manifest, widths = c(1.75, 1.25), heights = c(1, 1), labels = c("(a)", "(b)"), label.y = 1.01, font.label=list(color="black",size=8), hjust = c(0,0,0)) +
    theme(
      plot.margin = margin(10, 0, 0, 0)
    )
  metadata_age_manifest[[i]] <- factor(metadata_age_manifest[[i]], levels = c("S", "I", "R"))
  ### year plot
  if(!"I" %in% names(pivot_plot)){
    year <- metadata_age_manifest |>
      filter(!!((sym(i))) %in% c("S", "I", "R")) |>
      filter(!Year %in% "_" & !is.na(Year)) |>
      ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
      geom_bar() +
      scale_fill_manual(values = c("#0069EC", "#FF2722")) +
      theme(
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 5)
      ) +
      labs(
        y = "Number of isolates",
        x = "Year",
        fill = "Resistance"
      ) +
      geom_vline(xintercept = "2000", linetype = "dashed") +
      geom_vline(xintercept = "2010", linetype = "dashed") +
      geom_vline(xintercept = "2020", linetype = "dashed") 
    
  } else
    year <- metadata_age_manifest |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    filter(!Year %in% "_" & !is.na(Year)) |>
    ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
    geom_bar() +
    scale_fill_manual(values = c("#0069EC", "#F797B1", "#FF2722")) +
    theme(
      axis.text = element_text(size = 7),
      legend.position = "none",
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    ) +
    labs(
      y = "Number of isolates",
      x = "Year",
      fill = "Resistance"
    ) +
    geom_vline(xintercept = "2000", linetype = "dashed") +
    geom_vline(xintercept = "2010", linetype = "dashed") +
    geom_vline(xintercept = "2020", linetype = "dashed") 
  
  ### serotype GPSC plot
  
  metadata_age_manifest$GPSC <- sub(".*?;", "", metadata_age_manifest$GPSC)
  
  n50 <- metadata_age_manifest |>
    filter(!!((sym(i))) %in% "R" & !is.na(GPSC)) |>
    group_by(GPSC) |>
    tally() 
  
  tally_50 <- n50[order(n50$n, decreasing = TRUE),] |>
    slice_head(n = 30)
  
  GPSC_list <- unique(tally_50$GPSC) 
  
  tidy_modified <- metadata_age_manifest |>
    mutate(
      gpsc_tidy = ifelse(GPSC %in% GPSC_list, metadata_age_manifest$GPSC, "others")
    )
  
  serotype_colours_vector <- setNames(tidy_modified$serotype_colour, tidy_modified$serotype_clean)
  
  GPSC_serotype_1 <- tidy_modified |>
    filter(!!((sym(i))) %in% "R" & gpsc_tidy %in% "1") |>
    ggplot(aes(x = as.factor(gpsc_tidy), fill = serotype_clean)) +
    geom_bar() +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colours_vector) +
    labs(
      x = "GPSC",
      y = "Number of isolates"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(15, 0, 0, 0)
    )
  
  GPSC_serotype_others <- tidy_modified |>
    filter(!!((sym(i))) %in% "R" & !gpsc_tidy %in% "1") |>
    ggplot(aes(x = as.factor(gpsc_tidy), fill = serotype_clean)) +
    geom_bar() +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colours_vector) +
    labs(
      x = "GPSC",
      y = "Number of isolates"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(15, 0, 0, 0)
    )
  
  
  GPSC_serotype_plot <- ggarrange(GPSC_serotype_1, GPSC_serotype_others, widths = c(1, 5)) 
  
  GPSC_serotype_year <- ggarrange(year, GPSC_serotype_plot, ncol = 1, labels = c("(c)", "(d)"), font.label=list(color="black",size=8), hjust = c(0,0,0))
  
  full_plot <- ggarrange(map_manifest, GPSC_serotype_year, ncol = 1, heights = c(1.5, 2)) +
    theme(
      plot.margin = margin(0.5,0.1,0.5,0.1, "cm")
    )
  
  
  full_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/full plot.svg", sep = ""), fsep = "/")
  
  ggsave(full_plot, filename = full_plot_path, dpi = 320, width = 26.58, height = 16.91, units = "cm")

  
  # forest plots
  
  newdir <- paste0(i, " forest plot")
  
  dir.create(newdir)
  
  napath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/na forest.svg", sep = ""), fsep = "/")
  
  complete <- metadata_age_manifest |> filter(!!((sym(i))) %in% c("S", "I", "R"))
  

  
  #forest plot for resistance in north american countries
  NorthAmerica <- complete |>
    filter(Continent %in% "NORTH AMERICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  if(!"R" %in% names(NorthAmerica)) {NorthAmerica$R <- 0
  NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
  NorthAmerica[is.na(NorthAmerica)] <- 0
  nameta <- metaprop(R, total, data = NorthAmerica)
  }
  
  else NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
  NorthAmerica[is.na(NorthAmerica)] <- 0
  nameta <- metaprop(R, total, data = NorthAmerica)
  
  svg(file = napath, width = 12.5, height = 3)
  forest(nameta, comb.fixed = nameta$comb.fixed, comb.random = nameta$comb.random, 
         overall = TRUE, studlab = NorthAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  
  #latin america
  LatinAmerica <- complete |>
    filter(Continent %in% "LATIN AMERICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  lapath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/la forest.svg", sep = ""), fsep = "/")
  
  if(!"R" %in% names(LatinAmerica)) {LatinAmerica$R <- 0
  LatinAmerica$total <- rowSums(LatinAmerica[, !colnames(LatinAmerica) %in% "Country"], na.rm = TRUE)
  LatinAmerica[is.na(LatinAmerica)] <- 0
  lameta <- metaprop(R, total, data = LatinAmerica)
  }
  
  else LatinAmerica$total <- rowSums(LatinAmerica[, !colnames(LatinAmerica) %in% "Country"], na.rm = TRUE)
  LatinAmerica[is.na(LatinAmerica)] <- 0
  lameta <- metaprop(R, total, data = LatinAmerica)
  
  svg(file = lapath, width = 12.5, height = 3)
  forest(lameta, comb.fixed = lameta$comb.fixed, comb.random = lameta$comb.random, 
         overall = TRUE, studlab = LatinAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #europe
  Europe <- complete |>
    filter(Continent %in% "EUROPE") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  eupath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/eu forest.svg", sep = ""), fsep = "/")
  
  if(!"R" %in% names(Europe)){Europe$R <- 0
  Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
  Europe[is.na(Europe)] <- 0
  eumeta <- metaprop(R, total, data = Europe)
  }
  
  else Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
  Europe[is.na(Europe)] <- 0
  eumeta <- metaprop(R, total, data = Europe)
  
  svg(file = eupath, width = 12.5, height = 6)
  forest(eumeta, comb.fixed = eumeta$comb.fixed, comb.random = eumeta$comb.random, 
         overall = TRUE, studlab = Europe$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #asia
  Asia <- complete |>
    filter(Continent %in% "ASIA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  aspath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/as forest.svg", sep = ""), fsep = "/")
  
  if(!"R" %in% names(Asia)){Asia$R <- 0
  Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
  Asia[is.na(Asia)] <- 0
  asmeta <- metaprop(R, total, data = Asia)
  }
  
  else Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
  Asia[is.na(Asia)] <- 0
  asmeta <- metaprop(R, total, data = Asia)
  
  svg(file = aspath, width = 12.5, height = 6)
  forest(asmeta, comb.fixed = asmeta$comb.fixed, comb.random = asmeta$comb.random, 
         overall = TRUE, studlab = Asia$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #africa
  Africa <- complete |>
    filter(Continent %in% "AFRICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  afpath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/af forest.svg", sep = ""), fsep = "/")
  
  if(!"R" %in% names(Africa)){Africa$R <- 0
  Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
  Africa[is.na(Africa)] <- 0
  afmeta <- metaprop(R, total, data = Africa)
  }
  
  else Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
  Africa[is.na(Africa)] <- 0
  afmeta <- metaprop(R, total, data = Africa)
  
  svg(file = afpath, width = 12.5, height = 6)
  forest(afmeta, comb.fixed = afmeta$comb.fixed, comb.random = afmeta$comb.random, 
         overall = TRUE, studlab = Africa$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #oceania
  Oceania <- complete |>
    filter(Continent %in% "OCEANIA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  
  ocepath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/oce forest.svg", sep = ""), fsep = "/")
  
  if(!"R" %in% names(Oceania)){Oceania$R <- 0
  Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
  Oceania[is.na(Oceania)] <- 0
  ocemeta <- metaprop(R, total, data = Oceania) 
  }
  
  else Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
  Oceania[is.na(Oceania)] <- 0
  ocemeta <- metaprop(R, total, data = Oceania)
  
  svg(file = ocepath, width = 12.5, height = 2.5)
  forest(ocemeta, comb.fixed = ocemeta$comb.fixed, comb.random = ocemeta$comb.random, 
         overall = TRUE, studlab = Oceania$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  
  ### overview table of proportion resistance in each continent
  
  df_name <- paste(i, "_proportion", sep = "")
  
  # back conversion of proportion values and CI for each continent
  random_est_na <- c(nameta$TE.random, nameta$lower.random, nameta$upper.random)
  
  random_est_na_tr <- as.data.frame(t(backtransf(random_est_na, sm = "PLOGIT")))
  
  random_est_eu <- c(eumeta$TE.random, eumeta$lower.random, eumeta$upper.random)
  
  random_est_eu_tr <- as.data.frame(t(backtransf(random_est_eu, sm = "PLOGIT")))
  
  random_est_af <- c(afmeta$TE.random, afmeta$lower.random, afmeta$upper.random)
  
  random_est_af_tr <- as.data.frame(t(backtransf(random_est_af, sm = "PLOGIT")))
  
  random_est_as <- c(asmeta$TE.random, asmeta$lower.random, asmeta$upper.random)
  
  random_est_as_tr <- as.data.frame(t(backtransf(random_est_as, sm = "PLOGIT")))
  
  random_est_oce <- c(ocemeta$TE.random, ocemeta$lower.random, ocemeta$upper.random)
  
  random_est_oce_tr <- as.data.frame(t(backtransf(random_est_oce, sm = "PLOGIT")))
  
  random_est_la <- c(lameta$TE.random, lameta$lower.random, lameta$upper.random)
  
  random_est_la_tr <- as.data.frame(t(backtransf(random_est_la, sm = "PLOGIT")))
  
  
  # bind proportions for each continent into a singular dataframe
  proportion_table <- rbind(random_est_na_tr, random_est_eu_tr, random_est_af_tr,
                            random_est_la_tr, random_est_as_tr, random_est_oce_tr)
  
  # create column for continent
  proportion_table <- proportion_table |>
    mutate(
      Continent = c("North America", "Europe", "Africa",
                    "Latin America", "Asia", "Oceania")) |>
    relocate(Continent) |>
    head()
  
  # rename columns
  colnames(proportion_table) <- c("Continent", "Proportion", "Lower", "Upper")
  
  # save proportion table for each antibiotic into their respective folder
  proportion_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/proportion.csv", sep = ""), fsep = "/")
  write_csv(proportion_table, file = proportion_path)
  
  # round to 2 significant figures
  proportion_table <- proportion_table |>
    mutate(across(c("Proportion", "Lower", "Upper"), round, 2))
  
  # paste proportion and CIs togther
  proportion_table$proportion_ci <- paste0(proportion_table$Proportion, " (", proportion_table$Lower, ", ", proportion_table$Upper, ")")
  
  # select for continent and proportion/CI columns
  proportion_table <- proportion_table |>
    select(c(Continent, proportion_ci))
  
  # rename proportion/CI column based on looped variable
  colnames(proportion_table)[which(names(proportion_table) == "proportion_ci")] <- df_name
  
  # rename proportion table to looped variable
  assign(df_name, proportion_table)

} 

### compile random effects model proportion table
# create list of proportion table dataframes for each antibiotic
df_list <- list(WGS_AMO_SIR_proportion, WGS_PEN_SIR_Meningitis_proportion, WGS_CFX_SIR_proportion,
                                       WGS_TAX_SIR_Meningitis_proportion, WGS_CFT_SIR_Meningitis_proportion, WGS_MER_SIR_proportion,
                                       WGS_CHL_SIR_proportion, WGS_CLI_SIR_proportion, WGS_ERY_SIR_proportion,
                                       WGS_TET_SIR_proportion, WGS_COT_SIR_proportion, WGS_RIF_SIR_proportion)

# join all tables together
complete_proportion_table <- df_list |>
  reduce(left_join, by = "Continent")

colnames(complete_proportion_table) <- c("Continent", "Amoxicillin", "Penicillin (meningitis)", "Cefuroxime", "Cefotaxime (meningitis)", 
                                         "Ceftriaxone (meningitis)", "Meropenem", "Chloramphenicol", "Clindamycin", "Erythromycin", 
                                         "Tetracycline", "Co-trimoxazole", "Rifampin")

# save proportion table for each antibiotic into their respective folder
complete_proportion_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summarised data", "/complete proportion post pcv10 IPD.csv", fsep = "/")

write_csv(complete_proportion_table, file = complete_proportion_path)





#--------------------------- gpsc, serotypes and resistance of 10% most recent IPD samples -----------------------------

metadata_age_manifest |>
  filter(Manifest_category %in% "IPD" & !is.na(Year)) |>
  group_by(Continent) |>
  arrange(desc(Year)) |>
  slice_head(prop = 0.1) |>
  tally() |>
  view()



ipd_continent <- metadata_age_manifest |>
  filter(Manifest_category %in% "IPD" & !is.na(Year)) |>
  group_by(Continent) |>
  arrange(desc(Year)) |>
  slice_head(prop = 0.1)

ipd_p10 <- metadata_age_manifest |>
  filter(Manifest_category %in% "IPD" & !is.na(Year)) |>
  group_by(Continent) |>
  arrange(desc(Year)) |>
  slice_head(prop = 0.1)

# extract gpscs n > 4 for each continent
AF <- ipd_continent |>
  filter(Continent %in% "AFRICA") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4)|>
  arrange(desc(n)) |>
  slice_head(n = 20)

AF_gpsc <- AF$GPSC

AS <- ipd_continent |>
  filter(Continent %in% "ASIA") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4)

AS_gpsc <- AS$GPSC

NorthA <- ipd_continent |>
  filter(Continent %in% "NORTH AMERICA") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4)

NorthA_gpsc <- NorthA$GPSC

LA <- ipd_continent |>
  filter(Continent %in% "LATIN AMERICA") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4)

LA_gpsc <- LA$GPSC

OCE <- ipd_continent |>
  filter(Continent %in% "OCEANIA") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4)

OCE_gpsc <- OCE$GPSC

EU <- ipd_continent |>
  filter(Continent %in% "EUROPE") |>
  group_by(GPSC) |>
  tally() |>
  filter(n > 4) |>
  arrange(desc(n)) |>
  slice_head(n = 20)

EU_gpsc <- EU$GPSC


# classify gpscs with < 5 as others
AF_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% AF_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  ) |>
  filter(Continent %in% "AFRICA")

AS_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% AS_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  )|>
  filter(Continent %in% "ASIA")

NorthA_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% NorthA_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  )|>
  filter(Continent %in% "NORTH AMERICA")

LA_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% LA_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  )|>
  filter(Continent %in% "LATIN AMERICA")

OCE_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% OCE_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  )|>
  filter(Continent %in% "OCEANIA")

EU_complete <- ipd_p10 |>
  mutate(
    gpsc_clean = ifelse(GPSC %in% EU_gpsc & !GPSC %in% "unknown", GPSC, "Others")
  ) |>
  filter(Continent %in% "EUROPE")


# join all continent tables

complete_IPD <- rbind(AF_complete, NorthA_complete, LA_complete, EU_complete, AS_complete, OCE_complete)

recent_ipd_plot <- complete_IPD |>
  ggplot(aes(x = gpsc_clean, fill = serotype_clean)) +
  geom_bar() +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 2)) +
  facet_wrap(~Continent, scale = "free") +
  labs(
    title = "GPSC and serotype distribution of 10% most recent IPD cases in each continent",
    y = "Number of isolates",
    x = "GPSC"
  ) +
  theme(
    axis.title = element_text(size = 12),
    title = element_text(size = 16)
  )

recent_ipd_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "recent ipd plot.svg", fsep = "/")

ggsave(recent_ipd_plot, filename = recent_ipd_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")


##### think about other information that can be shown

write_csv(metadata_age_manifest, file = "metadata.csv")
#--------------------------- change in serotypes over time -----------------------------


# geom tile
# empty dataframe 


antibiotic <- c("WGS_PEN_SIR_Meningitis", "WGS_AMO_SIR", "WGS_TAX_SIR_Meningitis", "WGS_CFT_SIR_Meningitis",
                "WGS_MER_SIR", "WGS_CFX_SIR", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR",
                "WGS_CHL_SIR", "WGS_RIF_SIR")

continent <- c("NORTH AMERICA", "LATIN AMERICA", "ASIA", "OCEANIA")

serotype_colors_vector <- setNames(metadata_age_manifest$In_silico_serotype__colour.y, metadata_age_manifest$serotype_clean)


for(i in continent){

  
  forloop_output <- data.frame()
  csv_output = list()
  
  
  as_gpsc <- metadata_age_manifest |>
    group_by(GPSC) |>
    filter(Continent %in% i & Manifest_category %in% "IPD") |>
    tally() |>
    filter(n > 9) 
  
  as_gpsc_list <- as_gpsc$GPSC
  
  title_continent <- str_to_title(i)
  
  gpsc_serotype_tally <- metadata_age_manifest |>
    filter(Continent %in% i) |>
    mutate(
      GPSC_tidy_pcv = ifelse(GPSC %in% as_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
    ) |>
    group_by(Vaccine_period_clean, GPSC_tidy_pcv, serotype_clean) |>
    tally() 
  
  # path to write out file
  gpsc_reference_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", paste(i, "gpsc reference.csv", sep = " "), fsep = "/")
  
  # write file
  write_csv(gpsc_serotype_tally, file = gpsc_reference_path)
  
  
  # annotate vaccine types
  gpsc_serotype_tally$vaccine <- ifelse(gpsc_serotype_tally$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")
  
  
  continent_dynamics <- gpsc_serotype_tally |>
    ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv)) +
    geom_col(position = position_dodge(preserve = "single"), colour = "black") +
    facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
    scale_y_continuous(expand = c(0,0)) +
    scale_colour_manual(values = "black") +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colors_vector,
                      guide = guide_legend(ncol = 3)) +
    theme(title = element_text(size = 20),
          axis.title=element_text(size=18,face="bold"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x =element_text(size=10,color="black", angle = 90, hjust=1, vjust=0.2),
          axis.text.y = element_text(size = 12),
          strip.placement='outside',
          strip.background.x=element_blank(),
          strip.text=element_text(size=12,color="black",face="bold"),
          panel.spacing.x=unit(1,"pt"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.8, "cm"),
          plot.margin = margin(15, 0, 0, 0)) +
    labs(
      y = "Number of isolates",
      title = paste("Dynamics of GPSCs causing invasive disease in", title_continent, sep = " "),
      x = "GPSC"
    )
  
  dynamics_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", paste(i, " dynamics plot.svg", sep = " "), fsep = "/")
  
  ggsave(continent_dynamics, filename = dynamics_plot_path, dpi = 320, width = 60, height = 25.91, units = "cm")
  
  
  
  for (t in antibiotic) {
    
    antibiotic_tile <- metadata_age_manifest |>
      mutate(
        GPSC_tidy_pcv = ifelse(GPSC %in% as_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
      ) |>
      group_by(GPSC_tidy_pcv, !!sym((t))) |>
      tally() |>
      pivot_wider(names_from = !!sym((t)), values_from = n)
    
    antibiotic_tile[is.na(antibiotic_tile)] <- 0
    
    antibiotic_tile$total <- rowSums(antibiotic_tile[, !colnames(antibiotic_tile) %in% "GPSC_tidy_pcv"], na.rm = TRUE)
    
    
    row_for_dataframe <- antibiotic_tile |>
      mutate(
        antibiotic = t,
        proportion = R/total
      ) |>
      select(GPSC_tidy_pcv, proportion, antibiotic)
    
    row_for_dataframe[is.na(row_for_dataframe)] <- 0
    
    dataframe_row <- row_for_dataframe |>
      pivot_wider(names_from = GPSC_tidy_pcv, values_from = proportion)
    
    csv_output[[t]] <- dataframe_row
    
    forloop_output <- rbind(row_for_dataframe, forloop_output)
    
    
  }

  forloop_output <- left_join(forloop_output, antibiotic_tidying, by = "antibiotic") 
  
  forloop_output$antibiotic_tidy <- factor(forloop_output$antibiotic_tidy, levels = c("Rifampin", "Chloramphenicol", "Tetracycline",
                                                                                      "Co-trimoxazole", "Clindamycin", "Erythromycin",
                                                                                      "Meropenem", "Cefotaxime (meningitis)", "Ceftriaxone (meningitis)",
                                                                                      "Cefuroxime", "Amoxicillin", "Penicillin (meningitis)"))
  
  ### write out reference file for each continent
  resistance_reference <- forloop_output |>
    select(-antibiotic) |>
    pivot_wider(names_from = antibiotic_tidy, values_from = proportion)
  
  # path to write out file
  antibiotic_reference_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", paste(i, "reference.csv", sep = " "), fsep = "/")
  
  # write file
  write_csv(resistance_reference, file = antibiotic_reference_path)
  
  
  # tile plot
  tile <- forloop_output |>
    ggplot(aes(x = GPSC_tidy_pcv, y = antibiotic_tidy, fill = proportion)) +
    geom_tile(colour = "black") +
    scale_fill_gradientn(colours = c("white", "pink", "red"), values = c("0", "0.25", "0.5", "0.75", "1.00")) +
    coord_fixed() +
    theme(
      plot.margin = unit(c(-4,0,0,0), 'lines')
    ) +
    labs(
      y = "Antibiotic",
      fill = "Proportion",
      x = "GPSC"
    ) +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key.size = unit(0.8, "cm")
    )
  
  
  tile_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", paste(i, "tile plot.svg", sep = " "), fsep = "/")
  
  ggsave(tile, filename = tile_plot_path, dpi = 320, width = 65.58, height = 25.5, units = "cm")
  
  
}


gpsc_serotype_tally |>
  ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv)) +
  geom_col(position = position_dodge(preserve = "single"), colour = "black") +
  facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
  scale_y_continuous(expand = c(0,0)) +
  scale_colour_manual(values = "black") +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(title = element_text(size = 20),
        axis.title=element_text(size=18,face="bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(size=12,color="black", angle = 90, hjust=1, vjust=0.2),
        axis.text.y = element_text(size = 12),
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(size=12,color="black",face="bold"),
        panel.spacing.x=unit(1,"pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.8, "cm"),
        plot.margin = margin(15, 0, 0, 0)) +
  labs(
    y = "Number of isolates",
    title = paste("Dynamics of GPSCs causing invasive disease in", title_continent, sep = " "),
    x = "GPSC"
  )




# Africa
af_gpsc <- metadata_age_manifest |>
  group_by(GPSC) |>
  filter(Continent %in% "AFRICA" & Manifest_category %in% "IPD") |>
  tally() |>
  filter(n > 9) |>
  slice_head(n = 35)

af_gpsc_list <- af_gpsc$GPSC

af_gpsc_serotype_tally <- metadata_age_manifest |>
  filter(Continent %in% "AFRICA") |>
  mutate(
    GPSC_tidy_pcv = ifelse(GPSC %in% af_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
  ) |>
  group_by(Vaccine_period_clean, GPSC_tidy_pcv, serotype_clean) |>
  tally()

# path to write out file
af_gpsc_reference_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Africa gpsc reference.csv", fsep = "/")

# write file
write_csv(af_gpsc_serotype_tally, file = af_gpsc_reference_path)


# annotate vaccine types
gpsc_serotype_tally$vaccine <- ifelse(gpsc_serotype_tally$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")

africa_dynamics <- gpsc_serotype_tally |>
  ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv)) +
  geom_col(position = position_dodge(preserve = "single"), colour = "black") +
  facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(axis.title=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(size=8,color="black", angle = 90, hjust=1, vjust=0.2),
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(size=12,color="black",face="bold"),
        panel.spacing.x=unit(1,"pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 18)) +
  labs(
    y = "Number of isolates",
    title = "Dynamics of GPSCs causing invasive disease in Africa"
  )


africa_dynamics_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Africa dynamics plot.svg", fsep = "/")

ggsave(africa_dynamics, filename = africa_dynamics_plot_path, dpi = 320, width = 60, height = 25.91, units = "cm")


forloop_output <- data.frame()
csv_output = list()


for (t in antibiotic) {
  
  antibiotic_tile <- metadata_age_manifest |>
    filter(Continent %in% "AFRICA") |>
    mutate(
      GPSC_tidy_pcv = ifelse(GPSC %in% as_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
    ) |>
    group_by(GPSC_tidy_pcv, !!sym((t))) |>
    tally() |>
    pivot_wider(names_from = !!sym((t)), values_from = n)
  
  antibiotic_tile[is.na(antibiotic_tile)] <- 0
  
  antibiotic_tile$total <- rowSums(antibiotic_tile[, !colnames(antibiotic_tile) %in% "GPSC_tidy_pcv"], na.rm = TRUE)
  
  
  row_for_dataframe <- antibiotic_tile |>
    mutate(
      antibiotic = t,
      proportion = R/total
    ) |>
    select(GPSC_tidy_pcv, proportion, antibiotic)
  
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  
  dataframe_row <- row_for_dataframe |>
    pivot_wider(names_from = GPSC_tidy_pcv, values_from = proportion)
  
  csv_output[[t]] <- dataframe_row
  
  forloop_output <- rbind(row_for_dataframe, forloop_output)
  
  
}

forloop_output <- left_join(forloop_output, antibiotic_tidying, by = "antibiotic") 

forloop_output$antibiotic_tidy <- factor(forloop_output$antibiotic_tidy, levels = c("Rifampin", "Chloramphenicol", "Tetracycline",
                                                                                    "Co-trimoxazole", "Clindamycin", "Erythromycin",
                                                                                    "Meropenem", "Cefotaxime (meningitis)", "Ceftriaxone (meningitis)",
                                                                                    "Cefuroxime", "Amoxicillin", "Penicillin (meningitis)"))

### write out reference file for each continent
resistance_reference <- forloop_output |>
  select(-antibiotic) |>
  pivot_wider(names_from = antibiotic_tidy, values_from = proportion)

# path to write out file
antibiotic_reference_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Africa reference.csv", fsep = "/")

# write file
write_csv(resistance_reference, file = antibiotic_reference_path)


# tile plot
tile <- forloop_output |>
  ggplot(aes(x = GPSC_tidy_pcv, y = antibiotic_tidy, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = c("white", "pink", "red"), values = c("0", "0.25", "0.5", "0.75", "1.00")) +
  coord_fixed() +
  theme(
    plot.margin = unit(c(-4,0,0,0), 'lines')
  ) +
  labs(
    y = "Antibiotic",
    fill = "Proportion",
    x = "GPSC"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )


tile_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Africa tile plot.svg", fsep = "/")

ggsave(tile, filename = tile_plot_path, dpi = 320, width = 65.58, height = 25.5, units = "cm")









# Europe complete
eu_gpsc <- metadata_age_manifest |>
  group_by(GPSC) |>
  filter(Continent %in% "EUROPE" & Manifest_category %in% "IPD") |>
  tally() |>
  filter(n > 19) |>
  slice_head(n = 34)

eu_gpsc_list <- eu_gpsc$GPSC

# Europe tail
eu_gpsc <- metadata_age_manifest |>
  group_by(GPSC) |>
  filter(Continent %in% "EUROPE" & Manifest_category %in% "IPD") |>
  tally() |>
  filter(n > 19 & !GPSC %in% "3") |>
  slice_head(n = 33)

eu_gpsc_tail <- eu_gpsc$GPSC

gpsc_serotype_tally <- metadata_age_manifest |>
  filter(Continent %in% "EUROPE" & !GPSC %in% "3") |>
  mutate(
    GPSC_tidy_pcv = ifelse(GPSC %in% eu_gpsc_tail & !GPSC %in% "unknown", GPSC, "other")
  ) |>
  group_by(Vaccine_period_clean, GPSC_tidy_pcv, serotype_clean) |>
  tally()

# annotate vaccine types
gpsc_serotype_tally$vaccine <- ifelse(gpsc_serotype_tally$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")

europe_dynamics_tail <- gpsc_serotype_tally |>
  ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv)) +
  geom_col(position = position_dodge(preserve = "single"), colour = "black") +
  facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(axis.title=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(size=8,color="black", angle = 90, hjust=1, vjust=0.2),
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(size=12,color="black",face="bold"),
        panel.spacing.x=unit(1,"pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 18)) +
  labs(
    y = "Number of isolates",
    title = "Dynamics of GPSCs causing invasive disease in Europe"
  )




# Europe
eu_gpsc_head <- metadata_age_manifest |>
  group_by(GPSC) |>
  filter(Continent %in% "EUROPE" & Manifest_category %in% "IPD") |>
  tally() |>
  filter(GPSC %in% "3") 

eu_gpsc_head <- eu_gpsc_head$GPSC

eu_serotype_tally <- metadata_age_manifest |>
  filter(Continent %in% "EUROPE" & GPSC %in% "3") |>
  mutate(
    GPSC_tidy_pcv = ifelse(GPSC %in% eu_gpsc_head & !GPSC %in% "unknown", GPSC, "other")
  ) |>
  group_by(Vaccine_period_clean, GPSC_tidy_pcv, serotype_clean) |>
  tally()

# annotate vaccine types
eu_serotype_tally$vaccine <- ifelse(eu_serotype_tally$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")

europe_dynamics_head <- eu_serotype_tally |>
  ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 2)) +
  theme(axis.title=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(size=8,color="black", angle = 90, hjust=1, vjust=0.2),
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(size=12,color="black",face="bold"),
        panel.spacing.x=unit(1,"pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +
  labs(
    y = "Number of isolates",
    title = " "
  )

europe_dynamics <- ggarrange(europe_dynamics_head, europe_dynamics_tail, widths = c(1, 20)) 

europe_dynamics_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Europe dynamics plot.svg", fsep = "/")

ggsave(europe_dynamics, filename = europe_dynamics_plot_path, dpi = 320, width = 60, height = 25.91, units = "cm")


forloop_output <- data.frame()
csv_output = list()


for (t in antibiotic) {
  
  antibiotic_tile <- metadata_age_manifest |>
    filter(Continent %in% "EUROPE") |>
    mutate(
      GPSC_tidy_pcv = ifelse(GPSC %in% eu_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
    ) |>
    group_by(GPSC_tidy_pcv, !!sym((t))) |>
    tally() |>
    pivot_wider(names_from = !!sym((t)), values_from = n)
  
  antibiotic_tile[is.na(antibiotic_tile)] <- 0
  
  antibiotic_tile$total <- rowSums(antibiotic_tile[, !colnames(antibiotic_tile) %in% "GPSC_tidy_pcv"], na.rm = TRUE)
  
  
  row_for_dataframe <- antibiotic_tile |>
    mutate(
      antibiotic = t,
      proportion = R/total
    ) |>
    select(GPSC_tidy_pcv, proportion, antibiotic)
  
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  
  dataframe_row <- row_for_dataframe |>
    pivot_wider(names_from = GPSC_tidy_pcv, values_from = proportion)
  
  csv_output[[t]] <- dataframe_row
  
  forloop_output <- rbind(row_for_dataframe, forloop_output)
  
  
}

forloop_output <- left_join(forloop_output, antibiotic_tidying, by = "antibiotic") 

forloop_output$antibiotic_tidy <- factor(forloop_output$antibiotic_tidy, levels = c("Rifampin", "Chloramphenicol", "Tetracycline",
                                                                                    "Co-trimoxazole", "Clindamycin", "Erythromycin",
                                                                                    "Meropenem", "Cefotaxime (meningitis)", "Ceftriaxone (meningitis)",
                                                                                    "Cefuroxime", "Amoxicillin", "Penicillin (meningitis)"))


### write out reference file for each continent
resistance_reference <- forloop_output |>
  select(-antibiotic) |>
  pivot_wider(names_from = antibiotic_tidy, values_from = proportion)

# path to write out file
antibiotic_reference_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Europe reference.csv", fsep = "/")

# write file
write_csv(resistance_reference, file = antibiotic_reference_path)


# split into gpsc 3 and others
forloop_output_head <- forloop_output |>
  filter(GPSC_tidy_pcv %in% "3")

# tile plot
tile_head <- forloop_output_head |>
  ggplot(aes(x = GPSC_tidy_pcv, y = antibiotic_tidy, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = c("white", "pink", "red"), values = c("0", "0.25", "0.5", "0.75", "1.00"), limits = c(0, 1)) +
  coord_fixed() +
  labs(
    y = "Antibiotic",
    fill = "Proportion",
    x = "GPSC"
  ) +
  theme(axis.title.y = element_text(size=12,face="bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )

# split into gpsc 3 and others
forloop_output_tail <- forloop_output |>
  filter(!GPSC_tidy_pcv %in% "3")

# tile plot
tile_tail <- forloop_output_tail |>
  ggplot(aes(x = GPSC_tidy_pcv, y = antibiotic_tidy, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = c("white", "pink", "red"), values = c("0", "0.25", "0.5", "0.75", "1.00"), limits = c(0, 1)) +
  coord_fixed() +
  theme(
    plot.margin = unit(c(-4,0,0,0), 'lines')
  ) +
  labs(
    y = "Antibiotic",
    fill = "Proportion",
    x = "GPSC"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )


tile_head_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Europe tile head plot.svg", fsep = "/")

ggsave(tile_head, filename = tile_head_path, dpi = 320, width = 65.58, height = 25.5, units = "cm")

tile_tail_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", "Europe tile tail plot.svg", fsep = "/")

ggsave(tile_tail, filename = tile_tail_path, dpi = 320, width = 65.58, height = 25.5, units = "cm")



metadata_age_manifest |>
  group_by(Country, Vaccine_period_clean) |>
  tally() |>
  filter(n > 19) |>
  view()

### rotate bar chart horizontally, facet by gpsc like in pcv period plot

#------------------------------------ notes beyond this -------------------------------------


# with VT/NVT pattern
for(i in continent){
  
  as_gpsc <- metadata_age_manifest |>
    group_by(GPSC) |>
    filter(Continent %in% i & Manifest_category %in% "IPD") |>
    tally() |>
    filter(n > 9) 
  
  as_gpsc_list <- as_gpsc$GPSC
  
  title_continent <- str_to_title(i)
  
  gpsc_serotype_tally <- metadata_age_manifest |>
    filter(Continent %in% i) |>
    mutate(
      GPSC_tidy_pcv = ifelse(GPSC %in% as_gpsc_list & !GPSC %in% "unknown", GPSC, "other")
    ) |>
    group_by(Vaccine_period_clean, GPSC_tidy_pcv, serotype_clean) |>
    tally() 
  
  # annotate vaccine types
  gpsc_serotype_tally$vaccine <- ifelse(gpsc_serotype_tally$serotype_clean %in% c(pcv7_serotypes, pcv10_serotypes, pcv13_serotypes), "VT", "NVT")
  
  
  continent_dynamics <- gpsc_serotype_tally |>
    ggplot(aes(x = Vaccine_period_clean, y = n, fill = serotype_clean, group = GPSC_tidy_pcv, pattern = vaccine)) +
    geom_col_pattern(
      pattern_spacing = 0.01, 
      pattern_fill = "black", 
      pattern_size = 0.5, 
      pattern_linetype = 0.5,
      pattern_key_scale_factor = 0.5,
      pattern_orientation = "vertical", show.legend = TRUE,
      position = position_dodge(preserve = "single")) +
    facet_grid(cols = vars(GPSC_tidy_pcv), space="free_x", switch='x') +
    scale_y_continuous(expand = c(0,0)) +
    scale_colour_manual(values = "black") +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colors_vector,
                      guide = guide_legend(ncol = 2)) +
    theme(axis.title=element_text(size=12,face="bold"),
          axis.ticks.x = element_blank(),
          axis.text.x =element_text(size=8,color="black", angle = 90, hjust=1, vjust=0.2),
          strip.placement='outside',
          strip.background.x=element_blank(),
          strip.text=element_text(size=12,color="black",face="bold"),
          panel.spacing.x=unit(1,"pt"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.6, "cm"),
          plot.margin = margin(15, 0, 0, 0)) +
    labs(
      y = "Number of isolates",
      title = paste("Dynamics of GPSCs causing invasive disease in", title_continent, sep = " "),
      x = "GPSC"
    )
  
  dynamics_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", "dynamics plot", paste(i, "plot.svg", sep = " "), fsep = "/")
  
  ggsave(continent_dynamics, filename = dynamics_plot_path, dpi = 320, width = 40.58, height = 18.91, units = "cm")
  
}

