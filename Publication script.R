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

# read in file containing tidied serotypes
clean_serotype <- read_csv(file = "tidied_serotypes.csv")

# read in serotype colours
serotype_colours <- read_csv(file = "serotype_colours.csv")


#----------------------------------------- metadata cleanup -----------------------------------------

# fill in missing continent data
metadata$Continent <- ifelse(metadata$Country %in% "NIGERIA" & is.na(metadata$Continent), "AFRICA",
                             ifelse(metadata$Country %in% "THAILAND" & is.na(metadata$Continent), "ASIA", metadata$Continent))

# had one singular isolate from central america, grouped the isolate with latin america instead
metadata$Continent <- ifelse(metadata$Continent %in% "CENTRAL AMERICA", "LATIN AMERICA", metadata$Continent)

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

# set NA in GPSC to unknonwn
metadata$GPSC <- ifelse(is.na(metadata$GPSC), "unknown", metadata$GPSC)

# remove year after vaccine period
metadata$Vaccine_period_clean <- gsub("\\-.*", "", metadata$Vaccine_period)

# arrange PCV
metadata$Vaccine_period_clean <- factor(metadata$Vaccine_period_clean, 
                                            levels = c("PREPCV", "POSTPCV7", "POSTPCV10", "POSTPCV13"))

metadata[31:68] <- lapply(metadata[31:68], gsub, pattern = "FLAG", replacement = "S", fixed = TRUE)


# left join to match serotype with tidied serotype
metadata <- left_join(metadata, clean_serotype, by = "In_silico_serotype") |>
  relocate(serotype_clean, .after = In_silico_serotype)

# replace NAs in clean serotype column with in silico serotype
metadata$serotype_clean <- ifelse(is.na(metadata$serotype_clean), metadata$In_silico_serotype, metadata$serotype_clean)


# join in colours 
metadata <- left_join(metadata, serotype_colours, by = "serotype_clean") |>
  rename(
    serotype_colour = In_silico_serotype__colour.y 
  )

serotype_colors_vector <- setNames(metadata$In_silico_serotype__colour, metadata$serotype_clean)


#----------------------------------------- maps -----------------------------------------


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
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
  ) 

world_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "world plot.svg", fsep = "/")

ggsave(world_plot, filename = world_plot_path, dpi = 320, width =60.58, height = 28.91, units = "cm")


#----------------------------------------- vaccine reference table ----------------------------------------

# read in table of serotypes covered by each vaccine
vaccine_table <- read_csv(file = "vaccine_chart.csv")

# create a column with values as 1
vaccine_table$value <- 1

# reorder factors
vaccine_table$vaccine <- factor(vaccine_table$vaccine, levels = c("PCV31",
                                                                  "PCV25",
                                                                  "PCV24_GSK",
                                                                  "PCV21_SK_SANOFI",
                                                                  "PCV20", 
                                                                  "PCV15",
                                                                  "PCV13",
                                                                  "PCV10_PNEUMOSIL",
                                                                  "PCV10_GSK",
                                                                  "PCV7"))

# labels
vaccine_labels <- c("PCV31",
                    "PCV25",
                    "PCV24 (GSK)",
                    "PCV21 (SK/SANOFI)",
                    "PCV20", 
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
    axis.text = element_text(size = 24)
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



# tally table of serotypes in the 40 most abundant gpsc 
graph_dataset <- tidy |>
  group_by(GPSC_tidy, serotype_clean) |>
  tally()


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

#----------------------------------------- year/manifest plot ----------------------------------------

# isolates with >=15 or >5 will be labelled as NA 
# re-classify and tidying up age and clinical manifestation
metadata_age_manifest <- metadata |>
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
    x = "Year",
    y = "Number of isolates",
    fill = "Clinical \nManifestation",
    
  ) +
  geom_vline(xintercept = "2000", linetype = "dashed") +
  geom_vline(xintercept = "2010", linetype = "dashed") +
  geom_vline(xintercept = "2020", linetype = "dashed") +
  theme(
    title = element_text(size = 26),
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    strip.text.y = element_text(size = 18)
  )



continent_year_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "summary plots", "continent year plot.svg", fsep = "/")

ggsave(continent_year_plot, filename = continent_year_plot_path, dpi = 320, width =80.58, height = 36.91, units = "cm")

#----------------------------------------- panels ----------------------------------------

# vancomycin not included as R isolates were not found
# no CHL resistance found for GPSC1
# join serotype/gpsc table for CLI, COT, ERY, PEN, RIF, TET
antibiotic <- c("WGS_PEN_SIR_Meningitis", "WGS_AMO_SIR", "WGS_TAX_SIR_Meningitis", "WGS_CFT_SIR_Meningitis",
                "WGS_MER_SIR", "WGS_CFX_SIR", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR",
                "WGS_CHL_SIR", "WGS_RIF_SIR")

antibiotic_join <- c("WGS_PEN_SIR_Meningitis", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR", "WGS_CHL_SIR", "WGS_RIF_SIR")

antibiotic_diverse_serotype <- c("WGS_PEN_SIR_Meningitis", "WGS_ERY_SIR", "WGS_TET_SIR", "WGS_CLI_SIR", "WGS_COT_SIR")

Continent <- c("AFRICA", "NORTH AMERICA", "LATIN AMERICA", "EUROPE", "ASIA", "OCEANIA")

world <- map_data("world")
world <- subset(world, region != "Antarctica")

setwd("C:/Users/Max Choong/OneDrive/Desktop/GPS-AMR/figures")



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
    reframe(
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
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.6, "cm"),
      plot.margin = margin(0, 0, 0, 5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(size = 10)
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
    
    
    pivot_SIR <- tidy_dataset |>
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
    
    
    pivot_SIR <- tidy_dataset |>
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
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),
        legend.position.inside = c(.1,.35),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(-20, 0, 0, 10)
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
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.6, "cm"),
      legend.position.inside = c(.1,.35),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(-20, 0, 0, 10)
    ) 
  
  SIR_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/SIR table.csv", sep = ""), fsep = "/")
  write_csv(pivot_plot, file = SIR_path)
  
  map_manifest <- ggarrange(plot, age_manifest, widths = c(1.5, 1.25), heights = c(1, 1), labels = c("(a)", "(b)"), label.y = 1.01, font.label=list(color="black",size=12), hjust = c(0.5,0,0)) +
    theme(
      plot.margin = margin(10, 0, 10, 0)
    )
  tidy_dataset[[i]] <- factor(tidy_dataset[[i]], levels = c("S", "I", "R"))
 
  ### year plot
  if(!"I" %in% names(pivot_plot)){
    year <- tidy_dataset |>
      filter(!!((sym(i))) %in% c("S", "I", "R")) |>
      filter(!Year %in% "_" & !is.na(Year)) |>
      ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
      geom_bar() +
      scale_fill_manual(values = c("#0069EC", "#FF2722")) +
      theme(
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 8),
        legend.position = "none",
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1)
      ) +
      labs(
        y = "Number of\nisolates",
        x = "Year",
        fill = "Resistance"
      )
    
  } else
    year <- tidy_dataset |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    filter(!Year %in% "_" & !is.na(Year)) |>
    ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
    geom_bar() +
    scale_fill_manual(values = c("#0069EC", "#F797B1", "#FF2722")) +
    theme(
      axis.text = element_text(size = 10),
      legend.position = "none",
      axis.title = element_text(size = 11),
      axis.text.x = element_text(size = 10, angle = 45, vjust = 1)
    ) +
    labs(
      y = "Number of\nisolates",
      x = "Year",
      fill = "Resistance"
    ) 
  
  ### serotype GPSC plot
  
  tidy_dataset$GPSC <- sub(".*?;", "", tidy_dataset$GPSC)
  
  n50 <- tidy_dataset |>
    filter(!!((sym(i))) %in% "R" & !is.na(GPSC)) |>
    group_by(GPSC) |>
    tally() 
  
  tally_50 <- n50[order(n50$n, decreasing = TRUE),] |>
    slice_head(n = 25)
  
  GPSC_list <- sort(unique(as.numeric(tally_50$GPSC)))
  
  tidy_modified <- tidy_dataset |>
    mutate(
      gpsc_tidy = ifelse(GPSC %in% GPSC_list, tidy_dataset$GPSC, "others")
    ) |>
    filter(!gpsc_tidy %in% c("others", "unknown"))
  
  serotype_colours_vector <- setNames(tidy_modified$serotype_colour, tidy_modified$serotype_clean)

  if (!i %in% antibiotic_join) {
    GPSC_serotype_1 <- tidy_modified |>
      filter(!!((sym(i))) %in% "R" & gpsc_tidy %in% "1") |>
      ggplot(aes(x = factor(gpsc_tidy, levels = GPSC_list), fill = serotype_clean)) +
      geom_bar() +
      scale_fill_manual(name = "Serotype",
                        values = serotype_colours_vector,
                        guide = guide_legend(ncol = 1)) +
      labs(
        x = "GPSC",
        y = "Number of isolates"
      ) +
      theme(
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        plot.margin = margin(15, 0, 0, 0)
      )
    
    GPSC_serotype_others <- tidy_modified |>
      filter(!!((sym(i))) %in% "R" & !gpsc_tidy %in% "1") |>
      ggplot(aes(x = factor(gpsc_tidy, levels = GPSC_list), fill = serotype_clean)) +
      geom_bar() +
      scale_fill_manual(name = "Serotype",
                        values = serotype_colours_vector,
                        guide = guide_legend(ncol = 3)) +
      labs(
        x = "GPSC",
        y = "Number of isolates"
      ) +
      theme(
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(15, 0, 0, 0)
      )
    
    GPSC_serotype_plot <- ggarrange(GPSC_serotype_1, GPSC_serotype_others, widths = c(1.25, 5)) 
    
    GPSC_serotype_year <- ggarrange(year, GPSC_serotype_plot, ncol = 1, heights = c(1.5, 2), labels = c("(c)", "(d)"), font.label=list(color="black",size=12), hjust = c(0,0,0))
    
    full_plot <- ggarrange(map_manifest, GPSC_serotype_year, ncol = 1, heights = c(1, 2)) +
      theme(
        plot.margin = margin(0.5,0.3,0.5,0.3, "cm")
      )
    
    
  } else if (i %in% antibiotic_diverse_serotype) {
    GPSC_serotype_plot <- tidy_modified |>
      filter(!!((sym(i))) %in% "R") |>
      ggplot(aes(x = factor(gpsc_tidy, levels = GPSC_list), fill = serotype_clean)) +
      geom_bar() +
      scale_fill_manual(name = "Serotype",
                        values = serotype_colours_vector,
                        guide = guide_legend(ncol = 4)) +
      labs(
        x = "GPSC",
        y = "Number of isolates"
      ) +
      theme(
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(15, 0, 0, 0)
      )
    
    GPSC_serotype_year <- ggarrange(year, GPSC_serotype_plot, ncol = 1, heights = c(1.5, 2), labels = c("(c)", "(d)"), font.label=list(color="black",size=12), hjust = c(0,0,0))
    
    full_plot <- ggarrange(map_manifest, GPSC_serotype_year, ncol = 1, heights = c(1, 2)) +
      theme(
        plot.margin = margin(0.5,0.3,0.5,0.3, "cm")
      )
    
    
  } else {
    
    GPSC_serotype_plot <- tidy_modified |>
      filter(!!((sym(i))) %in% "R") |>
      ggplot(aes(x = factor(gpsc_tidy, levels = GPSC_list), fill = serotype_clean)) +
      geom_bar() +
      scale_fill_manual(name = "Serotype",
                        values = serotype_colours_vector,
                        guide = guide_legend(ncol = 3)) +
      labs(
        x = "GPSC",
        y = "Number of isolates"
      ) +
      theme(
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(15, 0, 0, 0)
      )
    
    GPSC_serotype_year <- ggarrange(year, GPSC_serotype_plot, ncol = 1, heights = c(1.5, 2), labels = c("(c)", "(d)"), font.label=list(color="black",size=12), hjust = c(0,0,0))
    
    full_plot <- ggarrange(map_manifest, GPSC_serotype_year, ncol = 1, heights = c(1, 2)) +
      theme(
        plot.margin = margin(0.5,0.3,0.5,0.3, "cm")
      )
    
  }
  
  
  
  full_plot_path <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/full plot.svg", sep = ""), fsep = "/")
  
  ggsave(full_plot, filename = full_plot_path, dpi = 320, width = 26.58, height = 16.91, units = "cm")

  
  # forest plots
  
  newdir <- paste0(i, " forest plot")
  
  dir.create(newdir)
  
  napath <- file.path("/Users", "Max Choong", "OneDrive", "Desktop", "GPS-AMR", "figures", paste(newdir, "/na forest.svg", sep = ""), fsep = "/")
  
  complete <- tidy_dataset |> filter(!!((sym(i))) %in% c("S", "I", "R"))
  

  
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
  
  svg(file = napath, width = 14.5, height = 3)
  forest(nameta, comb.fixed = nameta$comb.fixed, comb.random = nameta$comb.random, 
         overall = TRUE, studlab = NorthAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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
  
  svg(file = lapath, width = 14.5, height = 3)
  forest(lameta, comb.fixed = lameta$comb.fixed, comb.random = lameta$comb.random, 
         overall = TRUE, studlab = LatinAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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
  
  svg(file = eupath, width = 14.5, height = 6)
  forest(eumeta, comb.fixed = eumeta$comb.fixed, comb.random = eumeta$comb.random, 
         overall = TRUE, studlab = Europe$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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
  
  svg(file = aspath, width = 14.5, height = 6)
  forest(asmeta, comb.fixed = asmeta$comb.fixed, comb.random = asmeta$comb.random, 
         overall = TRUE, studlab = Asia$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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
  
  svg(file = afpath, width = 14.5, height = 6)
  forest(afmeta, comb.fixed = afmeta$comb.fixed, comb.random = afmeta$comb.random, 
         overall = TRUE, studlab = Africa$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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
  
  svg(file = ocepath, width = 14.5, height = 2.5)
  forest(ocemeta, comb.fixed = ocemeta$comb.fixed, comb.random = ocemeta$comb.random, 
         overall = TRUE, studlab = Oceania$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
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










#----------------------------------------- GPSC and serotypes driving resistance pre/post PCV10/13 ----------------

# classify vaccine period into pre/post pcv10

metadata_vaccine_period <- metadata |>
  mutate(
    vaccine_binom = ifelse(Vaccine_period_clean %in% c("PREPCV", "POSTPCV7", "POSTPCV10"), "PREPCV13",
                           ifelse(Vaccine_period_clean %in% "POSTPCV13", "POSTPCV13", "NA"))
  )

# metadata_vaccine_period$vaccine_binom <- factor(metadata_vaccine_period$vaccine_binom, levels = c("PREPCV10", "POSTPCV10"))

serotype_colors_vector <- setNames(metadata_vaccine_period$serotype_colour, metadata_vaccine_period$serotype_clean)


# obtain list of countries with 5 or more pre and post pcv10 isolates
continent_country_period <- metadata_vaccine_period |>
  group_by(Continent, Country, vaccine_binom) |>
  tally() |>
  filter(n > 10)

# filter for countries with both pre/post PCV10/13 isolates
continent_country_period_vaccine <- continent_country_period |>
  group_by(Country) |>
  filter(n() > 1) |>
  pivot_wider(names_from = vaccine_binom, values_from = n)

# obtain list of countries
vaccine_country_list <- continent_country_period_vaccine$Country




# test for EU

# obtain list of gpsc in this continent with 10 or more isolates pre and post pcv10
gpsc_period_vaccine <- metadata_vaccine_period |>
  filter(Country %in% vaccine_country_list & Continent %in% "EUROPE") |>
  group_by(GPSC, vaccine_binom) |>
  tally() |>
  pivot_wider(names_from = vaccine_binom, values_from = n) |>
  filter(PREPCV13 > 10 & POSTPCV13 > 10 & !GPSC %in% "unknown")

vaccine_gpsc_list <- gpsc_period_vaccine$GPSC

vaccine_gpsc_list_num <- sort(as.numeric(vaccine_gpsc_list))


eu_gpsc <- metadata_vaccine_period |>
  filter(Continent %in% "EUROPE" & Country %in% vaccine_country_list & GPSC %in% vaccine_gpsc_list) |>
  group_by(vaccine_binom, GPSC, serotype_clean) |>
  tally()


eu_gpsc_vector <- unique(sort(as.numeric(eu_gpsc$GPSC)))


eu_gpsc |>
  ggplot(aes(x = vaccine_binom, y = n, fill = serotype_clean, group = GPSC)) +
  geom_col(position = position_fill(), colour = "black") +
  facet_grid(~factor(GPSC, levels = eu_gpsc_vector), space="free_x", switch='x') +
  scale_y_continuous(expand = c(0,0)) +
  scale_colour_manual(values = "black") +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(title = element_text(size = 20),
        axis.title=element_text(size=18,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(size=10,color="black", angle = 90, hjust=1, vjust=0.2, face = "bold"),
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
    title = paste("Serotype distribution of GPSCs Europe pre- and post-PCV10", sep = " "),
    x = "GPSC"
  )



forloop_output_pre <- data.frame()
forloop_output_post <- data.frame()
forloop_output <- data.frame()
reference <- data.frame()

for (t in antibiotic) {
  
  # create SIR table for GPSCs pre-PCV10
  antibiotic_tile_pre <- metadata_vaccine_period |>
    filter(Continent %in% "EUROPE" & Country %in% vaccine_country_list & GPSC %in% vaccine_gpsc_list & vaccine_binom %in% "PREPCV13") |>
    group_by(GPSC, !!sym((t))) |>
    tally() |>
    pivot_wider(names_from = !!sym((t)), values_from = n)
  
  # set NAs to 0
  antibiotic_tile_pre[is.na(antibiotic_tile_pre)] <- 0
  
  # calculate total number of isolates
  antibiotic_tile_pre$total <- rowSums(antibiotic_tile_pre[, !colnames(antibiotic_tile_pre) %in% "GPSC"], na.rm = TRUE)
  
  # calculate proportion of resistant isolates for each antibiotic pre-PCV10
  row_for_dataframe_pre <- antibiotic_tile_pre |>
    mutate(
      antibiotic = t,
      proportion = R/total
    ) |>
    select(GPSC, proportion, antibiotic) |>
    arrange(factor(GPSC, levels = eu_gpsc_vector))
  
  # set NAs to 0
  row_for_dataframe_pre[is.na(row_for_dataframe_pre)] <- 0
  
  # pivot wider to allow rbind
  dataframe_row_pre <- row_for_dataframe_pre |>
    pivot_wider(names_from = GPSC, values_from = proportion)
  
  # rbind
  forloop_output_pre <- rbind(row_for_dataframe_pre, forloop_output_pre)
  
  
  
  # create SIR table for GPSCs pre-PCV10
  antibiotic_tile_post <- metadata_vaccine_period |>
    filter(Continent %in% "EUROPE" & Country %in% vaccine_country_list & GPSC %in% vaccine_gpsc_list & vaccine_binom %in% "POSTPCV13") |>
    group_by(GPSC, !!sym((t))) |>
    tally() |>
    pivot_wider(names_from = !!sym((t)), values_from = n)
  
  # set NAs to 0
  antibiotic_tile_post[is.na(antibiotic_tile_post)] <- 0
  
  # calculate total number of isolates
  antibiotic_tile_post$total <- rowSums(antibiotic_tile_post[, !colnames(antibiotic_tile_post) %in% "GPSC"], na.rm = TRUE)
  
  # calculate proportion of resistant isolates for each antibiotic pre-PCV10
  row_for_dataframe_post <- antibiotic_tile_post |>
    mutate(
      antibiotic = t,
      proportion = R/total
    ) |>
    select(GPSC, proportion, antibiotic) |>
    arrange(factor(GPSC, levels = eu_gpsc_vector))
  
  # set NAs to 0
  row_for_dataframe_post[is.na(row_for_dataframe_post)] <- 0
  
  # pivot wider to allow rbind
  dataframe_row_post <- row_for_dataframe_post |>
    pivot_wider(names_from = GPSC, values_from = proportion)


  
  # joining pre and post to calculate change
  
  row_for_dataframe <- left_join(row_for_dataframe_pre, row_for_dataframe_post, by = "GPSC", suffix = c("_pre", "_post"))
  
  
  # calculate change in percentage resistance in each GPSC for each antibiotic post-PCV10
  row_percentage_change <- row_for_dataframe |>
    mutate(
      proportion_change = proportion_pre - proportion_post
    ) |>
    select(-c(antibiotic_post, proportion_pre, proportion_post)) |>
    rename(
      antibiotic = antibiotic_pre
    )
  
  # rbind output
  
  forloop_output <- rbind(row_percentage_change, forloop_output)

  # rbind reference
  
  reference <- rbind(row_for_dataframe, reference)
  
}


# change names of antibiotic column and relevel
forloop_output <- left_join(forloop_output, antibiotic_tidying, by = "antibiotic") 

forloop_output$antibiotic_tidy <- factor(forloop_output$antibiotic_tidy, levels = c("Rifampin", "Chloramphenicol", "Tetracycline", 
                                                                                    "Co-trimoxazole", "Clindamycin", "Erythromycin", 
                                                                                    "Meropenem", "Cefotaxime", 
                                                                                    "Ceftriaxone", "Cefuroxime", 
                                                                                    "Amoxicillin", "Penicillin"))

reference_clean <- reference |>
  rename(
    antibiotic = antibiotic_pre
  ) |>
  select(-antibiotic_post)

reference_raw <- left_join(reference_clean, antibiotic_tidying, by = "antibiotic")

reference_raw$antibiotic_tidy <- factor(reference_raw$antibiotic_tidy, levels = c("Penicillin", "Amoxicillin", "Cefuroxime",
                                                                          "Ceftriaxone", "Cefotaxime", "Meropenem",
                                                                          "Erythromycin", "Clindamycin", "Co-trimoxazole",
                                                                          "Tetracycline", "Chloramphenicol", "Rifampin"))



# write out raw and calculated reference file for each continent
resistance_reference_raw <- reference_raw |>
  mutate(
    proportion_change = proportion_pre - proportion_post
  ) |>
  select(-antibiotic) |>
  relocate(antibiotic_tidy, .before = GPSC)



resistance_reference_change <- forloop_output |>
  select(-antibiotic) |>
  pivot_wider(names_from = GPSC, values_from = proportion_change)




# tile plot
forloop_output |>
  ggplot(aes(x = factor(GPSC, levels = eu_gpsc_vector), y = antibiotic_tidy, fill = proportion_change)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1)) +
  coord_fixed() +
  theme(
    plot.margin = unit(c(-4,0,0,0), 'lines')
  ) +
  labs(
    y = "Antibiotic",
    fill = "Proportion \nchange",
    x = "GPSC",
    title = "Change in prevalence of resistance pre- and post-PCV13 in Europe"
  ) +
  theme(
    title = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(size = 32),
    axis.text.x = element_text(size = 28),
    axis.ticks = element_blank(),
    legend.title = element_text(margin = margin(b = 20, unit = "pt"), size = 24),
    legend.text = element_text(size = 24),
    legend.key.size = unit(1.2, "cm")
  )







#----------------------------------------- Proportion of serotypes covered by vaccines ----------------------------


#metadata_res_any <- metadata |>
#  mutate(
#    resistance_c = paste(WGS_PEN_SIR_Meningitis, WGS_AMO_SIR, WGS_TAX_SIR_Meningitis, WGS_CFT_SIR_Meningitis,
#                          WGS_MER_SIR, WGS_CFX_SIR, WGS_ERY_SIR, WGS_CLI_SIR, WGS_COT_SIR, WGS_TET_SIR,
#                          WGS_CHL_SIR, WGS_RIF_SIR, sep = "_")
#  ) |>
#  filter(str_detect(resistance_c, "R") & Vaccine_period_clean %in% "POSTPCV13") 


#metadata_res_any |>
#  filter(WGS_PEN_SIR_Meningitis %in% "R") |>
#  group_by(Continent, serotype_clean) |>
#  tally()


# create empty data frame to store number of isolates of each serotype that are resistant to each antibiotic in each continent
serotype_coverage_df <- data.frame()
serotype_coverage_ref <- data.frame()

for(i in antibiotic){
  
  # tally out post pcv-10/13 isolates of each serotype in each continent that are resistant to the antibiotic
  row_for_serotype_coverage_df <- metadata |>
    filter(!!((sym(i))) %in% "R" & Vaccine_period_clean %in% "POSTPCV13") |>
    group_by(Continent, serotype_clean) |>
    tally() |>
    mutate(
      antibiotic = i
    )
  

  # pivot wider to allow rbind
  dataframe_row_serotype_coverage <- row_for_serotype_coverage_df |>
    pivot_wider(names_from = antibiotic, values_from = n)
  
  # rbind reference
  serotype_coverage_ref <- rbind(dataframe_row_serotype_coverage, serotype_coverage_ref)
  
  # set NAs to 0 
  serotype_coverage_ref[is.na(serotype_coverage_ref)] <- 0
  
  
  # rbind table for plotting
  serotype_coverage_df <- rbind(row_for_serotype_coverage_df, serotype_coverage_df)
  
  
}

write_csv(serotype_coverage_vaccine, file = "coverage reference.csv")



# if else statement to class serotype into vaccine coverage
# show colours only for vaccine types, grey for all NVTs


vaccine_table_unique_serotypes <- vaccine_table |>
  filter(!is.na(vaccine_unique) & !vaccine %in% "PCV25") |>
  select(serotype, vaccine_unique) |>
  rename(
    serotype_clean = serotype
  )

serotype_coverage_vaccine <- left_join(serotype_coverage_df, vaccine_table_unique_serotypes, by = "serotype_clean") 


serotype_coverage_vaccine_unique <- serotype_coverage_vaccine |>
  mutate(
    unique_serotype = ifelse(is.na(vaccine_unique), "NVT", serotype_clean),
    PCV13_vaccine = ifelse(serotype_clean %in% pcv13_serotypes, "VT", "NVT")
  )

serotype_coverage_vaccine_unique_clean <- left_join(serotype_coverage_vaccine_unique, antibiotic_tidying, by = "antibiotic")

serotype_coverage_vaccine_unique_clean$antibiotic_tidy <- factor(serotype_coverage_vaccine_unique_clean$antibiotic_tidy, 
                                                                 levels = c("Penicillin", "Amoxicillin", "Cefuroxime",
                                                                            "Ceftriaxone", "Cefotaxime", "Meropenem",
                                                                            "Erythromycin", "Clindamycin", "Co-trimoxazole",
                                                                            "Tetracycline", "Chloramphenicol", "Rifampin"))


# set colours
vaccine_colour <- serotype_colours |>
  rename(
    unique_serotype = serotype_clean
  )

serotype_coverage_vaccine_unique_clean_coloured <- left_join(serotype_coverage_vaccine_unique_clean, vaccine_colour, by = "unique_serotype")

vaccine_colors_vector <- setNames(serotype_coverage_vaccine_unique_clean_coloured$In_silico_serotype__colour, 
                                  serotype_coverage_vaccine_unique_clean_coloured$unique_serotype)





# plot PCV7 serotypes to obtain legend
pcv7 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV7") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV7 serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv7_legend <- get_legend(pcv7)


# plot PCV10 serotypes to obtain legend
pcv10 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV10_GSK") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV10 (GSK) unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv10_legend <- get_legend(pcv10)


# plot PCV13 serotypes to obtain legend
pcv13 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV13") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV13 unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv13_legend <- get_legend(pcv13)


# plot PCV15 serotypes to obtain legend
pcv15 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV15") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV15 unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv15_legend <- get_legend(pcv15)


# plot PCV20 serotypes to obtain legend
pcv20 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV20") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV20 unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv20_legend <- get_legend(pcv20)


# plot PCV21 serotypes to obtain legend
pcv21 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV21_SK_SANOFI") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV21 (SK/SANOFI) unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv21_legend <- get_legend(pcv21)


# plot PCV20 serotypes to obtain legend
pcv24 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV24_GSK") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV24 (GSK) unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv24_legend <- get_legend(pcv24)



# plot PCV20 serotypes to obtain legend
pcv31 <- serotype_coverage_vaccine_unique_clean |>
  filter(vaccine_unique %in% "PCV31") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(name = "PCV31 unique serotypes",
                    values = vaccine_colors_vector,
                    guide = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

pcv31_legend <- get_legend(pcv31)


# plot NVT serotypes to obtain legend
NVT <- serotype_coverage_vaccine_unique_clean |>
  filter(unique_serotype %in% "NVT") |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  scale_fill_manual(values = vaccine_colors_vector,
                    name = "Non-vaccine serotypes") +
  theme(
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm")
  )

NVT_legend <- get_legend(NVT)

legend <- ggarrange(pcv7_legend, pcv10_legend, pcv13_legend, pcv15_legend, pcv20_legend, 
                    pcv21_legend, pcv24_legend, pcv31_legend, NVT_legend,
          ncol = 1, align = "v", vjust = 20)

ggsave(legend, file = "sample figures/vaccine coverage legend.svg", , dpi = 320, width = 14, height = 29, units = "cm")





coverage_figure <- serotype_coverage_vaccine_unique_clean |>
  ggplot(aes(x = antibiotic_tidy, y = n, fill = unique_serotype, colour = unique_serotype, pattern = PCV13_vaccine)) +
  geom_bar(position = position_stack(), stat = "identity") +
  geom_col_pattern(
    pattern_spacing = 0.01, 
    pattern_fill = "black", 
    pattern_size = 0.5, 
    pattern_linetype = 0.5,
    pattern_key_scale_factor = 0.5,
    pattern_orientation = "vertical", show.legend = TRUE) +
  facet_wrap(~Continent, ncol = 3, scales = "free_y")+
  scale_fill_manual(name = "Serotype",
                    values = vaccine_colors_vector,
                    guide = "none") +
  scale_colour_manual(values = vaccine_colors_vector,
                      guide = "none") +
  scale_pattern_manual(
    name = "PCV13 serotypes",
    values = c("VT" = "stripe", "NVT" = "none")
  )  + 
  labs(
    x = "Antibiotic",
    y = "Number of isolates",
    title = "Serotypes driving AMR post-PCV13 implementation in each continent"
  ) +
  theme(
    title = element_text(size = 24),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 22), 
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(size = 16),
    legend.background = element_blank()
  )

ggsave(coverage_figure, file = "sample figures/vaccine coverage figure post pcv13 shaded.svg", , dpi = 320, width = 48, height = 36, units = "cm")














#------------------------------------ notes beyond this -------------------------------------

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

continent <- c("ASIA", "AFRICA", "LATIN AMERICA", "OCEANIA", "EUROPE", "NORTH AMERICA")

serotype_colors_vector <- setNames(metadata_vaccine_period$In_silico_serotype__colour, metadata_vaccine_period$serotype_clean)


for(i in continent){
  
  
  forloop_output <- data.frame()
  csv_output = list()
  
  
  as_gpsc <- metadata_vaccine_period |>
    group_by(GPSC) |>
    filter(Continent %in% i & Manifest_category %in% "IPD") |>
    tally() |>
    filter(n > 9) 
  
  as_gpsc_list <- as_gpsc$GPSC
  
  title_continent <- str_to_title(i)
  
  gpsc_serotype_tally <- metadata_vaccine_period |>
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
    
    antibiotic_tile <- metadata_vaccine_period |>
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
af_gpsc <- metadata_vaccine_period |>
  group_by(GPSC) |>
  filter(Continent %in% "AFRICA" & Manifest_category %in% "IPD") |>
  tally() |>
  filter(n > 9) |>
  slice_head(n = 35)

af_gpsc_list <- af_gpsc$GPSC

af_gpsc_serotype_tally <- metadata_vaccine_period |>
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
  
  antibiotic_tile <- metadata_vaccine_period |>
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
