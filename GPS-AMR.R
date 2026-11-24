library(tidyverse)
library(scatterpie)
library(ggpubr)
library(meta)
library(UpSetR)
library(rnaturalearth)

filter <- dplyr::filter
rename <- dplyr::rename
summarise <- dplyr::summarise



#-------------------------------------------- read in files -------------------------------------------
# set working directory
setwd("C:/Users/choong/OneDrive - EMBLEBI/Desktop/GPS-AMR")

# monocle metadata
full_metadata <- read_csv(file = "raw data/monocle-metadata-all-columns.csv")

# australian isolates
australia <- read_csv(file = "raw data/australia.csv")

# read in unpublished cambodian metadata
cambodia <- read.csv(file = "raw data/cambodia.csv", na.strings = c("","NA"))

cambodia <- cambodia |>
  filter(Published %in% "N")

# rbind Monocle data, australian and unpublished cambodian metadata
# requires the use of plyr package but has to unload to prevent conflicts and masking
library(plyr)

metadata <- rbind.fill(full_metadata, australia, cambodia)

detach("package:plyr", unload = TRUE)

# tidy antibiotics for plotting
antibiotic_tidying <- read.csv(file = "raw data/antibiotic.csv")

# read in file containing tidied serotypes
clean_serotype <- read_csv(file = "raw data/tidied_serotypes.csv")

# read in serotype colours
serotype_colours <- read_csv(file = "raw data/serotype_colours.csv")

# read in table of serotypes covered by each vaccine
vaccine_table <- read_csv(file = "vaccine_chart.csv")

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

# continent tidying
metadata$Continent <- str_to_title(metadata$Continent)

# country name tidying
metadata$Country <- str_to_title(metadata$Country)

metadata$Country_mapping <- ifelse(metadata$Country %in% "Russian Federation", "Russia",
                                   ifelse(metadata$Country %in% "Trinidad And Tobago", "Trinidad and Tobago",
                                          ifelse(metadata$Country %in% "United States", "United States of America",
                                                 ifelse(metadata$Country %in% "The Gambia", "Gambia", metadata$Country))))


# remove France and Spain with heavy selection bias of serogroup 24
metadata <- metadata |>
  filter(!Country %in% c("France", "Spain", "Italy", "Czech Republic"))

# change japan to post-pcv13 (from a bug in previous GPS db release)
metadata <- metadata |>
  mutate(
    Vaccine_period = ifelse(Country == "Japan", gsub("7", "13", metadata$Vaccine_period), Vaccine_period),
    PCV_type = ifelse(Country == "Japan", gsub("7", "13", metadata$PCV_type), PCV_type)
  )


# set NA in GPSC to unknonwn
metadata$GPSC <- ifelse(is.na(metadata$GPSC), "unknown", metadata$GPSC)

# remove year after vaccine period
metadata$Vaccine_period_clean <- gsub("\\-.*", "", metadata$Vaccine_period)

# arrange PCV
metadata$Vaccine_period_clean <- factor(metadata$Vaccine_period_clean, 
                                        levels = c("PREPCV", "POSTPCV7", "POSTPCV10", "POSTPCV13"))


# left join to match serotype with tidied serotype
metadata <- left_join(metadata, clean_serotype, by = "In_silico_serotype") |>
  relocate(serotype_clean, .after = In_silico_serotype)

# replace NAs in clean serotype column with in silico serotype
metadata$serotype_clean <- ifelse(is.na(metadata$serotype_clean), metadata$In_silico_serotype, metadata$serotype_clean)


# join in colours 
metadata <- left_join(metadata, serotype_colours, by = "serotype_clean")

colnames(metadata)[colnames(metadata) == 'In_silico_serotype__colour.y'] <- 'serotype_colour'

serotype_colors_vector <- setNames(metadata$In_silico_serotype__colour, metadata$serotype_clean)

### kanamycin resistance
metadata <- metadata |>
  mutate(
    WGS_KAN_SIR = ifelse(startsWith(Other, "APH3PRIMEIII"), "R", "S")
  ) |>
  relocate(WGS_KAN_SIR, .after = WGS_VAN_SIR)

# substitute FLAG for S
antibiotic <- c("WGS_PEN_SIR_Meningitis", "WGS_AMO_SIR", "WGS_TAX_SIR_Meningitis", "WGS_CFT_SIR_Meningitis",
                "WGS_MER_SIR", "WGS_CFX_SIR", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR",
                "WGS_CHL_SIR", "WGS_RIF_SIR", "WGS_KAN_SIR")

metadata[antibiotic] <- lapply(metadata[antibiotic], gsub, pattern = "FLAG", replacement = "S", fixed = TRUE)

# classify vaccine period into pre/post pcv13

metadata <- metadata |>
  mutate(
    vaccine_binom = ifelse(Vaccine_period_clean %in% c("PREPCV", "POSTPCV7", "POSTPCV10"), "PREPCV13",
                           ifelse(Vaccine_period_clean %in% "POSTPCV13", "POSTPCV13", "NA"))
  )

# isolates with >=15 or >5 will be labelled as NA 
# re-classify and tidying up age and clinical manifestation
metadata <- metadata |>
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

# fix continent name
metadata <- metadata |>
  mutate(
    Continent = str_to_title(Continent)
  )

#----------------------------------------- MDR -----------------------------------------

# add resistogram column for beta lactams
metadata$resistogram_beta_lactam <- paste(metadata$WGS_PEN_SIR_Meningitis, metadata$WGS_AMO_SIR, metadata$WGS_TAX_SIR_Meningitis, 
                                          metadata$WGS_CFT_SIR_Meningitis, metadata$WGS_MER_SIR, metadata$WGS_CFX_SIR)

# column to count beta lactam resistance status
metadata$resistance_beta_lactam <- ifelse(str_detect(metadata$resistogram_beta_lactam, 'R'), "R", "S")

# add resistogram column for beta lactams and other MDR antibiotics
metadata$resistogram <- paste(metadata$resistance_beta_lactam, metadata$WGS_CHL_SIR, metadata$WGS_ERY_SIR, 
                              metadata$WGS_COT_SIR, metadata$WGS_TET_SIR, metadata$WGS_KAN_SIR)

# column to count MDR status (if at least 3 or more R or I are found = MDR)
metadata$MDR <- ifelse(str_count(metadata$resistogram, 'R') > 2, "R", "S")


# pre-pcv13 MDR visualisation
pre_resistogram_tally <- metadata |>
  filter(MDR == "R" & vaccine_binom == "PREPCV13") |>
  group_by(resistogram) |>
  tally() |>
  arrange(desc(n)) |>
  slice_head(n = 11)

pre_resistogram <- pre_resistogram_tally$resistogram

mdr_isolates_pre <- metadata |>
  filter(MDR == "R" & vaccine_binom == "PREPCV13" & resistogram %in% pre_resistogram)

beta_lactam <- mdr_isolates_pre |>
  filter(resistance_beta_lactam == "R")
chloramphenicol <- mdr_isolates_pre |>
  filter(WGS_CHL_SIR == "R")
erythromycin <- mdr_isolates_pre |>
  filter(WGS_ERY_SIR == "R")
cotrimoxazole <- mdr_isolates_pre |>
  filter(WGS_COT_SIR == "R")
tetracycline <- mdr_isolates_pre |>
  filter(WGS_TET_SIR == "R")
kanamycin <- mdr_isolates_pre |>
  filter(WGS_KAN_SIR == "R")

beta_lactam_vector <- beta_lactam$Lane_id
chloramphenicol_vector <- chloramphenicol$Lane_id
erythromycin_vector <- erythromycin$Lane_id
cotrimoxazole_vector <- cotrimoxazole$Lane_id
tetracycline_vector <- tetracycline$Lane_id
kanamycin_vector <- kanamycin$Lane_id

listInput_pre <- list(BetaLactam = beta_lactam_vector, Chloramphenicol = chloramphenicol_vector,
                      Erythromycin = erythromycin_vector, Cotrimoxazole = cotrimoxazole_vector,
                      Tetracycline = tetracycline_vector, Kanamycin = kanamycin_vector)

svg(file = "pre_upset.svg", width = 18, height = 12)

upset(fromList(listInput_pre), nsets = 6, order.by = "freq", mainbar.y.label = "Number of isolates",
      sets.x.label = "Number of isolates", text.scale = c(3, 3, 3, 2, 3, 3), point.size = 4, line.size = 1)
dev.off()

# Post-pcv13 MDR visualisation
post_resistogram_tally <- metadata |>
  filter(MDR == "R" & vaccine_binom == "POSTPCV13") |>
  group_by(resistogram) |>
  tally() |>
  arrange(desc(n)) |>
  slice_head(n = 11)

post_resistogram <- post_resistogram_tally$resistogram

mdr_isolates_post <- metadata |>
  filter(MDR == "R" & vaccine_binom == "POSTPCV13" & resistogram %in% post_resistogram)

beta_lactam <- mdr_isolates_post |>
  filter(resistance_beta_lactam == "R")
chloramphenicol <- mdr_isolates_post |>
  filter(WGS_CHL_SIR == "R")
erythromycin <- mdr_isolates_post |>
  filter(WGS_ERY_SIR == "R")
cotrimoxazole <- mdr_isolates_post |>
  filter(WGS_COT_SIR == "R")
tetracycline <- mdr_isolates_post |>
  filter(WGS_TET_SIR == "R")
kanamycin <- mdr_isolates_post |>
  filter(WGS_KAN_SIR == "R")

beta_lactam_vector <- beta_lactam$Lane_id
chloramphenicol_vector <- chloramphenicol$Lane_id
erythromycin_vector <- erythromycin$Lane_id
cotrimoxazole_vector <- cotrimoxazole$Lane_id
tetracycline_vector <- tetracycline$Lane_id
kanamycin_vector <- kanamycin$Lane_id

listInput_post <- list(BetaLactam = beta_lactam_vector, Chloramphenicol = chloramphenicol_vector,
                       Erythromycin = erythromycin_vector, Cotrimoxazole = cotrimoxazole_vector,
                       Tetracycline = tetracycline_vector, Kanamycin = kanamycin_vector)

svg(file = "post_upset.svg", width = 18, height = 12)
upset(fromList(listInput_post), nsets = 6, order.by = "freq", mainbar.y.label = "Number of isolates",
      sets.x.label = "Number of isolates", text.scale = c(3, 3, 3, 2, 3, 3), point.size = 4, line.size = 1)
dev.off()

#----------------------------------------- Summary tables -----------------------------------------

dir.create("summarised data")
dir.create("summarised data/data overview")

setwd("C:/Users/choong/OneDrive - EMBLEBI/Desktop/GPS-AMR/")

continent_period_supplementary <- metadata |>
  group_by(Continent, vaccine_binom) |>
  tally() |>
  pivot_wider(names_from = vaccine_binom, values_from = n) |>
  relocate(PREPCV13, .after = Continent)


write_csv(continent_period_supplementary, file = "summarised data/data overview/continent period table.csv")

country_period_supplementary <- metadata |>
  group_by(Country, vaccine_binom) |>
  tally() |>
  pivot_wider(names_from = vaccine_binom, values_from = n, values_fill = 0) |>
  mutate(
    Total = rowSums(across(where(is.numeric)))
  ) |>
  relocate(PREPCV13, .after = Country)

write_csv(country_period_supplementary, file = "summarised data/data overview/country period table.csv")


continent_manifest_supplementary <- metadata |>
  group_by(Continent, Manifest_category) |>
  tally() |>
  pivot_wider(names_from = Manifest_category, values_from = n, values_fill = 0) 

write_csv(continent_manifest_supplementary, file = "summarised data/data overview/continent manifest table.csv")



#----------------------------------------- maps -----------------------------------------

dir.create("summary plots")

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
  geom_sf(aes(fill = log10(n + 1))) +
  scale_fill_gradientn(colours = c("blue", "steelblue", "pink", "red"), na.value = "white", limits = c(0, 4.2)) +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 23),
    legend.key.size = unit(1.2, "cm")
  ) 

world_plot_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "summary plots", "world plot.svg", fsep = "/")

ggsave(world_plot, filename = world_plot_path, dpi = 320, width =50.58, height = 18.91, units = "cm")


#----------------------------------------- vaccine reference table ----------------------------------------

# create a column with values as 1
vaccine_table$value <- 1

# reorder factors
vaccine_table$vaccine <- factor(vaccine_table$vaccine, levels = c("PCV31_VAXCYTE",
                                                                  "PCV25_IVT",
                                                                  "PCV24_GSK",
                                                                  "PCV21_SK_SANOFI",
                                                                  "PCV20_PFIZER", 
                                                                  "PCV15_MERCK",
                                                                  "PCV13_PFIZER",
                                                                  "PCV10_PNEUMOSIL",
                                                                  "PCV10_GSK",
                                                                  "PCV7_PFIZER"))

# labels
vaccine_labels <- c("PCV31 (Vaxcyte)",
                    "PCV25 (Inventprise)",
                    "PCV24 (GSK)",
                    "PCV21 (SK/SANOFI)",
                    "PCV20 (Pfizer)", 
                    "PCV15 (Merck)",
                    "PCV13 (Pfizer)",
                    "PCV10 (PNEUMOSIL)",
                    "PCV10 (GSK)",
                    "PCV7 (Pfizer)")

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
    axis.text = element_text(size = 20)
  )

pcv7 <- vaccine_table |>
  filter(vaccine %in% "PCV7_PFIZER")

pcv7_serotypes <- pcv7$serotype

pcv10 <- vaccine_table |>
  filter(vaccine %in% "PCV10_GSK")

pcv10_serotypes <- pcv10$serotype

pcv13 <- vaccine_table |>
  filter(vaccine %in% "PCV13_PFIZER")

pcv13_serotypes <- pcv13$serotype

pcv20 <- vaccine_table |>
  filter(vaccine %in% "PCV20_PFIZER")

pcv20_serotypes <- pcv20$serotype

PCV25 <- vaccine_table |>
  filter(vaccine %in% "PCV25_IVT")

pcv25_serotypes <- PCV25$serotype



#----------------------------------------- year/manifest plot ----------------------------------------



# table with continent manifestation for plotting
continent_year_tally <- metadata |>
  filter(!is.na(Year) & !is.na(Continent) & !Year %in% "_") |>
  group_by(Continent, Year, Manifest_category) |>
  tally() |>
  mutate(
    continent_tidy = str_to_title(Continent),
    join = paste0(continent_tidy, " ", Manifest_category)
  )

# relevel factors
continent_year_tally$Manifest_category <- factor(continent_year_tally$Manifest_category, levels = c("CARRIAGE", "DIS_CAR", "IPD", "NONINVASIVE_DISEASE", "UNKNOWN"))

# manifest, time of collection and continent plot
continent_year_plot <- continent_year_tally |>
  filter(as.integer(Year) > 1989 & !Continent %in% "CENTRAL AMERICA") |>
  ggplot(aes(x = as.factor(Year), y = n, fill = Manifest_category)) +
  geom_bar(stat = "identity") +
  facet_grid(continent_tidy ~ ., scales = "free_y") +
  scale_fill_manual(values = c("#BA55D3", "#2ECC71", "#00BFFF", "#FF8C00", "#808080"),
                    labels = c("Carriage", "Diseased carriage", "IPD", "Non IPD", "Unknown")) +
  labs(
    x = "Year",
    y = "Number of isolates",
    fill = "Clinical \nManifestation",
    
  ) +
  geom_vline(xintercept = "2000", linetype = "dashed") +
  geom_vline(xintercept = "2010", linetype = "dashed") +
  geom_vline(xintercept = "2020", linetype = "dashed") +
  theme(
    title = element_text(size = 28),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 30),
    legend.key.size = unit(1.5, "cm"),
    axis.text = element_text(size = 30),
    strip.text.y = element_text(size = 28),
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
    axis.title = element_text(size = 30),
    panel.spacing.y = unit(2, "lines")
  )


continent_year_plot_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "summary plots", "continent year plot.svg", fsep = "/")

ggsave(continent_year_plot, filename = continent_year_plot_path, dpi = 320, width =80.58, height = 50, units = "cm")


dir.create("summarised data/data overview")

continent_country_table <- metadata |>
  group_by(Continent, Country) |>
  tally()

continent_manifest_table <- metadata |>
  group_by(Continent, Manifest_category) |>
  tally() |>
  pivot_wider(names_from = Manifest_category, values_from = n)

write_csv(continent_country_table, file = "summarised data/data overview/continent country table.csv")

write_csv(continent_manifest_table, file = "summarised data/data overview/continent manifest table.csv")

#----------------------------------------- panels ----------------------------------------

dir.create("figures")

# vancomycin not included as R isolates were not found
# no CHL resistance found for GPSC1
# join serotype/gpsc table for CLI, COT, ERY, PEN, RIF, TET
antibiotic <- c("WGS_PEN_SIR_Meningitis", "WGS_AMO_SIR", "WGS_TAX_SIR_Meningitis", "WGS_CFT_SIR_Meningitis",
                "WGS_MER_SIR", "WGS_CFX_SIR", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR",
                "WGS_CHL_SIR", "WGS_RIF_SIR", "WGS_KAN_SIR", "MDR")

antibiotic_join <- c("WGS_PEN_SIR_Meningitis", "WGS_ERY_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "WGS_TET_SIR", "WGS_CHL_SIR", "WGS_RIF_SIR", "WGS_KAN_SIR", "MDR")

antibiotic_diverse_serotype <- c("WGS_PEN_SIR_Meningitis", "WGS_ERY_SIR", "WGS_TET_SIR", "WGS_CLI_SIR", "WGS_COT_SIR", "MDR")

continent <- c("Africa", "North America", "South America", "Oceania", "Asia", "Europe", "Overall")

period <- c("PREPCV13", "POSTPCV13")

world <- map_data("world")
world <- subset(world, region != "Antarctica")

setwd("C:/Users/choong/OneDrive - EMBLEBI/Desktop/GPS-AMR/figures")

dir.create("Proportion")


analysis_proportion_df <- data.frame()



### panel and forest plot for loops
for(p in period){
  
  metadata_period <- metadata |>
    filter(vaccine_binom %in% p)
  
  maindir <- paste0("/",p," plots")
  
  dir.create(maindir)
  
  for(i in antibiotic){
    
    #filter out isolates that failed resistance prediction
    tidy_dataset <- metadata_period |> filter(!!((sym(i))) %in% c("S", "I", "R") & !Country %in% c("France", "Spain"))
    
    ### create new directory for each looped variable
    newdir <- paste0(maindir,"/", i, " plots ", p)
    
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
    age_manifest_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", "Proportion", paste(newdir, "/age manifest table.csv", sep = ""), fsep = "/")
    
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
      
      pivot_SIR <- pivot_SIR |>
        mutate(
          R = ifelse("R" %in% names(pivot_SIR), R, NA)
        )
      
      pivot_SIR[is.na(pivot_SIR)] = 0
      
      pivot_plot <- pivot_SIR |>
        mutate(
          Longitude = as.numeric(case_when(
            Continent == "AFRICA" ~ 21.890909,
            Continent == "NORTH AMERICA" ~ -100.538118,
            Continent == "SOUTH AMERICA" ~ -53.604523,
            Continent == "EUROPE" ~ 15.738565,
            Continent == "ASIA" ~ 110.308879,
            Continent == "OCEANIA" ~ 150.562787
          )),
          Latitude = as.numeric(case_when(
            Continent == "AFRICA" ~ 6.208920,
            Continent == "NORTH AMERICA" ~ 49.249542,
            Continent == "SOUTH AMERICA" ~ -12.581874,
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
      
      pivot_SIR <- pivot_SIR |>
        mutate(
          R = ifelse("R" %in% names(pivot_SIR), R, NA)
        )
      
      pivot_SIR[is.na(pivot_SIR)] = 0
      
      pivot_plot <- pivot_SIR |>
        mutate(
          Longitude = as.numeric(case_when(
            Continent == "AFRICA" ~ 21.890909,
            Continent == "NORTH AMERICA" ~ -100.538118,
            Continent == "SOUTH AMERICA" ~ -53.604523,
            Continent == "EUROPE" ~ 15.738565,
            Continent == "ASIA" ~ 110.308879,
            Continent == "OCEANIA" ~ 150.562787
          )),
          Latitude = as.numeric(case_when(
            Continent == "AFRICA" ~ 6.208920,
            Continent == "NORTH AMERICA" ~ 49.249542,
            Continent == "SOUTH AMERICA" ~ -12.581874,
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
    
    SIR_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", "Proportion", paste(newdir, "/SIR table.csv", sep = ""), fsep = "/")
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
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
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
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
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
      filter(!gpsc_tidy %in% c("others"))
    
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
    
    
    
    full_plot_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", "Proportion", paste(newdir, "/full plot.svg", sep = ""), fsep = "/")
    
    ggsave(full_plot, filename = full_plot_path, dpi = 320, width = 26.58, height = 16.91, units = "cm")
    
    
    # forest plots
    
    newdir <- paste0("Proportion", maindir,"/", i, " forest plot ", p)
    
    dir.create(newdir)
    
    napath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/na forest.svg", sep = ""), fsep = "/")
    
    complete <- tidy_dataset |> filter(!!((sym(i))) %in% c("S", "I", "R"))
    
    
    
    #forest plot for resistance in north american countries
    NorthAmerica <- complete |>
      filter(Continent %in% "North America") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    
    if(nrow(NorthAmerica) > 0){
      if(!"R" %in% names(NorthAmerica)) {NorthAmerica$R <- 0
      NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
      NorthAmerica[is.na(NorthAmerica)] <- 0
      nameta <- metaprop(R, total, data = NorthAmerica)
      }
      
      else NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
      NorthAmerica[is.na(NorthAmerica)] <- 0
      nameta <- metaprop(R, total, data = NorthAmerica)
      
      svg(file = napath, width = 11, height = 3)
      forest(nameta, comb.fixed = nameta$comb.fixed, comb.random = nameta$comb.random, 
             overall = TRUE, studlab = NorthAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    #South america
    SouthAmerica <- complete |>
      filter(Continent %in% "South America") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    lapath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/la forest.svg", sep = ""), fsep = "/")
    
    
    if(nrow(SouthAmerica) > 0){
      if(!"R" %in% names(SouthAmerica)) {SouthAmerica$R <- 0
      SouthAmerica$total <- rowSums(SouthAmerica[, !colnames(SouthAmerica) %in% "Country"], na.rm = TRUE)
      SouthAmerica[is.na(SouthAmerica)] <- 0
      lameta <- metaprop(R, total, data = SouthAmerica)
      }
      
      else SouthAmerica$total <- rowSums(SouthAmerica[, !colnames(SouthAmerica) %in% "Country"], na.rm = TRUE)
      SouthAmerica[is.na(SouthAmerica)] <- 0
      lameta <- metaprop(R, total, data = SouthAmerica)
      
      svg(file = lapath, width = 11, height = 3)
      forest(lameta, comb.fixed = lameta$comb.fixed, comb.random = lameta$comb.random, 
             overall = TRUE, studlab = SouthAmerica$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    #europe
    Europe <- complete |>
      filter(Continent %in% "Europe") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    eupath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/eu forest.svg", sep = ""), fsep = "/")
    
    
    if(nrow(Europe) > 0){
      if(!"R" %in% names(Europe)){Europe$R <- 0
      Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
      Europe[is.na(Europe)] <- 0
      eumeta <- metaprop(R, total, data = Europe)
      }
      
      else Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
      Europe[is.na(Europe)] <- 0
      eumeta <- metaprop(R, total, data = Europe)
      
      svg(file = eupath, width = 11, height = 6)
      forest(eumeta, comb.fixed = eumeta$comb.fixed, comb.random = eumeta$comb.random, 
             overall = TRUE, studlab = Europe$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    #asia
    Asia <- complete |>
      filter(Continent %in% "Asia") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    aspath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/as forest.svg", sep = ""), fsep = "/")
    
    
    if(nrow(Asia) > 0){
      if(!"R" %in% names(Asia)){Asia$R <- 0
      Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
      Asia[is.na(Asia)] <- 0
      asmeta <- metaprop(R, total, data = Asia)
      }
      
      else Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
      Asia[is.na(Asia)] <- 0
      asmeta <- metaprop(R, total, data = Asia)
      
      svg(file = aspath, width = 11, height = 6)
      forest(asmeta, comb.fixed = asmeta$comb.fixed, comb.random = asmeta$comb.random, 
             overall = TRUE, studlab = Asia$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    #africa
    Africa <- complete |>
      filter(Continent %in% "Africa") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    afpath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/af forest.svg", sep = ""), fsep = "/")
    
    
    if(nrow(Africa) > 0){
      if(!"R" %in% names(Africa)){Africa$R <- 0
      Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
      Africa[is.na(Africa)] <- 0
      afmeta <- metaprop(R, total, data = Africa)
      }
      
      else Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
      Africa[is.na(Africa)] <- 0
      afmeta <- metaprop(R, total, data = Africa)
      
      svg(file = afpath, width = 11, height = 6)
      forest(afmeta, comb.fixed = afmeta$comb.fixed, comb.random = afmeta$comb.random, 
             overall = TRUE, studlab = Africa$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    #oceania
    Oceania <- complete |>
      filter(Continent %in% "Oceania") |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    ocepath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/oce forest.svg", sep = ""), fsep = "/")
    
    if(nrow(Oceania) > 0){
      if(!"R" %in% names(Oceania)){Oceania$R <- 0
      Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
      Oceania[is.na(Oceania)] <- 0
      ocemeta <- metaprop(R, total, data = Oceania) 
      }
      
      else Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
      Oceania[is.na(Oceania)] <- 0
      ocemeta <- metaprop(R, total, data = Oceania)
      
      svg(file = ocepath, width = 11, height = 2.5)
      forest(ocemeta, comb.fixed = ocemeta$comb.fixed, comb.random = ocemeta$comb.random, 
             overall = TRUE, studlab = Oceania$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
    # overall
    Overall <- complete |>
      group_by(Country, !!sym((i))) |>
      tally() |>
      pivot_wider(names_from = !!sym((i)), values_from = n)
    
    overallpath <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/overall forest.svg", sep = ""), fsep = "/")
    
    if(nrow(Overall) > 0){
      if(!"R" %in% names(Overall)){Overall$R <- 0
      Overall$total <- rowSums(Overall[, !colnames(Overall) %in% "Country"], na.rm = TRUE)
      Overall[is.na(Overall)] <- 0
      overallmeta <- metaprop(R, total, data = Overall) 
      }
      
      else Overall$total <- rowSums(Overall[, !colnames(Overall) %in% "Country"], na.rm = TRUE)
      Overall[is.na(Overall)] <- 0
      overallmeta <- metaprop(R, total, data = Overall)
      
      svg(file = overallpath, width = 14.5, height = 18)
      forest(overallmeta, comb.fixed = overallmeta$comb.fixed, comb.random = overallmeta$comb.random, 
             overall = TRUE, studlab = Overall$Country, plotwidth=unit(8, "cm"), xlim = c(0, 1), colgap.forest = "2.5cm", fontsize = 14)
      dev.off()
    }
    
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
    
    random_est_overall <- c(overallmeta$TE.random, overallmeta$lower.random, overallmeta$upper.random)
    
    random_est_overall_tr <- as.data.frame(t(backtransf(random_est_overall, sm = "PLOGIT")))
    
    
    # bind proportions for each continent into a singular dataframe
    proportion_table <- rbind(random_est_na_tr, random_est_eu_tr, random_est_af_tr,
                              random_est_la_tr, random_est_as_tr, random_est_oce_tr, random_est_overall_tr)
    
    # create column for continent
    proportion_table <- proportion_table |>
      mutate(
        Continent = c("North America", "Europe", "Africa",
                      "South America", "Asia", "Oceania", "Overall")) |>
      relocate(Continent, .before = V1)
    
    # rename columns
    colnames(proportion_table) <- c("Continent", "Percentage", "Lower", "Upper")
    
    # save proportion table for each antibiotic into their respective folder
    proportion_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", paste(newdir, "/proportion.csv", sep = ""), fsep = "/")
    write_csv(proportion_table, file = proportion_path)
    
    # round to 1 significant figure
    reference_proportion_table <- proportion_table |>
      mutate(across(c("Percentage", "Lower", "Upper"), ~ . * 100)) |>
      mutate(across(c("Percentage", "Lower", "Upper"), round, digits = 1))
    
    # paste proportion and CIs togther
    reference_proportion_table$proportion_ci <- paste0(reference_proportion_table$Percentage, " (", reference_proportion_table$Lower, ", ", reference_proportion_table$Upper, ")")
    
    # select for continent and proportion/CI columns
    reference_proportion_table <- reference_proportion_table |>
      select(c(Continent, proportion_ci))
    
    # rename proportion/CI column based on looped variable
    colnames(reference_proportion_table)[which(names(reference_proportion_table) == "reference_proportion_ci")] <- df_name
    
    # rename proportion table to looped variable
    assign(df_name, reference_proportion_table)
    
    analysis_proportion_table <- proportion_table |>
      mutate(
        period = p,
        antibiotic = i
      )
    
    total <- complete |>
      filter(!!sym((i)) %in% c("S", "I", "R")) |>
      group_by(Continent) |>
      tally() |>
      mutate(
        Continent = str_to_title(Continent)
      )
    
    analysis_proportion_table <- left_join(analysis_proportion_table, total)
    
    analysis_proportion_df <- rbind(analysis_proportion_table, analysis_proportion_df)
    
    
  }
  
  newdir <- paste0(maindir,"/complete proportion ",p)
  
  dir.create(newdir)
  
  ### compile random effects model proportion table
  # create list of proportion table dataframes for each antibiotic
  df_list <- list(WGS_PEN_SIR_Meningitis_proportion, WGS_AMO_SIR_proportion, WGS_CFX_SIR_proportion,
                  WGS_CFT_SIR_Meningitis_proportion, WGS_TAX_SIR_Meningitis_proportion, WGS_MER_SIR_proportion,
                  WGS_CHL_SIR_proportion, WGS_CLI_SIR_proportion, WGS_ERY_SIR_proportion,
                  WGS_TET_SIR_proportion, WGS_COT_SIR_proportion, WGS_RIF_SIR_proportion, WGS_KAN_SIR_proportion, MDR_proportion)
  
  # join all tables together
  complete_proportion_table <- df_list |>
    reduce(left_join, by = "Continent")
  
  colnames(complete_proportion_table) <- c("Continent", "Penicillin", "Amoxicillin", "Cefuroxime", 
                                           "Ceftriaxone", "Cefotaxime", "Meropenem", "Chloramphenicol", "Clindamycin", "Erythromycin", 
                                           "Tetracycline", "Co-trimoxazole", "Rifampin", "Kanamycin", "MDR")
  
  # save proportion table for each antibiotic into their respective folder
  complete_proportion_path <- file.path("/Users", "choong", "OneDrive - EMBLEBI", "Desktop", "GPS-AMR", "figures", "Proportion", paste(newdir,"/complete proportion.csv", sep = ""), fsep = "/")
  
  write_csv(complete_proportion_table, file = complete_proportion_path)
  
  
}

# relevel pcv period
analysis_proportion_df$period <- factor(analysis_proportion_df$period, levels = c("PREPCV13", "POSTPCV13"))

# reformat antibiotic names
analysis_proportion_df <- left_join(analysis_proportion_df, antibiotic_tidying, by = "antibiotic") |>
  select(-c(antibiotic, antibiotic_tidy))

# relevel antibiotic names
analysis_proportion_df$antibiotic_short <- factor(analysis_proportion_df$antibiotic_short, levels = c("PEN", "AMO", "CFX", "CFT", "TAX", "MER", "CHL", "CLI", "ERY", "TET", "COT", "RIF", "KAN", "MDR"))


for(t in continent){
  
  # summary figure for proportion analysis
  proportion_plot <- analysis_proportion_df |>
    filter(Continent %in% t) |>
    ggplot(aes(x = period, y = Percentage*100, group = antibiotic_short)) +
    geom_col(width = 0.8, stat = "identity", position = position_dodge(width = 0.8), aes(fill = period)) +
    geom_errorbar(aes(ymin = Lower*100, ymax = Upper*100), position = position_dodge(width = 0.8)) +
    facet_wrap(~antibiotic_short, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"),
                      labels = c("Pre-PCV13", "Post-PCV13")) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      x = "Antibiotic",
      y = "Percentage",
      fill = "Period",
      title = t
    ) +
    theme(
      title = element_text(size = 28),
      axis.title = element_text(size = 28),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 26),
      legend.title = element_text(size = 26),
      legend.text = element_text(size = 24),
      legend.key.size = unit(1, "cm"),
      legend.background = element_blank(),
      strip.placement='outside',
      strip.background.x=element_blank(),
      strip.text=element_text(size=28,color="black",face="bold"),
      panel.spacing.x=unit(1,"pt"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  ggsave(proportion_plot, filename = paste0("Proportion/overview figure/",t,".svg"), dpi = 320, width = 43.58, height = 18.91, units = "cm")
  
  
}

write_csv(analysis_proportion_df, file = "Proportion/overview figure/proportion reference.csv")


#----------------------------------------- GPSCs driving resistance post PCV10/13 ----------------

setwd("C:/Users/choong/OneDrive - EMBLEBI/Desktop/GPS-AMR/")

post_pcv13 <- metadata |>
  filter(Vaccine_period_clean == "POSTPCV13")
GPSCs <- unique(post_pcv13$GPSC)

continents <- unique(post_pcv13$Continent)

gpsc_colors_vector <- setNames(post_pcv13$GPSC__colour, post_pcv13$GPSC)

betalactam_results <- data.frame("GPSC" = character(), 
                                 "Antimicrobial" = character(), 
                                 "Continent" = character(), 
                                 "Proportion_R" = numeric())

for (GPSC in GPSCs) {
  for (continent in continents) {
    
    temp_post_pcv13 <- filter(post_pcv13,
                              GPSC == GPSC,
                              Continent == continent)
    
    betalactam_results <- rbind(
      betalactam_results,
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "PEN",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_PEN_SIR_Meningitis == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_PEN_SIR_Meningitis == "R", na.rm=TRUE)
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "AMO",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_AMO_SIR == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_AMO_SIR == "R", na.rm=TRUE)
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "CFX",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CFX_SIR == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CFX_SIR == "R", na.rm=TRUE)
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "CFT",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CFT_SIR_Meningitis == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CFT_SIR_Meningitis == "R", na.rm=TRUE)
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "TAX",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_TAX_SIR_Meningitis == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_TAX_SIR_Meningitis == "R", na.rm=TRUE)
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "MER",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_MER_SIR == "R", na.rm=TRUE) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_MER_SIR == "R", na.rm=TRUE)
      )
    )
  }
}

betalactam_gpscs <- betalactam_results %>%
  group_by(Continent, Antimicrobial) %>%
  top_n(3, Proportion_R) %>%
  ungroup()

betalactam_gpscs <- filter(betalactam_gpscs, Proportion_R != 0)

betalactam_gpscs <- betalactam_gpscs %>%
  mutate(Antimicrobial = fct_relevel(Antimicrobial, "PEN", "AMO", "CFX", "CFT", "TAX", "MER"))

beta_lactam_plot <- ggplot(betalactam_gpscs, aes(x = GPSC, y = Proportion_R, fill = GPSC)) + 
  geom_col() + 
  theme_minimal() + 
  scale_fill_manual(values = gpsc_colors_vector) +
  facet_grid(Continent ~ Antimicrobial, scales = "free_x", space = "free_x") +
  ylab("Proportion of Resistant Isolates per Continent") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))



## Non beta lactam plot
nonbl_results <- data.frame("GPSC" = character(), 
                            "Antimicrobial" = character(), 
                            "Continent" = character(), 
                            "Proportion_R" = numeric())


for (GPSC in GPSCs) {
  for (continent in continents) {
    
    temp_post_pcv13 <- filter(post_pcv13,
                              GPSC == GPSC,
                              Continent == continent)
    
    nonbl_results <- rbind(
      nonbl_results,
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "CLI",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CLI_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CLI_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "ERY",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_ERY_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_ERY_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "COT",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_COT_SIR %in% c("R")) / sum(post_pcv13$Continent == continent & post_pcv13$WGS_COT_SIR %in% c("R"))
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "TET",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_TET_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_TET_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "CHL",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CHL_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CHL_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "RIF",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CHL_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CHL_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "KAN",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$WGS_CHL_SIR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$WGS_CHL_SIR == "R")
      ),
      data.frame(
        GPSC = GPSC,
        Antimicrobial = "MDR",
        Continent = continent,
        Proportion_R = sum(post_pcv13$Continent == continent & post_pcv13$GPSC == GPSC & post_pcv13$MDR == "R") / sum(post_pcv13$Continent == continent & post_pcv13$MDR == "R")
      )
    )
  }
}

nonbl_gpscs <- nonbl_results %>%
  group_by(Continent, Antimicrobial) %>%
  top_n(3, Proportion_R) %>%
  ungroup()

nonbl_gpscs <- filter(nonbl_gpscs, Proportion_R != 0)

nonbl_gpscs <- nonbl_gpscs %>%
  mutate(Antimicrobial = fct_relevel(Antimicrobial, "CHL", "CLI", "ERY", "TET", "COT", "RIF", "KAN", "MDR"))

nonbl_plot <- ggplot(nonbl_gpscs, aes(x = GPSC, y = Proportion_R, fill = GPSC)) + 
  geom_col() + 
  theme_minimal() + 
  scale_fill_manual(values = gpsc_colors_vector) +
  facet_grid(Continent ~ Antimicrobial, scales = "free_x", space = "free_x") +
  ylab("Proportion of Resistant Isolates per Continent") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

## into one figure
gpsc_plot <- ggarrange(beta_lactam_plot, nonbl_plot, nrow = 2, labels=c("(a)", "(b)"))

ggsave(gpsc_plot, file ="figures/gpsc.svg", dpi = 320, width = 32, height = 36, units = "cm")


#------------------------------ Potential vaccine coverage and PCV13 serotype persistence post vaccine implementation ----------------------------

setwd("C:/Users/choong/OneDrive - EMBLEBI/Desktop/GPS-AMR/")

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

Continent = c("North America", "Europe", "Africa",
              "South America", "Asia", "Oceania")

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



# if else statement to class serotype into vaccine coverage
# show colours only for vaccine types, grey for all NVTs






vaccine_table_unique_serotypes <- vaccine_table |>
  filter(!is.na(vaccine_unique) & !vaccine %in% "PCV25_IVT") |>
  select(serotype, vaccine_unique) |>
  rename(
    serotype_clean = serotype
  )

serotype_coverage_vaccine <- left_join(serotype_coverage_df, vaccine_table_unique_serotypes, by = "serotype_clean") 

write_csv(serotype_coverage_vaccine, file = "coverage reference.csv")


# mutate columns for unique serotypes for each PCV and non-PCV13 serotypes
# setting PCV7 nad PCV10 (GSK) serotypes as PCV13 serotypes as PCV7/10 are not of interest
serotype_coverage_vaccine_unique <- serotype_coverage_vaccine |>
  mutate(
    unique_serotype = ifelse(is.na(vaccine_unique), "NVT", 
                             ifelse(vaccine_unique %in% c("PCV7_PFIZER", "PCV10_GSK"), "PCV13_PFIZER" , vaccine_unique)),
    PCV13_vaccine = ifelse(serotype_clean %in% pcv13_serotypes, serotype_clean, "NVT"),
    PCV25_vaccine = ifelse(serotype_clean %in% pcv25_serotypes, "PCV25", "NVT")
  )

serotype_coverage_vaccine_unique_clean <- left_join(serotype_coverage_vaccine_unique, antibiotic_tidying, by = "antibiotic")

serotype_coverage_vaccine_unique_clean$antibiotic_short <- factor(serotype_coverage_vaccine_unique_clean$antibiotic_short, 
                                                                  levels = c("PEN", "AMO", "CFX", "CFT", "TAX", 
                                                                             "MER", "CHL", "CLI", "ERY", "TET", 
                                                                             "COT", "RIF", "KAN", "MDR"))


# set colours
pcv13_colour <- serotype_colours |>
  rename(
    PCV13_vaccine = serotype_clean,
    PCV13_colour = In_silico_serotype__colour
  )


vaccine_colour <- serotype_colours |>
  rename(
    unique_serotype = serotype_clean,
    vaccine_colour = In_silico_serotype__colour
  )


serotype_coverage_vaccine_unique_clean_coloured <- left_join(serotype_coverage_vaccine_unique_clean, pcv13_colour, by = "PCV13_vaccine")

serotype_coverage_vaccine_unique_clean_coloured <- left_join(serotype_coverage_vaccine_unique_clean_coloured, vaccine_colour, by = "unique_serotype")

PCV13_colors_vector <- setNames(serotype_coverage_vaccine_unique_clean_coloured$PCV13_colour, 
                                serotype_coverage_vaccine_unique_clean_coloured$PCV13_vaccine)


vaccine_colors_vector <- setNames(serotype_coverage_vaccine_unique_clean_coloured$vaccine_colour,
                                  serotype_coverage_vaccine_unique_clean_coloured$unique_serotype)




serotype_coverage_vaccine_unique_clean_coloured$PCV13_vaccine <- factor(serotype_coverage_vaccine_unique_clean_coloured$PCV13_vaccine,
                                                                        levels = c("NVT", "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A", "19F", "23F"))




pcv13_coverage_figure <- serotype_coverage_vaccine_unique_clean_coloured |>
  ggplot(aes(x = antibiotic_short, y = n, fill = PCV13_vaccine, colour = PCV13_vaccine)) +
  geom_bar(position = position_stack(), stat = "identity") +
  facet_wrap(~Continent, ncol = 3, scales = "free_y")+
  scale_fill_manual(name = "PCV13 (Pfizer) \nserotypes",
                    values = PCV13_colors_vector) +
  scale_colour_manual(values = PCV13_colors_vector,
                      guide = "none") +
  labs(
    x = "Antibiotic",
    y = "Number of isolates",
    title = "PCV13 (Pfizer) serotypes driving AMR post-PCV13 implementation in each continent"
  ) +
  theme(
    title = element_text(size = 24),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 22), 
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    strip.text = element_text(size = 16),
    legend.background = element_blank()
  )

ggsave(pcv13_coverage_figure, file = "figures/vaccine coverage figure post pcv13.svg", dpi = 320, width = 48, height = 36, units = "cm")

serotype_coverage_vaccine_unique_clean_coloured$unique_serotype <- factor(serotype_coverage_vaccine_unique_clean_coloured$unique_serotype,
                                                                          levels = c("NVT", "PCV31_VAXCYTE", "PCV24_GSK", "PCV21_SK_SANOFI", "PCV20_PFIZER", "PCV15_MERCK", "PCV13_PFIZER"))

potential_coverage_plot <- serotype_coverage_vaccine_unique_clean_coloured |>
  ggplot(aes(x = antibiotic_short, y = n, fill = unique_serotype, colour = unique_serotype)) +
  geom_bar(position = position_stack(), stat = "identity") +
  facet_wrap(~Continent, ncol = 3, scales = "free_y")+
  scale_fill_manual(name = "Vaccine \nserotypes",
                    values = vaccine_colors_vector,
                    labels = c("NVT", "PCV31 (Vaxcyte)", "PCV24 (GSK)", "PCV21 (SK/Sanofi)", "PCV20 (Pfizer)", "PCV15 (Merck)", "PCV13 (Pfizer)")) +
  scale_colour_manual(values = vaccine_colors_vector,
                      guide = "none") +
  labs(
    x = "Antibiotic",
    y = "Number of isolates",
    title = "Potential vaccine coverage of AMR isolates post-PCV13 implementation in each continent"
  ) +
  theme(
    title = element_text(size = 24),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 22), 
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    strip.text = element_text(size = 16),
    legend.background = element_blank()
  )

ggsave(potential_coverage_plot, file = "figures/potential coverage figure.svg", dpi = 320, width = 48, height = 36, units = "cm")

### separate figure for PCV25

pcv25_potential_coverage_plot <- serotype_coverage_vaccine_unique_clean_coloured |>
  mutate(
    PCV25_colour = ifelse(PCV25_vaccine %in% "NVT", "#808080", "#5D3A9B")
  ) |>
  ggplot(aes(x = antibiotic_short, y = n, fill = PCV25_vaccine, colour = PCV25_vaccine)) +
  geom_bar(position = position_stack(), stat = "identity") +
  facet_wrap(~Continent, ncol = 3, scales = "free_y") +
  scale_fill_manual(name = "Vaccine \nserotypes",
                    values = c("#808080", "#5D3A9B"),
                    labels = c("NVT", "PCV25 (Inventprise)")) +
  scale_colour_manual(values = c("#808080", "#5D3A9B"),
                      guide = "none") +
  labs(
    x = "Antibiotic",
    y = "Number of isolates",
    title = "Potential vaccine coverage of AMR isolates post-PCV13 implementation in each continent"
  ) +
  theme(
    title = element_blank(),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 22), 
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    strip.text = element_text(size = 16),
    legend.background = element_blank()
  )

ggsave(pcv25_potential_coverage_plot, file = "figures/pcv25 potential coverage figure.svg", dpi = 320, width = 48, height = 36, units = "cm")


### potential coverage table

coverage_table <- data.frame()


cols <- c("PCV31_VAXCYTE", "PCV24_GSK", "PCV21_SK_SANOFI", "PCV20_PFIZER", "PCV15_MERCK", "PCV13_PFIZER") |>
  map_dfr( ~tibble(!!.x := logical() ) )

for(i in Continent){
  coverage_raw <- serotype_coverage_vaccine_unique_clean |>
    filter(Continent %in% i) |>
    select(-Continent) |>
    group_by(antibiotic_tidy, unique_serotype) |>
    summarise(
      n = sum(n)
    ) |>
    pivot_wider(names_from = unique_serotype, values_from = n) |>
    mutate(
      continent = i
    )
  
  coverage_raw <- bind_rows(coverage_raw, cols)
  
  
  coverage_raw[is.na(coverage_raw)] <- 0
  
  coverage_raw <- coverage_raw |>
    relocate(PCV21_SK_SANOFI, .after = NVT) |>
    relocate(PCV24_GSK, .after = PCV21_SK_SANOFI)
  
  coverage_clean <- coverage_raw |>
    mutate(
      total = rowSums(across(where(is.numeric))),
      pcv13_percentage = paste(round(100 * PCV13_PFIZER/total, 2), "%"),
      pcv13_clean = paste0(pcv13_percentage, " (", PCV13_PFIZER, ")"),
      total_pcv15 = PCV13_PFIZER + PCV15_MERCK,
      pcv15_percentage = paste(round(100 * total_pcv15/total, 2), "%"),
      pcv15_clean = paste0(pcv15_percentage, " (", total_pcv15, ")"), 
      total_pcv20 = total_pcv15 + PCV20_PFIZER,
      pcv20_percentage = paste(round(100 * total_pcv20/total, 2), "%"),
      pcv20_clean = paste0(pcv20_percentage, " (", total_pcv20, ")"),
      total_pcv21 = total_pcv20 + PCV21_SK_SANOFI,
      pcv21_percentage = paste(round(100 * total_pcv21/total, 2), "%"),
      pcv21_clean = paste0(pcv21_percentage, " (", total_pcv21, ")"),
      total_pcv24 = total_pcv21 + PCV24_GSK,
      pcv24_percentage = paste(round(100 * total_pcv24/total, 2), "%"),
      pcv24_clean = paste0(pcv24_percentage, " (", total_pcv24, ")"),
      total_pcv31 = total_pcv24 + PCV31_VAXCYTE,
      pcv31_percentage = paste(round(100 * total_pcv31/total, 2), "%"),
      pcv31_clean = paste0(pcv31_percentage, " (", total_pcv31, ")"),
    ) |>
    select(continent, antibiotic_tidy, pcv13_clean, pcv15_clean, pcv20_clean, pcv21_clean, pcv24_clean, pcv31_clean, total)
  
  coverage_table <- rbind(coverage_clean, coverage_table)
}


colnames(coverage_table) <- c("Continent", "Antibiotic", "PCV13 (Pfizer)", "PCV15 (Merck)", "PCV20 (Pfizer)", "PCV21 (SK/Sanofi)", "PCV24 (GSK)", "PCV31 (Vaxcyte)", "Total")


write_csv(coverage_table, file = "figures/coverage table.csv")


cols_pcv25 <- c("PCV25") |>
  map_dfr( ~tibble(!!.x := logical() ) )

coverage_table_pcv25 <- data.frame()

for(i in Continent){
  coverage_raw_pcv25 <- serotype_coverage_vaccine_unique_clean |>
    filter(Continent %in% i) |>
    select(-Continent) |>
    group_by(antibiotic_tidy, PCV25_vaccine) |>
    summarise(
      n = sum(n)
    ) |>
    pivot_wider(names_from = PCV25_vaccine, values_from = n) |>
    mutate(
      continent = i
    ) 
  
  coverage_raw_pcv25 <- bind_rows(coverage_raw_pcv25, cols_pcv25)
  
  
  coverage_raw_pcv25[is.na(coverage_raw_pcv25)] <- 0
  
  
  coverage_clean_pcv25 <- coverage_raw_pcv25 |>
    mutate(
      total = rowSums(across(where(is.numeric))),
      pcv25_percentage = paste(round(100 * PCV25/total, 2), "%"),
      pcv25_clean = paste0(pcv25_percentage, " (", PCV25, ")")
    ) |>
    select(continent, antibiotic_tidy, pcv25_clean, total)
  
  coverage_table_pcv25 <- rbind(coverage_clean_pcv25, coverage_table_pcv25)
}

colnames(coverage_table_pcv25) <- c("Continent", "Antibiotic", "PCV25 (Inventprise)", "Total")

write_csv(coverage_table_pcv25, file = "figures/coverage table pcv25.csv")


nvt_serotypes <- serotype_coverage_vaccine_unique_clean |>
  filter(unique_serotype %in% "NVT") |>
  group_by(antibiotic_short, Continent) |>
  select(Continent, antibiotic_short, serotype_clean, n) |> 
  arrange(desc(n)) 

write_csv(nvt_serotypes, file = "figures/nvt serotypes.csv")

