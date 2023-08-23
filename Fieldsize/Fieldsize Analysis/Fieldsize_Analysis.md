---
title: "German-wide Biodiversity Metrics: Analysis of Field sizes"
author: "Jannes Uhlott @JKI"
date: "2023-08-04"
output: pdf_document
---

# Introduction

Calculation of Biodiversity metrics based on crop type classification data (Preidl et al. 2023, *not published*) which are mapped by majority voting to geometry segments (Tetteh et al. 2023, *not published*). The calculation of the field size is done by four different methods (mean, mean_cut, w_mean and c_mean) per hexagon as reference area. 

```{r, packages, warning = FALSE, message = FALSE, echo=FALSE}
library(sp)
library(sf)
library(ggplot2)
library("GGally") 
library(tidyverse)
library(tidyr)
library(units)
library(vegan)
library(purrr)
library(reshape)
library(dplyr)
source("~/MonViA_Indikatoren/JU_MonViA/Skripte/Vektor/Functions_Indikatoren.R")
```


# Input
All used Packages are presented at the end of document. 
Load field size dataset and hexagon dataset.
```{r, data_definitions, warning = FALSE, message = FALSE, results='hide', echo = FALSE}
wd <- "~/Daten/Biodiversitätsmetriken/"
CTC <- "_PRE"
fd <- "D"
year <- "2019"
fd_year <- paste0(fd,"_", year)
short = "PRE" 

mask <- "_crop" # possible: Agrarmask and Cropmask [_Crop_, _Agrar]
savename <- paste0(fd_year, CTC, mask)

## input files
# shp_data <- st_read(paste0(wd, fd_year, "_PRE_Seg_crop.shp"))
hexagonfile <- paste0("~/Daten/Hexagon_DE_UTM32/", fd, "/", fd, "_Hexagone_filtered.shp")


fieldsizes_file <- paste0(wd, "Feldgröße/", savename, "_fieldsizes.shp")

```

## Data

```{r, data_input, warning = FALSE, message = FALSE, results='hide'}

hexagon_data <- st_read(hexagonfile)
fieldsizes <- st_read(fieldsizes_file)

```
The data set *fieldsize* contains 6 columns and contain 354252 hexagons. The first column describes the HexagonID followed by four columns as results of the different calculation methods (mean, mean_cut, w_mean and c-mean). The last columen contains the polygon geometry infomration of each hexagon.

```{r, data_details, warning = FALSE, message = FALSE, echo = FALSE}

knitr::kable(fieldsizes[1:3, ], caption = "Details of dataset field size")

```
\newpage
## Data Statistics 
Calculation of mean, standard deviation, maximum and minimum for every column (mean, mean_cut, w_mean, c_mean). In addition the coordinate reference system (crs) is given. 

```{r, data_statistics, warning = FALSE, message = FALSE, results='asis'}

# mean, sd, max and min of fieldsizes

mean_results <- sapply(X=st_drop_geometry(fieldsizes[2:6]), FUN = mean, na.rm = T)
nice_mean_results <- broom::tidy(mean_results)

sd_results <- sapply(X=st_drop_geometry(fieldsizes[2:6]), FUN = sd, na.rm = T)
nice_sd_results <- broom::tidy(sd_results)

max_results <- sapply(X=st_drop_geometry(fieldsizes[2:6]), FUN = max, na.rm = T)
nice_max_results <- broom::tidy(max_results)

min_results <- sapply(X=st_drop_geometry(fieldsizes[2:6]), FUN = min, na.rm = T)
nice_min_results <- broom::tidy(min_results)

statistic_results <- data_frame(name = nice_mean_results$names, 
                                mean = nice_mean_results$x, 
                                sd = nice_sd_results$x, 
                                max = nice_max_results$x, 
                                min = nice_min_results$x)
# Geographic extend
crs <- st_crs(fieldsizes)
print(paste0("crs fieldsize = ", crs$input))
```

```{r, result transformation, warning = FALSE, message = FALSE, results='hide', echo = FALSE}
results <- as_data_frame(t(statistic_results))

statistic_results <- data_frame(name = statistic_results$name,
                                mean = results$V1[2:5], 
                                sd = results$V2[2:5], 
                                max = results$V3[2:5], 
                                min = results$V4[2:5])
```

```{r, statistic_details, warning = FALSE, message = FALSE, results='asis', echo = FALSE}
knitr::kable(statistic_results, caption = "Statistics of fieldsize dataset")
```
\newpage

## Pearson Korrelation (r) of all Fieldsizes
```{r, korr, warning = FALSE, message = FALSE, results='hide'}

subset <- st_drop_geometry(fieldsizes) %>% 
  dplyr::select(mean, mean_cut, c_mean, w_mean)

corr_tab_pair <- cor(subset, use = "pairwise.complete.obs")
```

```{r, corr_details, warning = FALSE, message = FALSE, results='asis', echo = FALSE}

knitr::kable(corr_tab_pair, caption = "Pearson Korrelation (r) of all field sizes methods")
```
## Scatter and correlation plot
Because of the significant difference in the field size based on historical reasons the data set is separated in East and West Germany. The SN_L column from the hexagon data set contain the number of federal states. The numbers 11 to 16 represents the federal states in the East of Germany. 
```{r, subset, warning = FALSE, message = FALSE, results='hide'}

raw_hexagon_data <- st_drop_geometry(hexagon_data) %>% 
  dplyr::select(SN_L, HexagonID)

subset_hex <- st_drop_geometry(fieldsizes) 

subset_SN <- left_join(subset_hex, raw_hexagon_data)

subset_E <- subset_SN %>% # define East subset
    filter(SN_L %in% c(11, 12, 13, 14, 15, 16)) %>% 
    mutate(SN_Code = "E")

subset_W <- subset_SN %>% # define West subset
    filter(! SN_L %in% c(11, 12, 13, 14, 15, 16)) %>% 
    mutate(SN_Code = "W") 

subset_EW <- rbind(subset_E, # combine East and West subset
                   subset_W)

subset <- data_frame(mean = subset_EW$mean, 
                     mean_cut = subset_EW$mean_cut, 
                     w_mean = subset_EW$w_mean, 
                     c_mean = subset_EW$c_mean,
                     SN = subset_EW$SN_Code)
```

```{r, scatter, warning = FALSE, message = FALSE, results='hide', echo = FALSE, fig.width=6}
gg_D <- ggpairs(subset, columns = 1:4, aes(color = SN), title = "Fieldsize") #+
        # scale_x_continuous(limits = c(0,3))+
        # scale_y_continuous(limits = c(0,3))
gg_D
# name = "_subset_D"
# ggsave(gg_D, filename = paste0(wd, "Feldgröße/", savename, name ,"_corr_scaled.png"), dpi=300, width = 6, height = 5)
```
\newpage

# Packages

```{r, end, warning = FALSE, message = FALSE, results='asis', echo = FALSE}
# print(sessionInfo())
sessionInfo()
# installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
```
