---
title: "German-wide Biodiversity Metrics: Field size"
author: "Jannes Uhlott @JKI"
date: "2023-08-04"
output: pdf_document
---

# Introduction

Calculation of Biodiversity metrics based on crop type classification data (Preidl et al. 2023, *not published*) based on Preidl et al. (2020) which are mapped by majority voting to geometry segments (Tetteh et al. 2023, *not published*) based on Tetteh et al. (2021). The field size gets calculated by four different methods (mean, mean_cut, w_mean and c_mean) per hexagon as reference area. 

```{r, packages, warning = FALSE, message = FALSE, echo=FALSE}
library(sp)
library(sf)
library(ggplot2)
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
All used Packages and functions are presented at the end of document. 
Load data which is already intersected with the hexagon layer. Polygon areas should have calculated before intersection.

## Data

```{r, data_input, warning = FALSE, message = FALSE, results='hide'}

hexagon_data <- st_read(hexagonfile)
shp_data <- st_read(shp_data_file)

```

The data are filtered by crop classes before. The agrarmask (0-19) contain all crop types and grassland classes. The Cropmask only contain the crop type classes (0-18). The following lines check the data for the right codes and the polygon area gets calculated. 
```{r, Codes, warning = FALSE, message = FALSE, results='asis'}

codes_PRE <- sort(unique(shp_data$PRE))
print(codes_PRE)

polygon_cutarea <- get_polygon_area(shp_data)

polygon_details <- data_frame(HexagonID = polygon_cutarea$HexagonID, 
                   code = polygon_cutarea$PRE,      
                   area = polygon_cutarea$area,
                   cutarea = polygon_cutarea$area.1)

raw_polygon_details <- drop_units(polygon_details)
```

\newpage
# Field Size

There a different ways to calculate the fieldsize. There is the "normal" mean where is to decide, whether the whole polygon geometry (mean with polygon total area) should be included or only the polygon part which are included into the specific hexagon (mean with polygon cut, mean_cut). Furthermore the mean of the whole polygon area can be calculated but only in hexagon where the centroid of the polygone is located (centroid mean, c_mean). Also an area weighted mean (w_mean) can be calculated where the total polygon area is used in the calculation but its weight by its cutted area in the specific hexagon. Each method have advantages and disadvantages which are described in more details in the documentation. 


## mean polygon with whole area (mean)

```{r, fieldsize_mean_whole, warning = FALSE, message = FALSE, results='hide'}
mean_data <- get_meanarea_per_hexagon(hexagon_data, raw_polygon_details)

raw_mean_data <- mean_data %>% 
  dplyr::select(HexagonID, meanarea)

# st_write(raw_mean_data, paste0(wd, "Feldgröße/", savename, "_mean_total.shp"), delete_dsn = TRUE)

```

## mean polygon with cutted area (mean_cut)
For the cutted mean only the partial area of the polygons which intersects the hexagon (cutarea) get used for the mean calculation. 

```{r, fieldsize_mean_cut, warning = FALSE, message = FALSE, results='hide'}
mean_cut_polygon_details <- data_frame(HexagonID = raw_polygon_details$HexagonID, 
                   area = raw_polygon_details$cutarea) # CUT polygon area

mean_cut_data <- get_meanarea_per_hexagon(hexagon_data, mean_cut_polygon_details)

raw_mean_cut_data <- mean_cut_data %>% 
  dplyr::select(HexagonID, mean_cut = meanarea)

# st_write(raw_mean_cut_data, paste0(wd, "Feldgröße/", savename, "_mean_CUT.shp"), delete_dsn = TRUE)
```

## area weighted mean (weight total polygon area by area per hexagon, w_mean)

```{r, fieldsize_area_weighted, warning = FALSE, message = FALSE, results='hide'}
w_mean_data <- get_weightedmean_per_hexagon(hexagon_data, raw_polygon_details)

raw_w_mean_data <- w_mean_data %>% 
  dplyr::select(HexagonID, w_mean)

# st_write(raw_wmean_data, paste0(wd, "Feldgröße/", savename, "_wmean.shp"), delete_dsn = TRUE)

```

## centroid mean
For the centroid mean the centroid of all polygons got calculated and intersected to the hexagon layer. 

```{r, centroid_data, warning = FALSE, message = FALSE, results='hide', echo = FALSE}
centroid_file <- paste0(wd, fd_year, "_PRE_Seg" , mask, "_centroids.shp")
```

```{r, mean_centroid, results='hide', message = FALSE}
shp_centroids <- st_read(centroid_file)

raw_polygon_details <- data_frame(HexagonID = shp_centroids$HexagonID, 
                   code = shp_centroids$PRE,      
                   area = shp_centroids$area)

c_mean_data <- get_meanarea_per_hexagon(hexagon_data, raw_polygon_details)

raw_c_mean_data <- c_mean_data %>% 
  dplyr::select(HexagonID, meanarea)

# st_write(raw_c_mean_data, paste0(wd, "Feldgröße/", savename, "_centroid.shp"), delete_dsn = TRUE)
```

# Summary to one dataset

```{r, mean_summary_data, results='hide', message = FALSE,}
# definde dataframes (only one dataset are allowed to include geometries)

raw_mean_data <- data_frame(HexagonID = mean_data$HexagonID,
                            mean = mean_data$meanarea,
                            geometry = mean_data$geometry)
raw_mean_cut_data <- data_frame(HexagonID = mean_cut_data$HexagonID, 
                                mean_cut = mean_cut_data$meanarea)
raw_wmean_data <- data_frame(HexagonID = w_mean_data$HexagonID, 
                                w_mean = w_mean_data$w_mean)
raw_c_mean_data <- data_frame(HexagonID = c_mean_data$HexagonID, 
                                c_mean = c_mean_data$meanarea)

# Summarize all fieldsize data
data_mean_1 <- left_join(x= raw_mean_data, y=raw_mean_cut_data)
data_mean_2 <- left_join(data_mean_1, raw_wmean_data)
data_mean_3 <- left_join(data_mean_2, raw_c_mean_data)

fieldsizes <- data_mean_3
# st_write(fieldsizes, paste0(wd, "Feldgröße/", savename, "_fieldsizes.shp"))
```

```{r, field_size_details, warning = FALSE, message = FALSE, results='asis', echo = FALSE}
knitr::kable(fieldsizes[1:3, ], caption = "Details of dataset field size")
```

# Functions

## get polygon area

```{r, get_polygon_area, warning = FALSE, message = FALSE, results='hide'}
## Input: intersected with hexid, VEG
## Do: Calculate area for each polygon
## Output: intersected with area for each polygon
## Packages: sf, sp, units, lwgeom

get_polygon_area <- function(data) {
  
  ## calculate area for every patch ##
  area <- data %>%
    st_geometry(.) %>% 
    st_area(.) 
    
  area <- units::set_units(x = area, value = km^2)
  
  ## add calculations to intersected 
  data_area <- cbind(data, area)
  
  return(data_area)
}

```

## get mean area per hexagon

```{r, get_mean_area, warning = FALSE, message = FALSE, results='hide'}
### get_meanarea_per_hexagon
# Input: hexagon_data, raw_polygon_details [HexagonID, fid, area]
# Do: Calculates mean of fieldsize per Hexagon 
# Output: fieldsize_data [HexagonID, mean area, geometry]
# Packages: sf, tidyverse

get_meanarea_per_hexagon <- function(hexagon_data, raw_polygon_details){
  
  meanarea_output <- raw_polygon_details %>% 
    group_by(HexagonID) %>% 
    filter(!is.na(area)) %>% 
    summarise(meanarea = mean(area))
  
  meanarea_data <- left_join(x= hexagon_data, y=meanarea_output)
  
  return(meanarea_data)
}

```

## get weighted mean per hexagon

```{r, get_w-mean, warning = FALSE, message = FALSE, results='hide'}
# Input: polygon_cutarea [HexagonID, area, cutarea] 
  # with area as total area of polygon before intersection 
  # and cutarea as part o the area inside the hexagon. 
# Do: Calculates weighted mean fieldsize per Hexagon 
# Output: weightedmean per hexagon [HexagonID, weighted_mean, hexagongeometry]
# Packages: sf, tidyverse

get_weightedmean_per_hexagon <- function(hexagon_data, data) {
  
  weightedmean_output <- data %>% 
    group_by(HexagonID) %>% 
    filter(!is.na(area)) %>% 
    mutate(total_hex_sum = sum(cutarea)) %>% 
    mutate(weight = (cutarea/total_hex_sum)) %>% 
    summarise(w_mean = weighted.mean(area, weight))

  weightedmean_data <- left_join(x= hexagon_data, y=weightedmean_output)
  
  return (weightedmean_data)
}
```
# Liturature 

Preidl, S., Lange, M., & Doktor, D. (2020). Introducing apic for regionalised land cover mapping on the national scale using sentinel-2a imagery. Remote Sensing of Environment, 240, 111673. https://doi.org/10.1016/j.rse.2020.111673 

Preidl, S., Lange, M., & Doktor, D. (2022). Land cover classification map of germany’s agricultural area
based on sentinel-2a data from 2017-2019: Not published.

Tetteh, G. O., Gocht, A., Erasmi, S., Schwieder, M., & Conrad, C. (2021). Evaluation of sentinel-1 and
sentinel-2 feature sets for delineating agricultural fields in heterogeneous landscapes. IEEE
Access, 9, 116702–116719. https://doi.org/10.1109/ACCESS.2021.3105903  

Tetteh, G. O., Schwieder, M., Erasmi, S., & Gocht, A. (2023). Delineated agricultural fields from sentinel-
1 and sentinel-2 feature sets: Not published.

\newpage 
# Packages

```{r, packages_info, warning = FALSE, message = FALSE, results='asis', echo = FALSE}
sessionInfo()
```

