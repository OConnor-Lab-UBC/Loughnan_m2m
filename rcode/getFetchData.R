# January 7, 2026 by D. Loughnan
# Aim of this code is to calculate Fetch values for the m2m sites
# Code is based off of code provided by N. Knight

rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(tidyverse)
# library(janitor)
library(sf)
library(waver)
library(raster)

setwd("~/Documents/github/Loughnan_m2m")

latLong <- read.csv("input/m2mSiteLatLong.csv")

projcrs <- "epsg:32198" #NAD83 Québec Lambert recommended by Mel

sites_sf <- st_as_sf(x = latLong,                         
                     coords = c("long", "lat"),
                     crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sites_sf <- st_transform(sites_sf, crs = "epsg:2952") # Canada on shore and offshore https://epsg.io/2952

#-Shoreline downloaded from https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21

shoreline <- read_sf("input/fetchData/lpr_000b21a_e/lpr_000b21a_e.shp") # need all the files in the folder for this to work, not just .shp
shoreline <- st_transform(shoreline, crs = "epsg:2952")
shoreline <- shoreline[shoreline$PRNAME %in% c("British Columbia / Colombie-Britannique"),]

compareCRS(shoreline, sites_sf)

ggplot(shoreline) +
  geom_sf() 

#-Calculate fetch-#

JB_fetch <- fetch_len_multi(
  sites_sf,
  bearings = seq(from = 0, to = 350, by = 10),
  shoreline,
  dmax = 10000, # cut off for sites with open ocean
  spread = 0,
  method = "btree"
)

# convention = sum or average---> often ppl get wind frequency and speed to get relative exposure index = better metric bc it accounts for wind direction
# Bc wind data: weatherCan---all wind stations with weather data


#-Convert to km-#

fetch <- JB_fetch %>%
  as_tibble() %>%
  add_column(site = sites_sf$site) %>%
  relocate(site, .before = '30') %>%
  pivot_longer(!site, names_to = "direction", values_to = "fetch") %>%
  mutate(direction = as.numeric(direction)) %>%
  mutate(fetch = fetch/1000)

write.csv(fetch, "fetchDataCalvert.csv", row.names = FALSE)

fetch <- read.csv("fetchDataCalvert.csv")

ggplot(fetch, aes(x= fetch)) +
  geom_histogram()  + 
  facet_wrap(~ site)
