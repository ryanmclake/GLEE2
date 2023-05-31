
library(sp)
library(lattice)
library(magrittr)
library(tibble)
library(leaflet)
library(leaflet.extras2)
library(sf)
library(lubridate)
remotes::install_github("r-spatial/leafpop")
library(leafpop)
remotes::install_github("r-spatial/mapview")
library(mapview)

area_hexes_avg


mapview(area_hexes_avg, popup = popupGraph(p, type = "png", width = 300, height = 300))


mapview(area_hexes_avg, popup = popupGraph(p, type = "png", width = 300, height = 300))

