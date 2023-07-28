#read a few packages
library(ncdf4)
library(sf)
library(lattice)
library(vroom)
library(tidyr)

#specify netcdf file name
ncpath <- "/central/groups/carnegie_poc/rmcclure/project-GLEE/data/HWSD_1247/data/"
ncname <- "AWT_T_SOC"  
ncfname <- paste(ncpath, ncname, ".nc4", sep="")

# read the netcdf file
ncin <- nc_open(ncfname)

# extract the specific names and arrays of interest from the NC file
lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
T_OC_array <- ncvar_get(ncin,"SUM_t_c_12")

# create a matrix of the coordinates and change the array of TOC to a vector
lonlat <- as.matrix(expand.grid(lon,lat))
T_OC_vec <- as.vector(T_OC_array)

# make it is data frame
TOC_df <- data.frame(cbind(lonlat,T_OC_vec))

# change the header names to lat lon for SF package
names(TOC_df) <- c("lon","lat", paste("T_OC", sep="_"))

#update the coordinates from the DF to a geometric object
TOC_df <- TOC_df %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

# read in the glcp and then sf join it with the adapted netcdf soil organic carbon file
glcp <- vroom::vroom("/central/groups/carnegie_poc/rmcclure/glcp-analysis/glcp_extended_thin_GLEE_short.csv") %>%  
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  st_join(., TOC_df, join = st_nearest_feature) %>%
  extract(., geometry, into = c('Lat', 'Lon'), '\\((.*),(.*)\\)', conv = T) %>%
  write_csv(., "/central/groups/carnegie_poc/rmcclure/glcp-analysis/glcp_extended_thin_GLEE_short_w_SOIL.csv")

