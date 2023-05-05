library(sf)
library(nimble)
library(tidybayes)
library(tidyverse)
library(coda)

glee_diff <- read_csv("./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_sd_added.csv") %>% 
  select(lat, lon, mean_diff, mean_obs_wtemp_k) %>% na.omit(.) %>%
  filter(mean_diff < 5000) %>%
  mutate(mean_obs_wtemp_C = mean_obs_wtemp_k - 273.15) %>%
  select(-mean_obs_wtemp_k) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

world <-  ne_download(scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext")

grid_spacing <- 555000 # CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree)

grid <- st_make_grid(
  world,
  cellsize = c(grid_spacing, grid_spacing),
  #n = c(200, 200), # grid granularity
  crs = st_crs(world),
  what = "polygons",
  square = FALSE) %>%
  st_intersection(world)

grid <- st_sf(index = 1:length(lengths(grid)), grid)

area_hexes <- st_join(glee_diff, grid, join = st_intersects)

area_hexes_avg <- area_hexes %>%
  ungroup() %>%
  group_by(index) %>%
  summarise(mean_diffusion = mean(mean_diff, na.rm = TRUE),
            mean_temp = mean(mean_obs_wtemp_C, na.rm = TRUE)) %>%
  st_drop_geometry(.) %>%
  arrange(index)

diff_data_spatial <- merge(x = grid, y = area_hexes_avg, by = "index", all.x = FALSE)

diff_data_spatial.nb <- poly2nb(diff_data_spatial, row.names =  rownames(diff_data_spatial))

nbInfo <- nb2WB(diff_data_spatial.nb)

# A vector of indices indicating which regions are neighbors of which.
nbInfo$adj

head(nbInfo$weights)

nbInfo$num

nregions <- nrow(diff_data_spatial)

code <- nimbleCode({
  # priors
  beta ~ dnorm(0, sd = 1000)
  theta ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 1000)   # prior for variance components based on Gelman (2006)
  tau <- 1 / sigma^2
  # latent process
  s[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 0)
  # likelihood
  for(i in 1:N) {
    lambda[i] <- (theta * (beta*x[i]-20)) + s[i]
    y[i] ~ dpois(lambda[i])
  }
})

x <- diff_data_spatial$mean_temp

set.seed(1)

constants <- list(N = nregions, L = length(nbInfo$adj), 
                  adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num)

data <- list(y = diff_data_spatial$mean_diffusion,
             x = diff_data_spatial$mean_temp)

nchain <- 3

inits <- list()

for(i in 1:nchain){
  inits[[i]] <- list(theta = rnorm(1,100, 30),
                     beta = rnorm(1, 1.1, 0.5),
                     sigma = runif(1,0.05, 0.15))
}


nimble.out <- nimbleMCMC(code = code, 
                         data = data, 
                         inits = inits,
                         constants = constants,
                         monitors = c("theta", "beta", "sigma"),
                         niter = 10000,
                         nchains = 3,
                         samplesAsCodaMCMC = T)

burnin <- 1000
nimble.burn <- window(nimble.out, start = burnin)

traceplot(nimble.burn)
gelman.diag(nimble.burn)


inits <- list(beta = 1.1, theta = 100, sigma = 1, s = rnorm(nregions))

model <- nimbleModel(code, constants = constants, data = data, inits = inits)

cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('beta', 'theta', 'sigma'))

conf$printSamplers()

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

samples <- runMCMC(cMCMC, niter = 10000, nburnin = 1000)


plot(density(samples[,"sigma"]))

plot(density(samples[,"s[1]"]))

plot(density(samples[,"beta"]))

plot(density(samples[,"theta"]))




ggplot(area_hexes_avg) +
  geom_sf(data = world, color = "black", lwd = 1.5)+
  geom_sf(lwd = 0.4,
          aes(fill = mean_diffusion),color = "white")+
  scale_fill_viridis(option = "A", na.value = "white",
                     direction = -1, discrete = F, name = "**Basin-level Predictor** <br> Random Forest Model")+
  #scale_color_manual(values = c("blue", NA, "black"), na.value = "black",name = "**Predictive Skill** <br> NSE") +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-8600000, 8600000), expand = FALSE) +
  #guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.09, 0.2),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

