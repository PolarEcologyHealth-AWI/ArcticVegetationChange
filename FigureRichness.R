library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
library(rnaturalearth)

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### World Map
circl <- st_point(c(0,0)) %>% st_sfc(crs = proj) %>% st_buffer(5000000)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(proj) %>% st_intersection(circl) %>%
  st_geometry()

load("/Volumes/projects/bioing/user/slisovsk/ArcticSDM/data/grid_25km.rda")
rast <- st_as_stars(st_bbox(grid), dy = 25001, dx = 25001, crs = st_crs(grid)$input, values = NA_real_)

load("/Volumes/projects/bioing/user/slisovsk/ArcticSDM/finalRun//predArray.rda")
### dim: species, cells, years, scenarios, unrestricted|restricted


###########################
### Emerging habitats #####
###########################

out <- lapply(c("unrestricted", "restricted"), function(u) {
  lapply(1:3, function(s) {
    t <- ifelse(u == "unrestricted", 1, 2)
    tibble(type = u, scenario = s, year = 2040, new = apply((predArray[,,1,s,t]==0) & (predArray[,,2,s,t]==1), 1, sum)) %>%
      bind_rows(tibble(type = u, scenario = s, year = 2070, 
                       new = apply((predArray[,,1,s,t]==0 & predArray[,,2,s,t]==0) & (predArray[,,3,s,t]==1), 1, sum))) %>%
      bind_rows(tibble(type = u, scenario = s, year = 2100, 
                       new = apply((predArray[,,1,s,t]==0 & predArray[,,2,s,t]==0 & predArray[,,3,s,t]==0) & (predArray[,,4,s,t]==1), 1, sum)))
  }) %>% Reduce("rbind", .)
}) %>% Reduce("rbind", .)

outSM <- out %>% group_by(type, scenario, year) %>% summarise(median = quantile(new[new>0], probs = 0.5,  na.rm = T)*25*25, 
                                                              lower  = quantile(new[new>0], probs = 0.25, na.rm = T)*25*25, 
                                                              upper  = quantile(new[new>0], probs = 0.75, na.rm = T)*25*25)


ggplot(outSM, aes(x = factor(year), y = median, fill = as.factor(scenario))) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.9)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  facet_wrap(~factor(type, levels = c("unrestricted", "restricted"))) +
  coord_flip() +
  theme_light()



###########################
### Northward movement ####
###########################

dist <- grid %>% mutate(distance = as.numeric(st_distance(., st_point(c(0,0)) %>% st_sfc(crs = proj)))/1000) %>% pull(distance)
# plot(dist %>% dplyr::select(distance))

distOut <- lapply(c("unrestricted", "restricted"), function(u) {
  lapply(1:3, function(s) {
    t <- ifelse(u == "unrestricted", 1, 2)
    tibble(type = u, scenario = s, year = 2040, dist = dist[unlist(apply((predArray[,,1,s,t]==0) & (predArray[,,2,s,t]==1), 1, which))]) %>%
      bind_rows(
        tibble(type = u, scenario = s, year = 2070, dist = dist[unlist(apply((predArray[,,1,s,t]==0 & predArray[,,2,s,t]==0) & (predArray[,,3,s,t]==1), 1, which))])
      ) %>% bind_rows(
        tibble(type = u, scenario = s, year = 2100, dist = dist[unlist(apply((predArray[,,1,s,t]==0 & predArray[,,2,s,t]==0 & predArray[,,3,s,t]==0) & (predArray[,,4,s,t]==1), 1, which))])
      )
  }) %>% Reduce("rbind", .)
}) %>% Reduce("rbind", .)
    

distSM <- distOut %>% group_by(type, scenario, year) %>% summarise(median = quantile(dist, probs = 0.5,  na.rm = T), 
                                                                   lower  = quantile(dist, probs = 0.25, na.rm = T), 
                                                                   upper  = quantile(dist, probs = 0.75, na.rm = T)) %>%
  mutate(year = ifelse(scenario==1, year - 5, ifelse(scenario==3, year + 5, year)))


ggplot(distSM, aes(x = year, y = median, fill = as.factor(scenario))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,
                position=position_dodge(0.9), col = "grey60") +
  # geom_line(mapping = aes(color = as.factor(scenario)), linetype = 2) +
  geom_point(shape = 21, size = 8, col = "grey40") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  facet_wrap(~factor(type, levels = c("unrestricted", "restricted"))) +
  theme_light()









grid$richness_1 <- apply(predArray[,,1,1,1], 2, sum, na.rm = T)
grid$richness_2 <- apply(predArray[,,4,3,1], 2, sum, na.rm = T)
grid$richness_3 <- apply(predArray[,,4,3,2], 2, sum, na.rm = T)

rich_stars <- c(st_rasterize(grid %>% dplyr::select(richness_1), rast),
                st_rasterize(grid %>% dplyr::select(richness_2), rast),
                st_rasterize(grid %>% dplyr::select(richness_3), rast)) %>% merge(name = "layer")

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "transparent") +
  geom_stars(data = rich_stars) +
  facet_wrap(~layer) +
  coord_sf() +
  scale_fill_gradient(low = "yellow", high = "darkgreen", na.value = "transparent") +
  theme_light()


dat <- rich_stars %>% as_tibble() %>% setNames(c("x", "y", "layer", "value"))

ggplot(dat, aes(x=as.factor(layer), y=value)) + 
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, mult = 1, 
               geom="pointrange", color="darkgreen") +
  theme_linedraw()



#### Delta
grid$delta  <- apply(cbind(grid$richness_1, grid$richness_3), 1, diff)
delta_stars <- st_rasterize(grid %>% dplyr::select(delta), rast)

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "transparent") +
  geom_stars(data = delta_stars) +
  coord_sf() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0, na.value = "transparent") +
  theme_light()
