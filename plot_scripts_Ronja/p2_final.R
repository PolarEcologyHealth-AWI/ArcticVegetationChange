#### packages ####
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
library(rnaturalearth)
library(terra)
library(dplyr)
library(ggplot2)


### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### World Map
circl <- st_point(c(0,0)) %>% st_sfc(crs = proj) %>% st_buffer(5000000)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(proj) %>% st_intersection(circl) %>%
  st_geometry()

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{data}/grid_25km.rda")) #grid

rast <- st_as_stars(st_bbox(grid), dy = 25001, dx = 25001, crs = st_crs(grid)$input, values = NA_real_)

load(glue::glue("{data}/predArray.rda")) #array
### dim: species, cells, years, scenarios, unrestricted|restricted

#for gform plot
load(glue::glue("{data}/spTable.rda")) #spTable
traits <- read.csv(glue::glue("{data}/migration_meters_all.csv"))

###########################
### Emerging habitats #####
###########################

out <- lapply(c("unrestricted", "restricted"), function(u) {
  lapply(1:3, function(s) {
    t <- ifelse(u == "unrestricted", 1, 2)
    tibble(type = u, scenario = s, year = 2040, new = apply((predArray[,,1,s,t]==0) & (predArray[,,2,s,t]==1), 1, sum)) %>%
      bind_rows(tibble(type = u, scenario = s, year = 2070, 
                       new = apply((predArray[,,1,s,t]==0) & (predArray[,,3,s,t]==1), 1, sum))) %>%
      bind_rows(tibble(type = u, scenario = s, year = 2100, 
                       new = apply((predArray[,,1,s,t]==0) & (predArray[,,4,s,t]==1), 1, sum)))
  }) %>% Reduce("rbind", .)
}) %>% Reduce("rbind", .)

outSM <- out %>% group_by(type, scenario, year) %>% summarise(median = quantile(new[new>0], probs = 0.5,  na.rm = T)*25*25, 
                                                              lower  = quantile(new[new>0], probs = 0.25, na.rm = T)*25*25, 
                                                              upper  = quantile(new[new>0], probs = 0.75, na.rm = T)*25*25)


#### pl1 ####
pl1 <- ggplot(outSM, aes(x = factor(year), y = median, fill = as.factor(scenario))) +
  geom_bar(stat="identity", position=position_dodge(), show.legend = F) +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                width=0,  # Remove horizontal ticks
                position=position_dodge(0.9),
                color="grey50",  # Make error bars grey
                size=0.5) +  # Adjust line thickness if needed
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  facet_wrap(~factor(type, levels = c("unrestricted", "restricted"))) +
  theme(strip.background = element_rect(fill = "grey65"))+
  labs(x = "", y = "Emerging habitat") +
 # labs(x = "", y = expression(paste("Emerging habitat [", Delta, " 2010, km"^2, "]"))) +
  theme_light()
 


plot(pl1)

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


#### pl2 ####
pl2 <- ggplot(distSM, aes(x = year, y = median, fill = as.factor(scenario))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,
                position=position_dodge(0.9), col = "grey60") +
  geom_point(shape = 21, size = 3, col = "grey40", show.legend = F) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, 'YlOrRd')) +
  facet_wrap(~factor(type, levels = c("unrestricted", "restricted"))) +
  theme(strip.background = element_rect(fill = "grey65"))+
  ylab("Distance to 90 degrees north") + xlab("") +
  theme_light()+
  scale_x_continuous(breaks = c(2040, 2070, 2100), labels = c("2040", "2070", "2100"))
  

plot(pl2)

##### slopes ####
lapply(c("unrestricted", "restricted"), function(t) {
  mod <- lm(dist~year, data = distOut %>% filter(type == t, scenario == 1))
  nsim   <- 2000
  bsim   <- arm::sim(mod, n.sim=nsim)
  tibble(type = t, slope = median(bsim@coef[,2]), upper.0.975 = quantile(bsim@coef[,2], probs = 0.975), lower.0.025 = quantile(bsim@coef[,2], probs = 0.025))
}) %>% Reduce("rbind", .)

#### gform ####

#### res generate h0/result ####
result <- as.data.frame(matrix(nrow=nrow(spTable), ncol=11))
result[,1] <- spTable[,1]
colnames(result) = c("species", "pres", "y2s2", "y2s3", "y2s4", 
                     "y3s1", "y3s2", "y3s3",
                     "y4s1","y4s2","y4s3")


spArray <- predArray[,,,,2]

# Precompute the reference values
reference_values <- spArray[, , 1, 2]

# Calculate the sum for the second column
result[, 2] <- rowSums(reference_values)

# Create a function to calculate percentages with ifelse
calculate_newcells <- function(sp_data, reference, result_column, result_matrix) {
  comparison <- sp_data > reference
  sum_values <- rowSums(reference)
  result_matrix[, result_column] <- rowSums(comparison)
  result_matrix
}

# Calculate percentages for y = 2, 3, 4
for (y in 2:4) {
  start_col <- 3 * y - 4
  result <- calculate_newcells(spArray[, , y, 1], reference_values, start_col + 1, result)
  result <- calculate_newcells(spArray[, , y, 2], reference_values, start_col + 2, result)
  result <- calculate_newcells(spArray[, , y, 3], reference_values, start_col + 3, result)
}

h0 <- result
#write.csv2(result, glue::glue("{data}/sp_pres_newcells.csv"), row.names = FALSE)


#### unres generate h0/result ####

spArray <- predArray[,,,,1]

# Precompute the reference values
reference_values <- spArray[, , 1, 1] 

# Calculate the sum for the second column
result[, 2] <- rowSums(reference_values)

# Create a function to calculate percentages with ifelse
calculate_newcells <- function(sp_data, reference, result_column, result_matrix) {
  comparison <- sp_data > reference
  sum_values <- rowSums(reference)
  result_matrix[, result_column] <- rowSums(comparison)
  result_matrix
}

# Calculate percentages for y = 2, 3, 4
for (y in 2:4) {
  start_col <- 3 * y - 4
  result <- calculate_newcells(spArray[, , y, 1], reference_values, start_col + 1, result)
  result <- calculate_newcells(spArray[, , y, 2], reference_values, start_col + 2, result)
  result <- calculate_newcells(spArray[, , y, 3], reference_values, start_col + 3, result)
}

unres <- result

#write.csv2(result, glue::glue("{data}/sp_pres_newcells_unrespres.csv"), row.names = FALSE)

#h0 <- read.csv(glue::glue("{data}/sp_pres_newcells.csv"), sep=";")

h0t <- h0 %>%
  left_join(traits, by = "species")

h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"
h0t$gform[134] <- "Shrub"

#unres <- read.csv(glue::glue("{data}/sp_pres_newcells_unrespres.csv"), sep=";")

h0t$unres <- unres$y4s3
h0t <- subset(h0t, h0t$unres > 0)
h0t$percent <- h0t$y4s3/h0t$unres

result1g <- h0t

summary(as.factor(result1g$gform))

result1g <- result1g %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "Aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))

result1g$type <- "restricted"


p1g <- ggplot(result1g, aes(x = reorder(gform, -percent, FUN = median), y = percent)) +
  geom_boxplot(fill = "#009E73", outlier.size = 1, outlier.color = "grey50") +
  labs(title = "",
       x = "Growth form",
       y = "Reachable habitat share") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_rect(fill = "grey65"),  # Darker grey background for strip
    strip.text = element_text(color = "white")  # White bold text on strip
  ) +
  facet_wrap(~factor("restricted", levels = c("unrestricted", "restricted")))


print(p1g)


gform_summary <- result1m %>%
  group_by(gform) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE),
            median_percent = median(percent, na.rm = TRUE))

print(gform_summary)


#### together ####

library(ggpubr)

ggarrange(pl1, pl2, p1g, 
          nrow = 3, 
          labels = c("a)", "b)", "c)"),
          font.label = list(size = 16, color = "black", face="plain"),
          label.x = 0.05) 
