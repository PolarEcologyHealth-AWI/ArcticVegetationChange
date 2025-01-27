#p3ABC


#### packages ####
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(geodata)

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{data}/predArray.rda")) #array
load(glue::glue("{data}/grid_25km.rda")) #grid
load(glue::glue("{data}/spTable.rda")) #spTable

spArray <- predArray[,,,,2]

gridsf <- st_as_sf(grid)

com1 <- gridsf %>%
  mutate(sp = apply(spArray[,,1,1], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))


summary(com1$sp)
        
summary(apply(predArray[,,1,1,1],2, sum, na.rm = TRUE))

pres <- as.data.frame(predArray[,,1,1,1])

com2 <- gridsf %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum, na.rm = TRUE))

#unconstrained
spArray <- predArray[,,,,1]

com3 <- gridsf %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum, na.rm = TRUE))

summary(com1$sp)
plot(com1)

sp = apply(spArray[,,4,3], 2, sum)

summary(sp)

#### ABC ####

  
library(ggplot2)
library(viridis)
library(stars)

p11 <- ggplot() +
  geom_stars(data = com1, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +  # Apply Viridis color scale and handle NA values
  labs(title = "A) 2010") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )

print(p11)
p2s1 <- ggplot() +
  geom_stars(data = com2, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +
  labs(title = "B) 2100 constrained") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )

p2s5 <- ggplot() +
  geom_stars(data = com3, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +
  labs(title = "C) 2100 unconstrained") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )


summary(com2)

s1 = as.data.frame(spArray[,,1,1])
s1 = apply(spArray[,,1,1], 2, sum)
s2 = apply(spArray[,,4,3], 2, sum)
summary(s2)



###########
library(cowplot)

# Create a combined plot without a separate legend
combined_plot <- plot_grid(
  p11 + theme(legend.position = "none"),  # Remove legend from p1
  p2s1 + theme(legend.position = "none"),  # Remove legend from p2s1
  p2s5 + theme(legend.position = "none"),  # Remove legend from p2s5
  ncol = 3, 
  nrow = 1
)

# Adjust the legend to reduce space around it
legend <- get_legend(
  p11 + 
    theme(
      legend.position = "right", 
      legend.box.margin = margin(0, 0, 0, 0),  # Remove margin around the legend box
      legend.spacing.y = unit(0, "cm"),        # Reduce vertical spacing between legend items
      legend.spacing.x = unit(0, "cm")         # Reduce horizontal spacing between legend items
    )
)

# Add the adjusted legend directly to the right of the combined plot
combined_plot_with_legend <- plot_grid(
  combined_plot, 
  legend, 
  rel_widths = c(3, 0.2)  # Adjust relative widths to minimize space further
)

combined_plot_with_legend <- plot_grid(
  p11 + theme(legend.position = "none"),  # Remove legend from p1
  p2s1 + theme(legend.position = "none"),  # Remove legend from p2s1
  p2s5 + theme(legend.position = "none"),  # Remove legend from p2s5
  legend, 
  nrow = 1,
  rel_widths = c(1,1,1,0.2)  # Adjust relative widths to minimize space further
)
# Print the final combined plot with a compact legend
print(combined_plot_with_legend)
ABC <- combined_plot_with_legend

#### violin plot ####

#data
summary(apply(predArray[,,1,1,1],2, sum, na.rm = TRUE))
sp1 <- apply(predArray[,,1,1,2],2, sum, na.rm = TRUE)
sp2 <- apply(predArray[,,4,3,2],2, sum, na.rm = TRUE)
sp3 <- apply(predArray[,,4,3,1],2, sum, na.rm = TRUE)

spdata <- as.data.frame(cbind(sp1,sp2,sp3))
spdata$cell <- 1:38657

boxplot(data=spdata[,1:3], x=spdata$cell)

str(spdata)

library(ggplot2)
library(tidyr)


# Reshape the data from wide to long format
spdata_long <- pivot_longer(spdata, cols = c(sp1, sp2, sp3), names_to = "sp", values_to = "value")


spdata_long <- spdata_long %>%
  mutate(sp = case_when(
    sp == "sp1" ~ "2010",
    sp == "sp2" ~ "2100, cons",
    sp == "sp3" ~ "2100, uncons",
    TRUE ~ as.character(sp)  # Keep the original value as a string if none of the above conditions match
  ))

# Create the boxplot

dataplot <- ggplot(spdata_long, aes(x = sp, y = value)) +
  geom_violin(fill="dodgerblue3") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "white") +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  labs(x = "year and scenario", y = "species", title = "F) richness", fill = "data") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )


print(dataplot)

summary(spdata$sp1)
summary(spdata$sp2)
summary(spdata$sp3)

#### DE ####

gridsf <- st_as_sf(grid)

# Create a function to count the number of new 1 cells for each grid cell
count_new_ones <- function(present, future) {
  new_ones <- future - present
  new_ones[new_ones < 0] <- 0
  return(new_ones)
}


##### sig1 ####

spArray <- predArray[,,,,2]

present <- as.data.frame(spArray[, , 1, 1])
future <- as.data.frame(spArray[, , 4, 3])


new_array <- as.data.frame(count_new_ones(present, future))

com1d <- gridsf %>%
  mutate(sp = apply(new_array, 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))

summary(com1d$sp)
plot(com1d)

##### plot ####

p1 <- ggplot() +
  geom_stars(data = com1d, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(option = "C", limits = c(0, 108), na.value = "transparent", name = "new species") +  # Apply Viridis color scale and handle NA values
  labs(title = expression(paste("D) ", Delta, " 2100 constrained"))) +  # Update labels    
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "right" , # Remove the legend
    options(repr.plot.width = 5, repr.plot.height = 4) 
  )
print(p1)

##### sig5 ####

spArray <- predArray[,,,,1]

present <- as.data.frame(spArray[, , 1, 1])
future <- as.data.frame(spArray[, , 4, 3])


new_array <- as.data.frame(count_new_ones(present, future))


com5d <- gridsf %>%
  mutate(sp = apply(new_array, 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))


##### plot ####

p2 <- ggplot() +
  geom_stars(data = com5d, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(option = "C", limits = c(0, 472), na.value = "transparent", name = "new species") +  # Apply Viridis color scale and handle NA values
  labs(title = expression(paste("E) ", Delta, " 2100 unconstrained"))) + # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    #plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "right", # Remove the legend +
    options(repr.plot.width = 5, repr.plot.height = 4) 
  )


print(p2)

##### both ####
library(cowplot)

# Combine plots with specified relative widths
DE <- plot_grid(
  NULL,   # Empty space for the first column (0.5)
  p1,     # Plot 1 in the second column (1)
  p2,     # Plot 2 in the third column (1)
  NULL,   # Empty space for the fourth column (0.5)
  ncol = 4,   # Specify number of columns
  rel_widths = c(0.4, 1.2, 1.2, 0.4) # Adjust relative widths
)

# DE <- plot_grid(
#   NULL,
#   p1,     # Plot 1 in the second column (1)
#   p2,     # Plot 2 in the third column (1)
#   ncol = 3   # Specify number of columns
# )

print(DE)

ABCDE <- plot_grid(ABC, DE, nrow = 2)

print(ABCDE)

summary(com5d$sp)

#### elevation ####

# Download elevation data 
#elevation_data <- worldclim_global(var = "elev", res = 10, data = "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/elevation")

#input: elevation and gg
#elevation_data <- rast("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/elevation/wc2.1_10m/wc2.1_10m_elev.tif")
elevation_data <- rast(glue::glue("{data}/wc2.1_10m_elev.tif"))
plot(elevation_data)

# Convert SpatRaster to a data frame with coordinates
elevation_df <- as.data.frame(elevation_data, xy = TRUE)

# Rename columns for clarity
colnames(elevation_df) <- c("Longitude", "Latitude", "Elevation")

# Convert to sf object
elevation_sf <- st_as_sf(elevation_df, coords = c("Longitude", "Latitude"), crs = 4326)

#crs
el_transformed <- st_transform(elevation_sf, crs = st_crs(grid))

# Crop the elevation sf object using the grid's extent
elevation_c_sf <- st_crop(el_transformed, st_bbox(grid) )

# Define new resolution (e.g., 1 km)
new_resolution <- 25050  # in meters

# Create a grid over the extent of your elevation data
new_grid <- st_make_grid(elevation_c_sf, cellsize = new_resolution, what = "polygons")

com <- st_as_sf(com1d)
comc <- st_transform(com, crs = crs(elevation_c_sf)) 

com5 <- st_as_sf(com5d)
comc5 <- st_transform(com5, crs = crs(elevation_c_sf))


# Step 2: Perform the spatial join
combined_data <- st_join(elevation_c_sf, comc, join = st_intersects)
combined_data5 <- st_join(elevation_c_sf, comc5, join = st_intersects)

dat <- na.omit(combined_data)
dat5 <- na.omit(combined_data5)


str(dat)

df <- dat %>%
  mutate(coords = st_coordinates(geometry)) %>%
  mutate(x = coords[,1],
         y = coords[,2]) %>%
  select(-coords) %>%
  st_drop_geometry()

names(df) <- c("altitude", "sig1", "x", "y")


df$sig5 <- dat5$sp

df_long <- df[,1:5]%>%
  pivot_longer(
    cols = c(sig1, sig5),
    names_to = "scenario",
    values_to = "newspecies"
  )

# Create altitude classes
df_long <- df_long %>%
  mutate(altitude_class = cut(altitude, breaks = seq(0, max(altitude, na.rm = TRUE) + 200, by = 200),
                              labels = seq(0, max(altitude, na.rm = TRUE), by = 200)))

df_long <- na.omit(df_long)


df_summary <- df_long %>% 
  group_by(scenario, altitude_class) %>% 
  summarise(median = median(newspecies[newspecies > 0], na.rm = TRUE),
            lower = quantile(newspecies[newspecies > 0], probs = 0.25, na.rm = TRUE),
            upper = quantile(newspecies[newspecies > 0], probs = 0.75, na.rm = TRUE),
            .groups = 'drop')

df_summary[18,3] <- 0

elev_plot <- ggplot(df_summary %>% filter(scenario %in% c("sig1", "sig5")), 
       aes(x = altitude_class, y = median, fill = scenario)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("sig1" = "deepskyblue2", "sig5" = "dodgerblue4"),
                    labels = c("sig1" = "cons", "sig5" = "uncons"),
                    name = "Scenario") +
  labs(x = "Altitude [m]", y = "New species", 
       title = "G) Elevation impact",) +
  theme_minimal() +
  theme(
          plot.title = element_text(size = 16, hjust = 0),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14),
        legend.position = c(0.95, 0.95),  # Position legend inside top-right
        legend.justification = c(1, 1),   # Anchor point for legend
        legend.box.just = "right",        # Justify legend to the right
        legend.margin = margin(6, 6, 6, 6),  # Add some margin around legend
        legend.background = element_rect(fill = "white", color = NA)) +  # Optional: add white background to legend) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1, check.overlap = TRUE),
                   breaks = function(x) x[seq(1, length(x), 2)])+
  
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                width=0,  # Remove horizontal ticks
                position=position_dodge(0.9),
                color="grey50",  # Make error bars grey
                size=0.5)  # Adjust line thickness if needed


line3 <- plot_grid(
  dataplot, elev_plot,
  ncol = 2   # Specify number of columns
)
# Print the combined plot
print(line3)

ABCDEFG <- plot_grid(ABC, DE, line3, nrow = 3)

print(ABCDEFG)
