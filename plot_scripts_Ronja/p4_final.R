#plot 4

#source louvain_analysis
#https://stackoverflow.com/questions/49834827/louvain-community-detection-in-r-using-igraph-format-of-edges-and-vertices

#### packages & data ####
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(psych)
library(igraph)
library(ggplot2)
library(glue)
library(units)

path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{path}/grid_25km.rda")) #grid
gridsf <- st_as_sf(grid)

load(glue::glue("{path}/predArray.rda")) #array

#############################################################################
#### A ####
#############################################################################

spArray <- predArray[,,,,2]

#### clustering loop ####

for(y in 1:4){
  
  dat <- spArray[,,y,3] #year
  
  dat[is.na(dat)] <- 0
  cormatrix <- cor (dat)
  
  #save(cormatrix, file = "cormatrix_cells_y4s3.rda")
  
  distancematrix <- cor2dist(cormatrix)
  
  DM2 <- as.matrix(distancematrix)
  
  DM2[cormatrix < 0.3] = 0
  DM2[is.na(DM2)] <- 0
  
  G2a <- graph_from_adjacency_matrix(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)
  cl2 <- cluster_louvain(G2a) 
  
  membercl <- data.frame(group = cl2$membership, label = 1:38657)
  
  write.csv(membercl, glue::glue("{path}membercl_cells_y{y}s3_ne.csv"), row.names = FALSE)
}

##### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[1] #tundra
tundra2 <- hcl.colors(6, "Greens 3")[2] #tundra
tundra3 <- hcl.colors(6, "Greens 3")[3] #tundra
tundra4 <- hcl.colors(6, "Greens 3")[4] #tundra

taiga_ne <- hcl.colors(6, "Purples")[1] #nearctic
taiga_ne2 <- hcl.colors(6, "Purples")[2] #nearctic
taiga_ne3 <- hcl.colors(6, "Purples")[3] #nearctic
taiga_ne4 <- hcl.colors(6, "Purples")[4] #nearctic

taiga_eu <- hcl.colors(6, "Blues 3")[1] #europe
taiga_eu_W <- hcl.colors(6, "Blues 3")[2] #europe
taiga_eu_E <- hcl.colors(6, "Blues 3")[3] #europe
taiga_eu_E2 <- hcl.colors(6, "Blues 3")[4] #europe


col_all <- c(taiga_eu, taiga_eu_E, taiga_eu_E2, taiga_eu_W,
             taiga_ne, taiga_ne2, taiga_ne3, taiga_ne4,
             tundra, tundra2, tundra3)

##### y1 ####
membercl <- read.csv(glue::glue("{path}/membercl_cells_y1s3_sigma1.csv"), sep =",")


# Step 1: Identify groups with only one label
single_label_groups <- names(which(table(membercl$group) == 1))

# Step 2: Update the group values to 0 for these labels
membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0

unique(membercl$group)

membercl <- membercl %>%
  mutate(new_group = case_when(
    group == 0 ~ NA,
    group == 1 ~ "tundra",
    group == 9 ~ "taiga, nearctic",
    group == 10 ~ NA,
    group == 12 ~ "taiga, eurasian",
    TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
  ))
write.csv(membercl, file = glue::glue("{path}membercl1.csv"), row.names = F)
col1 = c(taiga_eu, taiga_ne, tundra)

##### join ####

sp_with_groups <- gridsf %>%
  left_join(membercl %>% select(label, new_group), by = join_by(Id == label))

# Rasterizing the data with new_group
com1 <- sp_with_groups %>%
  mutate(sp = as.factor(new_group)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))


##### plotting ####
plot1 <- ggplot() +
  geom_stars(data = com1, aes(fill = as.factor(sp))) +  # Use fill for raster colors
  scale_fill_manual(values = col1, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
  theme_minimal() +
  ggtitle("a)", subtitle = "2010")+
  labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
  theme(
    plot.title = element_text(size = 16),
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )
print(plot1)

##### y2 ####
membercl <- read.csv(glue::glue("{path}/membercl_cells_y2s3_sigma1.csv"), sep =",")

single_label_groups <- names(which(table(membercl$group) == 1))
membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0

membercl <- membercl %>%
  mutate(new_group = case_when(
    group == 0 ~ NA,
    group == 1 ~ NA,
    group == 5 ~ "tundra",
    group == 10 ~ "taiga, nearctic",
    group == 11 ~ "taiga, eurasian",
    group == 19 ~ "taiga, eurasian, E",
    TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
  ))

write.csv(membercl, file = glue::glue("{path}membercl2.csv"), row.names = F)

col2 = c(taiga_eu, taiga_eu_E, taiga_ne, tundra)
com2 <- gridsf %>% mutate(sp = membercl$new_group) %>% dplyr::select(sp) %>%
  st_rasterize(.,st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_)) 

unique(membercl$new_group)

##### join ####

sp_with_groups <- gridsf %>%
  left_join(membercl %>% select(label, new_group), by = join_by(Id == label))

# Rasterizing the data with new_group
com2 <- sp_with_groups %>%
  mutate(sp = as.factor(new_group)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))

##### plotting ####
plot2 <- ggplot() +
  geom_stars(data = com2, aes(fill = as.factor(sp))) +  # Use fill for raster colors
  scale_fill_manual(values = col2, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
  theme_minimal() +
  ggtitle("", subtitle = "2040")+
  labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )
print(plot2)
  


  ##### y3 ####
  membercl <- read.csv(glue::glue("{path}/membercl_cells_y3s3_sigma1.csv"), sep =",")

  single_label_groups <- names(which(table(membercl$group) == 1))
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ NA,
      group == 6 ~ "tundra'",
      group == 10 ~ "taiga, nearctic'",
      group == 11 ~ "taiga, eurasian, W",
      group == 19 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))

  write.csv(membercl, file = glue::glue("{path}membercl3.csv"), row.names = F)
  col3 = c(taiga_eu_E, taiga_eu_W, taiga_ne2, tundra2)
  
  ##### join ####

  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com3 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
  ##### plotting ####
  plot3 <- ggplot() +
    geom_stars(data = com3, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col3, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("", subtitle = "2070")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot3)

  ##### y4 ####

  membercl <- read.csv(glue::glue("{path}membercl_cells_y4s3_sigma1.csv"), sep =",")
  
  single_label_groups <- names(which(table(membercl$group) == 1))
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~  "taiga, nearctic'''",
      group == 6 ~ "tundra''",
      group == 13 ~ "taiga, eurasian, W",
      group == 17 ~ "taiga, nearctic''",
      group == 19 ~ "taiga, eurasian, E'",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  
  write.csv(membercl, file = glue::glue("{path}membercl4.csv"), row.names = F)
  
  col4 = c(taiga_eu_E2, taiga_eu_W, taiga_ne3, taiga_ne4, tundra3)
  
  #### join ####
  library(dplyr)
  library(sf)
  library(stars)
  
  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com4 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
  ##### plotting ####
  plot4 <- ggplot() +
    geom_stars(data = com4, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col4, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("", subtitle = "2100")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot4)
  
  #### combined ####
  library(cowplot)
  A <- plot_grid(plot1, plot2, plot3, plot4 , ncol = 4)
  print(A)

#############################################################################  
#### B ####
#############################################################################
  
  spArray <- predArray[,,,,2]
  sp1 <- as.data.frame(spArray[,,1,1])
  sp1 <- as.data.frame(t(sp1))
  
  membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
  membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
  membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
  membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))
  
  #### ndist ####
  
  ##### distance ####
  
  center_point <- st_sfc(st_point(c(0, 0)), crs = st_crs(grid))
  
  distg <- grid %>%
    mutate(dist_n = st_distance(geometry, center_point))
  
  colnames(distg)[5] <- "distN"
  
  distN <- as.data.frame(distg[,5])
  
  membercldat1$dist <- distg$distN
  membercldat4$dist <- distg$distN
  
  # Assign categories based on the presence of keywords
  membercldat4$Category[grepl("tundra", membercldat4$new_group)] <- "tundra"
  membercldat4$Category[grepl("taiga, eurasian", membercldat4$new_group)] <- "taiga, eurasian"
  membercldat4$Category[grepl("taiga, nearctic", membercldat4$new_group)] <- "taiga, nearctic"
  
  membercldat1$Category = membercldat1$new_group
  
  membercldat4$dist <- as.numeric(set_units(membercldat4$dist, NULL))
  membercldat1$dist <- as.numeric(set_units(membercldat1$dist, NULL))
  
  m1 = as.data.frame(cbind(membercldat1$label, membercldat1$Category, membercldat1$dist))
  m1 = na.omit(m1)
  colnames(m1) = c("label", "Category", "dist")
  m1$dist = as.numeric(m1$dist)
  
  m4 = as.data.frame(cbind(membercldat4$label, membercldat4$Category, membercldat4$dist))
  m4 = na.omit(m4)
  colnames(m4) = c("label", "Category", "dist")
  m4$dist = as.numeric(m4$dist)
  
  ##### ndist delta ####
  
  m4$ddist = m4$dist - median(m1$dist)
  summary(m1$dist)
  
  m4a <- m4 %>%
    filter(Category != "tundra")
  
  
  ##### new cells ####
  
  ####### 3 coms #####
  taiga_ne1 <- m1 %>%
    filter(Category %in% "taiga, nearctic")
  
  taiga_ne2 <- m4 %>%
    filter(Category %in% "taiga, nearctic")
  
  taiga_ne_diff <- anti_join(taiga_ne2, taiga_ne1, by = "label")
  
  taiga_ne_diff$delta <- -taiga_ne_diff$dist + median(taiga_ne1$dist)
  
  taiga_eu1 <- m1 %>%
    filter(Category %in% "taiga, eurasian")
  
  taiga_eu2 <- m4 %>%
    filter(Category %in% "taiga, eurasian")
  
  taiga_eu_diff <- anti_join(taiga_eu2, taiga_eu1, by = "label")
  
  taiga_eu_diff$delta <- -taiga_eu_diff$dist + median(taiga_eu1$dist)
  
  tundra1 <- m1 %>%
    filter(Category %in% "tundra")
  
  tundra2 <- m4 %>%
    filter(Category %in% "tundra")
  
  
  tundra_diff <- anti_join(tundra2, tundra1, by = "label")
  #category4 has the category with new cells
  
  tundra_diff$delta <- -tundra_diff$dist + median(tundra1$dist)
  
  dif <- as.data.frame(rbind(taiga_eu_diff, taiga_ne_diff, tundra_diff))
  dif <- as.data.frame(rbind(taiga_eu_diff, taiga_ne_diff))
  
  ##### colours ####
  library(scales)
  
  tundra <- hcl.colors(6, "Greens 3")[2] #tundra
  taiga_ne <- hcl.colors(6, "Purples")[2] #america
  taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe
  
  
  B <- ggplot(data = dif, aes(x = Category, y = delta * 0.001, fill = Category)) +
    geom_boxplot(color = "black", outlier.size = 1, outlier.shape = 16) +
    labs(x = "",  y = "migration distance [km]", title = "new cells") +
    theme_minimal() +
    labs(title = "b)") + 
    # labs(title = ") new cells") +  # Update labels
    theme(
      axis.text.x = element_blank(), #element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
      axis.text.y = element_text(size = 14),  # Increase y-axis text size
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      plot.title = element_text(size = 16),  # Increase title
      legend.position = "none")+  # Remove the legend
    ylim(0, 150) +
    scale_fill_manual(values = c("taiga, nearctic" = taiga_ne, 
                                 "taiga, eurasian" = taiga_eu)) 
  
  print(B)
  
  
  eu <- subset(dif, dif$Category == "taiga, eurasian")
  ne <- subset(dif, dif$Category == "taiga, nearctic")
  
  
  #Anderson-Darling Test:
  library(nortest)
  ad_test_result <- ad.test(ne$delta)
  print(ad_test_result)
  
  wilcox.test(eu$delta, ne$delta)

#############################################################################
#### C ####
#############################################################################
  
  path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"
  
  membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
  membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
  membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
  membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))
  
  membercldat1 = na.omit(membercldat1)
  membercldat2 = na.omit(membercldat2)
  membercldat3 = na.omit(membercldat3)
  membercldat4 = na.omit(membercldat4)
  
  year1_clusters <-  cbind(membercldat1$new_group, membercldat1$label)
  year2_clusters <-  cbind(membercldat2$new_group, membercldat2$label)
  year3_clusters <-  cbind(membercldat3$new_group, membercldat3$label)
  year4_clusters <-  cbind(membercldat4$new_group, membercldat4$label)
  
  colnames(year1_clusters) <- c("cluster_group", "cell_number")
  colnames(year2_clusters) <- c("cluster_group", "cell_number")
  colnames(year3_clusters) <- c("cluster_group", "cell_number")
  colnames(year4_clusters) <- c("cluster_group", "cell_number")
  
  # Merge the data frames on cell_number
  merged_data <- merge(year1_clusters, year2_clusters, by = "cell_number", suffixes = c("_year1", "_year2"))
  
  merged_data <- merge(merged_data, year3_clusters, by = "cell_number")
  colnames(merged_data)[4] <- "cluster_group_year3"
  merged_data <- merge(merged_data, year4_clusters, by = "cell_number")
  colnames(merged_data)[5] <- "cluster_group_year4"
  
  colnames(merged_data) <- c("cell", "year1", "year2", "year3", "year4")
  
  #write.csv(merged_data, "clusters_dispersal.csv", row.names =F)
  
  ##### Sankey ####
  #remotes::install_github("davidsjoberg/ggsankey")
  library(ggsankey)
  
  colnames(merged_data) <- c("cells", "1", "2", "3", "4")
  
  df2 <- merged_data %>%
    make_long("1", "2", "3", "4")
  
  df2 <- df2 %>%
    mutate(year = case_when(
      x == 1 ~ 2010,
      x == 2 ~ 2040,
      x == 3 ~ 2070,
      x == 4 ~ 2100
      # Keep the original value as a string if none of the above conditions match
    ))
  
  col_sankey2 <- c(
    taiga_ne, taiga_ne2, taiga_eu, taiga_eu_E, taiga_eu_W,
    tundra, tundra2,taiga_ne4, taiga_ne3, taiga_eu_E2,  tundra3)
  
  df2 <- df2 %>% mutate_all(~ gsub("taiga", "boreal forest", .))
  
  df2 <- df2 %>% mutate_all(~ gsub("eurasian", "palearctic", .))
  
  df2 <- df2 %>%
    mutate(label = ifelse(x == 1, node, NA))
  
  unique(membercldat1$new_group)
  
  C <- ggplot(df2, aes(x = as.factor(year), 
                       next_x = next_x, 
                       node = node, 
                       next_node = next_node,
                       fill = factor(node),
                       label = label)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    scale_fill_manual(values = col_sankey2,
                      guide = guide_legend(title = "clusters"),
                      na.value = "transparent") +
    theme_sankey() +
    theme(legend.position = "none",
          plot.title = element_text(size = 16)) +
    ggtitle("c)")+
    labs(x = "", y = "") 
  
  print(C)
  
##############################################################################  
#### D ####
##############################################################################
  
  ##### summary table ####
  
  membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
  membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
  membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
  membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))
  
  # Convert tables to data frames
  df1 <- as.data.frame(table(membercldat1$new_group))
  df2 <- as.data.frame(table(membercldat2$new_group))
  df3 <- as.data.frame(table(membercldat3$new_group))
  df4 <- as.data.frame(table(membercldat4$new_group))
  
  # Set the names of the data frames to the first column (group names)
  names(df1) <- c("Group", "Count1")
  names(df2) <- c("Group", "Count2")
  names(df3) <- c("Group", "Count3")
  names(df4) <- c("Group", "Count4")
  
  # Merge all data frames by Group
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE), list(df1, df2, df3, df4))
  
  # Sort the merged_df by the first column (Group) alphabetically
  sorted_df <- merged_df[order(merged_df$Group, na.last = TRUE), ]
  
  sorted_df$Category <- NA  # Initialize the new column with NA
  
  # Assign categories based on the presence of keywords
  sorted_df$Category[grepl("tundra", sorted_df$Group)] <- "tundra"
  sorted_df$Category[grepl("taiga, eurasian", sorted_df$Group)] <- "Palearctic boreal forest"
  sorted_df$Category[grepl("taiga, nearctic", sorted_df$Group)] <- "Nearctic boreal forest"
  sorted_df$Category[grepl("other", sorted_df$Group)] <- "other"
  
  sorted_df[is.na(sorted_df)] <- 0
  summary_table <- aggregate(cbind(Count1, Count2, Count3, Count4) ~ Category,
                             data = sorted_df,
                             FUN = sum)
  
  ##### colours ####
  library(scales)
  
  tundra <- hcl.colors(6, "Greens 3")[2] #tundra
  taiga_ne <- hcl.colors(6, "Purples")[2] #america
  taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe
  
  
  # Assuming percent_table is your data frame
  summary_table_long <- summary_table %>%
    pivot_longer(cols = starts_with("Count"), 
                 names_to = "Count_Type", 
                 values_to = "Value")
  
  summary_table_long$year <- c(2010, 2040, 2070, 2100,2010, 2040, 2070, 2100,2010, 2040, 2070, 2100)
  
  custom_colors <- c(taiga_ne, taiga_eu, tundra)
  
  ##### no extinction (ne) ####
  
  ##### clustering ne ####

  y1 <- as.data.frame(spArray[,,1,3])
  y2 <- as.data.frame(spArray[,,2,3])
  y3 <- as.data.frame(spArray[,,3,3])
  y4 <- as.data.frame(spArray[,,4,3])
  
  
  # Create a function to count the number of new 1 cells for each grid cell
  
  make_comparison<- function(present, future) {
    
    comparison = present + future
    comparison[comparison == 2] <- 1
    return(comparison)
  }
  
  present <- y1
  
  for(y in 2:4){
    
    future <- spArray[,,y,3] #year
    dat <- as.data.frame(make_comparison(present, future))
    
    dat[is.na(dat)] <- 0
    cormatrix <- cor (dat)
    
    distancematrix <- cor2dist(cormatrix)
    
    DM2 <- as.matrix(distancematrix)
    
    DM2[cormatrix < 0.3] = 0
    DM2[is.na(DM2)] <- 0
    
    G2a <- graph_from_adjacency_matrix(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)
    cl2 <- cluster_louvain(G2a) 
    
    membercl <- data.frame(group = cl2$membership, label = 1:38657)
    
    write.csv(membercl, glue::glue("{path}membercl_cells_y{y}s3_ne.csv"), row.names = FALSE)
    
  }
  
  ##### name clusters #####
  
  ##### y1 ####
  membercl <- read.csv(glue::glue("{path}membercl_cells_y1s3_sigma1.csv"), sep =",")
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "taiga, nearctic",
      group == 6 ~ "tundra",
      group == 10 ~ NA,
      group == 11 ~ "taiga, eurasian",
      group == 31 ~ NA,
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, file = glue::glue("{path}membercl1_ne.csv"), row.names = F)

  ##### y2 ####
  membercl <- read.csv(glue::glue("{path}membercl_cells_y2s3_ne.csv"), sep =",")

  single_label_groups <- names(which(table(membercl$group) == 1))
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 18 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, file = glue::glue("{path}membercl2_ne.csv"), row.names = F)
  
  ##### y3 ####
  membercl <- read.csv(glue::glue("{path}membercl_cells_y3s3_ne.csv"), sep =",")
  
  single_label_groups <- names(which(table(membercl$group) == 1))
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 14 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, file = glue::glue("{path}membercl3_ne.csv"), row.names = F)
 
  ##### y4 ####
  membercl <- read.csv(glue::glue("{path}membercl_cells_y4s3_ne.csv"), sep =",")
 
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 14 ~ "taiga, eurasian, E", #18
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, file = glue::glue("{path}membercl4_ne.csv"), row.names = F)
  
  ##### read data ####
  membercldat1 <- read.csv(glue::glue("{path}/membercl1_ne.csv"))
  membercldat2 <- read.csv(glue::glue("{path}/membercl2_ne.csv"))
  membercldat3 <- read.csv(glue::glue("{path}/membercl3_ne.csv"))
  membercldat4 <- read.csv(glue::glue("{path}/membercl4_ne.csv"))
  
  summary(as.factor(membercldat4$new_group))
  
  membercldat4 <- membercldat4 %>%
    mutate(new_group = case_when(
      new_group == "14" ~ "taiga, eurasian, E",
      TRUE ~ as.character(new_group)
    ))
  
  membercldat1 <- membercldat1 %>%
    mutate(new_group = case_when(
      new_group == "other" ~ "NA",
      TRUE ~ as.character(new_group)
    ))
  
  membercldat1 <-  na.omit(membercldat1)
  
  summary(as.factor(membercldat1$new_group))
  
  ##### summary table ####
  
  # Convert tables to data frames
  df1 <- as.data.frame(table(membercldat1$new_group))
  df2 <- as.data.frame(table(membercldat2$new_group))
  df3 <- as.data.frame(table(membercldat3$new_group))
  df4 <- as.data.frame(table(membercldat4$new_group))
  
  # Set the names of the data frames to the first column (group names)
  names(df1) <- c("Group", "Count1")
  names(df2) <- c("Group", "Count2")
  names(df3) <- c("Group", "Count3")
  names(df4) <- c("Group", "Count4")
  
  # Merge all data frames by Group
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE), list(df1, df2, df3, df4))
  
  
  # Print the final table
  print(merged_df)
  
  # Sort the merged_df by the first column (Group) alphabetically
  sorted_df <- merged_df[order(merged_df$Group, na.last = TRUE), ]
  
  sorted_df$Category <- NA  # Initialize the new column with NA
  
  # Assign categories based on the presence of keywords
  sorted_df$Category[grepl("tundra", sorted_df$Group)] <- "tundra"
  sorted_df$Category[grepl("taiga, eurasian", sorted_df$Group)] <- "Palearctic boreal forest"
  sorted_df$Category[grepl("taiga, nearctic", sorted_df$Group)] <- "Nearctic boreal forest"
  sorted_df$Category[grepl("other", sorted_df$Group)] <- "other"
  
  sorted_df[is.na(sorted_df)] <- 0
  summary_table <- aggregate(cbind(Count1, Count2, Count3, Count4) ~ Category,
                             data = sorted_df,
                             FUN = sum)
  
  summary_table_long_ne <- summary_table %>%
    pivot_longer(cols = starts_with("Count"), 
                 names_to = "Count_Type", 
                 values_to = "Value")
  
  summary_table_long_ne$year <- c(2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100)
  
  
  ##### comparison plot ####
  D <- ggplot() +
    geom_point(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", position = position_dodge()) +
    geom_line(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dotted", position = position_dodge()) +
    geom_point(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), shape = 17, stat = "identity", position = position_dodge()) +
    geom_line(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dashed", position = position_dodge()) +
    scale_fill_manual(values = c(taiga_eu, taiga_ne, tundra)) +
    labs(title = "d)",
         x = "",
         y = "colonized cells")+
    # color = "Communities") +  # Change legeend title to "Cluster"
    scale_color_manual(values = custom_colors) +
    scale_x_continuous(breaks=c(2010, 2040, 2070, 2100), 
                       labels=c("2010", "2040", "2070", "2100")) + # Custom x-axis breaks and labels
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 16))
  
  print(D)
  
  ##### estimations ####
  
  summary_table_long$Value_ne <- summary_table_long_ne$Value
  
  summary(as.factor(membercldat4$new_group))
  
  d4 <- subset(summary_table_long, summary_table_long$Count_Type == "Count4" | summary_table_long$Count_Type == "Count1")
  
  (d4[2,5] -d4[1,3])/ (d4[2,3] -d4[1,3])  #nearc value 4 - #nearc value 1 /nearc value_ne 4 - #nearc value 1
  # 0.21 of increase nearc 
  
  (d4[4,5] -d4[3,3])/ (d4[4,3] -d4[3,3])  #palearc value 4 - #palearc value 1 /palearc value_ne 4 - #palearc value 1
  # 0.68 of increase palearc
  
  (d4[6,5] -d4[5,3])/ (d4[6,3] -d4[5,3])  #tundra value 4 - #tundra value 1 /tundra value_ne 4 - #tundra value 1
  # 0.5 of decrease tundra
  
##################################################################################
#### plot 4 combined ####
##################################################################################
  
  library(cowplot)

  BCD <- plot_grid(B, C, D, ncol = 3, nrow=1, rel_widths = c(1,2,2))
  ABCD <- plot_grid(A, BCD, ncol= 1, nrow=2)
  
  print(ABCD)
  