#### altitude mini ####

#### here is the data ####

data <- data.frame(
  scenario = c("sig1", "sig5"),
  mean = c(487.7332, 366.1173),
  q25 = c(184.2325, 126.0122),
  q75 = c(693.5037, 493.4743)
)


#or run after p3_final.R
#input for this scriot is the df_long generated there
  
  s1 <- subset(df_long, df_long$scenario == "sig1")
  s5 <- subset(df_long, df_long$scenario == "sig5")
  s1 <- subset(s1, s1$newspecies > 0)
  s5 <- subset(s5, s5$newspecies > 0)
  
  library(cNORM)
  q1 <- weighted.quantile(s1$altitude, c(0.25, 0.75), weights = s1$newspecies, type = "Harrell-Davis")
  q2 <- weighted.quantile(s5$altitude, c(0.25, 0.75), weights = s5$newspecies, type = "Harrell-Davis")
  
  dot1 <- weighted.mean(s1$altitude, w = s1$newspecies)
  dot2 <- weighted.mean(s5$altitude, w = s5$newspecies)
  
  
  library(ggplot2)
  
  summary_data <- data.frame(
    scenario = c("sig1", "sig5"),
    mean = c(dot1, dot2),
    q25 = c(q1[1], q2[1]),
    q75 = c(q1[2], q2[2])
  )
  
  
#### plotting ####
  
  ggplot(data, aes(x = scenario, y = mean, color = scenario)) +
    geom_segment(aes(x = scenario, xend = scenario, y = q25, yend = q75), 
                 size = 1, color = "grey50") +
    geom_point(size = 6) + #maybe larger later when we add into the other plot
    scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 200)) +
    scale_color_manual(values = c("sig1" = "deepskyblue2", "sig5" = "dodgerblue4")) +
    scale_x_discrete(expand = c(0, 2)) +  # Adjust this value to bring lines closer
    labs(x = NULL, y = "Altitude [m]") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(), 
      axis.text.x = element_text(size = 12),
      legend.position = "none",
      panel.grid.major.y = element_blank(),  # Remove vertical grid lines (horizontal after flip)
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "lightgrey", linetype = "solid"),  # Keep horizontal grid lines (vertical after flip)
      panel.grid.minor.x = element_blank()
    ) +
    coord_flip()
    
    
    