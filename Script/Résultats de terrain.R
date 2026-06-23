#######################
# Résultats de terrain
#######################

# Detection Rates ####

#pour l'abondance
library(dplyr)
library(ggplot2)

#Chercher fichier dans qPCR aquarium
dataqPCRV2 <- read.csv("Résultat qPCR aquarium.csv", sep=";", header=T)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Supprimer valeurs invalides
dataqPCRV2 <- dataqPCRV2[!is.na(dataqPCRV2$Ct) & dataqPCRV2$Ct > 0, ]

dataqPCRV2$Detection <- ifelse(dataqPCRV2$Ct < 55, 1, 0)

DetectionRate <- dataqPCRV2 %>%
  group_by(Nb.ecrevisse) %>%
  summarise(
    N_replicates = n(),
    Positive_replicates = sum(Detection),
    Detection_rate = round(100 * mean(Detection), 1)
  )

DetectionRate

# pour les jours
dataqPCRV2 <- read.csv("Résultat qPCR aquarium.csv", sep=";", header=T)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Supprimer valeurs invalides
dataqPCRV2 <- dataqPCRV2[!is.na(dataqPCRV2$Ct) & dataqPCRV2$Ct > 0, ]

dataqPCRV2$Detection <- ifelse(dataqPCRV2$Ct < 55, 1, 0)
DetectionRate2 <- dataqPCRV2 %>%
  group_by(Jour) %>%
  summarise(
    N_replicates = n(),
    Positive_replicates = sum(Detection),
    Detection_rate = round(100 * mean(Detection), 1)
  )

DetectionRate2

# Plot sites : ####

# Import données


library(ggplot2)
library(dplyr)
library(scales)

MoySite <- read.csv(
  "MoyCtSite.csv",
  header = TRUE,
  sep = ";",
  dec = ","
)

MoySite$Ct_value <- as.numeric(gsub(",", ".", MoySite$Ct_value))

MoySite$Ct_value[
  is.na(MoySite$Ct_value) | MoySite$Ct_value == 0
] <- 55

site_names <- c(
  "Trepot",
  "Tarcenay (nord molinaie)",
  "Les Monts-Ronds",
  "Tarcenay ReZo - exterieur etang",
  "Tarcenay ReZo - interieur etang",
  "Le Gratteris",
  "Etalans - etang des Durgeons",
  "Etalans - exutoire etang des Durgeons",
  "Marais de Saone - ruisseau proche etang FDC25",
  "Marais de Saone - etang FDC25",
  "Marais de Saone - ruisseau marais",
  "Saone - etang bord de route"
)

site_codes <- c(
  "1", "2", "3", "4a", "4b", "5",
  "6a", "6b", "7", "8", "9", "10"
)

MoySite <- MoySite %>%
  filter(Site %in% site_names)

MoySite$Site_code <- factor(
  MoySite$Site,
  levels = site_names,
  labels = site_codes
)

p <- ggplot(MoySite,
            aes(x = Site_code,
                y = Ct_value)) +
  
  geom_boxplot(
    fill = "lightblue",
    alpha = 0.7,
    width = 0.3,
    outlier.color = "black"
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 2
  ) +
  
  scale_y_continuous(
    breaks = seq(35, 55, by = 5)
  ) +
  
  coord_cartesian(
    ylim = c(35, 56)
  ) +
  
  labs(
    x = "Sampling site",
    y = "Ct value"
  ) +
  
  theme_classic(base_family = "Times", base_size = 14) +
  
  theme(
    text = element_text(family = "Times", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    panel.grid.major = element_line(
      colour = alpha("grey50", 0.1),
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_line(
      colour = alpha("grey50", 0.1),
      linewidth = 0.25
    )
  )

p

#Moyenne :
MoySite <- read.csv(
  "MoyCtSite.csv",
  header = TRUE,
  sep = ";",
  dec = ","
)

MoySite$Ct_value <- as.numeric(gsub(",", ".", MoySite$Ct_value))

MoySite2 <- MoySite %>%
  filter(!is.na(Ct_value),
         Ct_value <= 42)


mean(MoySite2$Ct_value)
sd(MoySite$Ct_value, na.rm = TRUE)

resume_site <- MoySite2 %>%
  group_by(Site) %>%
  summarise(
    n = n(),
    mean_Ct = mean(True_pos),
    sd_Ct = sd(True_pos),
    min_Ct = min(True_pos),
    max_Ct = max(True_pos),
    .groups = "drop"
  )

resume_site

# Légende encadrée

site_legend <- paste(
  paste(site_codes, "=", site_names),
  collapse = "\n"
)

legend_box <- ggparagraph(
  text = site_legend,
  family = "Times",
  size = 10,
  color = "black"
) +
  theme(
    plot.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.6
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# Figure finale

ggarrange(
  p,
  legend_box,
  ncol = 2,
  widths = c(3, 1.4)
)

# Distance entre les points : ####

library(sf)
library(dplyr)

points <- st_read("points_adne.shp")

points <- st_transform(points, 2154)

if ("NOM" %in% names(points)) {
  site_names <- points$NOM
} else if ("ID" %in% names(points)) {
  site_names <- points$ID
} else {
  site_names <- paste0("Point_", seq_len(nrow(points)))
}

dist_mat <- st_distance(points)
dist_mat <- as.matrix(dist_mat)
dist_mat <- round(dist_mat, 2)

rownames(dist_mat) <- site_names
colnames(dist_mat) <- site_names

write.csv(
  dist_mat,
  "distance_matrix_sites.csv",
  row.names = TRUE
)

dist_table <- as.data.frame(as.table(dist_mat))

names(dist_table) <- c(
  "Site_1",
  "Site_2",
  "Distance_m"
)

dist_table <- dist_table %>%
  filter(Site_1 != Site_2)

write.csv(
  dist_table,
  "distance_table_sites.csv",
  row.names = FALSE
)

dist_mat_no_diag <- dist_mat
diag(dist_mat_no_diag) <- Inf

nearest_table <- data.frame(
  Site = rownames(dist_mat_no_diag),
  Nearest_site = colnames(dist_mat_no_diag)[max.col(-dist_mat_no_diag)],
  Nearest_distance_m = apply(dist_mat_no_diag, 1, min)
)

write.csv(
  nearest_table,
  "nearest_site_distance.csv",
  row.names = FALSE
)

dist_mat
dist_table
nearest_table

#Cor graph 
library(ggplot2)
library(dplyr)

dist_df <- as.data.frame(as.table(dist_mat))

names(dist_df) <- c("Site_1", "Site_2", "Distance_m")

dist_df$Distance_m <- as.numeric(dist_df$Distance_m)

ggplot(dist_df,
       aes(x = Site_1,
           y = Site_2,
           fill = Distance_m)) +
  
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = round(Distance_m, 0)),
    size = 2.5,
    family = "Times",
    colour = "black"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "white",
      "#e3f2fd",
      "#90caf9",
      "#42a5f5",
      "#1565c0"
    ),
    name = "Distance (m)"
  ) +
  
  coord_equal() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_family = "Times", base_size = 12) +
  
  theme(
    text = element_text(family = "Times", colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 9
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 9
    ),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black")
  )


