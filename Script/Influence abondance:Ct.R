##############################
# Influence abondance / [ADNe]
##############################

# Test infuence de l'abondance ####

hist(dataqPCRV2$Ct, breaks=30)
ggplot(dataqPCRV2,
       aes(x = "",
           y = Ct)) +
  
  geom_boxplot(
    fill = "lightblue",
    alpha = 0.7,
    width = 0.3,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.1,
    alpha = 0.6,
    size = 2
  ) +
  
  labs(
    x = NULL,
    y = "Ct value"
  ) +
  
  theme_classic(base_family = "Times", base_size = 14) +
  
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

#Chercher fichier dans qPCR aquarium
dataqPCRV2 <- read.csv("Résultat qPCR aquarium.csv", sep=";", header=T)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Supprimer valeurs invalides
dataqPCRV2 <- dataqPCRV2[!is.na(dataqPCRV2$Ct) & dataqPCRV2$Ct > 0, ]
# Supprimer les Ct > 42
dataqPCRV2 <- dataqPCRV2[dataqPCRV2$Ct <= 42, ]

dataqPCRV2$Nb.écrevisses <- as.factor(dataqPCRV2$Nb.écrevisses)


#Anova
anova_ct <- aov(
  Ct ~ as.factor(Nb.ecrevisse),
  data = dataqPCRV2
)

summary(anova_ct)

TukeyHSD(anova_ct)

library(effectsize)

eta_squared(anova_ct)

# Boxplot abondance/[ADNe] ####
library(dplyr)
library(ggplot2)

# Importation des données
dataqPCRV2 <- read.csv(
  "Résultat qPCR aquarium.csv",
  sep = ";",
  header = TRUE
)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Suppression des valeurs invalides
dataqPCRV2 <- dataqPCRV2 %>%
  filter(
    !is.na(Ct),
    Ct > 0
  )

# Conversion de l'abondance en facteur
dataqPCRV2 <- dataqPCRV2 %>%
  mutate(
    Nb.ecrevisse = as.factor(Nb.ecrevisse)
  )

# Calcul des moyennes
resume <- dataqPCRV2 %>%
  group_by(Nb.ecrevisse) %>%
  summarise(
    Mean_Ct = mean(Ct),
    .groups = "drop"
  )

# Graphique
ggplot(
  dataqPCRV2,
  aes(
    x = Nb.ecrevisse,
    y = Ct
  )
) +
  
  geom_boxplot(
    fill = "lightblue",
    alpha = 0.7,
    width = 0.4,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 2
  ) +
  
  geom_point(
    data = resume,
    aes(
      x = Nb.ecrevisse,
      y = Mean_Ct
    ),
    color = "red",
    size = 4
  ) +
  
  labs(
    x = "Crayfish abundance",
    y = "Ct value"
  ) +
  
  theme_classic(
    base_family = "Times",
    base_size = 14
  ) +
  
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
# Boxplot diff LN et Gus ####
# Packages

library(dplyr)
library(ggplot2)
library(scales)

# Importation des données
dataqPCRV2 <- read.csv("Résultat qPCR aquarium.csv", sep=";", header=T)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Suppression des valeurs invalides
dataqPCRV2 <- dataqPCRV2 %>%
  filter(
    !is.na(Ct),
    Ct > 0
  )

# Conversion des variables en facteur
dataqPCRV2 <- dataqPCRV2 %>%
  mutate(
    Jour = as.factor(Jour),
    Nb.ecrevisse = as.factor(Nb.ecrevisse),
    Qui = as.factor(Qui)
  )

# Graphique Ct ~ Jour
# Colonnes = abondance
# Lignes = expérimentateur

graph_ct_jour <- ggplot(
  dataqPCRV2,
  aes(
    x = Jour,
    y = Ct
  )
) +
  
  geom_boxplot(
    fill = "lightblue",
    alpha = 0.7,
    outlier.color = "black",
    width = 0.5
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 2
  ) +
  
  facet_grid(
    Qui ~ Nb.ecrevisse,
    labeller = labeller(
      Qui = function(x) paste("Experimenter =", x),
      Nb.ecrevisse = function(x) paste("Abundance =", x)
    )
  ) +
  
  labs(
    x = "Days",
    y = "qPCR Ct values"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    panel.grid.major = element_line(
      colour = alpha("grey50", 0.1),
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_line(
      colour = alpha("grey50", 0.1),
      linewidth = 0.25
    )
  )

# Affichage
graph_ct_jour

