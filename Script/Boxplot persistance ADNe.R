###########################################
# Analyses de la persistance de l'ADNe
###########################################

library(dplyr)

#Chercher fichier dans qPCR aquarium
dataqPCRV2 <- read.csv("Résultat qPCR aquarium.csv", sep=";", header=T)

# Nettoyage des Ct
dataqPCRV2$Ct <- gsub(",", ".", dataqPCRV2$Ct)
dataqPCRV2$Ct <- as.numeric(dataqPCRV2$Ct)

# Supprimer valeurs invalides
dataqPCRV2 <- dataqPCRV2[!is.na(dataqPCRV2$Ct) & dataqPCRV2$Ct > 0, ]
# Supprimer les Ct > 42
dataqPCRV2 <- dataqPCRV2[dataqPCRV2$Ct <= 42, ]

library(lmerTest)

mod.lin <- lmer(
  Ct ~ Jour + Nb.ecrevisse + (1|Qui),
  data = dataqPCRV2
)

summary(mod.lin)
anova(mod.lin)
library(performance)

r2(mod.lin)

mod.quad <- lmer(
  Ct ~ Jour + I(Jour^2) + Nb.ecrevisse + (1|Qui),
  data = dataqPCRV2
)

summary(mod.quad)
anova(mod.quad)

AIC(mod.lin, mod.quad)

mod.int <- lmer(
  Ct ~ Jour * Nb.ecrevisse + (1|Qui),
  data = dataqPCRV2
)

summary(mod.int)
anova(mod.int)

AIC(mod.lin,
    mod.quad,
    mod.int)

# Moyennes des graphs 
library(dplyr)

moyennes_jour <- dataqPCRV2 %>%
  group_by(Jour) %>%
  summarise(
    Mean_Ct = mean(Ct, na.rm = TRUE),
    SD_Ct = sd(Ct, na.rm = TRUE),
    N = n()
  )

moyennes_jour

moyennes_ecrevisses <- dataqPCRV2 %>%
  group_by(Nb.ecrevisse) %>%
  summarise(
    Mean_Ct = mean(Ct, na.rm = TRUE),
    SD_Ct = sd(Ct, na.rm = TRUE),
    N = n()
  )

moyennes_ecrevisses





##################################
# Représentation graphique 
##################################
library(dplyr)
library(ggplot2)
library(scales)

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
dataqPCRV2 <- dataqPCRV2[
  !is.na(dataqPCRV2$Ct) &
    dataqPCRV2$Ct > 0,
]

# Suppression des Ct > 42
dataqPCRV2 <- dataqPCRV2[
  dataqPCRV2$Ct <= 42,
]

# Conversion en facteurs
dataqPCRV2$Jour <- as.factor(dataqPCRV2$Jour)
dataqPCRV2$Nb.ecrevisse <- as.factor(dataqPCRV2$Nb.ecrevisse)

# Graphique de persistance
ggplot(
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
  
  facet_wrap(
    ~ Nb.ecrevisse,
    nrow = 1
  ) +
  
  labs(
    x = "Days",
    y = "qPCR Ct values"
  ) +
  
  theme_classic(base_size = 14) +
  
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

