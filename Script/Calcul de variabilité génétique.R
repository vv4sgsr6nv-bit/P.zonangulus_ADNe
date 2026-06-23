##################################################################################
# Similarité génétique entre les espèces proche génétiquement et/ou co-occurrentes
##################################################################################

library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

# Import
seqs <- read.dna("pairwise comparison.fasta", format = "fasta")

# Renommer les séquences avec des noms courts
rownames(seqs) <- c(
  "P. zonangulus",
  "P. virginalis",
  "P. clarkii",
  "P. acutus"
)

# Nombre de mutations
mutations <- dist.dna(
  seqs,
  model = "N",
  as.matrix = TRUE
)

# Pourcentage d'identité
dist_prop <- dist.dna(
  seqs,
  model = "raw",
  as.matrix = TRUE
)

identity <- (1 - dist_prop) * 100

# Créer un tableau pour le graphique
plot_df <- as.data.frame(as.table(identity)) %>%
  rename(
    Species_1 = Var1,
    Species_2 = Var2,
    Identity = Freq
  )

mut_df <- as.data.frame(as.table(mutations)) %>%
  rename(
    Species_1 = Var1,
    Species_2 = Var2,
    Mutations = Freq
  )

plot_df <- left_join(
  plot_df,
  mut_df,
  by = c("Species_1", "Species_2")
)

# Graphique matrice
ggplot(plot_df,
       aes(x = Species_1,
           y = Species_2,
           fill = Identity)) +
  
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = paste0(
      round(Identity, 2),
      "%\n(",
      Mutations,
      ")"
    )),
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "white",
      "#f7fbff",
      "#deebf7",
      "#c6dbef",
      "#9ecae1",
      "#6baed6"
    ),
    limits = c(89, 100),
    name = "Identity (%)"
  ) +
  
  coord_equal() +
  
  theme_minimal() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank()
  )