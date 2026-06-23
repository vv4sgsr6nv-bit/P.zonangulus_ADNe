# Analyses ADNe du complexe *Procambarus acutus*

Ce dépôt contient les scripts R et les jeux de données utilisés pour analyser les résultats qPCR obtenus au cours des expérimentations en aquarium et des échantillonnages de terrain réalisés dans le cadre du développement et de la validation d'un test qPCR ciblant le complexe *Procambarus acutus* (*P. acutus* et *P. zonangulus*).

## Objectifs

Les analyses réalisées ont pour objectifs de :

* évaluer la sensibilité de l'étude ADNe ;
* analyser l'effet de l'abondance en écrevisses sur les valeurs de Ct ;
* évaluer la persistance de l'ADNe au cours du temps ;
* estimer les taux de détection ;
* analyser les résultats obtenus sur les sites de terrain.
* comparer les séquences des espèces étudiées ;

## Fichiers de données

Les scripts utilisent notamment les fichiers suivants :

* `Résultat qPCR aquarium.csv`
* `LODCalc_Troth_PzonPacu(PZON).csv`
* `LODCalc_Troth_PzonPacu(PACU).csv`
* `pairwise comparison.fasta`
* `MoyCtSite.csv`
* `points_adne.shp`

## Packages nécessaires

Les analyses reposent principalement sur les packages suivants :

```r
library(ggplot2)
library(dplyr)
library(scales)
library(lmerTest)
library(performance)
library(effectsize)
library(ape)
library(tidyr)
library(sf)
```

## Organisation des analyses

### 1. Limite de détection (LOD)

Les relations entre concentration d'ADN et valeurs de Ct sont étudiées à partir de séries de dilution pour *P. zonangulus* et *P. acutus*. Des régressions linéaires sont utilisées afin d'évaluer les performances du test qPCR.

### 2. Visualisation de la spécificité

Les résultats sont réalisés sur différentes espèces sont visualisés afin de vérifier la spécificité de la méthode ADNe. Les non-détections sont codées par une valeur de Ct égale à 55.

### 3. Effets de l'abondance et du temps

Les valeurs de Ct sont analysées en fonction :

* du nombre d'écrevisses présentes ;
* du temps écoulé depuis le début des expérimentations.

Les distributions sont représentées sous forme de boxplots.

### 4. Modèles linéaires mixtes

Des modèles linéaires mixtes sont utilisés afin d'étudier les effets du jour d'échantillonnage et de l'abondance en écrevisses sur les valeurs de Ct, tout en prenant en compte les deux répétitions indépendantes de l'expérience.

Les modèles linéaire, quadratique et avec interaction sont comparés à l'aide des valeurs d'AIC.

### 5. Statistiques descriptives

Les moyennes, écarts-types et effectifs sont calculés pour :

* chaque jour d'échantillonnage ;
* chaque niveau d'abondance.

### 6. Comparaisons de séquences

Des comparaisons par paires entre espèces sont réalisées afin d'estimer :

* le pourcentage d'identité entre séquences ;
* le nombre de mutations observées.

Les résultats sont représentés sous forme de matrice.

### 7. Taux de détection

Les taux de détection sont calculés :

* selon le nombre d'écrevisses ;
* selon les jours d'échantillonnage.

### 8. Analyses des données de terrain

Les données issues des sites de terrain permettent :

* d'étudier la distribution des valeurs de Ct entre sites ;
* de calculer des statistiques descriptives ;
* de représenter les valeurs de Ct par site.

### 9. Distances entre sites

Les coordonnées géographiques des sites sont utilisées afin de :

* calculer les distances entre tous les sites ;
* identifier le site le plus proche de chaque point d'échantillonnage ;
* produire des matrices de distances.

## Prétraitement des données

Avant chaque analyse, les valeurs de Ct sont converties au format numérique et les valeurs invalides sont supprimées.

Pour certaines analyses, les valeurs de Ct supérieures à 42 sont exclues. Les non-détections sont codées par une valeur de Ct égale à 55.

## Reproductibilité

Pour reproduire les analyses :

1. placer tous les fichiers de données dans le répertoire de travail ;
2. installer les packages nécessaires ;
3. ouvrir les scripts dans RStudio ;
4. exécuter les différentes sections dans l'ordre.

## Auteur

**Gurvan Le Borgne**

Laboratoire Écologie et Biologie des Interactions (EBI, UMR CNRS 7267)
Université de Poitiers
2026
