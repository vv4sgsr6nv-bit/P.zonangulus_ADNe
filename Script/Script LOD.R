#### LoD / LoQ analysis - P. zonangulus ####
## Corrected version for LODCalc_Troth_PzonPacu(PZON).csv

## Required packages
library(ggplot2)
library(drc)

## Import data
Data <- read.table(
  "LODCalc_Troth_PzonPacu(PZON).csv",
  sep = ";",
  dec = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

## Remove empty rows instead of removing a fixed number of rows
Data <- Data[!is.na(Data$Target) & !is.na(Data$SQ), ]

DAT <- Data

## Parameters
LOQ.Threshold <- 0.7
LOD.FCT <- "Best"
LOQ.FCT <- "Best"

## Create analysis log file
write(paste0("Analysis started: ", date(), "\n\n"), file = "Analysis Log.txt")

## Check column names
if(sum(colnames(DAT) == "Target") != 1) {
  A <- grep("target", colnames(DAT), ignore.case = TRUE)
  if(length(A) == 1) colnames(DAT)[A] <- "Target"
  if(length(A) != 1) stop("Problem with the 'Target' column.")
}

if(sum(colnames(DAT) == "Cq") != 1) {
  A <- grep("cq|ct|cycle", colnames(DAT), ignore.case = TRUE)
  if(length(A) == 1) colnames(DAT)[A] <- "Cq"
  if(length(A) != 1) stop("Problem with the 'Cq' column.")
}

if(sum(colnames(DAT) == "SQ") != 1) {
  A <- grep("sq|copies|starting|quantity", colnames(DAT), ignore.case = TRUE)
  if(length(A) == 1) colnames(DAT)[A] <- "SQ"
  if(length(A) != 1) stop("Problem with the 'SQ' column.")
}

## Format data
DAT$Target <- as.factor(DAT$Target)
DAT$Cq <- suppressWarnings(as.numeric(as.character(DAT$Cq)))
DAT$SQ <- suppressWarnings(as.numeric(as.character(DAT$SQ)))

## Important correction: Cq = 0 corresponds to non-detection, not a true Cq value
DAT$Cq[DAT$Cq == 0] <- NA

## Remove rows without valid SQ
DAT <- DAT[!is.na(DAT$SQ), ]

Targets <- unique(DAT$Target)
Standards <- unlist(lapply(Targets, function(t) unique(DAT$SQ[DAT$Target == t & !is.na(DAT$SQ)])))
Target <- unlist(lapply(Targets, function(t) rep(as.character(t), length(unique(DAT$SQ[DAT$Target == t & !is.na(DAT$SQ)])))))

## Standard curves and copy estimates
curve.list <- character(0)
DAT$Copy.Estimate <- NA
DAT$Mod <- 0

for(i in seq_along(Targets)) {
  STDS <- data.frame(S = unique(DAT$SQ[DAT$Target == Targets[i]]), R = NA)
  
  for(j in seq_len(nrow(STDS))) {
    STDS$R[j] <- sum(!is.na(DAT$Cq) & DAT$SQ == STDS$S[j] & DAT$Target == Targets[i]) /
      sum(DAT$SQ == STDS$S[j] & DAT$Target == Targets[i])
  }
  
  if(sum(STDS$R >= 0.5, na.rm = TRUE) >= 3) {
    STDS2 <- STDS$S[STDS$R >= 0.5 & !is.na(STDS$R) & !is.na(STDS$S)]
  } else {
    STDS2 <- STDS$S[order(STDS$R, decreasing = TRUE)][1:3]
  }
  
  for(j in seq_along(STDS2)) {
    D <- DAT$Cq[DAT$Target == Targets[i] & DAT$SQ == STDS2[j]]
    q <- quantile(D, na.rm = TRUE)
    DAT$Mod[DAT$Target == Targets[i] & DAT$SQ == STDS2[j] &
              DAT$Cq >= q[2] & DAT$Cq <= q[4] & !is.na(DAT$Cq)] <- 1
  }
  
  curve.name <- paste0("curve", i)
  assign(curve.name, lm(Cq ~ log10(SQ), data = DAT[DAT$Target == Targets[i] & DAT$Mod == 1, ]))
  curve.list <- c(curve.list, curve.name)
  
  Intercept <- coef(get(curve.name))[1]
  Slope <- coef(get(curve.name))[2]
  DAT$Copy.Estimate[DAT$Target == Targets[i]] <- 10^((DAT$Cq[DAT$Target == Targets[i]] - Intercept) / Slope)
}

## Summarize data
DAT2 <- data.frame(
  Standards = Standards,
  Target = Target,
  Reps = NA,
  Detects = NA,
  Cq.mean = NA,
  Cq.sd = NA,
  Copy.CV = NA,
  Cq.CV = NA
)

for(i in seq_len(nrow(DAT2))) {
  idx <- DAT$SQ == DAT2$Standards[i] & DAT$Target == DAT2$Target[i]
  DAT2$Reps[i] <- sum(idx, na.rm = TRUE)
  DAT2$Detects[i] <- sum(!is.na(DAT$Cq) & idx, na.rm = TRUE)
  DAT2$Cq.mean[i] <- mean(DAT$Cq[idx], na.rm = TRUE)
  DAT2$Cq.sd[i] <- sd(DAT$Cq[idx], na.rm = TRUE)
  DAT2$Copy.CV[i] <- sd(DAT$Copy.Estimate[idx], na.rm = TRUE) / mean(DAT$Copy.Estimate[idx], na.rm = TRUE)
  DAT2$Cq.CV[i] <- sqrt(2^(DAT2$Cq.sd[i]^2 * log(2)) - 1)
}

DAT2$Rate <- DAT2$Detects / DAT2$Reps
write.csv(DAT2, file = "Data_summary.csv", row.names = FALSE)

## Assay summary
DAT$Detect <- as.numeric(!is.na(DAT$Cq))
DAT3 <- data.frame(
  Assay = Targets,
  R.squared = NA,
  Slope = NA,
  Intercept = NA,
  Low.95 = NA,
  LOD = NA,
  LOQ = NA,
  rep2.LOD = NA,
  rep3.LOD = NA,
  rep4.LOD = NA,
  rep5.LOD = NA,
  rep8.LOD = NA
)

LOD.FCTS <- list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), W1.2(), W1.3(), W1.4(), W2.2(), W2.3(), W2.4(), AR.2(), AR.3(), MM.2(), MM.3())
LOD.list2 <- rep(NA_character_, length(Targets))
LOD.list3 <- rep(NA_character_, length(Targets))
LOD.CI <- NULL

for(i in seq_along(Targets)) {
  DAT3$R.squared[i] <- summary(get(curve.list[i]))$r.squared
  DAT3$Slope[i] <- coef(get(curve.list[i]))[2]
  DAT3$Intercept[i] <- coef(get(curve.list[i]))[1]
  DAT3$Low.95[i] <- min(DAT2$Standards[DAT2$Rate >= 0.95 & DAT2$Target == Targets[i]], na.rm = TRUE)
  
  informative <- sum(DAT2$Rate[DAT2$Target == Targets[i]] != 1 & DAT2$Rate[DAT2$Target == Targets[i]] != 0)
  
  if(informative == 0) {
    msg <- paste0("WARNING: For ", Targets[i], ", all standards detected fully or failed fully. The LoD dose-response model cannot converge. LoD should be reported as < ", DAT3$Low.95[i], " copies/reaction or retested with lower concentrations.")
    message(msg)
    write(paste0(msg, "\n\n"), file = "Analysis Log.txt", append = TRUE)
  } else {
    tryCatch({
      TEMP.DAT <- DAT[DAT$Target == Targets[i], ]
      LOD.mod <- drm(Detect ~ SQ, data = TEMP.DAT, fct = W2.4())
      LOD.FCT2 <- row.names(mselect(LOD.mod, LOD.FCTS))[1]
      LOD.FCT3 <- getMeanFunctions(fname = LOD.FCT2)
      mod.name <- paste0("LOD.mod2", i)
      assign(mod.name, drm(Detect ~ SQ, data = TEMP.DAT, fct = LOD.FCT3[[1]]))
      LOD.list2[i] <- mod.name
      LOD.list3[i] <- LOD.FCT2
      
      DAT3$LOD[i] <- ED(get(mod.name), 0.95, type = "absolute")[1]
      DAT3$rep2.LOD[i] <- ED(get(mod.name), 1 - sqrt(0.05), type = "absolute")[1]
      DAT3$rep3.LOD[i] <- ED(get(mod.name), 1 - 0.05^(1/3), type = "absolute")[1]
      DAT3$rep4.LOD[i] <- ED(get(mod.name), 1 - 0.05^0.25, type = "absolute")[1]
      DAT3$rep5.LOD[i] <- ED(get(mod.name), 1 - 0.05^0.2, type = "absolute")[1]
      DAT3$rep8.LOD[i] <- ED(get(mod.name), 1 - 0.05^0.125, type = "absolute")[1]
      
      tmpCI <- rbind(
        ED(get(mod.name), 0.95, interval = "delta", type = "absolute"),
        ED(get(mod.name), 1 - sqrt(0.05), interval = "delta", type = "absolute"),
        ED(get(mod.name), 1 - 0.05^(1/3), interval = "delta", type = "absolute"),
        ED(get(mod.name), 1 - 0.05^0.25, interval = "delta", type = "absolute"),
        ED(get(mod.name), 1 - 0.05^0.2, interval = "delta", type = "absolute"),
        ED(get(mod.name), 1 - 0.05^0.125, interval = "delta", type = "absolute")
      )
      tmpCI <- data.frame(tmpCI, LoD = c("1rep.LOD", "2rep.LOD", "3rep.LOD", "4rep.LOD", "5rep.LOD", "8rep.LOD"), Assay = Targets[i])
      LOD.CI <<- rbind(LOD.CI, tmpCI)
    }, error = function(e) {
      msg <- paste0("ERROR: LOD model cannot be defined for ", Targets[i], ": ", e$message)
      message(msg)
      write(paste0(msg, "\n\n"), file = "Analysis Log.txt", append = TRUE)
    })
  }
}

## Export summaries
write.csv(DAT3, file = "Assay_summary.csv", row.names = FALSE)

## Export LOD confidence intervals only if the object exists
if(!is.null(LOD.CI)) {
  write.csv(LOD.CI, file = "LOD_confint.csv", row.names = FALSE)
} else {
  message("LOD_confint.csv was not created because no LoD model converged.")
}

## Standard curve plots
DAT$Mod <- ifelse(DAT$Mod == 1, "Modeled", "Excluded")

for(i in seq_along(Targets)) {
  ggOut <- ggplot(data = DAT[DAT$Target == Targets[i] & !is.na(DAT$SQ), ],
                  aes(x = SQ, y = Cq, color = factor(Mod), shape = factor(Mod), size = factor(Mod))) +
    geom_jitter(width = 0.1, alpha = 0.75, na.rm = TRUE) +
    scale_shape_manual("", values = c(3, 20)) +
    scale_size_manual("", values = c(1, 3)) +
    scale_x_log10() +
    xlab("Standard concentrations (copies/reaction)") +
    ylab("Ct value") +
    geom_abline(intercept = coef(get(curve.list[i]))[1], slope = coef(get(curve.list[i]))[2]) +
    theme_bw() +
    annotate("text",
             y = min(DAT$Cq[DAT$Target == Targets[i]], na.rm = TRUE) * 1.05,
             x = min(DAT$SQ[DAT$Target == Targets[i]], na.rm = TRUE) * 1.01,
             hjust = 0,
             label = paste0("R-squared: ", round(DAT3$R.squared[i], 3),
                            "\ny = ", round(DAT3$Slope[i], 3), "x + ", round(DAT3$Intercept[i], 3),
                            "\nLow.95 = ", DAT3$Low.95[i], " copies/reaction"))
  print(ggOut)
}

library(ggplot2)
library(dplyr)

# Import P. zonangulus
LODdata <- read.csv(
  "LODCalc_Troth_PzonPacu(PZON).csv",
  header = TRUE,
  dec = ",",
  sep = ";"
)

LODdata$Cq <- as.numeric(gsub(",", ".", LODdata$Cq))
LODdata$LOG_SQ <- as.numeric(gsub(",", ".", LODdata$LOG_SQ))
LODdata <- LODdata[!is.na(LODdata$Cq) & LODdata$Cq > 0, ]



# Import P. acutus

LODdata2 <- read.csv(
  "LODCalc_Troth_PzonPacu(PACU).csv",
  header = TRUE,
  dec = ",",
  sep = ";"
)

LODdata2$Cq <- as.numeric(gsub(",", ".", LODdata2$Cq))
LODdata2$LOG_SQ <- as.numeric(gsub(",", ".", LODdata2$LOG_SQ))
LODdata2 <- LODdata2[!is.na(LODdata2$Cq) & LODdata2$Cq > 0, ]

# Ajouter noms des espèces

LODdata$espece  <- "P. zonangulus"
LODdata2$espece <- "P. acutus"

LOD_all <- bind_rows(LODdata, LODdata2)

# Graphique final

ggplot(LOD_all, aes(x = LOG_SQ, y = Cq, color = espece)) +
  
  geom_point(size = 2) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  
  labs(
    x = expression(Log[10]~DNA~concentration),
    y = "Ct value",
    color = NULL
  ) +
  
  theme_bw(base_family = "Times") +
  
  theme(
    panel.background = element_rect(
      fill = "grey97",
      colour = NA
    ),
    
    panel.grid.major = element_line(
      colour = "grey85",
      linewidth = 0.4
    ),
    
    panel.grid.minor = element_line(
      colour = "grey92",
      linewidth = 0.2
    ),
    
    text = element_text(color = "black"),
    
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    
    axis.title = element_text(
      color = "black",
      size = 14
    ),
    
    legend.position = "right",
    
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.8
    ),
    
    legend.key = element_rect(
      fill = "white",
      colour = "white"
    ),
    
    legend.text = element_text(
      colour = "black",
      size = 11
    ),
    
    legend.title = element_blank()
  )
