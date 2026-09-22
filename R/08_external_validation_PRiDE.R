# External validation analysis
# GitHub-ready version: original analytical code retained.
# Only the local input-data path was changed to a repository-relative path,
# and the incomplete trailing plotting fragment from the supplied script was removed.

library(pROC)
library(ggplot2)
library(plotROC)
library(dplyr)
library(readxl)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)


d <- read_excel("data/external_validation/PRiDE/PRIDE FULL 128 _GDM_Methylation validation_DATA_June 23 2026_CAUCASIAN ONLY_KULDEEP.xlsx")


names(d)


#Models

Model1a <- glm(
  as.factor(GDMgroup_Code) ~
    cg04985016_M2_TRIP12 +
    cg11001216_M3_A1BG +
    cg04539775_M4_CLSTN3 +
    cg10139015_M5_PSMC3 +
    cg09518293_M6_TTN +
    cg17967426_M7_PANK3 +
    cg00078968_M8_WNT3A +
    cg15980170_M9_ZNF778 +
    cg09217522_M10_RECQL +
    cg24745753_T02_CTNNA2_M11 +
    cg16311883_T13_ZNF664_M12,
  data = d,
  family = binomial
)

Model1b <- glm(as.factor(GDMgroup_Code) ~ cg04985016_M2_TRIP12 +
                 cg11001216_M3_A1BG +
                 cg04539775_M4_CLSTN3 +
                 cg10139015_M5_PSMC3 +
                 cg09518293_M6_TTN +
                 cg17967426_M7_PANK3 +
                 cg00078968_M8_WNT3A +
                 cg15980170_M9_ZNF778 +
                 cg09217522_M10_RECQL +
                 cg24745753_T02_CTNNA2_M11 +
                 cg16311883_T13_ZNF664_M12 + age + BMI + as.factor(v1r_Inc_FH), data = d, family = binomial)




pred1a <- predict(Model1a, type = "response")
pred1b <- predict(Model1b, type = "response")


#ROC objects

library(pROC)

roc1a <- roc(Model1a$y, pred1a)
roc1b <- roc(Model1b$y, pred1b)

#AUC and 95% CI
auc1a <- auc(roc1a)
auc1b <- auc(roc1b)


ci1a <- ci.auc(roc1a)
ci1b <- ci.auc(roc1b)



#Sensitivity and specificity at threshold = 0.5
pt1a <- coords(roc1a,
               x=0.5,
               input="threshold",
               ret=c("specificity","sensitivity"))

pt1b <- coords(roc1b,
               x=0.5,
               input="threshold",
               ret=c("specificity","sensitivity"))


#ROC plot
plot(roc1a,
     col="#666666",
     lwd=7,
     legacy.axes=TRUE,
     print.auc=FALSE)

plot(roc1b,
     add=TRUE,
     col="#D98A9D",
     lwd=7)


# Set thicker axes and larger tick labels
par(
  lwd = 3,          # Axis line thickness
  lwd.axis = 3,     # Tick mark thickness
  cex.axis = 1.5,   # Tick label size
  cex.lab = 1.8,    # Axis title size
  mgp = c(3, 1, 0), # Axis label spacing
  tcl = -0.5        # Tick mark length
)

plot(
  roc1a,
  col = "#666666",
  lwd = 7,
  legacy.axes = TRUE,
  print.auc = FALSE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity"
)

plot(
  roc1b,
  add = TRUE,
  col = "#D98A9D",
  lwd = 7
)



#Automatic legend
legend(
  "bottomright",
  
  legend = c(
    
    sprintf(
      "Model 1a: AUC %.3f (%.3f-%.3f)\n   Sensitivity=%.1f%%   Specificity=%.1f%%",
      auc1a, ci1a[1], ci1a[3],
      100*pt1a["sensitivity"],
      100*pt1a["specificity"]
    ),
    
    sprintf(
      "Model 1b: AUC %.3f (%.3f-%.3f)\n   Sensitivity=%.1f%%   Specificity=%.1f%%",
      auc1b, ci1b[1], ci1b[3],
      100*pt1b["sensitivity"],
      100*pt1b["specificity"]
    )
  ),
  
  col = c("#666666", "#D98A9D"),
  
  lwd = 3,
  
  bty = "n",
  
  cex = 1.2,
  
  y.intersp = 2.0,      # vertical spacing
  x.intersp = 1.2,      # horizontal spacing
  seg.len = 3,          # longer color lines
  inset = c(0.02,0.02) # margin from border
)
