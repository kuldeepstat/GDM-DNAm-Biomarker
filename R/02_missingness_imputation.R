# ==============================================================================
# 02_missingness_imputation.R
#
# Purpose:
# Missingness filtering and imputation of DNA methylation beta values.
#
# Input:
# DataProcessed.Rdata
#   - Processed methylation beta-value dataset
#   - 258 samples
#   - 597,856 CpGs before missingness filtering
#
# Procedure:
#   1. Retain CpGs with >80% non-missing observations
#   2. Impute remaining missing beta values using methyLImp
#   3. Perform imputation in blocks to facilitate computation
#   4. Combine the imputed blocks into the final dataset
#
# Output:
# FDATAV1.Rdata
#   - 258 samples
#   - 597,670 CpG beta-value variables
#   - 1 binary outcome variable (Gr: 0 = Control, 1 = GDM)
#   - Used as input for subsequent EWAS
#
# ==============================================================================


# ==============================================================================
# 1. Load required packages
# ==============================================================================

library(doParallel)
library(methyLImp)


# ==============================================================================
# 2. Load processed methylation data and filter CpGs by missingness
# ==============================================================================

load("data/DataProcessed.Rdata")

# Retain CpGs with >80% non-missing observations
data <- DataProcessed[, which(colMeans(!is.na(DataProcessed)) > 0.8)]

DATANA <- data


# ==============================================================================
# 3. Imputation: CpGs 1-100,000
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,1:10000] <- methyLImp(
  as.matrix(DATANA[,1:10000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()

DATANA[,10001:20000] <- methyLImp(
  as.matrix(DATANA[,10001:20000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,20001:30000] <- methyLImp(
  as.matrix(DATANA[,20001:30000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,30001:40000] <- methyLImp(
  as.matrix(DATANA[,30001:40000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,40001:50000] <- methyLImp(
  as.matrix(DATANA[,40001:50000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,50001:60000] <- methyLImp(
  as.matrix(DATANA[,50001:60000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()

DATANA[,60001:70000] <- methyLImp(
  as.matrix(DATANA[,60001:70000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,70001:80000] <- methyLImp(
  as.matrix(DATANA[,70001:80000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,80001:90000] <- methyLImp(
  as.matrix(DATANA[,80001:90000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,90001:100000] <- methyLImp(
  as.matrix(DATANA[,90001:100000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d1 <- DATANA[,1:100000]

save(d1, file = "data/D1.Rdata")


# ==============================================================================
# 4. Imputation: CpGs 100,001-200,000
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,100001:110000] <- methyLImp(
  as.matrix(DATANA[,100001:110000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,110001:120000] <- methyLImp(
  as.matrix(DATANA[,110001:120000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,120001:130000] <- methyLImp(
  as.matrix(DATANA[,120001:130000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,130001:140000] <- methyLImp(
  as.matrix(DATANA[,130001:140000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,140001:150000] <- methyLImp(
  as.matrix(DATANA[,140001:150000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,150001:160000] <- methyLImp(
  as.matrix(DATANA[,150001:160000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,160001:170000] <- methyLImp(
  as.matrix(DATANA[,160001:170000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,170001:180000] <- methyLImp(
  as.matrix(DATANA[,170001:180000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,180001:190000] <- methyLImp(
  as.matrix(DATANA[,180001:190000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,190001:200000] <- methyLImp(
  as.matrix(DATANA[,190001:200000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d2 <- DATANA[,100001:200000]

save(d2, file = "data/D2.Rdata")


# ==============================================================================
# 5. Imputation: CpGs 200,001-300,000
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,200001:210000] <- methyLImp(
  as.matrix(DATANA[,200001:210000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,210001:220000] <- methyLImp(
  as.matrix(DATANA[,210001:220000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,220001:230000] <- methyLImp(
  as.matrix(DATANA[,220001:230000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,230001:240000] <- methyLImp(
  as.matrix(DATANA[,230001:240000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,240001:250000] <- methyLImp(
  as.matrix(DATANA[,240001:250000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,250001:260000] <- methyLImp(
  as.matrix(DATANA[,250001:260000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,260001:270000] <- methyLImp(
  as.matrix(DATANA[,260001:270000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,270001:280000] <- methyLImp(
  as.matrix(DATANA[,270001:280000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,280001:290000] <- methyLImp(
  as.matrix(DATANA[,280001:290000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,290001:300000] <- methyLImp(
  as.matrix(DATANA[,290001:300000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d3 <- DATANA[,200001:300000]

save(d3, file = "data/D3.Rdata")


# ==============================================================================
# 6. Imputation: CpGs 300,001-400,000
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,300001:310000] <- methyLImp(
  as.matrix(DATANA[,300001:310000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,310001:320000] <- methyLImp(
  as.matrix(DATANA[,310001:320000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,320001:330000] <- methyLImp(
  as.matrix(DATANA[,320001:330000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,330001:340000] <- methyLImp(
  as.matrix(DATANA[,330001:340000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,340001:350000] <- methyLImp(
  as.matrix(DATANA[,340001:350000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,350001:360000] <- methyLImp(
  as.matrix(DATANA[,350001:360000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,360001:370000] <- methyLImp(
  as.matrix(DATANA[,360001:370000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,370001:380000] <- methyLImp(
  as.matrix(DATANA[,370001:380000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,380001:390000] <- methyLImp(
  as.matrix(DATANA[,380001:390000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,390001:400000] <- methyLImp(
  as.matrix(DATANA[,390001:400000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d4 <- DATANA[,300001:400000]

save(d4, file = "data/D4.Rdata")


# ==============================================================================
# 7. Imputation: CpGs 400,001-500,000
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,400001:410000] <- methyLImp(
  as.matrix(DATANA[,400001:410000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,410001:420000] <- methyLImp(
  as.matrix(DATANA[,410001:420000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,420001:430000] <- methyLImp(
  as.matrix(DATANA[,420001:430000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,430001:440000] <- methyLImp(
  as.matrix(DATANA[,430001:440000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,440001:450000] <- methyLImp(
  as.matrix(DATANA[,440001:450000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,450001:460000] <- methyLImp(
  as.matrix(DATANA[,450001:460000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,460001:470000] <- methyLImp(
  as.matrix(DATANA[,460001:470000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,470001:480000] <- methyLImp(
  as.matrix(DATANA[,470001:480000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,480001:490000] <- methyLImp(
  as.matrix(DATANA[,480001:490000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,490001:500000] <- methyLImp(
  as.matrix(DATANA[,490001:500000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d5 <- DATANA[,400001:500000]

save(d5, file = "data/D5.Rdata")


# ==============================================================================
# 8. Imputation: CpGs 500,001-597,670
# ==============================================================================

startTime <- Sys.time()
print(startTime)

DATANA[,500001:510000] <- methyLImp(
  as.matrix(DATANA[,500001:510000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,510001:520000] <- methyLImp(
  as.matrix(DATANA[,510001:520000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,520001:530000] <- methyLImp(
  as.matrix(DATANA[,520001:530000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,530001:540000] <- methyLImp(
  as.matrix(DATANA[,530001:540000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,540001:550000] <- methyLImp(
  as.matrix(DATANA[,540001:550000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,550001:560000] <- methyLImp(
  as.matrix(DATANA[,550001:560000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,560001:570000] <- methyLImp(
  as.matrix(DATANA[,560001:570000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,570001:580000] <- methyLImp(
  as.matrix(DATANA[,570001:580000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,580001:590000] <- methyLImp(
  as.matrix(DATANA[,580001:590000]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)


startTime <- Sys.time()
print(startTime)

DATANA[,590001:597670] <- methyLImp(
  as.matrix(DATANA[,590001:597670]),
  min = 0, max = 1, max.sv = NULL, col.list = NULL
)

endTime <- Sys.time()
print(endTime - startTime)

d6 <- DATANA[,500001:597670]

save(d6, file = "data/D6.Rdata")


# ==============================================================================
# 9. Combine all imputed blocks
# ==============================================================================

FDATAV1 <- as.data.frame(
  cbind(d1, d2, d3, d4, d5, d6)
)


# ==============================================================================
# 10. Create GDM outcome variable
# ==============================================================================

# Extract group membership from the second character of the sample identifier
FDATAV1$Gr <- substr(rownames(FDATAV1), 2, 2)

# Recode group:
#   1 = Control
#   2 = GDM
# as:
#   0 = Control
#   1 = GDM

FDATAV1$Gr <- factor(
  FDATAV1$Gr,
  levels = c("1", "2"),
  labels = c("0", "1")
)

# Check outcome distribution
table(FDATAV1$Gr)


# ==============================================================================
# 11. Save final imputed methylation dataset
# ==============================================================================

save(
  FDATAV1,
  file = "data/FDATAV1.Rdata"
)
