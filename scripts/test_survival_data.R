# Testing survival
library(survminer)
library(survival)

PLOT_DF = matrix_clin_df
PLOT_DF <- select(PLOT_DF, OVERALL_SURVIVAL_DAYS, OVERALL_SURVIVAL_EVENT, PAM50_SUBTYPE)
PLOT_DF$OVERALL_SURVIVAL_DAYS <- PLOT_DF$OVERALL_SURVIVAL_DAYS %>% as.character() %>% as.numeric()
PLOT_DF$OVERALL_SURVIVAL_EVENT <- as.integer(PLOT_DF$OVERALL_SURVIVAL_EVENT)


GSE96058_Clinical_DF_flipsurv <- GSE96058_Clinical_DF
GSE96058_Clinical_DF_flipsurv$OVERALL_SURVIVAL_EVENT <- ifelse(GSE96058_Clinical_DF_flipsurv$OVERALL_SURVIVAL_EVENT == 1, 1, ifelse(GSE96058_Clinical_DF_flipsurv$OVERALL_SURVIVAL_EVENT == 0, 0, NA))
PLOT_DF = GSE96058_Clinical_DF

SURV <- Surv(PLOT_DF$OVERALL_SURVIVAL_DAYS, PLOT_DF$OVERALL_SURVIVAL_EVENT) ~ PLOT_DF$PAM50_SUBTYPE
FIT <- survival::survfit(Surv(OVERALL_SURVIVAL_DAYS, OVERALL_SURVIVAL_EVENT) ~ PAM50_SUBTYPE, data = PLOT_DF)

# KM Plot
KMplot <- ggsurvplot(
  FIT,
  data = PLOT_DF,
  pval = T,
  pval.size = 5,
  pval.method = T,
  pval.method.coord = c(0.1,0.05),
  pval.coord = c(2,0.05, hjust = 0),
  legend.title = "",
  risk.table = F,
  risk.table.height = 0.30,
  risk.table.col = "strata",
  censor = T,
  xlab = "Time (Days)",
  ylab = "Overall Survival", 
  conf.int = F
)

KMplot

plot(FIT)
