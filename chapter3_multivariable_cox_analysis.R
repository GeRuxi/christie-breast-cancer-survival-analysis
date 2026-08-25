### Chapter 3: Multivariable Cox regression and diagnostics

rm(list = ls())
library(survival)

## Output setup

output_root <- path.expand("~/Downloads/chapter3_outputs")
fig_dir <- file.path(output_root, "figures")
tab_dir <- file.path(output_root, "tables")
report_dir <- file.path(output_root, "reports")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

## Data preparation

load("/Users/geruxi/Documents/硕士所有东西/硕士毕业论文/老师给的/OneDrive_1_2026-6-26/The breast cancer survival data/Data to use for project/bc survival data - to be used.RData")
bc <- bcdata.use

predictors <- setdiff(names(bc), c("stime", "cens"))
bc[predictors] <- lapply(bc[predictors], factor)

stopifnot(all(bc$stime > 0), all(bc$cens %in% c(0, 1)))
stopifnot(!anyNA(bc[, c("stime", "cens", predictors)]))

cat("Sample size and censoring status:\n")
print(table(cens = bc$cens))
cat("\nNumber of levels in each covariate:\n")
print(sapply(bc[predictors], nlevels))

## Full Cox model

full_formula <- as.formula(
  paste("Surv(stime, cens) ~", paste(predictors, collapse = " + "))
)

cox_null <- coxph(
  Surv(stime, cens) ~ 1,
  data = bc,
  ties = "exact"
)

cox_full <- coxph(
  full_formula,
  data = bc,
  ties = "exact",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

summary(cox_full)

full_term_tests <- drop1(cox_full, test = "Chisq")
full_term_tests

## AIC reduction

cox_step <- step(
  cox_full,
  direction = "both",
  trace = 1
)

cox_reduced <- update(
  cox_step,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

cat("\nSelected reduced model:\n")
print(formula(cox_reduced))
summary(cox_reduced)

reduced_term_tests <- drop1(cox_reduced, test = "Chisq")
reduced_term_tests

# Model comparison

aic_comparison <- AIC(cox_null, cox_full, cox_reduced)
aic_comparison

lr_comparison <- anova(cox_reduced, cox_full, test = "Chisq")
lr_comparison

# Coefficient tables and model output

cox_table <- function(fit) {
  s <- summary(fit)
  data.frame(
    term = rownames(s$coefficients),
    coefficient = s$coefficients[, "coef"],
    hazard_ratio = s$conf.int[, "exp(coef)"],
    lower_95_CI = s$conf.int[, "lower .95"],
    upper_95_CI = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

full_results <- cox_table(cox_full)
reduced_results <- cox_table(cox_reduced)

write.csv(full_results, file.path(tab_dir, "chapter3_full_cox_results.csv"), row.names = FALSE)
write.csv(reduced_results, file.path(tab_dir, "chapter3_reduced_cox_results.csv"), row.names = FALSE)
write.csv(aic_comparison, file.path(tab_dir, "chapter3_aic_comparison.csv"))
write.csv(full_term_tests, file.path(tab_dir, "chapter3_full_term_tests.csv"))
write.csv(reduced_term_tests, file.path(tab_dir, "chapter3_reduced_term_tests.csv"))

capture.output({
  cat("FULL COX MODEL\n")
  print(summary(cox_full))
  cat("\nFULL-MODEL WHOLE-TERM TESTS\n")
  print(full_term_tests)

  cat("\nSTEPWISE AIC PATH\n")
  print(cox_step$anova)

  cat("\nREDUCED COX MODEL\n")
  print(summary(cox_reduced))
  cat("\nREDUCED-MODEL WHOLE-TERM TESTS\n")
  print(reduced_term_tests)

  cat("\nAIC COMPARISON\n")
  print(aic_comparison)
  cat("\nLIKELIHOOD-RATIO COMPARISON: REDUCED VS FULL\n")
  print(lr_comparison)
}, file = file.path(report_dir, "chapter3_model_output.txt"))

## Reduced-model proportional-hazards diagnostics

ph_test <- cox.zph(cox_reduced)
print(ph_test)

capture.output(ph_test, file = file.path(report_dir, "chapter3_ph_test.txt"))

pdf(file.path(fig_dir, "chapter3_schoenfeld_residuals.pdf"), width = 7, height = 5)
plot(ph_test)
dev.off()

# Residual diagnostics

martingale_residuals <- residuals(cox_reduced, type = "martingale")
deviance_residuals <- residuals(cox_reduced, type = "deviance")
linear_predictor <- predict(cox_reduced, type = "lp")

cox_snell_residuals <- bc$cens - martingale_residuals
cox_snell_fit <- survfit(Surv(cox_snell_residuals, bc$cens) ~ 1)
cox_snell_hazard <- -log(cox_snell_fit$surv)
keep <- is.finite(cox_snell_hazard)

pdf(file.path(fig_dir, "chapter3_cox_snell_residuals.pdf"), width = 6, height = 5)
plot(
  cox_snell_fit$time[keep],
  cox_snell_hazard[keep],
  type = "s",
  xlab = "Cox--Snell residual",
  ylab = "Estimated cumulative hazard",
  main = "Cox--Snell residual diagnostic"
)
abline(0, 1, lty = 2)
dev.off()

pdf(file.path(fig_dir, "chapter3_martingale_deviance_residuals.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(
  linear_predictor,
  martingale_residuals,
  xlab = "Linear predictor",
  ylab = "Martingale residual",
  main = "Martingale residuals"
)
abline(h = 0, lty = 2)
lines(lowess(linear_predictor, martingale_residuals), lwd = 2)

plot(
  linear_predictor,
  deviance_residuals,
  xlab = "Linear predictor",
  ylab = "Deviance residual",
  main = "Deviance residuals"
)
abline(h = 0, lty = 2)
lines(lowess(linear_predictor, deviance_residuals), lwd = 2)

dev.off()

# Efron ties are used only for DFBETA because DFBETA is unavailable for the exact fit.

# Influence diagnostics for reduced model

cox_influence <- coxph(
  formula(cox_reduced),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

influence_coef_comparison <- data.frame(
  term = names(coef(cox_reduced)),
  exact = unname(coef(cox_reduced)),
  efron = unname(coef(cox_influence))
)

influence_coef_comparison$difference <-
  influence_coef_comparison$efron -
  influence_coef_comparison$exact

print(influence_coef_comparison)

write.csv(
  influence_coef_comparison,
  file.path(tab_dir, "chapter3_exact_efron_coefficient_comparison.csv"),
  row.names = FALSE
)

dfbeta_residuals <- residuals(
  cox_influence,
  type = "dfbeta"
)

max_abs_dfbeta <- if (is.null(dim(dfbeta_residuals))) {
  abs(dfbeta_residuals)
} else {
  apply(abs(dfbeta_residuals), 1, max)
}

pdf(file.path(fig_dir, "chapter3_influence_dfbeta.pdf"), width = 7, height = 5)

plot(
  max_abs_dfbeta,
  type = "h",
  xlab = "Observation number",
  ylab = "Maximum absolute DFBETA",
  main = "Influence diagnostic"
)

dev.off()


plot(ph_test, var = "presite.new")
abline(h = 0, lty = 2)

plot(ph_test, var = "nodes.axilla")
abline(h = 0, lty = 2)


# Selected Schoenfeld residual plots

ph_coef <- cox.zph(
  cox_reduced,
  terms = FALSE
)

print(ph_coef)

pdf(
  file.path(fig_dir, "chapter3_reduced_schoenfeld_selected_coefficients.pdf"),
  width = 11,
  height = 4.5
)

par(mfrow = c(1, 3))

plot(
  ph_coef,
  var = "presite.new2",
  resid = TRUE,
  se = TRUE
)
abline(
  h = unname(coef(cox_reduced)["presite.new2"]),
  lty = 3
)

plot(
  ph_coef,
  var = "presite.new9",
  resid = TRUE,
  se = TRUE
)
abline(
  h = unname(coef(cox_reduced)["presite.new9"]),
  lty = 3
)

plot(
  ph_coef,
  var = "nodes.axilla2",
  resid = TRUE,
  se = TRUE
)
abline(
  h = unname(coef(cox_reduced)["nodes.axilla2"]),
  lty = 3
)

par(mfrow = c(1, 1))

dev.off()


## Stratification by predominant tumour site

cox_stratified <- update(
  cox_reduced,
  . ~ . - presite.new + strata(presite.new),
  x = TRUE,
  y = TRUE,
  model = TRUE
)

formula(cox_stratified)
summary(cox_stratified)

ph_stratified <- cox.zph(
  cox_stratified,
  terms = TRUE
)

print(ph_stratified)

pdf(
  file.path(fig_dir, "chapter3_schoenfeld_residuals_stratified.pdf"),
  width = 7,
  height = 5
)

plot(ph_stratified)

dev.off()


out_pdf <- file.path(
  fig_dir,
  "chapter3_schoenfeld_residuals_stratified.pdf"
)

pdf(
  file = out_pdf,
  width = 7,
  height = 5,
  onefile = TRUE
)

plot(ph_stratified)

dev.off()

cat("PDF saved to:\n", out_pdf, "\n")
cat("File successfully created:", file.exists(out_pdf), "\n")


ph_stratified_coef <- cox.zph(
  cox_stratified,
  terms = FALSE
)

print(ph_stratified_coef)

write.table(
  ph_stratified_coef$table,
  file = file.path(
    report_dir,
    "chapter3_ph_test_stratified_coefficients.txt"
  ),
  quote = FALSE,
  sep = "\t"
)


## Final stratified model

final_model <- cox_stratified


# Final proportional-hazards diagnostics

ph_final <- cox.zph(
  final_model,
  terms = TRUE
)

print(ph_final)

write.table(
  ph_final$table,
  file = file.path(report_dir, "chapter3_final_ph_test.txt"),
  quote = FALSE,
  sep = "\t"
)

pdf(
  file.path(fig_dir, "chapter3_final_schoenfeld_residuals.pdf"),
  width = 7,
  height = 5,
  onefile = TRUE
)

plot(ph_final)

dev.off()


# Martingale and deviance residuals

martingale_residuals <- residuals(
  final_model,
  type = "martingale"
)

deviance_residuals <- residuals(
  final_model,
  type = "deviance"
)

linear_predictor <- predict(
  final_model,
  type = "lp"
)

pdf(
  file.path(fig_dir, "chapter3_final_martingale_deviance.pdf"),
  width = 10,
  height = 5
)

par(mfrow = c(1, 2))

plot(
  linear_predictor,
  martingale_residuals,
  xlab = "Linear predictor",
  ylab = "Martingale residual",
  main = "Martingale residuals"
)

abline(h = 0, lty = 2)

lines(
  lowess(linear_predictor, martingale_residuals),
  lwd = 2
)

plot(
  linear_predictor,
  deviance_residuals,
  xlab = "Linear predictor",
  ylab = "Deviance residual",
  main = "Deviance residuals"
)

abline(h = 0, lty = 2)

lines(
  lowess(linear_predictor, deviance_residuals),
  lwd = 2
)

dev.off()

par(mfrow = c(1, 1))


# Cox-Snell residuals

cox_snell_residuals <-
  bc$cens - martingale_residuals

cox_snell_fit <- survfit(
  Surv(cox_snell_residuals, bc$cens) ~ 1
)

cox_snell_hazard <-
  -log(cox_snell_fit$surv)

keep <- is.finite(cox_snell_hazard)

pdf(
  file.path(fig_dir, "chapter3_final_cox_snell.pdf"),
  width = 6,
  height = 5
)

plot(
  cox_snell_fit$time[keep],
  cox_snell_hazard[keep],
  type = "s",
  xlab = "Cox--Snell residual",
  ylab = "Estimated cumulative hazard",
  main = "Cox--Snell residual diagnostic"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

dev.off()


# Efron ties are used only for this DFBETA diagnostic; the reported final model remains exact.

# DFBETA influence diagnostic

influence_model <- coxph(
  formula(final_model),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

dfbeta_residuals <- residuals(
  influence_model,
  type = "dfbeta"
)

max_abs_dfbeta <- apply(
  abs(dfbeta_residuals),
  1,
  max
)

pdf(
  file.path(fig_dir, "chapter3_final_dfbeta.pdf"),
  width = 7,
  height = 5
)

plot(
  max_abs_dfbeta,
  type = "h",
  xlab = "Observation number",
  ylab = "Maximum absolute DFBETA",
  main = "Influence diagnostic"
)

dev.off()


list.files(
  output_root,
  pattern = "chapter3_final",
  recursive = TRUE,
  full.names = TRUE
)


## Influence and sensitivity analysis


dfbeta_matrix <- as.matrix(dfbeta_residuals)

colnames(dfbeta_matrix) <- names(coef(final_model))


max_abs_dfbeta <- apply(
  abs(dfbeta_matrix),
  1,
  max
)

most_influential <- which.max(max_abs_dfbeta)

affected_coefficient <- colnames(dfbeta_matrix)[
  which.max(abs(dfbeta_matrix[most_influential, ]))
]

cat(
  "Most influential observation:",
  most_influential,
  "\n"
)

cat(
  "Maximum absolute DFBETA:",
  max_abs_dfbeta[most_influential],
  "\n"
)

cat(
  "Coefficient most affected:",
  affected_coefficient,
  "\n"
)


influential_patient <- bc[
  most_influential,
  ,
  drop = FALSE
]

print(influential_patient)

write.csv(
  influential_patient,
  file.path(
    tab_dir,
    "chapter3_most_influential_patient.csv"
  ),
  row.names = TRUE
)


# The sensitivity refit continues to use the exact tie method.

bc_sensitivity <- bc[
  -most_influential,
  ,
  drop = FALSE
]

sensitivity_model <- coxph(
  formula(final_model),
  data = bc_sensitivity,
  ties = "exact",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

summary(sensitivity_model)


original_summary <- summary(final_model)
sensitivity_summary <- summary(sensitivity_model)

terms <- names(coef(final_model))

coefficient_comparison <- data.frame(
  term = terms,
  
  original_coefficient =
    unname(coef(final_model)[terms]),
  
  sensitivity_coefficient =
    unname(coef(sensitivity_model)[terms]),
  
  coefficient_change =
    unname(
      coef(sensitivity_model)[terms] -
        coef(final_model)[terms]
    ),
  
  original_HR =
    exp(unname(coef(final_model)[terms])),
  
  sensitivity_HR =
    exp(unname(coef(sensitivity_model)[terms])),
  
  HR_ratio =
    exp(
      unname(
        coef(sensitivity_model)[terms] -
          coef(final_model)[terms]
      )
    ),
  
  original_p =
    original_summary$coefficients[
      terms,
      "Pr(>|z|)"
    ],
  
  sensitivity_p =
    sensitivity_summary$coefficients[
      terms,
      "Pr(>|z|)"
    ],
  
  sign_changed =
    sign(coef(final_model)[terms]) !=
    sign(coef(sensitivity_model)[terms]),
  
  row.names = NULL
)

coefficient_comparison <- coefficient_comparison[
  order(
    abs(coefficient_comparison$coefficient_change),
    decreasing = TRUE
  ),
]

print(coefficient_comparison)

write.csv(
  coefficient_comparison,
  file.path(
    tab_dir,
    "chapter3_sensitivity_coefficient_comparison.csv"
  ),
  row.names = FALSE
)


ph_sensitivity <- cox.zph(
  sensitivity_model,
  terms = TRUE
)

print(ph_sensitivity)

write.table(
  ph_sensitivity$table,
  file = file.path(
    report_dir,
    "chapter3_sensitivity_ph_test.txt"
  ),
  quote = FALSE,
  sep = "\t"
)


head(
  coefficient_comparison,
  10
)


# Sensitivity model: remove menopausal status

final_model2 <- update(
  final_model,
  . ~ . - menop.new,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

summary(final_model2)

final_model2_term_tests <- drop1(
  final_model2,
  test = "Chisq"
)

print(final_model2_term_tests)

AIC(
  final_model,
  final_model2
)

anova(
  final_model2,
  final_model,
  test = "Chisq"
)

ph_final2 <- cox.zph(
  final_model2,
  terms = TRUE
)

print(ph_final2)

influence_model2 <- coxph(
  formula(final_model2),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

dfbeta2 <- as.matrix(
  residuals(
    influence_model2,
    type = "dfbeta"
  )
)

max_abs_dfbeta2 <- apply(
  abs(dfbeta2),
  1,
  max
)

cat(
  "Most influential observation:",
  which.max(max_abs_dfbeta2),
  "\n"
)

cat(
  "Maximum absolute DFBETA:",
  max(max_abs_dfbeta2),
  "\n"
)


# Sensitivity model: remove age group

final_model3 <- update(
  final_model,
  . ~ . - age.gp,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

summary(final_model3)

final_model3_term_tests <- drop1(
  final_model3,
  test = "Chisq"
)

print(final_model3_term_tests)

AIC(
  final_model,
  final_model2,
  final_model3
)

anova(
  final_model3,
  final_model,
  test = "Chisq"
)

ph_final3 <- cox.zph(
  final_model3,
  terms = TRUE
)

print(ph_final3)

influence_model3 <- coxph(
  formula(final_model3),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

dfbeta3 <- as.matrix(
  residuals(
    influence_model3,
    type = "dfbeta"
  )
)

max_abs_dfbeta3 <- apply(
  abs(dfbeta3),
  1,
  max
)

cat(
  "Most influential observation:",
  which.max(max_abs_dfbeta3),
  "\n"
)

cat(
  "Maximum absolute DFBETA:",
  max(max_abs_dfbeta3),
  "\n"
)


## Final model outputs

chapter3_final_model <- final_model

s <- summary(chapter3_final_model)

chapter3_final_results <- data.frame(
  term = rownames(s$coefficients),
  coefficient = s$coefficients[, "coef"],
  hazard_ratio = s$conf.int[, "exp(coef)"],
  lower_95_CI = s$conf.int[, "lower .95"],
  upper_95_CI = s$conf.int[, "upper .95"],
  p_value = s$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)

write.csv(
  chapter3_final_results,
  file.path(tab_dir, "chapter3_final_model_results.csv"),
  row.names = FALSE
)

chapter3_final_term_tests <- drop1(
  chapter3_final_model,
  test = "Chisq"
)

write.csv(
  chapter3_final_term_tests,
  file.path(tab_dir, "chapter3_final_model_term_tests.csv")
)

capture.output({
  cat("FINAL STRATIFIED COX MODEL\n\n")
  print(formula(chapter3_final_model))
  
  cat("\nMODEL SUMMARY\n\n")
  print(summary(chapter3_final_model))
  
  cat("\nWHOLE-VARIABLE TESTS\n\n")
  print(chapter3_final_term_tests)
  
  cat("\nPROPORTIONAL-HAZARDS TEST\n\n")
  print(cox.zph(chapter3_final_model, terms = TRUE))
  
  cat("\nSENSITIVITY MODEL COMPARISON\n\n")
  print(AIC(final_model, final_model2, final_model3))
  
  cat("\nREMOVE MENOPAUSAL STATUS\n\n")
  print(anova(final_model2, final_model, test = "Chisq"))
  
  cat("\nREMOVE AGE GROUP\n\n")
  print(anova(final_model3, final_model, test = "Chisq"))
  
  cat("\nINFLUENTIAL OBSERVATION\n\n")
  cat("Observation:", most_influential, "\n")
  cat("Maximum absolute DFBETA:",
      max_abs_dfbeta[most_influential], "\n")
  print(influential_patient)
  
}, file = file.path(report_dir, "chapter3_final_complete_output.txt"))


## Schoenfeld residual audit

ph_reduced_term <- cox.zph(
  cox_reduced,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

ph_reduced_coef <- cox.zph(
  cox_reduced,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

ph_final_term <- cox.zph(
  final_model,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

ph_final_coef <- cox.zph(
  final_model,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

cox_reduced_efron <- coxph(
  formula(cox_reduced),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

final_model_efron <- coxph(
  formula(final_model),
  data = bc,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

ph_reduced_term_efron <- cox.zph(
  cox_reduced_efron,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

ph_reduced_coef_efron <- cox.zph(
  cox_reduced_efron,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

ph_final_term_efron <- cox.zph(
  final_model_efron,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

capture.output({
  cat("SURVIVAL PACKAGE VERSION\n")
  print(packageVersion("survival"))
  
  cat("\nREDUCED MODEL: MAIN TERM-LEVEL TEST\n")
  print(ph_reduced_term)
  
  cat("\nREDUCED MODEL: SELECTED COEFFICIENT-LEVEL TESTS\n")
  print(
    ph_reduced_coef$table[
      c("presite.new2", "presite.new9", "nodes.axilla2", "GLOBAL"),
      ,
      drop = FALSE
    ]
  )
  
  cat("\nREDUCED MODEL: EFRON TERM-LEVEL SENSITIVITY TEST\n")
  print(ph_reduced_term_efron)
  
  cat("\nREDUCED MODEL: EFRON SELECTED COEFFICIENT-LEVEL TESTS\n")
  print(
    ph_reduced_coef_efron$table[
      c("presite.new2", "presite.new9", "nodes.axilla2", "GLOBAL"),
      ,
      drop = FALSE
    ]
  )
  
  cat("\nFINAL STRATIFIED MODEL: MAIN TERM-LEVEL TEST\n")
  print(ph_final_term)
  
  cat("\nFINAL STRATIFIED MODEL: ALL COEFFICIENT-LEVEL TESTS\n")
  print(ph_final_coef)
  
  cat("\nFINAL STRATIFIED MODEL: EFRON TERM-LEVEL SENSITIVITY TEST\n")
  print(ph_final_term_efron)
  
  cat("\nMAXIMUM ABSOLUTE COEFFICIENT DIFFERENCE: REDUCED EXACT VS EFRON\n")
  print(
    max(
      abs(
        coef(cox_reduced) -
          coef(cox_reduced_efron)[names(coef(cox_reduced))]
      )
    )
  )
  
  cat("\nMAXIMUM ABSOLUTE COEFFICIENT DIFFERENCE: FINAL EXACT VS EFRON\n")
  print(
    max(
      abs(
        coef(final_model) -
          coef(final_model_efron)[names(coef(final_model))]
      )
    )
  )
}, file = file.path(report_dir, "chapter3_schoenfeld_audit.txt"))


# Selected final-model Schoenfeld plots

ph_final_coef <- cox.zph(
  final_model,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

selected_ph_coefficients <- c(
  "menop.new9",
  "typeps.new8",
  "path.size2",
  "recept.new4"
)

pdf(
  file.path(fig_dir, "chapter3_final_schoenfeld_selected_coefficients.pdf"),
  width = 7,
  height = 5,
  onefile = TRUE
)

for (v in selected_ph_coefficients) {
  plot(
    ph_final_coef,
    var = v,
    resid = TRUE,
    se = TRUE
  )
  abline(
    h = unname(coef(final_model)[v]),
    lty = 3
  )
}

dev.off()
