### Chapter 2: Exploratory analysis of baseline covariates

rm(list = ls())


## Packages and data setup

packages <- c("survival", "survminer", "ggplot2", "gridExtra")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)


load("/Users/geruxi/Documents/硕士所有东西/硕士毕业论文/老师给的/OneDrive_1_2026-6-26/The breast cancer survival data/Data to use for project/bc survival data - to be used.RData")

stopifnot(exists("bcdata.use"))
stopifnot(all(c("stime", "cens") %in% names(bcdata.use)))


## Output folders

fig_dir <- "figures"
tab_dir <- "tables"

dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)


## Baseline covariate recoding

bcdata.use$age.gp.f <- factor(
  bcdata.use$age.gp,
  levels = c(1, 2, 3),
  labels = c("Age < 50 years", "50 <= age < 65 years", "Age >= 65 years")
)

bcdata.use$menop.new.f <- factor(
  bcdata.use$menop.new,
  levels = c(1, 2, 9),
  labels = c(
    "Last menstrual period within previous 3 years",
    "Last menstrual period more than 3 years previously",
    "Unknown"
  )
)

bcdata.use$presite.new.f <- factor(
  bcdata.use$presite.new,
  levels = c(1, 2, 9),
  labels = c("Non-central site", "Central or subareolar site", "Unknown")
)

bcdata.use$side.f <- factor(
  bcdata.use$side,
  levels = c(1, 2, 9),
  labels = c("Right breast", "Left breast", "Unknown")
)

bcdata.use$max.cdg.f <- factor(
  bcdata.use$max.cdg,
  levels = c(1, 2, 3, 9),
  labels = c(
    "Clinical tumour diameter < 2 cm",
    "Clinical tumour diameter 2 to < 5 cm",
    "Clinical tumour diameter >= 5 cm",
    "Unknown"
  )
)

bcdata.use$nodes.axilla.f <- factor(
  bcdata.use$nodes.axilla,
  levels = c(1, 2, 9),
  labels = c("Negative axillary nodes", "Positive axillary nodes", "Unknown")
)

bcdata.use$mcr.stage.f <- factor(
  bcdata.use$mcr.stage,
  levels = c(1, 2, 3, 4, 9),
  labels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Unknown")
)

bcdata.use$typeps.new.f <- factor(
  bcdata.use$typeps.new,
  levels = c(3, 4, 6, 7, 8),
  labels = c(
    "Excision biopsy or surgery unknown",
    "Simple mastectomy",
    "Wide local excision with axillary nodal clearance",
    "Mastectomy with axillary nodal clearance",
    "Quadrantectomy with axillary nodal clearance"
  )
)

bcdata.use$adj.radio.f <- factor(
  bcdata.use$adj.radio,
  levels = c(1, 2, 9),
  labels = c("No", "Yes", "Unknown")
)

bcdata.use$histo.new.f <- factor(
  bcdata.use$histo.new,
  levels = c(1, 2, 3, 99),
  labels = c(
    "Infiltrating ductal carcinoma",
    "Infiltrating lobular carcinoma",
    "Other histological type",
    "Unknown"
  )
)

bcdata.use$nr.cat.f <- factor(
  bcdata.use$nr.cat,
  levels = c(1, 2, 3, 9),
  labels = c(
    "Node ratio = 0",
    "0 < node ratio <= 0.4",
    "0.4 < node ratio <= 1",
    "Unknown"
  )
)

bcdata.use$path.size.f <- factor(
  bcdata.use$path.size,
  levels = c(1, 2, 3, 9),
  labels = c(
    "Pathological tumour size < 2 cm",
    "Pathological tumour size 2 to < 5 cm",
    "Pathological tumour size >= 5 cm",
    "Unknown"
  )
)

bcdata.use$tgrade.new.f <- factor(
  bcdata.use$tgrade.new,
  levels = c(2, 3, 9),
  labels = c("Grades 1 and 2 combined", "Grade 3", "Unknown")
)

bcdata.use$recept.new.f <- factor(
  bcdata.use$recept.new,
  levels = c(1, 2, 3, 4, 9),
  labels = c(
    "ER negative / PR negative",
    "ER negative / PR positive",
    "ER positive / PR negative",
    "ER positive / PR positive",
    "Both receptor values unknown"
  )
)


## Descriptive cross-tabulations

make_crosstab <- function(data, var, var_label) {

  x <- droplevels(data[[var]])

  tab <- table(
    x,
    factor(data$cens, levels = c(0, 1))
  )

  out <- data.frame(
    Variable = var_label,
    Category = rownames(tab),
    Censored = as.integer(tab[, "0"]),
    Event_observed = as.integer(tab[, "1"]),
    Total = rowSums(tab),
    Event_percent = round(100 * as.integer(tab[, "1"]) / rowSums(tab), 1),
    row.names = NULL
  )

  out
}


# Table 2.2: all 14 baseline covariates

ct_age <- make_crosstab(bcdata.use, "age.gp.f", "Age group")
ct_menop <- make_crosstab(bcdata.use, "menop.new.f", "Menopausal status")
ct_presite <- make_crosstab(bcdata.use, "presite.new.f", "Predominant tumour site")
ct_side <- make_crosstab(bcdata.use, "side.f", "Side of affected breast")
ct_maxcdg <- make_crosstab(bcdata.use, "max.cdg.f", "Maximum clinical tumour diameter")
ct_nodes <- make_crosstab(bcdata.use, "nodes.axilla.f", "Axillary node status")
ct_stage <- make_crosstab(bcdata.use, "mcr.stage.f", "Manchester stage")
ct_surgery <- make_crosstab(bcdata.use, "typeps.new.f", "Type of primary surgery")
ct_radio <- make_crosstab(bcdata.use, "adj.radio.f", "Adjuvant radiotherapy")
ct_histo <- make_crosstab(bcdata.use, "histo.new.f", "Histological type")
ct_ratio <- make_crosstab(bcdata.use, "nr.cat.f", "Node-ratio category")
ct_pathsize <- make_crosstab(bcdata.use, "path.size.f", "Pathological tumour size")
ct_grade <- make_crosstab(bcdata.use, "tgrade.new.f", "Tumour grade")
ct_receptor <- make_crosstab(bcdata.use, "recept.new.f", "Receptor status")

ct_all <- rbind(
  ct_age,
  ct_menop,
  ct_presite,
  ct_side,
  ct_maxcdg,
  ct_nodes,
  ct_stage,
  ct_surgery,
  ct_radio,
  ct_histo,
  ct_ratio,
  ct_pathsize,
  ct_grade,
  ct_receptor
)

write.csv(
  ct_all,
  file = file.path(tab_dir, "all14_covariate_crosstabs.csv"),
  row.names = FALSE
)


## Log-rank tests

make_logrank <- function(data, var, var_label) {
  f <- as.formula(paste("Surv(stime, cens) ~", var))
  lr <- survdiff(f, data = data)

  chisq <- unname(lr$chisq)
  df <- length(lr$n) - 1
  pval <- pchisq(chisq, df = df, lower.tail = FALSE)

  data.frame(
    Variable = var_label,
    Chi_square = round(chisq, 1),
    df = df,
    p_value = pval,
    p_value_print = ifelse(
      pval < 2e-16,
      "<2e-16",
      formatC(pval, format = "e", digits = 2)
    ),
    row.names = NULL
  )
}

lr_stage <- make_logrank(
  bcdata.use,
  "mcr.stage.f",
  "Manchester stage"
)

lr_nodes <- make_logrank(
  bcdata.use,
  "nodes.axilla.f",
  "Axillary node status"
)

lr_ratio <- make_logrank(
  bcdata.use,
  "nr.cat.f",
  "Node ratio category"
)

lr_all <- rbind(lr_stage, lr_nodes, lr_ratio)

write.csv(
  lr_all,
  file = file.path(tab_dir, "logrank_tests.csv"),
  row.names = FALSE
)


## LaTeX and table exports

latex_escape <- function(x) {
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

write_crosstab_latex <- function(df, file_name) {

  con <- file(file_name, open = "w")

  writeLines("\\begingroup", con)
  writeLines("\\small", con)
  writeLines("\\setlength{\\tabcolsep}{3pt}", con)
  writeLines("\\renewcommand{\\arraystretch}{1.08}", con)
  writeLines("\\setlength{\\LTleft}{0pt}", con)
  writeLines("\\setlength{\\LTright}{0pt}", con)

  writeLines(
    "\\begin{longtable}{@{}L{0.24\\textwidth}L{0.34\\textwidth}C{0.09\\textwidth}C{0.11\\textwidth}C{0.07\\textwidth}C{0.08\\textwidth}@{}}",
    con
  )

  writeLines(
    "\\caption{Cross-tabulations of all baseline covariates by censoring status}",
    con
  )
  writeLines("\\label{tab:selected-covariate-crosstabs}\\\\", con)

  writeLines("\\toprule", con)
  writeLines("Variable & Category & Censored & Event observed & Total & Event (\\%) \\\\", con)
  writeLines("\\midrule", con)
  writeLines("\\endfirsthead", con)

  writeLines(
    "\\caption[]{Cross-tabulations of all baseline covariates by censoring status (continued)}\\\\",
    con
  )
  writeLines("\\toprule", con)
  writeLines("Variable & Category & Censored & Event observed & Total & Event (\\%) \\\\", con)
  writeLines("\\midrule", con)
  writeLines("\\endhead", con)

  writeLines("\\midrule", con)
  writeLines("\\multicolumn{6}{r}{Continued on next page}\\\\", con)
  writeLines("\\endfoot", con)

  writeLines("\\bottomrule", con)
  writeLines("\\endlastfoot", con)

  previous_variable <- ""

  for (i in seq_len(nrow(df))) {

    current_variable <- df$Variable[i]

    variable_print <- ifelse(
      current_variable == previous_variable,
      "",
      latex_escape(current_variable)
    )

    line <- paste0(
      variable_print, " & ",
      latex_escape(df$Category[i]), " & ",
      df$Censored[i], " & ",
      df$Event_observed[i], " & ",
      df$Total[i], " & ",
      sprintf("%.1f", df$Event_percent[i]), " \\\\"
    )

    writeLines(line, con)
    previous_variable <- current_variable
  }

  writeLines("\\end{longtable}", con)
  writeLines("\\endgroup", con)

  close(con)
}

write_crosstab_latex(
  ct_all,
  file.path(tab_dir, "all14_covariate_crosstabs.tex")
)


# Table 2.2 PDF and PNG export

make_table2_2_grob <- function(df) {

  display_df <- df

  display_df$Variable[duplicated(display_df$Variable)] <- ""

  names(display_df) <- c(
    "Variable",
    "Category",
    "Censored",
    "Event observed",
    "Total",
    "Event (%)"
  )

  table_theme <- ttheme_minimal(
    base_size = 8,
    core = list(
      fg_params = list(hjust = 0, x = 0.02)
    ),
    colhead = list(
      fg_params = list(fontface = "bold")
    )
  )

  tableGrob(
    display_df,
    rows = NULL,
    theme = table_theme
  )
}

table2_2_grob <- make_table2_2_grob(ct_all)

pdf(
  file = file.path(tab_dir, "table2_2_all14_covariates.pdf"),
  width = 11.7,
  height = 16.5,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(table2_2_grob)
dev.off()

png(
  filename = file.path(tab_dir, "table2_2_all14_covariates.png"),
  width = 2400,
  height = 3400,
  res = 220
)
grid::grid.newpage()
grid::grid.draw(table2_2_grob)
dev.off()


# Log-rank table export

write_logrank_latex <- function(df, file_name) {
  con <- file(file_name, open = "w")

  writeLines("\\begin{table}[htbp]", con)
  writeLines("\\centering", con)
  writeLines("\\caption{Log-rank tests for selected grouped survival comparisons}", con)
  writeLines("\\label{tab:logrank-tests}", con)
  writeLines("\\begin{tabular}{lrrr}", con)
  writeLines("\\hline", con)
  writeLines("Variable & $\\chi^2$ & df & p-value \\\\", con)
  writeLines("\\hline", con)

  for (i in seq_len(nrow(df))) {
    p_latex <- ifelse(
      df$p_value_print[i] == "<2e-16",
      "$<2\\times 10^{-16}$",
      paste0("$", df$p_value_print[i], "$")
    )

    line <- paste0(
      df$Variable[i], " & ",
      df$Chi_square[i], " & ",
      df$df[i], " & ",
      p_latex, " \\\\"
    )

    writeLines(line, con)
  }

  writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\end{table}", con)

  close(con)
}

write_logrank_latex(
  lr_all,
  file.path(tab_dir, "logrank_tests.tex")
)


## Kaplan-Meier figures

make_km_plot <- function(data, var, file_name, risk_table_height = 0.25) {

  plot_data <- as.data.frame(data)

  stopifnot(all(c("stime", "cens", var) %in% names(plot_data)))

  plot_data$group <- plot_data[[var]]
  plot_data$group <- factor(plot_data$group)

  group_labels <- levels(plot_data$group)
  n_groups <- length(group_labels)

  fit <- survfit(Surv(stime, cens) ~ group, data = plot_data)

  km_plot <- ggsurvplot(
    fit = fit,
    data = plot_data,

    conf.int = FALSE,
    censor = FALSE,

    risk.table = TRUE,
    risk.table.height = risk_table_height,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = FALSE,

    break.time.by = 500,
    xlim = c(0, 2600),
    ylim = c(0, 1),

    xlab = "Time in days",
    ylab = "Survival probability",

    legend.title = "",
    legend.labs = group_labels,
    legend = "bottom",

    palette = rep("black", n_groups),
    linetype = "strata",
    size = 0.8,

    ggtheme = theme_classic(base_size = 12),
    tables.theme = theme_cleantable(base_size = 10)
  )

  km_plot$plot <- km_plot$plot +
    theme(
      plot.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(1.2, "cm"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  pdf(
    file = file.path(fig_dir, file_name),
    width = 6.8,
    height = 5.8,
    useDingbats = FALSE
  )

  print(km_plot)
  dev.off()
}


make_km_plot(
  data = bcdata.use,
  var = "mcr.stage.f",
  file_name = "km_manchester_stage_final.pdf",
  risk_table_height = 0.25
)

make_km_plot(
  data = bcdata.use,
  var = "nodes.axilla.f",
  file_name = "km_axillary_nodes_final.pdf",
  risk_table_height = 0.23
)

make_km_plot(
  data = bcdata.use,
  var = "nr.cat.f",
  file_name = "km_node_ratio_final.pdf",
  risk_table_height = 0.30
)


## Console outputs

cat("\nTable 2.2: cross-tabulations of all 14 baseline covariates:\n")
print(ct_all)

cat("\nLog-rank tests:\n")
print(lr_all)

cat("\nFiles generated:\n")
cat("- tables/all14_covariate_crosstabs.csv\n")
cat("- tables/all14_covariate_crosstabs.tex\n")
cat("- tables/table2_2_all14_covariates.pdf\n")
cat("- tables/table2_2_all14_covariates.png\n")
cat("- tables/logrank_tests.csv\n")
cat("- tables/logrank_tests.tex\n")
cat("- figures/km_manchester_stage_final.pdf\n")
cat("- figures/km_axillary_nodes_final.pdf\n")
cat("- figures/km_node_ratio_final.pdf\n")
