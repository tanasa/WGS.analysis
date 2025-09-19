## ==========================================
## Robust per-gene survival (parallel, tolerant)
## Requires in env: df (cols: time, vital_status_bin, genes...)
## Reads genes from: A2780_targets_few_genes.txt  -> sets gene_list
## Outputs: per-gene plots/tables under survival_by_gene/ + summary CSV
## ==========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(future)
  library(future.apply)
})

options(warn = 1)  # show warnings immediately

## ---- Load gene list from file ----
# Check if gene list file exists
if (!file.exists("A2780_targets.txt")) {
  cat("Gene list file 'A2780_targets_few_genes.txt' not found. Creating a test list...\n")
  # Create a small test list with some common genes
  test_genes <- c("ENSG00000149948", "ENSG00000163605", "ENSG00000123456", "ENSG00000178901", "ENSG00000134567")
  writeLines(test_genes, "A2780_targets_few_genes.txt")
  genes <- test_genes
} else {
  # Read, trim, drop empties, deduplicate, strip version suffixes like ".12"
  genes <- readr::read_lines("A2780_targets.txt", progress = FALSE)
  genes <- unique(trimws(genes))
  genes <- genes[nzchar(genes)]
  genes <- sub("\\.[0-9]+$", "", genes)  # ENSG00000149948.12 -> ENSG00000149948
}

# Check which are present in your data frame 'df' (genes are columns)
stopifnot(exists("df"), is.data.frame(df))
genes_present <- intersect(genes, colnames(df))
genes_missing <- setdiff(genes, colnames(df))
cat("Present:", length(genes_present), " | Missing:", length(genes_missing), "\n")
if (length(genes_missing)) writeLines(genes_missing, "missing_genes.txt")

gene_list <- genes_present
print(paste("Number of genes:", length(gene_list)))

## ---- Basic checks on df ----
stopifnot(all(c("time", "vital_status_bin") %in% colnames(df)))
stopifnot(is.numeric(df$time))
stopifnot(all(na.omit(df$vital_status_bin) %in% c(0,1)))

## ---- Parallel plan ----
workers <- max(1, future::availableCores() - 1)
future::plan(multisession, workers = workers)
cat("Using workers:", workers, "\n")

## ---- One-gene analysis ----
analyze_gene_one <- function(df, gene_id, out_dir = "survival_by_gene",
                             min_group_n = 5, min_total_events = 10) {
  if (!gene_id %in% colnames(df)) return(NULL)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  base <- file.path(out_dir, gsub("[^[:alnum:]_.-]", "_", gene_id))

  # analysis subset
  d <- df %>%
    filter(!is.na(time), time > 0,
           !is.na(vital_status_bin),
           !is.na(.data[[gene_id]]))
  if (nrow(d) < (2 * min_group_n)) return(NULL)
  stopifnot(all(na.omit(d$vital_status_bin) %in% c(0,1)))

  # gene column numeric
  if (!is.numeric(d[[gene_id]])) {
    d[[gene_id]] <- suppressWarnings(as.numeric(as.character(d[[gene_id]])))
  }
  if (all(is.na(d[[gene_id]]))) return(NULL)

  # tertiles; guard collapsed cutpoints
  qs <- suppressWarnings(quantile(d[[gene_id]], c(0.33, 0.67), na.rm = TRUE, type = 1))
  if (!all(is.finite(qs)) || qs[1] >= qs[2]) return(NULL)

  d <- d %>%
    mutate(expr_group = dplyr::case_when(
      .data[[gene_id]] >= qs[2] ~ "High",
      .data[[gene_id]] <= qs[1] ~ "Low",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(expr_group)) %>%
    mutate(expr_group = factor(expr_group, levels = c("Low","High")))

  # size / events checks
  gsz <- table(d$expr_group)
  if (length(gsz) < 2 || any(gsz < min_group_n)) return(NULL)
  ev_by_group <- tapply(d$vital_status_bin, d$expr_group, sum, na.rm = TRUE)
  if (any(ev_by_group < 1) || sum(ev_by_group) < min_total_events) return(NULL)

  # ===== KM + log-rank (Surv in formula) =====
  fit <- survfit(Surv(time, vital_status_bin) ~ expr_group, data = d)
  lr  <- tryCatch(survdiff(Surv(time, vital_status_bin) ~ expr_group, data = d), error = function(e) NULL)
  if (is.null(lr) || is.na(lr$chisq)) return(NULL)
  p_logrank <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

  # p-value shown on plot (redundant with log-rank) for record
  p_plot <- tryCatch({
    as.numeric(survminer::surv_pvalue(fit, data = d)$pval)
  }, error = function(e) NA_real_)

  # ===== KM plot (forced save) =====
  km_plot <- ggsurvplot(
    fit, data = d, pval = TRUE, conf.int = TRUE,
    risk.table = TRUE, risk.table.height = 0.22,
    surv.median.line = "hv",
    censor.shape = 124, censor.size = 3,
    xlab = "Time (days)", ylab = "Overall survival probability",
    title = paste("Survival by expression of", gene_id),
    legend.title = "Expression", legend.labs = c("Low","High"),
    ggtheme = theme_minimal(base_size = 12)
  )

  km_expr_file <- paste0(base, "_KM_plot_expression.jpg")
  message("Saving: ", km_expr_file)
  ggsave(km_expr_file, plot = km_plot$plot, width = 8, height = 6, dpi = 300)

  # ===== Proper risk table as ggplot (avoids grob/gList) =====
  rt <- tryCatch(
    survminer::ggsurvplot_risktable(
      fit, data = d,
      ggtheme = theme_minimal(base_size = 12),
      legend.title = "Expression",
      legend.labs = c("Low","High")
    ),
    error = function(e) NULL
  )
  km_risk_file <- paste0(base, "_KM_plot_risktable.jpg")
  if (!is.null(rt) && inherits(rt, "ggplot")) {
    message("Saving: ", km_risk_file)
    ggsave(km_risk_file, plot = rt, width = 8, height = 3, dpi = 300)
  } else {
    warning(sprintf("Risk table save skipped for %s (could not build ggplot)", gene_id))
  }

  # ===== Cox PH (High vs Low) =====
  cox_fit <- tryCatch(coxph(Surv(time, vital_status_bin) ~ expr_group, data = d), error = function(e) NULL)
  HR <- lo <- hi <- p_cox <- c_index <- NA_real_
  ph_global <- NA_real_

  if (!is.null(cox_fit)) {
    cs <- summary(cox_fit)
    HR      <- unname(cs$coefficients[1, "exp(coef)"])
    p_cox   <- unname(cs$coefficients[1, "Pr(>|z|)"])
    lo      <- unname(cs$conf.int[1, "lower .95"])
    hi      <- unname(cs$conf.int[1, "upper .95"])
    c_index <- unname(cs$concordance[1])

    # Forest plot: try ggforest, else manual
    forest <- tryCatch({
      ggforest(cox_fit, data = d, main = paste("Hazard Ratio for", gene_id))
    }, error = function(e) {
      coef_names <- rownames(cs$coefficients)
      HRv  <- unname(cs$coefficients[, "exp(coef)"])
      pval <- unname(cs$coefficients[, "Pr(>|z|)"])
      lo95 <- unname(cs$conf.int[, "lower .95"])
      hi95 <- unname(cs$conf.int[, "upper .95"])
      stopifnot(length(coef_names) == length(HRv),
                length(HRv) == length(pval),
                length(lo95) == length(hi95))
      coef_data <- data.frame(
        Variable = coef_names,
        HR       = HRv,
        CI_lower = lo95,
        CI_upper = hi95,
        p_value  = pval,
        check.names = FALSE
      )
      ggplot(coef_data, aes(x = Variable, y = HR)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        coord_flip() +
        labs(title = paste("Hazard Ratio for", gene_id),
             x = "Variable", y = "Hazard Ratio") +
        theme_minimal()
    })

    tryCatch({
      ggsave(paste0(base, "_cox_forest.png"), plot = forest, width = 7, height = 5, dpi = 300)
    }, error = function(e) warning(sprintf("Forest save failed for %s: %s", gene_id, conditionMessage(e))))

    # PH diagnostics
    ph <- tryCatch(cox.zph(cox_fit), error = function(e) NULL)
    if (!is.null(ph)) {
      ph_global <- suppressWarnings({
        pr <- tryCatch(as.numeric(ph$table["GLOBAL","p"]), error = function(e) NA_real_)
        pr
      })

      # png(paste0(base, "_schoenfeld_residuals.png"), width = 8, height = 6, units = "in", res = 300)
      # op <- par(no.readonly = TRUE); par(oma = c(0,0,2.2,0))
      # plot(ph); mtext("Schoenfeld Residuals Test", outer = TRUE, line = 0.8, cex = 1.1)
      # par(op); dev.off()

      # png(paste0(base, "_schoenfeld_residuals_scaled.png"), width = 8, height = 6, units = "in", res = 300)
      # op <- par(no.readonly = TRUE); par(oma = c(0,0,2.2,0))
      # plot(cox.zph(cox_fit, transform = "identity"))
      # mtext("Scaled Schoenfeld Residuals (identity)", outer = TRUE, line = 0.8, cex = 1.1)
      # par(op); dev.off()
    }

    # Influential obs (dfbeta)
    # dfbeta_vals <- tryCatch(residuals(cox_fit, type = "dfbeta"), error = function(e) NULL)
    # if (!is.null(dfbeta_vals)) {
    #   png(paste0(base, "_dfbeta_influential.png"), width = 8, height = 6, units = "in", res = 300)
    #   matplot(dfbeta_vals, type = "l", lty = 1,
    #           main = paste("Influential Observations (dfbeta) -", gene_id),
    #           xlab = "Observation index", ylab = "dfbeta values")
    #   abline(h = 0, lty = 2)
    #   dev.off()
    # }
  }

  # ---- Per-gene p-values (append) ----
  pvals_df <- tibble::tibble(
    gene = gene_id,
    p_logrank = p_logrank,
    p_cox_wald = p_cox,
    p_plot_logrank = p_plot,
    p_ph_global = ph_global
  )
  out_file <- paste0(base, "_pvalues.csv")
  append_mode <- file.exists(out_file) && file.size(out_file) > 0
  write.table(
    pvals_df, file = out_file, sep = ",", quote = FALSE,
    row.names = FALSE,
    col.names = !append_mode,
    append = append_mode
  )

  # ---- Return summary row ----
  tibble::tibble(
    gene = gene_id,
    n_low  = as.integer(gsz["Low"]),
    n_high = as.integer(gsz["High"]),
    events_low  = unname(ev_by_group["Low"]),
    events_high = unname(ev_by_group["High"]),
    p_logrank = p_logrank,
    p_cox_wald = p_cox,
    p_plot_logrank = p_plot,
    p_ph_global = ph_global,
    cox_HR = HR, cox_CI_low = lo, cox_CI_high = hi,
    c_index = c_index
  )
}

## ---- Use the file-derived gene_list ----
present <- gene_list
missing <- character(0)  # already computed above
cat("Genes present:", length(present), " | missing:", length(missing), "\n")
if (length(missing)) writeLines(missing, "genes_missing_from_df.txt")

## ---- Minimal df for workers ----
keep_cols <- unique(c("time", "vital_status_bin", present))
df_min_all <- df[, keep_cols, drop = FALSE]

## ---- Error log ----
if (!dir.exists("survival_by_gene")) dir.create("survival_by_gene", recursive = TRUE)
error_log_file <- file.path("survival_by_gene", "errors.log")
if (file.exists(error_log_file)) file.remove(error_log_file)

safe_analyze <- function(df_min, g) {
  tryCatch({
    res <- analyze_gene_one(df_min, g, out_dir = "survival_by_gene")
    list(ok = TRUE, gene = g, res = res)
  }, error = function(e) {
    msg <- sprintf("[%s] %s", g, conditionMessage(e))
    cat(msg, "\n", file = error_log_file, append = TRUE)
    warning(msg)
    list(ok = FALSE, gene = g, res = NULL, err = conditionMessage(e))
  })
}

## ---- Parallel run over genes (error-tolerant) ----
res_list <- future_lapply(
  present,
  function(g) {
    # load libs inside worker (some clusters require this)
    library(dplyr); library(survival); library(survminer); library(ggplot2); library(ggpubr)
    df_min <- df_min_all[, c("time","vital_status_bin", g), drop = FALSE]
    df_min <- df_min[complete.cases(df_min), ]
    if (nrow(df_min) < 10) {
      msg <- sprintf("[%s] Skipped (rows < 10 after filtering)", g)
      cat(msg, "\n", file = error_log_file, append = TRUE)
      return(list(ok = FALSE, gene = g, res = NULL, err = "rows<10"))
    }
    safe_analyze(df_min, g)
  },
  future.seed = TRUE
)

## ---- Collect & save summary (with BH FDR) ----
ok_results <- Filter(function(x) is.list(x) && isTRUE(x$ok) && !is.null(x$res), res_list)
if (length(ok_results)) {
  summary_tbl <- dplyr::bind_rows(lapply(ok_results, `[[`, "res"))
  # BH-adjust p-values (per family across processed genes)
  summary_tbl <- summary_tbl %>%
    mutate(
      p_logrank_adj_BH = p.adjust(p_logrank, method = "BH"),
      p_cox_wald_adj_BH = p.adjust(p_cox_wald, method = "BH")
    ) %>%
    arrange(p_logrank_adj_BH, p_logrank)

  out_csv <- file.path("survival_by_gene", "summary_results_with_BH.csv")
  readr::write_csv(summary_tbl, out_csv)
  print(head(summary_tbl, 10))
  message("Summary written: ", out_csv)
} else {
  message("No genes passed the data-quality checks (or all failed). See survival_by_gene/errors.log")
}
