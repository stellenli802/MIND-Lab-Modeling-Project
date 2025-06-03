
#' predictSurvProb.survreg
#'
#' allows the function pec::predictSurvProb to work on survreg objects
#'
#' @param object 
#' @param newdata 
#' @param times 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
predictSurvProb.survreg <- function(object, newdata, times, ...) {
  stopifnot(requireNamespace("pec", quietly = TRUE))
  lp <- predict(object, newdata = newdata, type = "lp")
  scale <- object$scale
  dist <- object$dist
  
  if (dist != "loglogistic") {
    stop("This function currently supports only 'loglogistic' distributions.")
  }
  
  surv_probs <- outer(times, lp, function(t, mu) {
    1 / (1 + (t / exp(mu)) ^ (1 / scale))
  })
  
  return(t(surv_probs))  # rows = newdata, cols = times
}

#' Compute Harrell's C-index for a Survival Model
#'
#' @param fit A fitted Cox model (from `coxph()`).
#' @param data A data frame with the original survival data.
#' @param time_col Column name for time-to-event.
#' @param status_col Column name for event indicator (1 = event, 0 = censored).
#'
#' @return A numeric value for Harrell’s concordance index.
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' getHarrellC(fit, data = lung, time_col = "time", status_col = "status")
getHarrellC <- function(fit, data, time_col = "W", status_col = "delta") {
  stopifnot(requireNamespace("survival", quietly = TRUE))
  
  if (inherits(fit, "survreg")) {
    # Negatve linear predictor for survreg models: higher = higher risk
    lp <- -predict(fit, newdata = data, type = "lp")  
  } else {
    lp <- predict(fit, newdata = data, type = "lp")
  }
  
  concord <- survival::concordance(Surv(data[[time_col]], data[[status_col]]) ~ lp, data = data)
  
  c_stat <- concord$concordance
  n_pairs <- sum(unlist(concord$count[c("concordant", "discordant", "tied.risk")]), na.rm = TRUE)
  
  message("C = ", round(c_stat, 3), "; usable pairs = ", n_pairs)
  
  return(c_stat)
}

#' Compute Uno's C-index for a Survival Model
#'
#' @param fit A fitted Cox model (from `coxph()`).
#' @param data A data frame with the original survival data.
#' @param time_col Column name for time-to-event.
#' @param status_col Column name for event indicator (1 = event, 0 = censored).
#'-
#' @return A numeric value for Uno’s concordance index.
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' getUnoC(fit, data = lung, time_col = "time", status_col = "status")
getUnoC_multiple_times <- function(fit, data, time_col, status_col, model_name, t_star_vec) {
  stopifnot(requireNamespace("survAUC", quietly = TRUE))
  
  # Adjust for AFT models (e.g., survreg)
  if (inherits(fit, "survreg")) {
    lp <- -predict(fit, newdata = data, type = "lp")
  } else {
    lp <- predict(fit, newdata = data, type = "lp")
  }
  
  surv_obj <- Surv(data[[time_col]], data[[status_col]])
  
  results <- lapply(t_star_vec, function(t_star) {
    tryCatch({
      result <- survAUC::AUC.uno(
        Surv.rsp = surv_obj,
        Surv.rsp.new = surv_obj,
        lpnew = lp,
        times = t_star
      )
      data.frame(
        model = model_name,
        t_star = t_star,
        UnoC = result$iauc
      )
    }, error = function(e) {
      warning(sprintf("Failed at t = %s for model %s: %s", t_star, model_name, conditionMessage(e)))
      data.frame(model = model_name, t_star = t_star, UnoC = NA)
    })
  })
  
  do.call(rbind, results)
}


#' Get ROC Data and AUC for Survival Models at Multiple Timepoints
#'
#' Computes time-dependent ROC curves and AUC values using the `survivalROC` package.
#'
#' @param fit A fitted Cox or survival model.
#' @param times A vector of timepoints at which to compute ROC.
#' @param data A data frame containing the event/censoring info and covariates.
#' @param time_col Name of the time-to-event variable.
#' @param status_col Name of the event indicator (1 = event, 0 = censored).
#' @param model_name Optional name to label the model in outputs.
#'
#' @return A data frame containing FPR, TPR, AUC, time, and model info.
#' @export
getSurvivalROCdata <- function(fit, times, data, time_col, status_col, model_name = "Model") {
  stopifnot(requireNamespace("survivalROC", quietly = TRUE))
  
  # Adjust for AFT models (e.g., survreg)
  if (inherits(fit, "survreg")) {
    lp <- -predict(fit, newdata = data, type = "lp")
  } else {
    lp <- predict(fit, newdata = data, type = "lp")
  }
  
  all_data <- list()
  
  for (t in times) {
    roc_obj <- survivalROC::survivalROC(
      Stime = data[[time_col]],
      status = data[[status_col]],
      marker = lp,
      predict.time = t,
      method = "NNE",
      span = 0.25 * nrow(data)^(-0.20)
    )
    
    df <- data.frame(
      FPR = roc_obj$FP,
      TPR = roc_obj$TP,
      AUC = roc_obj$AUC,
      time = t,
      model = model_name
    )
    all_data[[as.character(t)]] <- df
  }
  
  dplyr::bind_rows(all_data)
}

#' Plot Survival ROC Curves with AUC
#'
#' @param roc_df Data frame returned by `getSurvivalROCdata()`.
#' @param color_by Column used for color grouping (default: "time").
#' @param facet_by Column used for faceting (default: "model").
#' @param show_auc Logical; if TRUE, include AUC in the legend.
#'
#' @return A ggplot object.
#' @export
plotSurvivalROCgg <- function(roc_df, color_by = "model", show_auc = TRUE) {
  library(ggplot2)
  library(dplyr)
  
  if (show_auc) {
    roc_df <- roc_df %>%
      group_by(model, time) %>%
      mutate(label_auc = paste0(.data[[color_by]], " (AUC=", round(first(AUC), 3), ")")) %>%
      ungroup()
    
    ggplot(roc_df, aes(x = FPR, y = TPR, color = label_auc, group = interaction(model, time))) +
      geom_line(size = 1.1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(model), cols = vars(time)) +
      labs(
        title = "Time-dependent ROC Curves",
        x = "False Positive Rate",
        y = "True Positive Rate",
        color = tools::toTitleCase(color_by)
      ) +
      theme_minimal() +
      theme(strip.text = element_text(size = 12))
  } else {
    ggplot(roc_df, aes(x = FPR, y = TPR, color = .data[[color_by]], group = interaction(model, time))) +
      geom_line(size = 1.1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(model), cols = vars(time)) +
      labs(
        title = "Time-dependent ROC Curves",
        x = "False Positive Rate",
        y = "True Positive Rate",
        color = tools::toTitleCase(color_by)
      ) +
      theme_minimal() +
      theme(strip.text = element_text(size = 12))
  }
}

#' Collect Calibration Curve Data for ggplot Visualization
#'
#' Runs `calPlot()` across multiple timepoints and exports observed vs. predicted calibration
#' data, including quantile or NNE method-specific labels for faceting or coloring in ggplot.
#'
#' @param fit A fitted model object (used by `pec::calPlot()`).
#' @param times A numeric vector of timepoints to plot.
#' @param data A data frame used for the calibration evaluation.
#' @param method Either `"quantile"` or `"nne"` for calibration method.
#' @param q Quantile value (used only if `method == "quantile"`).
#' @param model_name A string used to label the model.
#'
#' @return A data frame containing observed, predicted, time, method, and model labels.
#' @export
getCalibrationData <- function(fit, times, data,
                               method,
                               q = NULL,
                               model_name = "Model") {
  all_data <- list()
  
  # if (inherits(fit, "survreg")) {
  #   fit <- pec::as.predictionModel(fit, data = data, type = "survreg")
  # }
  
  for (t in times) {
    cal_obj <- if (method == "quantile") {
      pec::calPlot(object = fit, time = t, data = data, method = "quantile", q = q)
    } else {
      pec::calPlot(object = fit, time = t, data = data, method = "nne")
    }
    
    cal_df <- extract_calPlot_data(cal_obj, model_id = "Model.1", model_name = model_name, time = t)
    cal_df$method_value <- if (method == "quantile") q else "NNE"
    all_data[[as.character(t)]] <- cal_df
  }
  
  dplyr::bind_rows(all_data)
}

#' extract_calPlot_data
#'
#' @param cal_obj 
#' @param model_id 
#' @param model_name 
#' @param time 
#'
#' @return
#' @export
#'
#' @examples
extract_calPlot_data <- function(cal_obj,
                                 model_id = "Model.1",
                                 model_name = "Model",
                                 time = NULL) {
  stopifnot(inherits(cal_obj, "calibrationPlot"))
  
  if (!model_id %in% names(cal_obj$plotFrames)) {
    stop("Model name '", model_id, "' not found in cal_obj$plotFrames.")
  }
  
  df <- cal_obj$plotFrames[[model_id]]
  qvals <- attr(df, "quantiles")
  if (!is.null(qvals) && length(qvals) == nrow(df) + 1) {
    qnames <- names(qvals)
    quantile_labels <- paste(head(qnames, -1), tail(qnames, -1), sep = "–")
    df$quantile_group <- quantile_labels
  } else {
    df$quantile_group <- paste0("Q", seq_len(nrow(df)))
  }
  df$model <- model_name
  df$time <- if (!is.null(time)) time else cal_obj$time
  
  return(df)
}

#' Plot Calibration Curves Using ggplot2
#'
#' @param cal_df Data frame returned by `getCalibrationData()`.
#' @param color_by Column to color lines by (e.g., `"time"` or `"model"`).
#' @param facet_by Optional column to facet by (e.g., `"model"` or `"method_value"`).
#' @param reference_line Logical; add 45-degree reference line (default TRUE).
#'
#' @return A `ggplot` object.
#' @export
plotCalibrationGG <- function(cal_df,
                              color_by = "time",
                              facet_by = NULL,
                              reference_line = TRUE) {
  library(ggplot2)
  
  p <- ggplot(cal_df, aes(x = Pred, y = Obs, color = as.factor(.data[[color_by]]))) +
    geom_line(size = 1.1) +
    labs(x = "Predicted Probability", 
         y = "Observed Probability", 
         title = "Calibration Plot",
         color = "Time") +
    theme_minimal()
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  if (reference_line) {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  }
  
  return(p)
}

#' Compute Brier Score and Integrated Brier Score (IBS) for a Survival Model
#'
#' Uses the `pec` package to evaluate the Brier score at specified timepoints, accounting for censoring.
#'
#' @param fit A fitted Cox or survival model (e.g., from `coxph()` or `rfsrc()`).
#' @param data The data frame used for evaluation.
#' @param formula A survival formula, e.g., `Surv(time, status) ~ .`
#' @param times A numeric vector of timepoints at which to compute the Brier score.
#'
#' @return An object of class `pec` containing Brier scores and IBS.
#' @export
#'
#' @examples
#' library(survival)
#' cox_fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' brier_obj <- getBrierScore(cox_fit, data = lung, formula = Surv(time, status) ~ age + sex, times = c(90, 180, 365))
#' print(brier_obj)
getBrierScore <- function(fit, data, formula, times) {
  stopifnot(requireNamespace("pec", quietly = TRUE))
  
  pec::pec(object = fit,
           formula = formula,
           data = data,
           times = times,
           exact = FALSE,
           cens.model = "cox",
           splitMethod = "none")
}

#' Extract Brier Scores and IBS from pec Object
#'
#' @param pec_obj An object returned by `pec::pec()`.
#'
#' @return A list with `brier_scores` (data.frame) and `ibs` (numeric).
#' @export
extractBrierSummary <- function(pec_obj) {
  scores <- as.data.frame(pec_obj$AppErr)
  times <- pec_obj$time
  model_names <- colnames(scores)
  
  brier_df <- tibble::tibble(
    time = rep(times, each = length(model_names)),
    model = rep(model_names, times = length(times)),
    brier = as.vector(as.matrix(scores))
  )
  
  ibs_vals <- pec::crps(pec_obj)
  
  list(
    brier_scores = brier_df,
    ibs = ibs_vals
  )
}

#' Plot Brier Scores Over Time with IBS Labels
#'
#' Takes a named list of model Brier score summaries (output from `extractBrierSummary()`)
#' and returns a ggplot object with one line per model, labeled with its Integrated Brier Score (IBS).
#'
#' @param summaries A named list of `extractBrierSummary()` outputs (one per model).
#'
#' @return A ggplot object showing Brier score curves and IBS-labeled legend.
#' @export
#'
#' @examples
#' summaries <- list(
#'   MRS = extractBrierSummary(getBrierScore(...)),
#'   CAP = extractBrierSummary(getBrierScore(...))
#' )
#' plotBrierSummaries(summaries)
plotBrierSummaries <- function(summaries) {
  # Step 1: Extract IBS values and build labels
  ibs_df <- purrr::map2_dfr(summaries, names(summaries), function(model_list, model_name) {
    ibs_mat <- model_list$ibs
    model_labels <- rownames(ibs_mat)
    
    method_row <- which(model_labels != "Reference")[1]
    method_name <- model_labels[method_row]
    ibs_value <- ibs_mat[method_row, 1]
    
    tibble(
      model_name = model_name,
      method = method_name,
      ibs = ibs_value,
      label = paste0(model_name, " (IBS = ", round(ibs_value, 3), ")")
    )
  })
  
  # Step 2: Combine all brier_scores and join IBS labels
  all_scores <- purrr::map2_dfr(summaries, names(summaries), function(model_list, model_name) {
    bs <- model_list$brier_scores
    bs$model_name <- model_name
    bs
  }) %>%
    dplyr::filter(model != "Reference") %>%
    dplyr::left_join(ibs_df, by = "model_name")
  
  # Step 3: Plot with IBS in the legend
  ggplot(all_scores, aes(x = time, y = brier, color = label, linetype = model)) +
    geom_line(size = 1.2) +
    labs(
      title = "Brier Score Over Time",
      y = "Brier Score",
      x = "Time (Years)",
      color = "Model (with IBS)"
    ) +
    theme_minimal()
}