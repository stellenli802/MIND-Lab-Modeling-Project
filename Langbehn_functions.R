library(survival)
library(dplyr)
library(nloptr)
library(purrr)
library(tibble)
library(glue)

Langbehn_LogLik <- function(params, age, cag, status, verbose = FALSE) {
  if (params[3] < 0.01) return(-1e10)
  
  mu <- params[1] + exp(params[2] - params[3] * cag)
  sigma <- sqrt(params[4] + exp(params[5] - params[6] * cag))
  
  if (any(!is.finite(mu)) || any(!is.finite(sigma)) || any(sigma <= 1e-6)) return(-1e10)
  
  z <- (age - mu) / sigma
  z_scaled <- pi * z / sqrt(3)
  if (verbose && any(abs(z_scaled) > 100)) cat("Extreme z_scaled\n")
  if (any(abs(z_scaled) > 100)) return(-1e10)
  
  loglik <- -status * log(sigma) - z_scaled - (1 + status) * log1p(exp(-z_scaled))
  ll <- sum(loglik)
  if (!is.finite(ll)) return(-1e10)
  
  return(ll)
}

Langbehn_survival <- function(params, age, cag) {
  mu <- params[1] + exp(params[2] - params[3] * cag)
  sigma <- sqrt(params[4] + exp(params[5] - params[6] * cag))
  z <- (age - mu) / sigma
  1 - plogis(pi * z / sqrt(3))
}

Langbehn_cond_CDF <- function(params, current_age, time_to_diag, cag) {
  S_now <- Langbehn_survival(params, current_age, cag)
  S_future <- Langbehn_survival(params, current_age + time_to_diag, cag)
  (S_now - S_future) / S_now
}

grid_search_mu_decay <- function(decay_grid,
                                 age,
                                 cag,
                                 status,
                                 init = c(21.54, 9.56, 35.55, 17.72, 0.327),  # no decay
                                 lb = c(0, 0, 0.001, 0, 0.01),
                                 ub = c(100, 100, 50, 50, 5),
                                 verbose = TRUE) {
  
  results <- list()
  optim_algs <- c("NLOPT_LN_BOBYQA", "NLOPT_LN_NELDERMEAD")
  
  for (decay in decay_grid) {
    if (verbose) cat("===== Decay =", decay, "=====\n")
    
    obj_fn <- function(p) {
      full_params <- append(p, decay, after = 2)
      val <- -Langbehn_LogLik(full_params, age, cag, status)
      if (!is.finite(val)) val <- 1e6
      return(val)
    }
    
    best_fit <- NULL
    best_obj <- Inf
    best_alg <- NA
    for (alg in optim_algs) {
      if (verbose) cat("Trying optimizer:", alg, "...\n")
      fit_try <- tryCatch({
        nloptr(
          x0 = init,
          eval_f = obj_fn,
          lb = lb,
          ub = ub,
          opts = list(
            algorithm = alg,
            maxeval = 1000,
            xtol_rel = 1e-6
          )
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit_try) && fit_try$objective < best_obj) {
        best_fit <- fit_try
        best_obj <- fit_try$objective
        best_alg <- alg
      }
      
      if (verbose && !is.null(fit_try)) {
        cat(glue("Status: {fit_try$status}, Obj: {round(fit_try$objective, 2)}\n"))
      }
    }
    
    results[[as.character(decay)]] <- list(
      decay = decay,
      value = if (!is.null(best_fit)) best_fit$objective else NA,
      par = if (!is.null(best_fit)) append(best_fit$solution, decay, after = 2) else rep(NA, 6),
      converged = if (!is.null(best_fit)) best_fit$status %in% c(1, 4, 5) else FALSE,
      best_alg = best_alg
    )
  }
  
  tibble(
    mu_decay = decay_grid,
    negloglik = sapply(results, function(r) r$value),
    converged = sapply(results, function(r) r$converged),
    params = lapply(results, function(r) r$par),
    best_optimizer = sapply(results, function(r) r$best_alg)
  )
}

fit_langbehn_model <- function(df,
                               decay_vals = seq(0.10, 0.20, by = 0.01),
                               published_params = c(21.54, 9.56, 0.146, 35.55, 17.72, 0.327),
                               verbose = FALSE) {
  grid_results <- grid_search_mu_decay(
    decay_grid = decay_vals,
    age = df$W_age,
    cag = df$CAG,
    status = df$delta,
    init = published_params[-3],
    verbose = verbose
  )
  
  converged_fits <- grid_results %>%
    filter(converged) %>%
    arrange(negloglik)
  
  if (nrow(converged_fits) == 0) {
    warning("Langbehn model failed to converge.")
    return(rep(NA, 6))
  }
  
  best_fit <- converged_fits %>%
    slice(1) %>%
    pull(params) %>%
    .[[1]]
  
  if (verbose) {
    print(grid_results %>% select(mu_decay, negloglik, best_optimizer, converged))
  }
  
  print(glue::glue("Best negloglik = {round(min(grid_results$negloglik, na.rm = TRUE), 3)}"))
  return(best_fit)
}