sdmTMB:::get_index
function (obj, bias_correct = FALSE, level = 0.95, area = 1, 
          silent = TRUE, ...) 
{
  d <- get_generic(obj, value_name = "link_total", bias_correct = bias_correct, 
                   level = level, trans = exp, area = area, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d
}
<bytecode: 0x7f836bf26ce0>
  <environment: namespace:sdmTMB>
  > 
  > sdmTMB:::get_generic
function (obj, value_name, bias_correct = FALSE, level = 0.95, 
          trans = I, area = 1, silent = TRUE, ...) 
{
  if ((!isTRUE(obj$do_index) && value_name[1] == "link_total") || 
      value_name[1] == "cog_x") {
    if (is.null(obj[["obj"]])) {
      cli_abort(paste0("`obj` needs to be created with ", 
                       "`predict(..., return_tmb_object = TRUE).`"))
    }
    test <- suppressWarnings(tryCatch(obj$obj$report(obj$obj$env$last.par), 
                                      error = function(e) NA))
    if (all(is.na(test))) 
      cli_abort(c("It looks like the model was built with an older version of sdmTMB. ", 
                  "Please refit with the current version."))
    if (bias_correct && obj$fit_obj$control$parallel > 1) {
      cli_warn("Bias correction can be slower with multiple cores; using 1 core.")
      obj$fit_obj$control$parallel <- 1L
    }
    predicted_time <- sort(unique(obj$data[[obj$fit_obj$time]]))
    fitted_time <- sort(unique(obj$fit_obj$data[[obj$fit_obj$time]]))
    if (!all(fitted_time %in% predicted_time)) {
      cli_abort(paste0("Some of the fitted time elements were not predicted ", 
                       "on with `predict.sdmTMB()`. Please include all time elements."))
    }
    if (length(area) != nrow(obj$pred_tmb_data$proj_X_ij[[1]]) && 
        length(area) != 1L) {
      cli_abort("`area` should be of the same length as `nrow(newdata)` or of length 1.")
    }
    if (length(area) == 1L) 
      area <- rep(area, nrow(obj$pred_tmb_data$proj_X_ij[[1]]))
    tmb_data <- obj$pred_tmb_data
    tmb_data$area_i <- area
    if (value_name[1] == "link_total") 
      tmb_data$calc_index_totals <- 1L
    if (value_name[1] == "cog_x") 
      tmb_data$calc_cog <- 1L
    pars <- get_pars(obj$fit_obj)
    eps_name <- "eps_index"
    pars[[eps_name]] <- numeric(0)
    new_obj <- TMB::MakeADFun(data = tmb_data, parameters = pars, 
                              map = obj$fit_obj$tmb_map, random = obj$fit_obj$tmb_random, 
                              DLL = "sdmTMB", silent = silent)
    old_par <- obj$fit_obj$model$par
    new_obj$fn(old_par)
    sr <- TMB::sdreport(new_obj, bias.correct = FALSE, ...)
  }
  else {
    sr <- obj$sd_report
    pars <- get_pars(obj)
    tmb_data <- obj$tmb_data
    obj <- list(fit_obj = obj)
    eps_name <- "eps_index"
  }
  sr_est <- as.list(sr, "Estimate", report = TRUE)
  if (bias_correct && value_name[1] == "link_total") {
    pars[[eps_name]] <- rep(0, length(sr_est$total))
    new_values <- rep(0, length(sr_est$total))
    names(new_values) <- rep(eps_name, length(new_values))
    fixed <- c(obj$fit_obj$model$par, new_values)
    new_obj2 <- TMB::MakeADFun(data = tmb_data, parameters = pars, 
                               map = obj$fit_obj$tmb_map, random = obj$fit_obj$tmb_random, 
                               DLL = "sdmTMB", silent = silent)
    gradient <- new_obj2$gr(fixed)
    corrected_vals <- gradient[names(fixed) == eps_name]
  }
  else {
    if (value_name[1] == "link_total") 
      cli_inform(c("Bias correction is turned off.", "\n        It is recommended to turn this on for final inference."))
  }
  conv <- get_convergence_diagnostics(sr)
  ssr <- summary(sr, "report")
  log_total <- ssr[row.names(ssr) %in% value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  time_name <- obj$fit_obj$time
  names(d) <- c("trans_est", "se")
  if (bias_correct) {
    d$trans_est <- log(corrected_vals)
    d$est <- corrected_vals
  }
  else {
    d$est <- as.numeric(trans(d$trans_est))
  }
  d$lwr <- as.numeric(trans(d$trans_est + stats::qnorm((1 - 
                                                          level)/2) * d$se))
  d$upr <- as.numeric(trans(d$trans_est + stats::qnorm(1 - 
                                                         (1 - level)/2) * d$se))
  d[[time_name]] <- sort(unique(obj$fit_obj$data[[time_name]]))
  d[, c(time_name, "est", "lwr", "upr", "trans_est", "se"), 
    drop = FALSE]
}
<bytecode: 0x7f836bf4d4f8>
  <environment: namespace:sdmTMB>