> sdmTMB:::predict.sdmTMB
function (object, newdata = NULL, type = c("link", "response"), 
          se_fit = FALSE, re_form = NULL, re_form_iid = NULL, nsim = 0, 
          sims_var = "est", model = c(NA, 1, 2), offset = NULL, mcmc_samples = NULL, 
          return_tmb_object = FALSE, return_tmb_report = FALSE, return_tmb_data = FALSE, 
          tmbstan_model = deprecated(), sims = deprecated(), area = deprecated(), 
          ...) 
{
  if ("version" %in% names(object)) {
    check_sdmTMB_version(object$version)
  }
  else {
    cli_abort(c("This looks like a very old version of a model fit.", 
                "Re-fit the model before predicting with it."))
  }
  if (!"xy_cols" %in% names(object$spde)) {
    cli_warn(c("It looks like this model was fit with make_spde().", 
               "Using `xy_cols`, but future versions of sdmTMB may not be compatible with this.", 
               "Please replace make_spde() with make_mesh()."))
  }
  else {
    xy_cols <- object$spde$xy_cols
  }
  if (is_present(tmbstan_model)) {
    deprecate_stop("0.2.2", "predict.sdmTMB(tmbstan_model)", 
                   "predict.sdmTMB(mcmc_samples)")
  }
  if (is_present(area)) {
    deprecate_stop("0.0.22", "predict.sdmTMB(area)", "get_index(area)")
  }
  else {
    area <- 1
  }
  if (is_present(sims)) {
    deprecate_warn("0.0.21", "predict.sdmTMB(sims)", "predict.sdmTMB(nsim)")
  }
  else {
    sims <- nsim
  }
  assert_that(model[[1]] %in% c(NA, 1, 2), msg = "`model` argument not valid; should be one of NA, 1, 2")
  if (missing(model)) {
    if (.has_delta_attr(object)) 
      model <- attr(object, "delta_model_predict")
  }
  model <- model[[1]]
  type <- match.arg(type)
  if (is.null(re_form) && isTRUE(se_fit)) {
    msg <- paste0("Prediction can be slow when `se_fit = TRUE` and random fields ", 
                  "are included (i.e., `re_form = NA`). Consider using the `nsim` argument ", 
                  "to take draws from the joint precision matrix and summarizing the standard ", 
                  "devation of those draws.")
    cli_inform(msg)
  }
  nd_arg_was_null <- FALSE
  if (is.null(newdata)) {
    if (is_delta(object) || nsim > 0 || type == "response" || 
        !is.null(mcmc_samples) || se_fit || !is.null(re_form) || 
        !is.null(re_form_iid)) {
      newdata <- object$data
      nd_arg_was_null <- TRUE
    }
  }
  if (any(grepl("t2\\(", object$formula[[1]])) && !is.null(newdata)) {
    cli_abort("There are unresolved issues with predicting on newdata when the formula includes t2() terms. Either predict with `newdata = NULL` or use s(). Post an issue if you'd like us to prioritize fixing this.")
  }
  sys_calls <- unlist(lapply(sys.calls(), deparse))
  vr <- check_visreg(sys_calls)
  visreg_df <- vr$visreg_df
  if (visreg_df) {
    re_form <- vr$re_form
    se_fit <- vr$se_fit
  }
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, 
                                                                  NA)))
  pop_pred_iid <- (!is.null(re_form_iid) && ((re_form_iid == 
                                                ~0) || identical(re_form_iid, NA)))
  if (pop_pred_iid) {
    exclude_RE <- rep(1L, length(object$tmb_data$exclude_RE))
  }
  else {
    exclude_RE <- object$tmb_data$exclude_RE
  }
  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L
  no_spatial <- as.logical(object$tmb_data$no_spatial)
  if (!is.null(newdata)) {
    if (any(!xy_cols %in% names(newdata)) && isFALSE(pop_pred) && 
        !no_spatial) 
      cli_abort(c("`xy_cols` (the column names for the x and y coordinates) are not in `newdata`.", 
                  "Did you miss specifying the argument `xy_cols` to match your data?", 
                  "The newer `make_mesh()` (vs. `make_spde()`) takes care of this for you."))
    if (object$time == "_sdmTMB_time") 
      newdata[[object$time]] <- 0L
    if (visreg_df) {
      if (!object$time %in% names(newdata)) {
        newdata[[object$time]] <- max(object$data[[object$time]], 
                                      na.rm = TRUE)
      }
    }
    check_time_class(object, newdata)
    original_time <- as.numeric(sort(unique(object$data[[object$time]])))
    new_data_time <- as.numeric(sort(unique(newdata[[object$time]])))
    if (!all(new_data_time %in% original_time)) 
      cli_abort(c("Some new time elements were found in `newdata`. ", 
                  "For now, make sure only time elements from the original dataset are present.", 
                  "If you would like to predict on new time elements,", 
                  "see the `extra_time` argument in `?sdmTMB`."))
    if (!identical(new_data_time, original_time) & isFALSE(pop_pred)) {
      if (isTRUE(return_tmb_object) || nsim > 0) {
        cli_warn(c("The time elements in `newdata` are not identical to those in the original dataset.", 
                   "This is normally fine, but may create problems for index standardization."))
      }
      missing_time <- original_time[!original_time %in% 
                                      new_data_time]
      fake_nd_list <- list()
      fake_nd <- newdata[1L, , drop = FALSE]
      for (.t in seq_along(missing_time)) {
        fake_nd[[object$time]] <- missing_time[.t]
        fake_nd_list[[.t]] <- fake_nd
      }
      fake_nd <- do.call("rbind", fake_nd_list)
      newdata[["_sdmTMB_fake_nd_"]] <- FALSE
      fake_nd[["_sdmTMB_fake_nd_"]] <- TRUE
      newdata <- rbind(newdata, fake_nd)
      if (!is.null(offset)) 
        offset <- c(offset, rep(0, nrow(fake_nd)))
    }
    fake_spatial_added <- FALSE
    if (pop_pred) {
      for (i in c(1, 2)) {
        if (!xy_cols[[i]] %in% names(newdata)) {
          suppressWarnings({
            newdata[[xy_cols[[i]]]] <- mean(object$data[[xy_cols[[i]]]], 
                                            na.rm = TRUE)
            fake_spatial_added <- TRUE
          })
        }
      }
    }
    if (sum(is.na(new_data_time)) > 0) 
      cli_abort(c("There is at least one NA value in the time column.", 
                  "Please remove it."))
    newdata$sdm_orig_id <- seq(1L, nrow(newdata))
    if (!no_spatial) {
      if (requireNamespace("dplyr", quietly = TRUE)) {
        unique_newdata <- dplyr::distinct(newdata[, xy_cols, 
                                                  drop = FALSE])
      }
      else {
        unique_newdata <- unique(newdata[, xy_cols, drop = FALSE])
      }
      unique_newdata[["sdm_spatial_id"]] <- seq(1, nrow(unique_newdata)) - 
        1L
      if (requireNamespace("dplyr", quietly = TRUE)) {
        newdata <- dplyr::left_join(newdata, unique_newdata, 
                                    by = xy_cols)
      }
      else {
        newdata <- base::merge(newdata, unique_newdata, 
                               by = xy_cols, all.x = TRUE, all.y = FALSE)
        newdata <- newdata[order(newdata$sdm_orig_id), 
                           , drop = FALSE]
      }
      proj_mesh <- fmesher::fm_basis(object$spde$mesh, 
                                     loc = as.matrix(unique_newdata[, xy_cols, drop = FALSE]))
    }
    else {
      proj_mesh <- object$spde$A_st
      newdata[[xy_cols[1]]] <- NA_real_
      newdata[[xy_cols[2]]] <- NA_real_
      newdata[["sdm_spatial_id"]] <- rep(0L, nrow(newdata))
    }
    if (length(object$formula) == 1L) {
      thresh <- list(check_and_parse_thresh_params(object$formula[[1]], 
                                                   newdata))
      formula <- list(thresh[[1]]$formula)
    }
    else {
      thresh <- list(check_and_parse_thresh_params(object$formula[[1]], 
                                                   newdata), check_and_parse_thresh_params(object$formula[[2]], 
                                                                                           newdata))
      formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
    }
    nd <- newdata
    response <- get_response(object$formula[[1]])
    sdmTMB_fake_response <- FALSE
    if (!response %in% names(nd)) {
      nd[[response]] <- 0
      sdmTMB_fake_response <- TRUE
    }
    if (!"mgcv" %in% names(object)) 
      object[["mgcv"]] <- FALSE
    RE_names <- object$split_formula[[1]]$barnames
    proj_RE_indexes <- vapply(RE_names, function(x) as.integer(nd[[x]]) - 
                                1L, rep(1L, nrow(nd)))
    if (isFALSE(pop_pred_iid)) {
      for (i in seq_along(RE_names)) {
        levels_fit <- levels(object$data[[RE_names[i]]])
        levels_nd <- levels(newdata[[RE_names[i]]])
        if (sum(!levels_nd %in% levels_fit)) {
          msg <- paste0("Extra levels found in random intercept factor levels for `", 
                        RE_names[i], "`. Please remove them.")
          cli_abort(msg)
        }
      }
    }
    proj_X_ij <- list()
    for (i in seq_along(object$formula)) {
      f2 <- remove_s_and_t2(object$split_formula[[i]]$form_no_bars)
      tt <- stats::terms(f2)
      attr(tt, "predvars") <- attr(object$terms[[i]], "predvars")
      Terms <- stats::delete.response(tt)
      mf <- model.frame(Terms, newdata, xlev = object$xlevels[[i]])
      proj_X_ij[[i]] <- model.matrix(Terms, mf, contrasts.arg = object$contrasts[[i]])
    }
    sm <- parse_smoothers(object$formula[[1]], data = object$data, 
                          newdata = nd, basis_prev = object$smoothers$basis_out)
    if (!is.null(object$time_varying)) 
      proj_X_rw_ik <- model.matrix(object$time_varying, 
                                   data = nd)
    else proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1)
    if (length(area) != nrow(proj_X_ij[[1]]) && length(area) != 
        1L) {
      cli_abort("`area` should be of the same length as `nrow(newdata)` or of length 1.")
    }
    if (!is.null(offset)) {
      if (nrow(proj_X_ij[[1]]) != length(offset)) 
        cli_abort("Prediction offset vector does not equal number of rows in prediction dataset.")
    }
    tmb_data$proj_offset_i <- if (!is.null(offset)) 
      offset
    else rep(0, nrow(proj_X_ij[[1]]))
    if (nd_arg_was_null) 
      tmb_data$proj_offset_i <- tmb_data$offset_i
    tmb_data$proj_X_threshold <- thresh[[1]]$X_threshold
    tmb_data$area_i <- if (length(area) == 1L) 
      rep(area, nrow(proj_X_ij[[1]]))
    else area
    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_X_rw_ik <- proj_X_rw_ik
    tmb_data$proj_RE_indexes <- proj_RE_indexes
    tmb_data$proj_year <- make_year_i(nd[[object$time]])
    tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
    tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$pop_pred <- as.integer(pop_pred)
    tmb_data$exclude_RE <- exclude_RE
    tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
    tmb_data$proj_Zs <- sm$Zs
    tmb_data$proj_Xs <- sm$Xs
    if (!is.null(object$spatial_varying)) {
      z_i <- model.matrix(object$spatial_varying_formula, 
                          newdata)
      .int <- grep("(Intercept)", colnames(z_i))
      if (sum(.int) > 0) 
        z_i <- z_i[, -.int, drop = FALSE]
    }
    else {
      z_i <- matrix(0, nrow(newdata), 0L)
    }
    tmb_data$proj_z_i <- z_i
    epsilon_covariate <- rep(0, length(unique(newdata[[object$time]])))
    if (tmb_data$est_epsilon_model) {
      time_steps <- unique(newdata[[object$time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- newdata[newdata[[object$time]] == 
                                          time_steps[i], object$epsilon_predictor, drop = TRUE][[1]]
      }
    }
    tmb_data$epsilon_predictor <- epsilon_covariate
    if (return_tmb_data) {
      return(tmb_data)
    }
    new_tmb_obj <- TMB::MakeADFun(data = tmb_data, parameters = get_pars(object), 
                                  map = object$tmb_map, random = object$tmb_random, 
                                  DLL = "sdmTMB", silent = TRUE)
    old_par <- object$model$par
    new_tmb_obj$fn(old_par)
    if (sims > 0 && is.null(mcmc_samples)) {
      if (!"jointPrecision" %in% names(object$sd_report) && 
          !has_no_random_effects(object)) {
        message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
        sd_report <- TMB::sdreport(object$tmb_obj, getJointPrecision = TRUE)
      }
      else {
        sd_report <- object$sd_report
      }
      if (has_no_random_effects(object)) {
        t_draws <- t(mvtnorm::rmvnorm(n = sims, mean = sd_report$par.fixed, 
                                      sigma = sd_report$cov.fixed))
        row.names(t_draws) <- NULL
      }
      else {
        t_draws <- rmvnorm_prec(mu = new_tmb_obj$env$last.par.best, 
                                tmb_sd = sd_report, n_sims = sims)
      }
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(mcmc_samples)) {
      t_draws <- mcmc_samples
      if (nsim > 0) {
        if (nsim > ncol(t_draws)) {
          cli_abort("`nsim` must be <= number of MCMC samples.")
        }
        else {
          t_draws <- t_draws[, seq(ncol(t_draws) - nsim + 
                                     1, ncol(t_draws)), drop = FALSE]
        }
      }
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(mcmc_samples) || sims > 0) {
      if (return_tmb_report) 
        return(r)
      .var <- switch(sims_var, est = "proj_eta", est_rf = "proj_rf", 
                     omega_s = "proj_omega_s_A", zeta_s = "proj_zeta_s_A", 
                     epsilon_st = "proj_epsilon_st_A_vec", sims_var)
      out <- lapply(r, `[[`, .var)
      predtype <- as.integer(model[[1]])
      if (isTRUE(object$family$delta) && sims_var == "est") {
        if (predtype %in% c(1L, NA)) {
          out1 <- lapply(out, function(x) x[, 1L, drop = TRUE])
          out1 <- do.call("cbind", out1)
        }
        if (predtype %in% c(2L, NA)) {
          out2 <- lapply(out, function(x) x[, 2L, drop = TRUE])
          out2 <- do.call("cbind", out2)
        }
        if (is.na(predtype)) {
          out <- object$family[[1]]$linkinv(out1) * object$family[[2]]$linkinv(out2)
          if (type != "response") 
            out <- object$family[[2]]$linkfun(out)
        }
        else if (predtype == 1L) {
          out <- out1
          if (type == "response") 
            out <- object$family[[1]]$linkinv(out)
        }
        else if (predtype == 2L) {
          out <- out2
          if (type == "response") 
            out <- object$family[[2]]$linkinv(out)
        }
        else {
          cli_abort("`model` type not valid.")
        }
      }
      else {
        if (isTRUE(object$family$delta) && sims_var != 
            "est" && is.na(model[[1]])) {
          cli_warn("`model` argument was left as NA; defaulting to 1st model component.")
          model <- 1L
        }
        else {
          model <- as.integer(model)
        }
        if (!isTRUE(object$family$delta)) {
          model <- 1L
        }
        if (length(dim(out[[1]])) == 2L) {
          out <- lapply(out, function(.x) .x[, model])
          out <- do.call("cbind", out)
        }
        else if (length(dim(out[[1]])) == 3L) {
          xx <- list()
          for (i in seq_len(dim(out[[1]])[2])) {
            xx[[i]] <- lapply(out, function(.x) .x[, 
                                                   i, model])
            xx[[i]] <- do.call("cbind", xx[[i]])
          }
          out <- xx
          if (sims_var == "zeta_s") 
            names(out) <- object$spatial_varying
          if (length(out) == 1L) 
            out <- out[[1]]
        }
        else {
          cli_abort("Too many dimensions returned from model. Try `return_tmb_report = TRUE` and parse the output yourself.")
        }
        if (type == "response") 
          out <- object$family$linkinv(out)
      }
      if (sims_var == "est") {
        rownames(out) <- nd[[object$time]]
        attr(out, "time") <- object$time
        if (type == "response") {
          attr(out, "link") <- "response"
        }
        else {
          if (isTRUE(object$family$delta)) {
            if (is.na(predtype)) {
              attr(out, "link") <- object$family[[2]]$link
            }
            else if (predtype == 1L) {
              attr(out, "link") <- object$family[[1]]$link
            }
            else if (predtype == 2L) {
              attr(out, "link") <- object$family[[2]]$link
            }
            else {
              cli_abort("`model` type not valid.")
            }
          }
          else {
            attr(out, "link") <- object$family$link
          }
        }
      }
      return(out)
    }
    lp <- new_tmb_obj$env$last.par.best
    r <- new_tmb_obj$report(lp)
    if (return_tmb_report) 
      return(r)
    if (isFALSE(pop_pred)) {
      if (isTRUE(object$family$delta)) {
        nd$est1 <- r$proj_eta[, 1]
        nd$est2 <- r$proj_eta[, 2]
        nd$est_non_rf1 <- r$proj_fe[, 1]
        nd$est_non_rf2 <- r$proj_fe[, 2]
        nd$est_rf1 <- r$proj_rf[, 1]
        nd$est_rf2 <- r$proj_rf[, 2]
        nd$omega_s1 <- r$proj_omega_s_A[, 1]
        nd$omega_s2 <- r$proj_omega_s_A[, 2]
        for (z in seq_len(dim(r$proj_zeta_s_A)[2])) {
          nd[[paste0("zeta_s_", object$spatial_varying[z], 
                     "1")]] <- r$proj_zeta_s_A[, z, 1]
          nd[[paste0("zeta_s_", object$spatial_varying[z], 
                     "2")]] <- r$proj_zeta_s_A[, z, 2]
        }
        nd$epsilon_st1 <- r$proj_epsilon_st_A_vec[, 1]
        nd$epsilon_st2 <- r$proj_epsilon_st_A_vec[, 2]
        if (type == "response" && !se_fit) {
          nd$est1 <- object$family[[1]]$linkinv(nd$est1)
          nd$est2 <- object$family[[2]]$linkinv(nd$est2)
          if (object$tmb_data$poisson_link_delta) {
            .n <- nd$est1
            .p <- 1 - exp(-.n)
            .w <- nd$est2
            .r <- (.n * .w)/.p
            nd$est1 <- .p
            nd$est2 <- .r
            nd$est <- .n * .w
          }
          else {
            nd$est <- nd$est1 * nd$est2
          }
        }
      }
      else {
        nd$est <- r$proj_eta[, 1]
        nd$est_non_rf <- r$proj_fe[, 1]
        nd$est_rf <- r$proj_rf[, 1]
        nd$omega_s <- r$proj_omega_s_A[, 1]
        for (z in seq_len(dim(r$proj_zeta_s_A)[2])) {
          nd[[paste0("zeta_s_", object$spatial_varying[z])]] <- r$proj_zeta_s_A[, 
                                                                                z, 1]
        }
        nd$epsilon_st <- r$proj_epsilon_st_A_vec[, 1]
        if (type == "response") {
          nd$est <- object$family$linkinv(nd$est)
        }
      }
    }
    nd$sdm_spatial_id <- NULL
    nd$sdm_orig_id <- NULL
    obj <- new_tmb_obj
    if ("visreg_model" %in% names(object)) {
      model <- object$visreg_model
    }
    else {
      if (visreg_df) 
        model <- 1L
    }
    if (se_fit) {
      sr <- TMB::sdreport(new_tmb_obj, bias.correct = FALSE)
      sr_est_rep <- as.list(sr, "Estimate", report = TRUE)
      sr_se_rep <- as.list(sr, "Std. Error", report = TRUE)
      if (pop_pred) {
        proj_eta <- sr_est_rep[["proj_fe"]]
        se <- sr_se_rep[["proj_fe"]]
      }
      else {
        proj_eta <- sr_est_rep[["proj_eta"]]
        se <- sr_se_rep[["proj_eta"]]
      }
      if (is.na(model)) 
        model_temp <- 1L
      else model_temp <- model
      proj_eta <- proj_eta[, model_temp, drop = TRUE]
      se <- se[, model_temp, drop = TRUE]
      nd$est <- proj_eta
      nd$est_se <- se
    }
    if (type == "response" && se_fit) {
      est_name <- if (isTRUE(object$family$delta)) 
        "'est1' and 'est2'"
      else "'est'"
      msg <- paste0("predict(..., type = 'response', se_fit = TRUE) detected; ", 
                    "returning the prediction ", est_name, " in link space because the standard errors ", 
                    "are calculated in link space.")
      cli_warn(msg)
      type <- "link"
    }
    if (pop_pred) {
      if (isTRUE(object$family$delta)) {
        if (type == "response") {
          nd$est1 <- object$family[[1]]$linkinv(r$proj_fe[, 
                                                          1])
          nd$est2 <- object$family[[2]]$linkinv(r$proj_fe[, 
                                                          2])
          nd$est <- nd$est1 * nd$est2
        }
        else {
          nd$est1 <- r$proj_fe[, 1]
          nd$est2 <- r$proj_fe[, 2]
          if (is.na(model)) {
            p1 <- object$family[[1]]$linkinv(r$proj_fe[, 
                                                       1])
            p2 <- object$family[[2]]$linkinv(r$proj_fe[, 
                                                       2])
            nd$est <- object$family[[2]]$linkfun(p1 * 
                                                   p1)
            if (se_fit) {
              nd$est <- sr_est_rep$proj_rf_delta
              nd$est_se <- sr_se_rep$proj_rf_delta
            }
          }
        }
      }
      else {
        if (type == "response") {
          nd$est <- object$family$linkinv(r$proj_fe[, 
                                                    1])
        }
        else {
          nd$est <- r$proj_fe[, 1]
        }
      }
    }
    if (pop_pred && visreg_df) {
      nd$est <- r$proj_fe[, model, drop = TRUE]
    }
    orig_dat <- object$tmb_data$y_i
    if (model == 2L && nrow(nd) == nrow(orig_dat) && visreg_df) {
      nd <- nd[!is.na(orig_dat[, 2]), , drop = FALSE]
    }
    if ("sdmTMB_fake_year" %in% names(nd)) {
      nd <- nd[!nd$sdmTMB_fake_year, , drop = FALSE]
      nd$sdmTMB_fake_year <- NULL
    }
    if (fake_spatial_added) {
      for (i in 1:2) nd[[xy_cols[[i]]]] <- NULL
    }
    if (sdmTMB_fake_response) {
      nd[[response]] <- NULL
    }
  }
  else {
    if (se_fit) {
      cli_warn(paste0("Standard errors have not been implemented yet unless you ", 
                      "supply `newdata`. In the meantime you could supply your original data frame ", 
                      "to the `newdata` argument."))
    }
    if (isTRUE(object$family$delta)) {
      cli_abort(c("Delta model prediction not implemented for `newdata = NULL` yet.", 
                  "Please provide your data to `newdata` and include the `offset` vector if needed."))
    }
    nd <- object$data
    lp <- object$tmb_obj$env$last.par.best
    r <- object$tmb_obj$report(lp)
    nd$est <- r$eta_i[, 1]
    nd$est_non_rf <- r$eta_fixed_i[, 1] + r$eta_rw_i[, 1] + 
      r$eta_iid_re_i[, 1]
    nd$est_rf <- r$omega_s_A[, 1] + r$epsilon_st_A_vec[, 
                                                       1]
    if (!is.null(object$spatial_varying_formula)) 
      cli_abort(c("Prediction with `newdata = NULL` is not supported with spatially varying coefficients yet.", 
                  "Please provide your data to `newdata`."))
    nd$omega_s <- r$omega_s_A[, 1]
    nd$epsilon_st <- r$epsilon_st_A_vec[, 1]
    nd <- nd[!nd[[object$time]] %in% object$extra_time, , 
             drop = FALSE]
    obj <- object
  }
  if (!object$tmb_data$include_spatial[1]) {
    nd$omega_s1 <- NULL
    nd$omega_s <- NULL
  }
  if (isTRUE(object$family$delta)) {
    if (!object$tmb_data$include_spatial[2]) {
      nd$omega_s2 <- NULL
    }
  }
  if (as.logical(object$tmb_data$spatial_only)[1]) {
    nd$epsilon_st1 <- NULL
    nd$epsilon_st <- NULL
  }
  if (isTRUE(object$family$delta)) {
    if (as.logical(object$tmb_data$spatial_only)[2]) {
      nd$epsilon_st2 <- NULL
    }
  }
  if (!object$tmb_data$spatial_covariate) {
    nd$zeta_s1 <- NULL
    nd$zeta_s1 <- NULL
    nd$zeta_s <- NULL
  }
  if (no_spatial) 
    nd[, xy_cols] <- NULL
  nd[["_sdmTMB_time"]] <- NULL
  if (no_spatial) 
    nd[["est_rf"]] <- NULL
  if (no_spatial) 
    nd[["est_non_rf"]] <- NULL
  if ("_sdmTMB_fake_nd_" %in% names(nd)) {
    nd <- nd[!nd[["_sdmTMB_fake_nd_"]], , drop = FALSE]
  }
  nd[["_sdmTMB_fake_nd_"]] <- NULL
  row.names(nd) <- NULL
  if (return_tmb_object) {
    return(list(data = nd, report = r, obj = obj, fit_obj = object, 
                pred_tmb_data = tmb_data))
  }
  else {
    if (visreg_df) {
      if (isTRUE(se_fit)) {
        return(list(fit = nd$est, se.fit = nd$est_se))
      }
      else {
        return(nd$est)
      }
    }
    else {
      return(nd)
    }
  }
}
<bytecode: 0x7f836bea1898>
  <environment: namespace:sdmTMB>
  > 