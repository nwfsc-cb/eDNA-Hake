> sdmTMB::sdmTMB
function (formula, data, mesh, time = NULL, family = gaussian(link = "identity"), 
          spatial = c("on", "off"), spatiotemporal = c("iid", "ar1", 
                                                       "rw", "off"), share_range = TRUE, time_varying = NULL, 
          time_varying_type = c("rw", "rw0", "ar1"), spatial_varying = NULL, 
          weights = NULL, offset = NULL, extra_time = NULL, reml = FALSE, 
          silent = TRUE, anisotropy = FALSE, control = sdmTMBcontrol(), 
          priors = sdmTMBpriors(), knots = NULL, bayesian = FALSE, 
          previous_fit = NULL, do_fit = TRUE, do_index = FALSE, predict_args = NULL, 
          index_args = NULL, experimental = NULL) 
{
  data <- droplevels(data)
  delta <- isTRUE(family$delta)
  n_m <- if (delta) 
    2L
  else 1L
  if (!missing(spatial)) {
    if (length(spatial) > 1 && !is.list(spatial)) {
      cli_abort("`spatial` should be a single value or a list")
    }
  }
  if (!missing(spatiotemporal)) {
    if (length(spatiotemporal) > 1 && !is.list(spatiotemporal)) {
      cli_abort("`spatiotemporal` should be a single value or a list")
    }
    if (delta && !is.list(spatiotemporal)) {
      spatiotemporal <- rep(spatiotemporal[[1]], 2L)
    }
    spatiotemporal <- vapply(seq_along(spatiotemporal), function(i) check_spatiotemporal_arg(spatiotemporal, 
                                                                                             time = time, .which = i), FUN.VALUE = character(1L))
  }
  else {
    if (is.null(time)) 
      spatiotemporal <- rep("off", n_m)
    else spatiotemporal <- rep("iid", n_m)
  }
  if (is.null(time)) {
    spatial_only <- rep(TRUE, n_m)
  }
  else {
    spatial_only <- ifelse(spatiotemporal == "off", TRUE, 
                           FALSE)
  }
  if (is.list(spatial)) {
    spatial <- vapply(spatial, parse_spatial_arg, FUN.VALUE = character(1L))
  }
  else {
    spatial <- rep(parse_spatial_arg(spatial), n_m)
  }
  include_spatial <- "on" %in% spatial
  if (!include_spatial && !is.null(spatial_varying)) {
    omit_spatial_intercept <- TRUE
    include_spatial <- TRUE
    spatial <- rep("on", length(spatial))
  }
  else {
    omit_spatial_intercept <- FALSE
  }
  if (!include_spatial && all(spatiotemporal == "off") || !include_spatial && 
      all(spatial_only)) {
    no_spatial <- TRUE
    if (missing(mesh)) {
      mesh <- sdmTMB::pcod_mesh_2011
    }
  }
  else {
    no_spatial <- FALSE
  }
  share_range <- unlist(share_range)
  if (length(share_range) == 1L) 
    share_range <- rep(share_range, n_m)
  share_range[spatiotemporal == "off"] <- TRUE
  share_range[spatial == "off"] <- TRUE
  spde <- mesh
  epsilon_model <- NULL
  epsilon_predictor <- NULL
  if (!is.null(experimental)) {
    if ("epsilon_predictor" %in% names(experimental)) {
      epsilon_predictor <- experimental$epsilon_predictor
    }
    else {
      epsilon_predictor <- NULL
    }
    if ("epsilon_model" %in% names(experimental)) {
      epsilon_model <- experimental$epsilon_model
    }
    else {
      epsilon_model <- NULL
    }
  }
  normalize <- control$normalize
  nlminb_loops <- control$nlminb_loops
  newton_loops <- control$newton_loops
  quadratic_roots <- control$quadratic_roots
  start <- control$start
  multiphase <- control$multiphase
  map <- control$map
  lower <- control$lower
  upper <- control$upper
  get_joint_precision <- control$get_joint_precision
  upr <- control$censored_upper
  dot_checks <- c("lower", "upper", "profile", "parallel", 
                  "censored_upper", "nlminb_loops", "newton_steps", "mgcv", 
                  "quadratic_roots", "multiphase", "newton_loops", "start", 
                  "map", "get_joint_precision", "normalize")
  .control <- control
  for (i in dot_checks) .control[[i]] <- NULL
  ar1_fields <- spatiotemporal == "ar1"
  rw_fields <- spatiotemporal == "rw"
  assert_that(is.logical(reml), is.logical(anisotropy), is.logical(share_range), 
              is.logical(silent), is.logical(multiphase), is.logical(normalize))
  if (!is.null(spatial_varying)) 
    assert_that(class(spatial_varying) %in% c("formula", 
                                              "list"))
  if (!is.null(time_varying)) 
    assert_that(class(time_varying) %in% c("formula", "list"))
  if (!is.null(previous_fit)) 
    assert_that(identical(class(previous_fit), "sdmTMB"))
  assert_that(is.list(priors))
  assert_that(is.list(.control))
  if (!is.null(time)) 
    assert_that(is.character(time))
  assert_that(inherits(spde, "sdmTMBmesh"))
  assert_that(class(formula) %in% c("formula", "list"))
  assert_that(inherits(data, "data.frame"))
  time_varying_type <- match.arg(time_varying_type)
  if (!is.null(map) && length(map) != length(start)) {
    cli_warn(c("`length(map) != length(start)`.", "You likely want to specify `start` values if you are setting the `map` argument."))
  }
  if (!is.null(time)) {
    assert_that(time %in% names(data), msg = "Specified `time` column is missing from `data`.")
  }
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }
  else {
    if (sum(is.na(data[[time]])) > 1) 
      cli_abort("There is at least one NA value in the time column. Please remove it.")
  }
  if (is.factor(data[[time]])) {
    if (length(levels(data[[time]])) > length(unique(data[[time]]))) {
      cli_abort("The time column is a factor and there are extra factor levels. ", 
                "Please remove these or turn your time column into an integer.")
    }
  }
  if (!no_spatial) {
    if (!identical(nrow(spde$loc_xy), nrow(data))) {
      msg <- c("Number of x-y coordinates in `mesh` does not match `nrow(data)`.", 
               "Is it possible you passed a different data frame to `make_mesh()` and `sdmTMB()`?")
      cli::cli_abort(msg)
    }
  }
  if (family$family[1] == "censored_poisson") {
    if ("lwr" %in% names(experimental) || "upr" %in% names(experimental)) {
      cli_abort("Detected `lwr` or `upr` in `experimental`. `lwr` is no longer needed and `upr` is now specified as `control = sdmTMBcontrol(censored_upper = ...)`.")
    }
    if (is.null(upr)) 
      cli_abort("`censored_upper` must be defined in `control = sdmTMBcontrol()` to use the censored Poisson distribution.")
    assert_that(length(upr) == nrow(data))
  }
  if (is.null(upr)) 
    upr <- Inf
  if (inherits(formula, "formula")) {
    original_formula <- replicate(n_m, list(formula))
    thresh <- list(check_and_parse_thresh_params(formula, 
                                                 data))
    if (delta) {
      formula <- list(thresh[[1]]$formula, thresh[[1]]$formula)
    }
    else {
      formula <- list(thresh[[1]]$formula)
    }
  }
  else {
    original_formula <- formula
    thresh <- list(check_and_parse_thresh_params(formula[[1]], 
                                                 data), check_and_parse_thresh_params(formula[[2]], 
                                                                                      data))
    formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
  }
  if (is.character(offset)) {
    offset <- data[[offset]]
  }
  if (!is.null(extra_time)) {
    data <- expand_time(df = data, time_slices = extra_time, 
                        time_column = time, weights = weights, offset = offset, 
                        upr = upr)
    if (!is.null(offset)) 
      offset <- data[["__sdmTMB_offset__"]]
    if (!is.null(weights)) 
      weights <- data[["__weight_sdmTMB__"]]
    if (!is.null(upr)) 
      upr <- data[["__dcens_upr__"]]
    data[["__dcens_upr__"]] <- NULL
    spde$loc_xy <- as.matrix(data[, spde$xy_cols, drop = FALSE])
    spde$A_st <- fmesher::fm_basis(spde$mesh, loc = spde$loc_xy)
    spde$sdm_spatial_id <- seq(1, nrow(data))
  }
  check_irregalar_time(data, time, spatiotemporal, time_varying)
  spatial_varying_formula <- spatial_varying
  if (!is.null(spatial_varying)) {
    mf1 <- model.frame(spatial_varying, data)
    for (i in seq_len(ncol(mf1))) {
      if (is.character(mf1[[i]])) {
        cli_warn(paste0("Detected '{colnames(mf1)[i]}' as a character term in the ", 
                        "'spatial_varying' formula. We suggest you make this a factor if you plan ", 
                        "to predict with only some factor levels. `as.factor({colnames(mf1)[i]})`."))
      }
    }
    z_i <- model.matrix(spatial_varying, data)
    .int <- sum(grep("(Intercept)", colnames(z_i)) > 0)
    if (length(attr(z_i, "contrasts")) && !.int && !omit_spatial_intercept) {
      msg <- c("Detected predictors with factor levels in `spatial_varying` with the intercept omitted from the `spatial_varying` formula.", 
               "You likely want to set `spatial = 'off'` since the constant spatial field (`omega_s`) also represents a spatial intercept.`")
      cli_inform(paste(msg, collapse = " "))
    }
    .int <- grep("(Intercept)", colnames(z_i))
    if (sum(.int) > 0) 
      z_i <- z_i[, -.int, drop = FALSE]
    spatial_varying <- colnames(z_i)
  }
  else {
    z_i <- matrix(0, nrow(data), 0L)
  }
  n_z <- ncol(z_i)
  if (any(grepl("offset\\(", formula))) 
    cli_abort("Detected `offset()` in formula. Offsets in sdmTMB must be specified via the `offset` argument.")
  contains_offset <- check_offset(formula[[1]])
  split_formula <- list()
  RE_indexes <- list()
  nobs_RE <- list()
  ln_tau_G_index <- list()
  X_ij = list()
  mf <- list()
  mt <- list()
  sm <- list()
  for (ii in seq_along(formula)) {
    contains_offset <- check_offset(formula[[ii]])
    split_formula[[ii]] <- split_form(formula[ii][[1]])
    RE_names <- split_formula[[ii]]$barnames
    fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], 
                                                                        .name = x), TRUE)
    RE_indexes[[ii]] <- vapply(RE_names, function(x) as.integer(data[[x]]) - 
                                 1L, rep(1L, nrow(data)))
    nobs_RE[[ii]] <- unname(apply(RE_indexes[[ii]], 2L, max)) + 
      1L
    if (length(nobs_RE[[ii]]) == 0L) 
      nobs_RE[[ii]] <- 0L
    formula[[ii]] <- split_formula[[ii]]$form_no_bars
    ln_tau_G_index[[ii]] <- unlist(lapply(seq_along(nobs_RE[[ii]]), 
                                          function(i) rep(i, each = nobs_RE[[ii]][i]))) - 1L
    formula_no_sm <- remove_s_and_t2(formula[[ii]])
    X_ij[[ii]] <- model.matrix(formula_no_sm, data)
    mf[[ii]] <- model.frame(formula_no_sm, data)
    if (length(split_formula[[ii]]$bars)) {
      termsfun <- function(x) {
        ff <- eval(substitute(~foo, list(foo = x[[2]])))
        tt <- try(terms(ff, data = mf[[ii]]), silent = TRUE)
        tt
      }
      reXterms <- lapply(split_formula[[ii]]$bars, termsfun)
      if (length(attr(reXterms[[1]], "term.labels"))) 
        cli_abort("This model appears to have a random slope specified (e.g., y ~ (1 + b | group)). sdmTMB currently can only do random intercepts (e.g., y ~ (1 | group)).")
    }
    mt[[ii]] <- attr(mf[[ii]], "terms")
    sm[[ii]] <- parse_smoothers(formula = formula[[ii]], 
                                data = data, knots = knots)
  }
  if (delta) {
    if (any(unlist(lapply(nobs_RE, function(.x) .x > 0)))) {
      if (original_formula[[1]] != original_formula[[2]]) {
        msg <- paste0("For now, if delta models contain random intercepts, both ", 
                      "components must have the same main-effects formula.")
        cli_abort(msg)
      }
    }
    if (any(unlist(lapply(sm, `[[`, "has_smooths")))) {
      if (original_formula[[1]] != original_formula[[2]]) {
        msg <- paste0("For now, if delta models contain smoothers, both components ", 
                      "must have the same main-effects formula.")
        cli_abort(msg)
      }
    }
  }
  RE_indexes <- RE_indexes[[1]]
  nobs_RE <- nobs_RE[[1]]
  ln_tau_G_index <- ln_tau_G_index[[1]]
  sm <- sm[[1]]
  y_i <- model.response(mf[[1]], "numeric")
  if (delta) {
    y_i2 <- model.response(mf[[2]], "numeric")
    if (!identical(y_i, y_i2)) 
      cli_abort("Response variable should be the same in both parts of the delta formula.")
  }
  if (family$family[1] %in% c("Gamma", "lognormal") && min(y_i) <= 
      0 && !delta) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }
  if (family$family[1] == "censored_poisson") {
    assert_that(mean(upr - y_i, na.rm = TRUE) >= 0)
  }
  size <- rep(1, nrow(X_ij[[1]]))
  if (identical(family$family[1], "binomial") && !delta) {
    y_i <- model.response(mf[[1]], type = "any")
    if (is.character(y_i)) {
      y_i <- model.response(mf[[1]], type = "factor")
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
    }
    if (is.factor(y_i)) {
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
      y_i <- pmin(as.numeric(y_i) - 1, 1)
      size <- rep(1, length(y_i))
    }
    else {
      if (is.matrix(y_i)) {
        size <- y_i[, 1] + y_i[, 2]
        yobs <- y_i[, 1]
        y_i <- yobs
      }
      else {
        if (all(y_i %in% c(0, 1))) {
          size <- rep(1, length(y_i))
        }
        else {
          y_i <- weights * y_i
          size <- weights
          weights <- rep(1, length(y_i))
        }
      }
    }
    if (is.logical(y_i)) {
      msg <- paste0("We recommend against using `TRUE`/`FALSE` ", 
                    "response values if you are going to use the `visreg::visreg()` ", 
                    "function after. Consider converting to integer with `as.integer()`.")
      cli_warn(msg)
    }
  }
  if (identical(family$link[1], "log") && min(y_i, na.rm = TRUE) < 
      0 && !delta) {
    cli_abort("`link = 'log'` but the reponse data include values < 0.")
  }
  if (is.null(offset)) 
    offset <- rep(0, length(y_i))
  assert_that(length(offset) == length(y_i), msg = "Offset doesn't match length of data")
  if (!is.null(time_varying)) {
    X_rw_ik <- model.matrix(time_varying, data)
  }
  else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
  }
  n_s <- nrow(spde$mesh$loc)
  barrier <- "spde_barrier" %in% names(spde)
  if (barrier && anisotropy) {
    cli_warn("Using a barrier mesh; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }
  if (any(c(!is.na(priors$matern_s[1:2]), !is.na(priors$matern_st[1:2]))) && 
      anisotropy) {
    cli_warn("Using PC Matern priors; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }
  df <- if (family$family[1] == "student" && "df" %in% names(family)) 
    family$df
  else 3
  est_epsilon_model <- 0L
  epsilon_covariate <- rep(0, length(unique(data[[time]])))
  if (!is.null(epsilon_predictor) & !is.null(epsilon_model)) {
    if (epsilon_model %in% c("trend", "trend-re")) {
      time_steps <- unique(data[[time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- data[data[[time]] == 
                                       time_steps[i], epsilon_predictor, drop = TRUE][[1]]
      }
      est_epsilon_model <- 1L
    }
  }
  est_epsilon_slope <- 0
  if (!is.null(epsilon_model)) {
    if (epsilon_model %in% c("trend", "trend-re")) {
      est_epsilon_slope <- 1L
      est_epsilon_model <- 1L
    }
  }
  est_epsilon_re <- 0
  if (!is.null(epsilon_model)) {
    if (epsilon_model[1] %in% c("re", "trend-re")) {
      est_epsilon_re <- 1L
      est_epsilon_model <- 1L
    }
  }
  priors_b <- priors$b
  .priors <- priors
  .priors$b <- NULL
  if (nrow(priors_b) == 1L && ncol(X_ij[[1]]) > 1L) {
    if (!is.na(priors_b[[1]])) {
      message("Expanding `b` priors to match model matrix.")
    }
    priors_b <- mvnormal(rep(NA, ncol(X_ij[[1]])))
  }
  if (ncol(X_ij[[1]]) > 0 & !identical(nrow(priors_b), ncol(X_ij[[1]]))) 
    cli_abort("The number of 'b' priors does not match the model matrix.")
  if (ncol(priors_b) == 2 && attributes(priors_b)$dist == "normal") {
    if (length(priors_b[, 2]) == 1L) {
      if (is.na(priors_b[, 2])) 
        priors_b[, 2] <- 1
    }
    priors_b <- mvnormal(location = priors_b[, 1], scale = diag(as.numeric(priors_b[, 
                                                                                    2]), ncol = nrow(priors_b)))
  }
  not_na <- which(!is.na(priors_b[, 1]))
  if (length(not_na) == 0L) {
    priors_b_Sigma <- diag(2)
  }
  else {
    Sigma <- as.matrix(priors_b[, -1])
    priors_b_Sigma <- as.matrix(Sigma[not_na, not_na])
  }
  priors_sigma_G <- tidy_sigma_G_priors(.priors$sigma_G, ln_tau_G_index)
  .priors$sigma_G <- NULL
  if (!"A_st" %in% names(spde)) 
    cli_abort("`mesh` was created with an old version of `make_mesh()`.")
  if (delta) 
    y_i <- cbind(ifelse(y_i > 0, 1, 0), ifelse(y_i > 0, y_i, 
                                               NA_real_))
  if (!delta) 
    y_i <- matrix(y_i, ncol = 1L)
  X_ij_list <- list()
  for (i in seq_len(n_m)) X_ij_list[[i]] <- X_ij[[i]]
  n_t <- length(unique(data[[time]]))
  random_walk <- if (!is.null(time_varying)) 
    switch(time_varying_type, rw = 1L, rw0 = 2L, ar1 = 0L)
  else 0L
  tmb_data <- list(y_i = y_i, n_t = n_t, z_i = z_i, offset_i = offset, 
                   proj_offset_i = 0, A_st = spde$A_st, sim_re = if ("sim_re" %in% 
                                                                     names(experimental)) as.integer(experimental$sim_re) else rep(0L, 
                                                                                                                                   6), A_spatial_index = spde$sdm_spatial_id - 1L, year_i = make_year_i(data[[time]]), 
                   ar1_fields = ar1_fields, simulate_t = rep(1L, n_t), rw_fields = rw_fields, 
                   X_ij = X_ij_list, X_rw_ik = X_rw_ik, Zs = sm$Zs, Xs = sm$Xs, 
                   proj_Zs = list(), proj_Xs = matrix(nrow = 0L, ncol = 0L), 
                   b_smooth_start = sm$b_smooth_start, proj_lon = 0, proj_lat = 0, 
                   do_predict = 0L, calc_se = 0L, pop_pred = 0L, short_newdata = 0L, 
                   exclude_RE = rep(0L, ncol(RE_indexes)), weights_i = if (!is.null(weights)) weights else rep(1, 
                                                                                                               length(y_i)), area_i = rep(1, length(y_i)), normalize_in_r = 0L, 
                   flag = 1L, calc_index_totals = 0L, calc_cog = 0L, random_walk = random_walk, 
                   ar1_time = as.integer(!is.null(time_varying) && time_varying_type == 
                                           "ar1"), priors_b_n = length(not_na), priors_b_index = not_na - 
                     1L, priors_b_mean = priors_b[not_na, 1], priors_b_Sigma = priors_b_Sigma, 
                   priors_sigma_G = priors_sigma_G, priors = as.numeric(unlist(.priors)), 
                   share_range = as.integer(if (length(share_range) == 1L) rep(share_range, 
                                                                               2L) else share_range), include_spatial = as.integer(include_spatial), 
                   omit_spatial_intercept = as.integer(omit_spatial_intercept), 
                   proj_mesh = Matrix::Matrix(c(0, 0, 2:0), 3, 5), proj_X_ij = list(matrix(0, 
                                                                                           ncol = 1, nrow = 1)), proj_X_rw_ik = matrix(0, ncol = 1, 
                                                                                                                                       nrow = 1), proj_year = 0, proj_spatial_index = 0, 
                   proj_z_i = matrix(0, nrow = 1, ncol = n_m), spde_aniso = make_anisotropy_spde(spde, 
                                                                                                 anisotropy), spde = get_spde_matrices(spde), barrier = as.integer(barrier), 
                   spde_barrier = make_barrier_spde(spde), barrier_scaling = if (barrier) spde$barrier_scaling else c(1, 
                                                                                                                      1), anisotropy = as.integer(anisotropy), family = .valid_family[family$family], 
                   size = c(size), link = .valid_link[family$link], df = df, 
                   spatial_only = as.integer(spatial_only), spatial_covariate = as.integer(!is.null(spatial_varying)), 
                   calc_quadratic_range = as.integer(quadratic_roots), X_threshold = thresh[[1]]$X_threshold, 
                   proj_X_threshold = 0, threshold_func = thresh[[1]]$threshold_func, 
                   RE_indexes = RE_indexes, proj_RE_indexes = matrix(0, 
                                                                     ncol = 0, nrow = 1), nobs_RE = nobs_RE, ln_tau_G_index = ln_tau_G_index, 
                   n_g = length(unique(ln_tau_G_index)), est_epsilon_model = as.integer(est_epsilon_model), 
                   epsilon_predictor = epsilon_covariate, est_epsilon_slope = as.integer(est_epsilon_slope), 
                   est_epsilon_re = as.integer(est_epsilon_re), has_smooths = as.integer(sm$has_smooths), 
                   upr = upr, lwr = 0L, poisson_link_delta = as.integer(isTRUE(family$type == 
                                                                                 "poisson_link_delta")), stan_flag = as.integer(bayesian), 
                   no_spatial = no_spatial)
  if (is.na(sum(nobs_RE))) {
    cli_abort("One of the groups used in the factor levels is NA - please remove")
  }
  b_thresh <- matrix(0, 2L, n_m)
  if (thresh[[1]]$threshold_func == 2L) 
    b_thresh <- matrix(0, 3L, n_m)
  tmb_params <- list(ln_H_input = matrix(0, nrow = 2L, ncol = n_m), 
                     b_j = rep(0, ncol(X_ij[[1]])), b_j2 = if (delta) rep(0,ncol(X_ij[[2]])) else numeric(0), 
                     bs = if (sm$has_smooths) matrix(0,nrow = ncol(sm$Xs), ncol = n_m) else array(0), ln_tau_O = rep(0, n_m), 
                     ln_tau_Z = matrix(0, n_z, n_m), ln_tau_E = rep(0, n_m), ln_kappa = matrix(0, 2L, n_m), thetaf = 0, 
                     logit_p_mix = 0, log_ratio_mix = 0, ln_phi = rep(0, n_m), 
                     ln_tau_V = matrix(0, ncol(X_rw_ik), n_m), rho_time_unscaled = matrix(0, 
                                                                                          ncol(X_rw_ik), n_m), ar1_phi = rep(0, n_m), ln_tau_G = matrix(0, 
                                                                                                                                                        ncol(RE_indexes), n_m), RE = matrix(0, sum(nobs_RE), 
                                                                                                                                                                                            n_m), b_rw_t = array(0, dim = c(tmb_data$n_t, ncol(X_rw_ik), 
                                                                                                                                                                                                                            n_m)), omega_s = matrix(0, if (!omit_spatial_intercept) n_s else 0L, 
                                                                                                                                                                                                                                                    n_m), zeta_s = array(0, dim = c(n_s, n_z, n_m)), 
                     epsilon_st = array(0, dim = c(n_s, tmb_data$n_t, n_m)), 
                     b_threshold = if (thresh[[1]]$threshold_func == 2L) matrix(0, 
                                                                                3L, n_m) else matrix(0, 2L, n_m), b_epsilon = rep(0, 
                                                                                                                                  n_m), ln_epsilon_re_sigma = rep(0, n_m), epsilon_re = matrix(0, 
                                                                                                                                                                                               tmb_data$n_t, n_m), b_smooth = if (sm$has_smooths) matrix(0, 
                                                                                                                                                                                                                                                         sum(sm$sm_dims), n_m) else array(0), ln_smooth_sigma = if (sm$has_smooths) matrix(0, 
                                                                                                                                                                                                                                                                                                                                           length(sm$sm_dims), n_m) else array(0))
  if (identical(family$link, "inverse") && family$family[1] %in% 
      c("Gamma", "gaussian", "student") && !delta) {
    fam <- family
    if (family$family == "student") 
      fam$family <- "gaussian"
    temp <- mgcv::gam(formula = formula[[1]], data = data, 
                      family = fam)
    tmb_params$b_j <- stats::coef(temp)
  }
  tmb_map <- map_all_params(tmb_params)
  tmb_map$b_j <- NULL
  if (delta) 
    tmb_map$b_j2 <- NULL
  if (family$family[[1]] == "tweedie") 
    tmb_map$thetaf <- NULL
  if (family$family[[1]] %in% c("gamma_mix", "lognormal_mix", 
                                "nbinom2_mix")) {
    tmb_map$log_ratio_mix <- NULL
    tmb_map$logit_p_mix <- NULL
  }
  if (delta) {
    if (family$family[[2]] %in% c("gamma_mix", "lognormal_mix", 
                                  "nbinom2_mix")) {
      tmb_map$log_ratio_mix <- NULL
      tmb_map$logit_p_mix <- NULL
    }
  }
  tmb_map$ln_phi <- rep(1, n_m)
  if (family$family[[1]] %in% c("binomial", "poisson", "censored_poisson")) 
    tmb_map$ln_phi[1] <- factor(NA)
  if (delta) {
    if (family$family[[2]] %in% c("binomial", "poisson", 
                                  "censored_poisson")) 
      tmb_map$ln_phi[2] <- factor(NA)
    else tmb_map$ln_phi[2] <- 2
  }
  tmb_map$ln_phi <- as.factor(tmb_map$ln_phi)
  if (!is.null(thresh[[1]]$threshold_parameter)) 
    tmb_map$b_threshold <- NULL
  if (est_epsilon_re == 1L) {
    tmb_map <- unmap(tmb_map, c("ln_epsilon_re_sigma", "epsilon_re"))
  }
  if (est_epsilon_slope == 1L) {
    tmb_map <- unmap(tmb_map, "b_epsilon")
  }
  if (multiphase && is.null(previous_fit) && do_fit) {
    original_tmb_data <- tmb_data
    tmb_data$no_spatial <- 1L
    tmb_data$include_spatial <- rep(0L, length(spatial))
    if (family$family[[1]] == "censored_poisson") 
      tmb_data$family <- .valid_family["poisson"]
    tmb_obj1 <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params, 
                               map = tmb_map, DLL = "sdmTMB", silent = silent)
    lim <- set_limits(tmb_obj1, lower = lower, upper = upper, 
                      silent = TRUE)
    tmb_opt1 <- stats::nlminb(start = tmb_obj1$par, objective = tmb_obj1$fn, 
                              lower = lim$lower, upper = lim$upper, gradient = tmb_obj1$gr, 
                              control = .control)
    tmb_data <- original_tmb_data
    tmb_params <- tmb_obj1$env$parList()
    tmb_params$b_threshold <- if (thresh[[1]]$threshold_func == 
                                  2L) 
      matrix(0, 3L, n_m)
    else matrix(0, 2L, n_m)
  }
  tmb_random <- c()
  if (any(spatial == "on") && !omit_spatial_intercept) {
    tmb_random <- c(tmb_random, "omega_s")
    tmb_map <- unmap(tmb_map, c("omega_s", "ln_tau_O"))
  }
  if (!all(spatiotemporal == "off")) {
    tmb_random <- c(tmb_random, "epsilon_st")
    tmb_map <- unmap(tmb_map, c("ln_tau_E", "epsilon_st"))
  }
  if (!is.null(spatial_varying)) {
    tmb_random <- c(tmb_random, "zeta_s")
    tmb_map <- unmap(tmb_map, c("zeta_s", "ln_tau_Z"))
  }
  if (anisotropy) 
    tmb_map <- unmap(tmb_map, "ln_H_input")
  if (!is.null(time_varying)) {
    tmb_random <- c(tmb_random, "b_rw_t")
    tmb_map <- unmap(tmb_map, c("b_rw_t", "ln_tau_V"))
    if (time_varying_type == "ar1") 
      tmb_map <- unmap(tmb_map, "rho_time_unscaled")
  }
  if (est_epsilon_re) {
    tmb_random <- c(tmb_random, "epsilon_re")
    tmb_map <- unmap(tmb_map, c("epsilon_re"))
  }
  tmb_map$ar1_phi <- as.numeric(tmb_map$ar1_phi)
  for (i in seq_along(spatiotemporal)) {
    if (spatiotemporal[i] == "ar1") 
      tmb_map$ar1_phi[i] <- i
  }
  tmb_map$ar1_phi <- as.factor(as.integer(as.factor(tmb_map$ar1_phi)))
  if (nobs_RE[[1]] > 0) {
    tmb_random <- c(tmb_random, "RE")
    tmb_map <- unmap(tmb_map, c("ln_tau_G", "RE"))
  }
  if (reml) 
    tmb_random <- c(tmb_random, "b_j")
  if (reml && delta) 
    tmb_random <- c(tmb_random, "b_j2")
  if (sm$has_smooths) {
    if (reml) 
      tmb_random <- c(tmb_random, "bs")
    tmb_random <- c(tmb_random, "b_smooth")
    tmb_map <- unmap(tmb_map, c("b_smooth", "ln_smooth_sigma", 
                                "bs"))
  }
  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }
  tmb_data$normalize_in_r <- as.integer(normalize)
  tmb_data$include_spatial <- as.integer(spatial == "on")
  if (!is.null(previous_fit)) 
    tmb_map <- previous_fit$tmb_map
  tmb_map$ln_kappa <- get_kappa_map(n_m = n_m, spatial = spatial, 
                                    spatiotemporal = spatiotemporal, share_range = share_range)
  for (i in seq_along(start)) {
    cli_inform(c(i = paste0("Initiating `", names(start)[i], 
                            "` at specified starting value(s) of:"), paste0("  ", 
                                                                            paste(round(start[[i]], 3), collapse = ", "))))
    tmb_params[[names(start)[i]]] <- start[[i]]
  }
  if (!is.matrix(tmb_params[["ln_kappa"]]) && "ln_kappa" %in% 
      names(start)) {
    msg <- c("Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).", 
             "It should be the same value in each row if `share_range = TRUE`.")
    cli_abort(msg)
  }
  if (nrow(tmb_params[["ln_kappa"]]) != 2L && "ln_kappa" %in% 
      names(start)) {
    msg <- c("Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).", 
             "It should be the same value in each row if `share_range = TRUE`.")
    cli_abort(msg)
  }
  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL
  data$sdmTMB_X_ <- data$sdmTMB_Y_ <- NULL
  if (delta && "off" %in% spatiotemporal) {
    tmb_map$epsilon_st <- array(seq_len(length(tmb_params$epsilon_st)), 
                                dim = dim(tmb_params$epsilon_st))
    tmb_map$ln_tau_E <- seq_len(length(tmb_params$ln_tau_E))
    for (i in which(spatiotemporal == "off")) {
      tmb_map$epsilon_st[, , i] <- NA
      tmb_map$ln_tau_E[i] <- NA
    }
    tmb_map$epsilon_st <- as.factor(tmb_map$epsilon_st)
    tmb_map$ln_tau_E <- as.factor(tmb_map$ln_tau_E)
  }
  if (delta && "off" %in% spatial) {
    tmb_map$omega_s <- array(seq_len(length(tmb_params$omega_s)), 
                             dim = dim(tmb_params$omega_s))
    tmb_map$ln_tau_O <- seq_len(length(tmb_params$ln_tau_O))
    for (i in which(spatial == "off")) {
      tmb_map$omega_s[, i] <- NA
      tmb_map$ln_tau_O[i] <- NA
    }
    tmb_map$omega_s <- as.factor(tmb_map$omega_s)
    tmb_map$ln_tau_O <- as.factor(tmb_map$ln_tau_O)
  }
  if (anisotropy && delta && !"ln_H_input" %in% names(map)) {
    tmb_map$ln_H_input <- factor(c(1, 2, 1, 2))
  }
  if (tmb_data$threshold_func > 0) 
    tmb_map$b_threshold <- NULL
  if (control$profile && delta) 
    cli_abort("Profile not yet working with delta models.")
  for (i in seq_along(map)) {
    cli_inform(c(i = paste0("Fixing (mapping) `", names(map)[i], 
                            "` at specified starting value(s) of:"), paste0("  ", 
                                                                            paste(round(tmb_params[[names(map)[i]]], 3), collapse = ", "))))
  }
  tmb_map <- c(map, tmb_map)
  prof <- c("b_j")
  if (delta) 
    prof <- c(prof, "b_j2")
  out_structure <- structure(list(data = data, spde = spde, 
                                  formula = original_formula, split_formula = split_formula, 
                                  time_varying = time_varying, threshold_parameter = thresh[[1]]$threshold_parameter, 
                                  threshold_function = thresh[[1]]$threshold_func, epsilon_predictor = epsilon_predictor, 
                                  time = time, family = family, smoothers = sm, response = y_i, 
                                  tmb_data = tmb_data, tmb_params = tmb_params, tmb_map = tmb_map, 
                                  tmb_random = tmb_random, spatial_varying = spatial_varying, 
                                  spatial = spatial, spatiotemporal = spatiotemporal, spatial_varying_formula = spatial_varying_formula, 
                                  reml = reml, priors = priors, nlminb_control = .control, 
                                  control = control, contrasts = lapply(X_ij, attr, which = "contrasts"), 
                                  terms = lapply(mf, attr, which = "terms"), extra_time = extra_time, 
                                  xlevels = lapply(seq_along(mf), function(i) stats::.getXlevels(mt[[i]], 
                                                                                                 mf[[i]])), call = match.call(expand.dots = TRUE), 
                                  version = utils::packageVersion("sdmTMB")), class = "sdmTMB")
  if (do_index) {
    args <- list(object = out_structure, return_tmb_data = TRUE)
    args <- c(args, predict_args)
    tmb_data <- do.call(predict.sdmTMB, args)
    if (!"newdata" %in% names(predict_args)) {
      cli_warn("`newdata` must be supplied if `do_index = TRUE`.")
    }
    if ("bias_correct" %in% names(index_args)) {
      cli_warn("`bias_correct` must be done later with `get_index(..., bias_correct = TRUE)`.")
      index_args$bias_correct <- NULL
    }
    if (!"area" %in% names(index_args)) {
      cli_warn("`area` not supplied to `index_args` but `do_index = TRUE`. Using `area = 1`.")
      if (is.null(index_args)) 
        index_args <- list()
      index_args[["area"]] <- 1
    }
    if (length(index_args$area) == 1L) {
      tmb_data$area_i <- rep(index_args[["area"]], nrow(predict_args[["newdata"]]))
    }
    else {
      if (length(index_args$area) != nrow(predict_args[["newdata"]])) 
        cli_abort("`area` length does not match `nrow(newdata)`.")
      tmb_data$area_i <- index_args[["area"]]
    }
    tmb_data$calc_index_totals <- 1L
    tmb_params[["eps_index"]] <- numeric(0)
    out_structure$do_index <- TRUE
  }
  else {
    out_structure$do_index <- FALSE
  }
  tmb_obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params, 
                            map = tmb_map, profile = if (control$profile) 
                              prof
                            else NULL, random = tmb_random, DLL = "sdmTMB", silent = silent)
  lim <- set_limits(tmb_obj, lower = lower, upper = upper, 
                    loc = spde$mesh$loc, silent = FALSE)
  out_structure$tmb_obj <- tmb_obj
  out_structure$tmb_data <- tmb_data
  out_structure$tmb_params <- tmb_params
  out_structure$lower <- lim$lower
  out_structure$upper <- lim$upper
  if (!do_fit) 
    return(out_structure)
  if (normalize) 
    tmb_obj <- TMB::normalize(tmb_obj, flag = "flag", value = 0)
  if (length(tmb_obj$par)) {
    tmb_opt <- stats::nlminb(start = tmb_obj$par, objective = tmb_obj$fn, 
                             gradient = tmb_obj$gr, lower = lim$lower, upper = lim$upper, 
                             control = .control)
  }
  else {
    tmb_opt <- list(par = tmb_obj$par, objective = tmb_obj$fn(tmb_obj$par))
  }
  if (nlminb_loops > 1) {
    if (!silent) 
      cli_inform("running extra nlminb optimization\n")
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 
                                                1))) {
      temp <- tmb_opt[c("iterations", "evaluations")]
      tmb_opt <- stats::nlminb(start = tmb_opt$par, objective = tmb_obj$fn, 
                               gradient = tmb_obj$gr, control = .control, lower = lim$lower, 
                               upper = lim$upper)
      tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + 
        temp[["iterations"]]
      tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + 
        temp[["evaluations"]]
    }
  }
  if (newton_loops > 0) {
    if (!silent) 
      cli_inform("attempting to improve convergence with optimHess\n")
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(tmb_obj$gr(tmb_opt$par))
      h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, 
                            gr = tmb_obj$gr)
      tmb_opt$par <- tmb_opt$par - solve(h, g)
      tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
    }
  }
  check_bounds(tmb_opt$par, lim$lower, lim$upper)
  sd_report <- TMB::sdreport(tmb_obj, getJointPrecision = get_joint_precision)
  conv <- get_convergence_diagnostics(sd_report)
  out_structure$tmb_obj <- tmb_obj
  out <- c(out_structure, list(model = tmb_opt, sd_report = sd_report, 
                               gradients = conv$final_grads, bad_eig = conv$bad_eig, 
                               pos_def_hessian = sd_report$pdHess))
  `class<-`(out, "sdmTMB")
}
<bytecode: 0x7f835307d0b0>
  <environment: namespace:sdmTMB>