fit_model <-
  function(data, dv, iv, id = NULL, full_iv = NULL, predict_ancillary = NULL,
           family, priors, chains = 4, warmup = 5000, iter = 10000, cores = 2, refresh = 2,
             seed = runif(n = 1, min = 1, max = 99999), control = list(
               adapt_delta = .99, stepsize = .1, max_treedepth = 20), stancode_only = F) {
  require(rstan)
  # model matrices with treatment contrasts
  if(is.null(full_iv)) {
    full_iv <- iv
  }
  for(i in 1:length(full_iv)) {
    attr(data[[full_iv[i]]], "contrasts") <- contr.treatment(levels(data[[full_iv[i]]]))
  }
  full_needed <- !identical(iv, full_iv) & !is.null(full_iv)
  design <- model.matrix(formula(paste0(
    dv, "~", "1+", paste0(iv, collapse = "*"))), data)
  if(full_needed) {
    design_full <- model.matrix(formula(paste0(
      dv, "~", "1+", paste0(full_iv, collapse = "*"))), data)
  }

  # generate code
  family_list <- family(predict_ancillary = predict_ancillary,
                        repeated_measures = !is.null(id),
                        multiple_matrices = full_needed,
                        priors = priors)
  stan_code <- family_list[["stan_code"]]
  if(stancode_only) {
    return(cat(stan_code))
  }

  stan_data <- list(
    N = length(data[[dv]]),
    Y = data[[dv]],
    K = ncol(design),
    X = design,
    I = as.integer(data[[id]]),
    N_I = length(levels(data[[id]])),
    J = family_list[["npar"]]
  )
  if(full_needed) {
    stan_data$K_full <- ncol(design_full)
    stan_data$X_full <- design_full
  }

  # fit stan model
  rstan::sampling(
    object = rstan::stan_model(
      model_code = stan_code),
    data = stan_data,
    chains = chains,
    warmup = warmup,
    iter = iter,
    cores = cores,
    refresh = refresh,
    seed = seed,
    control = control
  )
  }
