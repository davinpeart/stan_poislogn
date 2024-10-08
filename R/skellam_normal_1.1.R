# skellam-normal mixture with random intercepts
fit_rm_skellam_normal <- function(data, dv, iv, id,
                                  standardize_beta = FALSE,
                                  prior_beta,
                                  prior_log_phi,
                                  prior_log_delta,
                                  stanargs = list(
  chains = 4, warmup = 5000, iter = 10000, cores = 4, refresh = 1,
  seed = runif(n = 1, min = 1, max = 99999), control = list(
    adapt_delta = .99, stepsize = .5))) {
  require(rstan)
  require(magrittr)
  for(i in 1:length(iv)) {
    data[[iv[i]]] <- magrittr::set_attr(data[[iv[i]]],  "contrasts",  contr.treatment(levels(data[[iv[i]]])))
  }
  design <- model.matrix(formula(paste0(dv, "~ 1 + ", paste0(iv, collapse = "*"))), data)
  rstan::sampling(
    object = rstan::stan_model(
      model_code = "
functions {
  // mean & variance parameterization of skellam lpmf
    real skellam_lpmf(int y, real mu, real sigma) {
      if(mu == 0) {
        return(- sigma + log_modified_bessel_first_kind(abs(y), sqrt(sigma^2)));
      } else {
        return(- sigma + ((log(sigma + mu) - log(sigma - mu)) * y / 2) +
          log_modified_bessel_first_kind(abs(y), sqrt(sigma^2 - mu^2)));
      }
    }

  // skellam rng
    int skellam_rng(real mu, real sigma) {
      return(poisson_rng((sigma + mu) / 2) - poisson_rng((sigma - mu) / 2));
    }
}
data {
  // observations and fixed effects
    int<lower=1> N;  // number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of fixed effects
    matrix[N, K] X;  // fixed effect design matrix

  // random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters with random intercepts
}
parameters {
  // fixed effects
    vector[K] beta;  // non-standardized regression coefficients
    real phi;  // log positive real dispersion parameter = standard deviation of normal
    real delta;  // log positive real difference between mean and variance of skellam

  // random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes
    vector[N] mu;  // mixture mean
}
transformed parameters {
  // random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;
}
model {
  // priors
    //  fixed effects
      beta ~ normal(0, 3);  // positive prior for generalization coefficients
      phi ~ normal(0, 5);  // weakly informative prior for ancillary dispersion
      delta ~ normal(0, 3);  // weakly informative prior for ancillary difference

    // random effects
      L_I ~ lkj_corr_cholesky(1);  // uniform lkj prior on cholesky factors
      sigma_I ~ normal(0, 2.5);  // half-normal on subject sd
      to_vector(z_I) ~ std_normal();  // standard normal on standardized effects

  // mixture
    for(n in 1:N) {  // normal hyperprior on means
      target += normal_lpdf(mu[n] | X[n, ] * beta + z[1, I[n]], exp(phi + z[2, I[n]]));
    }
    for(n in 1:N) {  // skellam likelihood on means
      target += skellam_lpmf(Y[n] | mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}
generated quantities {
  // recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

  // store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
      z_array[j] = to_vector(z[j, ]);
    }

  // replications for posterior predictive checks
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));
    array[N] real mu_rep; // mu replications
    for(n in 1:N) {
      mu_rep[n] = normal_rng(X[n, ] * beta + z_rep[1, I[n]], exp(phi + z_rep[2, I[n]]));
    }
    array[N] int y_rep; // Y replications
    for(n in 1:N) {
      y_rep[n] = skellam_rng(mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}"),
    data = list(
      N = length(data[[dv]]),
      Y = data[[dv]],
      K = ncol(design),
      X = design,
      I = as.integer(data[[id]]),
      N_I = length(levels(data[[id]])),
      J = 3
    ),                      # named list of data
    chains = stanargs[["chain"]],             # number of Markov chains
    warmup = stanargs[["warmup"]],          # number of warm up iterations per chain
    iter = stanargs[["iter"]],            # total number of iterations per chain
    cores = stanargs[["cores"]],              # number of cores (one per chain)
    refresh = stanargs[["refresh"]],            # progress shown
    seed = stanargs[["seed"]],       # set seed
    control = stanargs[["control"]]
    )}

enframe_prop_integer <- function(y) {
  if(is.vector(y)) {
    x <- table(y)
    z <- as.vector(unname(x))
    return(data.frame(integer = as.integer(names(x)), freq = z, prop = z/length(y)))
  }
  if(is.array(y)) {
    I <- dim(y)[1]
    mx <- vector("numeric", I)
    mn <- vector("numeric", I)
    x <- vector("list", I)
    nm <- vector("list", I)
    for(i in 1:I) {
      x[[i]] <- table(y[i, ])/dim(y)[2]
      nm[[i]] <- as.integer(names(x[[i]]))
      mx[i] <- max(nm[[i]])
      mn[i] <- min(nm[[i]])
    }
    z <- seq.int(min(mn), max(mx))
    J <- length(z)
    w <- vector("list", J)
    names(w) <- as.character(z)
    for(j in 1:J) {
      for(i in 1:I) {
        w[[j]][i] <- x[[i]][names(w)[j]]
      }
    }
    op <- data.frame(integer = z)
    op$mean <- 0
    op$lower <- 0
    op$upper <- 0
    for(k in 1:J) {
      op[k, "mean"] <- mean(w[[as.character(op[k, "integer"])]], na.rm = T)
      op[k, "lower"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .05, na.rm = T))
      op[k, "upper"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .95, na.rm = T))
    }
    return(op)
  }
}

enframe_descriptives <- function(y, stat1 = mean, stat2 = var) {
  if(is.vector(y)) {
    return(data.frame(mean = stat1(y), var = stat2(y)))
  }
  if(is.matrix(y)) {
    I <- dim(y)[1]
    s1 <- vector("numeric", I)
    s2 <- vector("numeric", I)
    for(i in 1:I) {
      s1[i] <- stat1(y[i, ])
      s2[i] <- stat2(y[i, ])
    }
    op <- data.frame(s1, s2)
    colnames(op) <- c(as.character(substitute(stat1)), as.character(substitute(stat2)))
    return(op)
  }
}

set_theme <- function() {
  require(tidyverse)
  require(ggdist)
  ggdist::theme_ggdist() +
    ggplot2::theme(text = ggplot2::element_text(
      family = "Times New Roman", colour = "gray30", size = 8),
          legend.text = ggplot2::element_text(
            family = "Times New Roman", colour = "gray30", size = 8),
          legend.title = ggplot2::element_text(
            family = "Times New Roman", colour = "gray30", size = 8),
          strip.background = ggplot2::element_rect(fill="white"),
          plot.title = ggplot2::element_text(hjust = 0.5))
}

pp_int <- function(y, yrep, ynew = NULL, ynew2 = NULL) {
  require(tidyverse)
  graph <-
  ggplot2::ggplot(enframe_prop_integer(y), ggplot2::aes(x = integer)) +
  ggplot2::geom_bar(
    stat = "identity", mapping = ggplot2::aes(y = prop), fill = "#FF88A5") +
  ggplot2::geom_errorbar(
    data = enframe_prop_integer(yrep), mapping = ggplot2::aes(ymin = lower, ymax = upper),
    colour = "#687EC9", width = 0) +
  ggplot2::geom_point(
      data = enframe_prop_integer(yrep), mapping = ggplot2::aes(y = mean), colour = "#687EC9")

  if(is.null(ynew) & is.null(ynew2)) {
    return(graph)
  }
    +
      ggplot2::geom_line(
        data = enframe_prop_integer(ynew), ggplot2::aes(y = prop), colour = "gray55") +
    ggplot2::geom_point(
      data = enframe_prop_integer(ynew2), mapping = ggplot2::aes(y = prop),
      shape = 21, fill = NA, colour = "gray55")
}

pp_stat_dens <- function(yrep, y, ynu1, ynu2) {
  require(tidyverse)
  ramp_pallette <- grDevices::colorRampPalette(c("white", "#687EC9"))
  ggplot2::ggplot(enframe_descriptives(yrep), ggplot2::aes(x = mean, y = var)) +
    ggplot2::geom_density2d_filled(contour_var = "ndensity", bins = 5,
                                   linewidth = 0, alpha = .8, adjust = 1) +
    ggplot2::scale_fill_manual(values = c(ramp_pallette(5))) +
    ggplot2::geom_segment(
      data = enframe_descriptives(y),
      ggplot2::aes(x = mean, xend = mean, y = -Inf,  yend = var),
      colour = "#D86A83", linetype = 2) +
    ggplot2::geom_segment(
      data = enframe_descriptives(y),
      ggplot2::aes(x = -Inf, xend = mean, y = var,  yend = var),
      colour = "#D86A83", linetype = 2) +
    ggplot2::geom_point(data = enframe_descriptives(y), colour = "#D86A83", fill = "#FF88A5",
               alpha = 1, shape = 21, size = 2.5)
  if(is.null(ynu1) & is.null(ynu2)) {
    return(graph)
  } else {
    graph +
      ggplot2::geom_segment(data = enframe_descriptives(ynu1),
                            ggplot2::aes(x = mean, xend = mean, y = -Inf,  yend = var),
                            colour = "gray55", linetype = 3) +
      ggplot2::geom_segment(data = enframe_descriptives(ynu1),
                            ggplot2::aes(x = -Inf, xend = mean, y = var,  yend = var),
                            colour = "gray55", linetype = 3) +
      ggplot2::geom_point(data = enframe_descriptives(ynu1), colour = "gray55", fill = NA,
                          alpha = 1, shape = 21, size = 2.5) +
      ggplot2::geom_segment(data = enframe_descriptives(ynu2),
                            ggplot2::aes(x = mean, xend = mean, y = -Inf,  yend = var),
                            colour = "gray55", linetype = 3) +
      ggplot2::geom_segment(data = enframe_descriptives(ynu2),
                            ggplot2::aes(x = -Inf, xend = mean, y = var,  yend = var),
                            colour = "gray55", linetype = 3) +
      ggplot2::geom_point(data = enframe_descriptives(ynu2),
                          colour = "gray55", fill = NA,
                          alpha = 1, shape = 21, size = 2.5)
  }
}


