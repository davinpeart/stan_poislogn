poisson_lognormal <- function(predict_ancillary, repeated_measures, priors,
                              multiple_matrices, formatted_stancode = F) {
  # functions for printing
  add_fixef <- function(ancillaries, parname, full, dpar = F) {
    if(parname %in% ancillaries) {
      paste0(
        if(!dpar) {
          "X_full[n, ] * "
        } else { paste0(
          "vector[K", if(full) {
            "_full"
            }, "] "
        )
        }, "beta_", parname)
    } else {
      paste0(
        if(dpar) {
          "real "
        },
        parname
      )
    }
  }
  add_ranef <- function(logical, index, rep = F) {
    if(logical) {
      paste0(" + z", if(rep) {
        "_rep"
      }, "[", index, ", I[n]]")
    }
  }

  # functions block
  functions <- ""

  # data block
  data <- paste0(
"// observations and fixed effects
    int<lower=1> N;  // number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of fixed effects
    matrix[N, K] X;  // fixed effect design matrix
", if(!is.null(predict_ancillary)) {
"    int<lower=1> K_full;  // number of effects in full design
    matrix[N, K_full] X_full;  // full design matrix"
},
if(repeated_measures) {
  "

// random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters"
}
  )

  # parameters block
  parameters <- paste0(
"// fixed effects
    vector[K] beta;  // regression coefficients
    ", add_fixef(predict_ancillary, "phi", T),
    ";  // log positive real dispersion parameter = standard deviation of normal
    vector[N] mu;  // mixture mean",
    if(repeated_measures) {
      "

// random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes"
    }
  )

  #print transformed parameters
  transformed_parameters <-
    if(repeated_measures) {
"// random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;"
    } else { "" }

  # print priors
  priors <- paste0(
"// priors
  // fixed effects
    beta ~ ", priors[["beta"]], ";
    ", if("phi" %in% predict_ancillary) {
    "beta_"
  },
  "phi ~ ", priors[["phi"]], ";
  "    ,
  if(repeated_measures) {
    paste0("
  // random effects
    L_I ~ ", priors[["cor"]], ";
    sigma_I ~ ", priors[["sd"]], ";
    to_vector(z_I) ~ std_normal();  // standard normal on standardized effects"
    )
  }
  )

  # print likelihood
  likelihood <- paste0("
// likelihood
    for(n in 1:N) {  // normal hyperprior on means
      target += normal_lpdf(mu[n] | X[n, ] * beta", add_ranef(repeated_measures, 1),
                       ", exp(", add_fixef(predict_ancillary, "phi"),
                       add_ranef(repeated_measures, 2), "));
    }
    for(n in 1:N) {  // poisson likelihood on means
      target += poisson_lpmf(Y[n] | exp(mu[n]));
    }"
  )

  # print generated quantities
  generated_quantities <- paste0(
    if(repeated_measures) {
"  // recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

  // store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
        z_array[j] = to_vector(z[j, ]);
    }

  // correlated replications of random intercepts
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));
      "
    }, "
  // replications for posterior predictive checks",
    paste0("
    array[N] real mu_rep; // mu replications
    for(n in 1:N) {
        mu_rep[n] = normal_rng(X[n, ] * beta", add_ranef(repeated_measures, 1, T),
           ", exp(", add_fixef(predict_ancillary, "phi"),
           add_ranef(repeated_measures, 2, T), "));
    }
    array[N] int y_rep; // Y replications
    for(n in 1:N) {
        y_rep[n] = poisson_rng(exp(mu_rep[n]));
    }"
    )
  )

  stan_code <- paste0(
"functions {", functions, "
}
data {
", data, "
}
parameters {
", parameters, "
}
transformed parameters {
", transformed_parameters, "
}
model {
", priors, "
", likelihood, "
}
generated quantities {
", generated_quantities, "
}")

if(formatted_stancode) {
  return(cat(stan_code))
} else {
  return(list(stan_code = stan_code, npar = 2))
}
}
