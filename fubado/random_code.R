rdata <- "
  // random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters with random intercepts
"

rpar <- "
  // random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes
"

rtpar <- "
  // random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;
"

rprior <- "
    // random effects
      L_I ~ lkj_corr_cholesky(1);  // uniform lkj prior on cholesky factors
      sigma_I ~ normal(0, 2.5);  // half-normal on subject sd
      to_vector(z_I) ~ std_normal();  // standard normal on standardized effects
"

rgen <- "
  // recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

  // store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
      z_array[j] = to_vector(z[j, ]);
    }

  // random intercept replictations
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));
"




"
data {

}
parameters {

}
transformed parameters {

}
model {
  // priors


  // mixture
    for(n in 1:N) {  // normal hyperprior on means
      target += normal_lpdf(mu[n] | X[n, ] * beta + z[1, I[n]], exp(phi + z[2, I[n]]));
    }
    for(n in 1:N) {  // skellam likelihood on means
      target += skellam_lpmf(Y[n] | mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}
generated quantities {


  // replications for posterior predictive checks
    array[N] real mu_rep; // mu replications
    for(n in 1:N) {
      mu_rep[n] = normal_rng(X[n, ] * beta + z_rep[1, I[n]], exp(phi + z_rep[2, I[n]]));
    }
    array[N] int y_rep; // Y replications
    for(n in 1:N) {
      y_rep[n] = skellam_rng(mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}"
