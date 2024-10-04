library(tidyverse)
library(readxl)
options(mc.cores = 2)

# import data
morestrous <- read_excel("morestrous_poisson.xlsx") %>%
  select(-c(Date, Box)) %>%
  mutate(across(.cols = c(Phase, Dose, Subject, Sex), .fns = as_factor)) %>%
  mutate(cPhase = as_factor(if_else(Phase == "di", "MD", if_else(Phase == "met",
                                                                 "MD", if_else(Phase == "M", "M",
                                                                               "PE"))))) %>%
  rename("DS" = "Elev", "CS" = "FirstCS", "PR" = "PreCS", "LM" = "Total Activity")

test_model <-
  fit_model(data = morestrous, dv = "DS", iv = c("Dose", "cPhase"), id = "Subject", family = skellam_normal,
            priors = list(beta = "normal(0, 2.5)", phi = "normal(0, 5)", delta = "normal(0, 3)",
                          cor = "lkj_corr_cholesky(1)", sd = "normal(0, 2.5)"))

skellam_normal("phi", TRUE, priors = list(beta = "normal(0,10)",
                                          phi = "normal(0,2)",
                                          delta = "cauchy(0, 0.7)",
                                          cor = "lkj_corr_cholesky(1)",
                                          sd = "normal(0,5)"))

poisson_lognormal("phi", TRUE, priors = list(beta = "normal(0,5)",
                                          phi = "normal(0,5)",
                                          cor = "lkj_corr_cholesky(1)",
                                          sd = "normal(0,5)"), T)

poislogn <-
  fit_model(data = morestrous, dv = "CS", iv = c("Dose", "cPhase"), id = "Subject",
            full_iv = c("Dose", "cPhase"), predict_ancillary = "phi",
            family = poisson_lognormal2,
            priors = list(beta = "normal(0, 3)", phi = "normal(0, 5)",
                          cor = "lkj_corr_cholesky(1)", sd = "normal(0, 5)"),
            chains = 2, warmup = 5000, iter = 10000, cores = 2, refresh = 1,
              seed = 5555, control = list(
                adapt_delta = .99, stepsize = .5, max_treedepth = 15), stancode_only = F)

extract(poislogn, "y_rep")
y <- morestrous$CS
y_rep <- extract(poislogn, "y_rep")$y_rep

pp_int(y, na.omit(y_rep)[, ], NULL, NULL) + set_theme() + xlim(-.5, 18.5)
na.omit(y_rep)[1:2, ] %>% enframe_prop_integer()
enframe_prop_integer(y_rep)

y_rep[, -is.nan(y_rep)] %>% dim
enframe_prop_integer(y)
test_poislogn
