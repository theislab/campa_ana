
# Function to extract fold changes and their corresponding
# p-values from measurements made across multiple cells in 
# different wells (nested design).

fit_mixed_model <- function(
  dat,
  var,
  channel_name,
  transform = NULL,
  random_effect = well_name,
  contrast_var = treatment,
  contrast_var_reference = "Control") {
  
  require(tidyverse)
  require(nlme)
  require(emmeans)
  
  var = enquo(var)
  random_effect = enquo(random_effect)
  contrast_var = enquo(contrast_var)
  
  contrast_nm = as_label(contrast_var)
  
  # variable types and factor levels must be adjusted for the model
  dat_local <- filter(dat,channel==channel_name) %>%
    mutate(var = !!var,
           random_effect = !!random_effect,
           contrast_var = relevel(factor(!!contrast_var), contrast_var_reference)) 
  
  # transform input data
  if (is.null(transform)) {
  } else if (transform == "log") {
    dat_local <- mutate(dat_local, var = log(var))
  } else {
    stop(paste0("unknown value ", transform," for transform parameter"))
  }
  
  mod <- nlme::lme(var ~ contrast_var, 
                   data = dat_local,
                   random = ~ 1 | random_effect)
  
  # Use emmeans to look at the overall changes
  results <- contrast(emmeans(mod, ~ contrast_var),
                      interaction = "trt.vs.ctrl1") %>%
    summary() %>% 
    as_tibble() %>%
    rename(contrast_var = contrast_var_trt.vs.ctrl1) %>%
    mutate(comparison = "unnormalised",
           channel = channel_name) %>%
    separate(contrast_var,into=c("contrast_var","reference_level"),sep=" - ") %>%
    rename(!!contrast_nm := contrast_var)
           
  # Transform the results if necessary
  if (is.null(transform)) {
    results <- results %>%
      mutate(mean_difference = estimate,
             mean_difference_95lcl = estimate - qnorm(0.975) * SE,
             mean_difference_95ucl = estimate + qnorm(0.975) * SE) %>%
      select(-estimate,-SE,-t.ratio)
  } else if (transform == "log") {
    results <- results %>%
      mutate(fold_change = exp(estimate),
             fold_change_95lcl = exp(estimate - qnorm(0.975) * SE),
             fold_change_95ucl = exp(estimate + qnorm(0.975) * SE),
             log2_fold_change = log2(fold_change),
             log2_fold_change_95lcl = log2(fold_change_95lcl),
             log2_fold_change_95ucl = log2(fold_change_95ucl)) %>%
      select(-estimate,-SE,-t.ratio)
  } else {
    stop(paste0("unknown value ", transform," for transform parameter"))
  }
  
  return(results)
}

# function to extract compartment-specific fold changes and their corresponding
# p-values as well as the p-values after normalisation to the overall
# fold-change using nlme::lme. Measurements made across multiple cells in 
# different wells (nested design).
fit_mixed_model_per_CSL <- function(
  dat,
  var,
  channel_name,
  transform = NULL,
  object_id = mapobject_id,
  random_effect = well_name,
  contrast_var = treatment,
  contrast_var_reference = "Control",
  group_var = cluster_annotation,
  group_var_reference = "All",
  unnormalised_only = FALSE) {
  
  require(tidyverse)
  require(nlme)
  require(emmeans)
  
  var = enquo(var)
  object_id = enquo(object_id)
  random_effect = enquo(random_effect)
  contrast_var = enquo(contrast_var)
  group_var = enquo(group_var)
  
  contrast_nm = as_label(contrast_var)
  group_nm = as_label(group_var)
  
  # variable types and factor levels must be adjusted for the model
  dat_local <- filter(dat,channel==channel_name) %>%
    mutate(var = !!var,
           object_id = factor(!!object_id),
           random_effect = !!random_effect,
           contrast_var = relevel(factor(!!contrast_var), contrast_var_reference),
           group_var = relevel(factor(!!group_var), group_var_reference),
           group_var_id = as.integer(group_var))
  
  # transform input data
  if (is.null(transform)) {
  } else if (transform == "log") {
    dat_local <- mutate(dat_local, var = log(var))
  } else {
    stop(paste0("unknown value ", transform," for transform parameter"))
  }
  
  mod <- nlme::lme(var ~ contrast_var * group_var, 
                   data = dat_local,
                   # Uncorrelated compartment-specific random intercepts per well
                   random = list(random_effect = nlme::pdDiag(~ 0 + group_var)),
                   # Correlation between compartments within a cell
                   correlation = nlme::corSymm(form = ~ group_var_id | random_effect / object_id),   
                   # Different variances for each compartment
                   weights = nlme::varIdent(form = ~ 1 | contrast_var * group_var),
                   # TODO: check sing.tol with Mark, this was required to get convergence for PABPC1
                   control = list(maxIter = 1000, msMaxIter = 1000,
                                  msMaxEval = 1000, sing.tol=1e-20))
  
  # Use emmeans to look at the overall (unnormalised) terms
  comparison_raw <- contrast(emmeans(mod, ~ contrast_var | group_var),
                             interaction = "trt.vs.ctrl1") %>%
    summary() %>% 
    as_tibble() %>%
    rename(contrast_var = contrast_var_trt.vs.ctrl1) %>%
    mutate(comparison = "unnormalised",
           channel = channel_name) %>%
    rename(!!group_nm := group_var,
           !!contrast_nm := contrast_var)
  
  if (!unnormalised_only) {
    # Provide the option to exclude this step
    
    # Use emmeans to look at the interaction terms
    # Note that we have difficulties estimating the degrees of freedom using 
    # "satterthwaite" method so have set this manually as df = 5
    comparison_normalised <- contrast(emmeans(mod, ~ contrast_var * group_var, df = 5),
                                      interaction = "trt.vs.ctrl1") %>%
      summary() %>% 
      as_tibble() %>%
      rename(contrast_var = contrast_var_trt.vs.ctrl1,
             group_var = group_var_trt.vs.ctrl1) %>%
      mutate(comparison = paste0("relative_to_",group_var_reference),
             channel = channel_name) %>% 
      rename(!!group_nm := group_var,
             !!contrast_nm := contrast_var) 
    results <- bind_rows(comparison_raw,comparison_normalised)
  } else {
    results <- comparison_raw
  }
  
  # Combine and transform the results if necessary
  
  if (is.null(transform)) {
    results <- results %>%
      mutate(mean_difference = estimate,
             mean_difference_95lcl = estimate - qnorm(0.975) * SE,
             mean_difference_95ucl = estimate + qnorm(0.975) * SE) %>%
      select(-estimate,-SE,-t.ratio)
  } else if (transform == "log") {
    results <- results %>%
      mutate(fold_change = exp(estimate),
             fold_change_95lcl = exp(estimate - qnorm(0.975) * SE),
             fold_change_95ucl = exp(estimate + qnorm(0.975) * SE),
             log2_fold_change = log2(fold_change),
             log2_fold_change_95lcl = log2(fold_change_95lcl),
             log2_fold_change_95ucl = log2(fold_change_95ucl)) %>%
      select(-estimate,-SE,-t.ratio)
  } else {
    stop(paste0("unknown value ", transform," for transform parameter"))
  }
  
  return(results)
}

