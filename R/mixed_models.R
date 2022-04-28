
# Function to extract fold changes and their corresponding
# p-values from measurements made across multiple cells in 
# different wells (nested design).

fit_mixed_model <- function(
  dat,
  var,
  channel_name = NULL,
  transform = NULL,
  random_effect = well_name,
  contrast_var = treatment,
  contrast_var_reference = "Unperturbed",
  normalisation = "unnormalised") {
  
  require(tidyverse)
  require(nlme)
  require(emmeans)
  
  var = enquo(var)
  random_effect = enquo(random_effect)
  contrast_var = enquo(contrast_var)
  
  contrast_nm = as_label(contrast_var)
  
  # variable types and factor levels must be adjusted for the model
  if (!is.null(channel_name)) {
    dat_local <- filter(dat,channel==channel_name)
  } else {
    dat_local <- dat
  }
  dat_local <- dat_local %>%
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
                      interaction = "trt.vs.ctrl1",
                      mode = "appx-satterthwaite") %>%
    summary() %>% 
    as_tibble() %>%
    rename(contrast_var = contrast_var_trt.vs.ctrl1) %>%
    mutate(comparison = normalisation,
           channel = if_else(is.null(channel_name),NA_character_,channel_name)) %>%
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

# function to extract CSL-specific fold changes and their corresponding
# p-values as well as the p-values after normalisation to the overall
# fold-change using nlme::lme. Measurements made across multiple cells in 
# different wells (nested design).
fit_mixed_model_per_CSL <- function(
  dat,
  var,
  channel_name = NULL,
  transform = NULL,
  object_id = mapobject_id,
  random_effect = well_name,
  contrast_var = treatment,
  contrast_var_reference = "Control",
  group_var = cluster_annotation,
  group_var_reference = "all",
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
  if (!is.null(channel_name)) {
    dat_local <- filter(dat,channel==channel_name)
  } else {
    dat_local <- dat
  }
  dat_local <- dat_local %>%
    mutate(var = !!var,
           object_id = factor(!!object_id),
           random_effect = !!random_effect,
           # set contrast reference level
           contrast_var = relevel(factor(!!contrast_var), contrast_var_reference),
           # set grouping variable reference level
           group_var = relevel(factor(!!group_var), group_var_reference),
           # convert grouping variable to integer
           group_var_id = as.integer(group_var)) %>%
    select(object_id,random_effect,var,contrast_var,group_var,group_var_id)
  
  # check that only two treatments are being compared
  if (length(levels(dat_local$contrast_var)) != 2) {
    stop(paste0(
    "Intended use of fit_mixed_model_per_CSL is to perform pairwise comparisons of contrast_var.
    contrast_var has ",length(levels(dat_local$contrast_var))," levels."))
  }
  
  contrast_var_condition <- setdiff(levels(dat_local$contrast_var),contrast_var_reference)[1]
    
  # transform input data
  if (is.null(transform)) {
  } else if (transform == "log") {
    dat_local <- mutate(dat_local, var = log(var))
  } else {
    stop(paste0("unknown value ", transform," for transform parameter"))
  }
  
  # fit this inside a function to allow error catching with repeats
  do_fit <- function() {
    model_fit <- nlme::lme(var ~ contrast_var * group_var, 
                           data = dat_local,
                           # Correlated compartment-specific random intercepts per well
                           random = list(random_effect = nlme::pdSymm(~ 0 + group_var)),
                           # Correlation between CSLs (groups) within a cell (object_id)
                           correlation = nlme::corSymm(form = ~ group_var_id | random_effect / object_id),   
                           # Different variances for each CSL (group)
                           weights = nlme::varIdent(form = ~ 1 | contrast_var * group_var),
                           control = list(maxIter = 1000, msMaxIter = 1000,
                                          msMaxEval = 1000, sing.tol=1e-20))
    return(model_fit)
  }
  
  # Try to fit the model 3 times
  mod <- NULL
  attempt <- 1
  while( is.null(mod) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      mod <- do_fit()
    )
  } 
  
  # if model cannot be fit with 3 retries, skip and return a tibble of NAs
  if (is.null(mod)) {
    return(tibble(!!group_nm := levels(dat_local$group_var)[2],
                  !!contrast_nm := levels(dat_local$contrast_var)[2],
                  df = NA,
                  p.value = NA,
                  comparison = "unnormalised",
                  channel = if_else(is.null(channel_name),NA_character_,channel_name),
                  fold_change = NA))
  }
  
  # Use emmeans to look at the overall (unnormalised) terms
  comparison_raw <- contrast(emmeans(mod, ~ contrast_var | group_var),
                             interaction = "trt.vs.ctrl1") %>%
    summary() %>% 
    as_tibble() %>%
    rename(contrast_var = contrast_var_trt.vs.ctrl1) %>%
    mutate(comparison = "unnormalised",
           channel = if_else(is.null(channel_name),NA_character_,channel_name)) %>%
    rename(!!group_nm := group_var,
           !!contrast_nm := contrast_var)
  
  if (!unnormalised_only) {
    # Provide the option to exclude this step

    # Use emmeans to look at the interaction terms
    # Note that we have difficulties estimating the degrees of freedom using 
    # "satterthwaite" method so have set this manually as using the containment method
    containment_df <- summary(mod)$tTable[paste0("contrast_var",contrast_var_condition), "DF"]
    comparison_normalised <- contrast(emmeans(mod, ~ contrast_var * group_var, df = containment_df),
                                      interaction = "trt.vs.ctrl1") %>%
      summary() %>% 
      as_tibble() %>%
      rename(contrast_var = contrast_var_trt.vs.ctrl1,
             group_var = group_var_trt.vs.ctrl1) %>%
      mutate(comparison = paste0("relative_to_",group_var_reference),
             channel = if_else(is.null(channel_name),NA_character_,channel_name),
             containment_df = containment_df) %>%
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

