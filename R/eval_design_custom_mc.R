#'@title Monte Carlo power evaluation for experimental designs with user-supplied libraries
#'
#'@description Evaluates the power of an experimental design, given its run matrix and the
#'statistical model to be fit to the data, using monte carlo simulation. Simulated data is fit using a
#'user-supplied fitting library and power is estimated by the fraction of times a parameter is significant. Returns
#'a data frame of parameter powers.
#'
#'@param design The experimental design. Internally, \code{eval_design_custom_mc} rescales each numeric column
#'to the range [-1, 1].
#'@param model The model used in evaluating the design. If this is missing and the design
#'was generated with skpr, the generating model will be used. It can be a subset of the model used to
#'generate the design, or include higher order effects not in the original design generation. It cannot include
#'factors that are not present in the experimental design.
#'@param alpha Default `0.05`. The type-I error. p-values less than this will be counted as significant.
#'@param nsim The number of simulations.
#'@param rfunction Random number generator function. Should be a function of the form f(X, b), where X is the
#'model matrix and b are the anticipated coefficients.
#'@param fitfunction Function used to fit the data. Should be of the form f(formula, X, contrasts)
#'where X is the model matrix. If contrasts do not need to be specified for the user supplied
#'library, that argument can be ignored.
#'@param pvalfunction Function that returns a vector of p-values from the object returned from the fitfunction.
#'@param coef_function Function that, when applied to a fitfunction return object, returns the estimated coefficients.
#'@param calceffect  Default `FALSE`. Calculates effect power for a Type-III Anova (using the car package) using a Wald test.
#'this ratio can be a vector specifying the variance ratio for each subplot. Otherwise, it will use a single value for all strata. To work, the
#'fit returned by `fitfunction` must have a method compatable with the car package.
#'@param parameternames Vector of parameter names if the coefficients do not correspond simply to the columns in the model matrix
#'(e.g. coefficients from an MLE fit).
#'@param detailedoutput Default `FALSE`. If `TRUE`, return additional information about evaluation in results.
#'@param advancedoptions Default `NULL`. Named list of advanced options. `advancedoptions$anovatype` specifies the Anova type in the car package (default type `III`),
#'user can change to type `II`). `advancedoptions$anovatest` specifies the test statistic if the user does not want a `Wald` test--other options are likelyhood-ratio `LR` and F-test `F`.
#'`advancedoptions$progressBarUpdater` is a function called in non-parallel simulations that can be used to update external progress bar.`advancedoptions$GUI` turns off some warning messages when in the GUI.
#'If `advancedoptions$save_simulated_responses = TRUE`, the dataframe will have an attribute `simulated_responses` that contains the simulated responses from the power evaluation. `advancedoptions$ci_error_conf` will
#'set the confidence level for power intervals, which are printed when `detailedoutput = TRUE`.
#'@param anticoef The anticipated coefficients for calculating the power. If missing, coefficients will be
#'automatically generated based on \code{effectsize}.
#'@param effectsize The signal-to-noise ratio. Default 2. For a gaussian model, and for
#'continuous factors, this specifies the difference in response between the highest
#'and lowest levels of a factor (which are +1 and -1 after normalization).
#'More precisely: If you do not specify \code{anticoef}, the anticipated coefficients will be
#'half of \code{effectsize}. If you do specify \code{anticoef}, \code{effectsize} will be ignored.
#'@param contrasts  Default \code{contr.sum}. Function used to generate the contrasts encoding for categorical variables. If the user has specified their own contrasts
#'for a categorical factor using the contrasts function, those will be used. Otherwise, skpr will use contr.sum.
#'@param parallel Default `FALSE`. If `TRUE`, the power simulation will use all but one of the available cores.
#' If the user wants to set the number of cores manually, they can do this by setting `options("cores")` to the desired number (e.g. `options("cores" = parallel::detectCores())`).
#' NOTE: If you have installed BLAS libraries that include multicore support (e.g. Intel MKL that comes with Microsoft R Open), turning on parallel could result in reduced performance.
#'@param progress Default `TRUE`. Whether to include a progress bar.
#'@param parallel_pkgs A vector of strings listing the external packages to be included in each parallel worker.
#'@param ... Additional arguments.
#'@return A data frame consisting of the parameters and their powers. The parameter estimates from the simulations are
#'stored in the 'estimates' attribute.
#'@import foreach doParallel stats doRNG
#'@export
#'@examples #To demonstrate how a user can use their own libraries for Monte Carlo power generation,
#'#We will recreate eval_design_survival_mc using the eval_design_custom_mc framework.
#'
#'#To begin, first let us generate the same design and random generation function shown in the
#'#eval_design_survival_mc examples:
#'
#'basicdesign = expand.grid(a = c(-1, 1), b = c("a","b","c"))
#'design = gen_design(candidateset = basicdesign, model = ~a + b + a:b, trials = 100,
#'                          optimality = "D", repeats = 100)
#'
#'#Random number generating function
#'
#'rsurvival = function(X, b) {
#'  Y = rexp(n = nrow(X), rate = exp(-(X %*% b)))
#'  censored = Y > 1
#'  Y[censored] = 1
#'  return(survival::Surv(time = Y, event = !censored, type = "right"))
#'}
#'
#'#We now need to tell the package how we want to fit our data,
#'#given the formula and the model matrix X (and, if needed, the list of contrasts).
#'#If the contrasts aren't required, "contrastslist" should be set to NULL.
#'#This should return some type of fit object.
#'
#'fitsurv = function(formula, X, contrastslist = NULL) {
#'  return(survival::survreg(formula, data = X, dist = "exponential"))
#'}
#'
#'
#'#We now need to tell the package how to extract the p-values from the fit object returned
#'#from the fit function. This is how to extract the p-values from the survreg fit object:
#'
#'pvalsurv = function(fit) {
#'  return(summary(fit)$table[, 4])
#'}
#'
#'#And now we evaluate the design, passing the fitting function and p-value extracting function
#'#in along with the standard inputs for eval_design_mc.
#'#This has the exact same behavior as eval_design_survival_mc for the exponential distribution.
#'eval_design_custom_mc(design = design, model = ~a + b + a:b,
#'                      alpha = 0.05, nsim = 100,
#'                      fitfunction = fitsurv, pvalfunction = pvalsurv,
#'                      rfunction = rsurvival, effectsize = 1)
#'#We can also use skpr's framework for parallel computation to automatically parallelize this
#'#to speed up computation
#'\dontrun{eval_design_custom_mc(design = design, model = ~a + b + a:b,
#'                          alpha = 0.05, nsim = 1000,
#'                          fitfunction = fitsurv, pvalfunction = pvalsurv,
#'                          rfunction = rsurvival, effectsize = 1,
#'                          parallel = TRUE, parallel_pkgs = "survival")
#'}
eval_design_custom_mc = function(
  design,
  model = NULL,
  alpha = 0.05,
  nsim,
  rfunction,
  fitfunction,
  pvalfunction,
  anticoef,
  effectsize = 2,
  contrasts = contr.sum,
  coef_function = coef,
  calceffect = FALSE,
  detailedoutput = FALSE,
  parameternames = NULL,
  advancedoptions = NULL,
  progress = TRUE,
  parallel = FALSE,
  parallel_pkgs = NULL,
  ...
) {
  if (!is.null(advancedoptions)) {
    if (is.null(advancedoptions$save_simulated_responses)) {
      advancedoptions$save_simulated_responses = FALSE
    }
    if (is.null(advancedoptions$GUI)) {
      advancedoptions$GUI = FALSE
    }
    if (!is.null(advancedoptions$progressBarUpdater)) {
      progressBarUpdater = advancedoptions$progressBarUpdater
    } else {
      progressBarUpdater = NULL
    }
    if (is.null(advancedoptions$alphacorrection)) {
      advancedoptions$alphacorrection = TRUE
    } else {
      if (!advancedoptions$alphacorrection) {
        advancedoptions$alphacorrection = FALSE
      }
    }
  } else {
    advancedoptions = list()
    advancedoptions$GUI = FALSE
    advancedoptions$alphacorrection = TRUE
    progressBarUpdater = NULL
    advancedoptions$save_simulated_responses = FALSE
  }
  if (!is.null(getOption("skpr_progress"))) {
    progress = getOption("skpr_progress")
  }

  if (is.null(advancedoptions$ci_error_conf)) {
    advancedoptions$ci_error_conf = 0.95
  }
  if (!is.null(advancedoptions$anovatype)) {
    anovatype = advancedoptions$anovatype
  } else {
    anovatype = "III"
  }
  if (missing(design)) {
    stop("skpr: No design detected in arguments.")
  }
  if (missing(model) || (is.numeric(model) && missing(alpha))) {
    if (is.numeric(model) && missing(alpha)) {
      alpha = model
    }
    if (is.null(attr(design, "generating.model"))) {
      stop("skpr: No model detected in arguments or in design attributes.")
    } else {
      model = attr(design, "generating.model")
    }
  }
  args = list(...)
  if ("RunMatrix" %in% names(args)) {
    stop("skpr: RunMatrix argument deprecated. Use `design` instead.")
  }
  #detect pre-set contrasts
  presetcontrasts = list()
  for (x in names(design)[
    lapply(design, class) %in% c("character", "factor")
  ]) {
    if (!is.null(attr(design[[x]], "contrasts"))) {
      presetcontrasts[[x]] = attr(design[[x]], "contrasts")
    }
  }
  if (attr(terms.formula(model, data = design), "intercept") == 1) {
    nointercept = FALSE
  } else {
    nointercept = TRUE
  }

  #covert tibbles
  run_matrix_processed = as.data.frame(design)

  #----- Convert dots in formula to terms -----#
  model = convert_model_dots(run_matrix_processed, model)

  #----- Rearrange formula terms by order -----#
  model = rearrange_formula_by_order(model, data = run_matrix_processed)

  #------Normalize/Center numeric columns ------#
  run_matrix_processed = normalize_design(run_matrix_processed)

  #Remove skpr-generated REML blocking indicators if present
  run_matrix_processed = remove_skpr_blockcols(run_matrix_processed)

  #---------- Generating model matrix ----------#
  #remove columns from variables not used in the model
  RunMatrixReduced = reduceRunMatrix(run_matrix_processed, model)

  contrastslist = list()
  for (x in names(RunMatrixReduced[
    lapply(RunMatrixReduced, class) %in% c("character", "factor")
  ])) {
    if (!(x %in% names(presetcontrasts))) {
      contrastslist[[x]] = contrasts
      stats::contrasts(RunMatrixReduced[[x]]) = contrasts
    } else {
      contrastslist[[x]] = presetcontrasts[[x]]
    }
  }

  if (length(contrastslist) < 1) {
    contrastslist = NULL
  }

  ModelMatrix = model.matrix(
    model,
    RunMatrixReduced,
    contrasts.arg = contrastslist
  )
  #We'll need the parameter and effect names for output

  #saving model for return attribute
  generatingmodel = model

  if (is.null(parameternames)) {
    parameter_names = colnames(ModelMatrix)
  } else {
    parameter_names = parameternames
  }

  # autogenerate anticipated coefficients
  if (!missing(effectsize) && !missing(anticoef)) {
    warning(
      "User defined anticipated coefficnets (anticoef) detected; ignoring effectsize argument."
    )
  }
  if (missing(anticoef)) {
    anticoef = gen_anticoef(RunMatrixReduced, model, nointercept) *
      effectsize /
      2
    if (!("(Intercept)" %in% colnames(ModelMatrix))) {
      anticoef = anticoef[-1]
    }
  }
  if (length(anticoef) != dim(ModelMatrix)[2]) {
    stop("skpr: Wrong number of anticipated coefficients")
  }

  num_updates = min(c(nsim, 500))
  progressbarupdates = floor(seq(1, nsim, length.out = num_updates))
  progresscurrent = 1

  model_formula = update.formula(model, Y ~ .)
  nparam = ncol(ModelMatrix)
  RunMatrixReduced$Y = 1
  issued_non_convergence_warning = FALSE
  if (!parallel) {
    power_values = rep(0, length(parameter_names))
    effect_pvals_list = list()
    effect_power_values = list()
    estimates = list()
    if (interactive() && progress) {
      pb = progress::progress_bar$new(
        format = sprintf(
          "  Evaluating [:bar] (:current/:total, :tick_rate sim/s) ETA: :eta"
        ),
        total = nsim,
        clear = TRUE,
        width = 100
      )
    }
    for (j in seq_len(nsim)) {
      #simulate the data.
      RunMatrixReduced$Y = rfunction(ModelMatrix, anticoef)

      #fit a model to the simulated data.
      fit = fitfunction(model_formula, RunMatrixReduced, contrastslist)

      #determine whether beta[i] is significant. If so, increment nsignificant
      pvals = pvalfunction(fit)
      if (any(is.na(pvals))) {
        pvals[is.na(pvals)] = 1
        if (!issued_non_convergence_warning) {
          warning(
            "skpr: NaN or NA values found in calculating p values, it is likely the design does not support the model. ",
            "Reduce the model or increase the number of runs to resolve."
          )
          issued_non_convergence_warning = TRUE
        }
      }
      power_values[pvals < alpha] = power_values[pvals < alpha] + 1
      if (calceffect) {
        effect_pvals = effectpowermc(fit, type = anovatype, test = "Pr(>Chisq)")
        effect_pvals_list[[j]] = effect_pvals
      }
      estimates[[j]] = coef_function(fit)
      if (interactive() && progress && !advancedoptions$GUI) {
        pb$tick()
      }
    }
    if (calceffect) {
      effect_results = do.call(rbind, effect_pvals_list)
      effect_power_names = colnames(effect_results)
      effect_power_matrix = matrix(
        0,
        nrow(effect_results),
        ncol(effect_results)
      )
      effect_power_matrix[effect_results < alpha] = 1
      effect_power_results = apply(effect_power_matrix, 2, sum) / nsim
    }
    power_values = power_values / nsim
  } else {
    if (!getOption("skpr_progress", TRUE)) {
      progressbarupdates = c()
    }
    if (!advancedoptions$GUI && progress) {
      set_up_progressr_handler("Evaluating", "sims")
    }
    run_search = function(iterations) {
      prog = progressr::progressor(steps = nsim)
      foreach::foreach(
        i = iterations,
        .errorhandling = "remove",
        .options.future = list(
          packages = parallel_pkgs,
          globals = c(
            "extractPvalues",
            "pvalfunction",
            "parameter_names",
            "progress",
            "progressbarupdates",
            "model_formula",
            "rfunction",
            "RunMatrixReduced",
            "model.matrix",
            "anticoef",
            "prog",
            "fitfunction",
            "contrastslist",
            "effectpowermc",
            "anovatype",
            "calceffect",
            "alpha",
            "coef_function",
            "nsim",
            "num_updates"
          ),
          seed = TRUE
        )
      ) %dofuture%
        {
          if (i %in% progressbarupdates) {
            prog(amount = nsim / num_updates)
          }
          power_values = rep(0, ncol(ModelMatrix))
          #simulate the data.
          RunMatrixReduced$Y = rfunction(ModelMatrix, anticoef)

          #fit a model to the simulated data.
          fit = fitfunction(model_formula, RunMatrixReduced, contrastslist)

          #determine whether beta[i] is significant. If so, increment nsignificant
          pvals = pvalfunction(fit)
          if (any(is.na(pvals))) {
            pvals[is.na(pvals)] = 1
            if (!issued_non_convergence_warning) {
              warning(
                "skpr: NaN or NA values found in calculating p values, it is likely the design does not support the model. ",
                "Reduce the model or increase the number of runs to resolve."
              )
              issued_non_convergence_warning = TRUE
            }
          }
          pvals = pvals[order(factor(names(pvals), levels = parameter_names))]
          pvals[is.na(pvals)] = 1
          stopifnot(all(names(pvals) == parameter_names))

          if (calceffect) {
            effect_pvals = effectpowermc(
              fit,
              type = anovatype,
              test = "Pr(>Chisq)"
            )
            effect_pvals[is.na(effect_pvals)] = 1
          }

          power_values[pvals < alpha] = 1
          estimates = coef_function(fit)

          if (!calceffect) {
            list(
              "parameter_power" = power_values,
              "estimates" = estimates,
              "pvals" = pvals
            )
          } else {
            list(
              "parameter_power" = power_values,
              "effect_power" = effect_pvals,
              "estimates" = estimates,
              "pvals" = pvals
            )
          }
        }
    }
    power_estimates = run_search(seq_len(nsim))
    power_values = apply(
      do.call("rbind", lapply(power_estimates, \(x) x$parameter_power)),
      2,
      sum
    ) /
      nsim
    pvals = do.call("rbind", lapply(power_estimates, \(x) x$pvals))
    estimates = do.call("rbind", lapply(power_estimates, \(x) x$estimates))
    if (calceffect) {
      effect_power_results = apply(
        do.call("rbind", lapply(power_estimates, \(x) x$effect_power)),
        2,
        sum
      ) /
        nsim
      effect_power_names = colnames(effect_power_values)
    }
  }
  #output the results (tidy data format)
  power_final = data.frame(
    parameter = parameter_names,
    type = "custom.parameter.power.mc",
    power = power_values
  )
  if (calceffect) {
    effect_power_final = data.frame(
      parameter = effect_power_names,
      type = "custom.effect.power.mc",
      power = effect_power_results
    )
    power_final = rbind(effect_power_final, power_final)
  }
  attr(power_final, "generating.model") = generatingmodel
  attr(power_final, "estimatesnames") = parameter_names
  attr(power_final, "estimates") = estimates

  attr(power_final, "alpha") = alpha
  attr(power_final, "runmatrix") = RunMatrixReduced
  attr(power_final, "anticoef") = anticoef

  if (detailedoutput) {
    if (nrow(power_final) != length(anticoef)) {
      power_final$anticoef = c(
        rep(NA, nrow(power_final) - length(anticoef)),
        anticoef
      )
    } else {
      power_final$anticoef = anticoef
    }
    power_final$alpha = alpha
    power_final$trials = nrow(run_matrix_processed)
    power_final$nsim = nsim
    power_final = add_ci_bounds_mc_power(
      power_final,
      nsim = nsim,
      conf = advancedoptions$ci_error_conf
    )
    attr(power_final, "mc.conf.int") = advancedoptions$ci_error_conf
  }

  if (!inherits(power_final, "skpr_eval_output")) {
    class(power_final) = c("skpr_eval_output", class(power_final))
  }
  return(power_final)
}
globalVariables("i")
