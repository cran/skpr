#'@title Calculate G Efficiency
#'
#'Either calculates G-Efficiency by Monte Carlo sampling from the design space (ignoring constraints),
#'searching for the maximum point (slower but higher quality), or by using a user-specified candidate
#'set for the design space (fastest).
#'
#'@description Normalizes the numeric columns in the runmatrix
#'@param RunMatrix The run matrix
#'@return Normalized run matrix
#'@keywords internal
calculate_gefficiency = function(design, calculation_type = "random", randsearches = 1000,  design_space_mm = NULL) {
  variables = all.vars(get_attribute(design,"model"))
  designmm = get_attribute(design,"model.matrix")
  modelentries = names(calculate_level_vector(design,get_attribute(design,"model"), FALSE))
  model_entries_mul = gsub(":", "*", modelentries)
  model_entries_mul[1] = 1
  fulllist = list()
  factorvars = attr(design,"contrastslist")
  variables = variables[!variables %in% names(factorvars)]
  rand_vector = function(factorlist, model_entries_mul) {
    fulloutput = unlist(lapply(model_entries_mul, lazyeval::lazy_eval, data = factorlist))
    matrix(fulloutput,1,length(fulloutput))
  }
  calculate_optimality_for_point = function(x, infval = FALSE) {
    for(i in seq_len(length(x))) {
      fulllist[[i]] = x[i]
    }
    names(fulllist) = variables
    for(i in seq_len(length(names(factorvars)))) {
      numberlevels = length(levels(design[[names(factorvars)[i] ]]))
      fulllist[[names(factorvars)[i] ]] = factorvars[[names(factorvars)[i] ]](numberlevels)[sample(seq_len(numberlevels),1),,drop=TRUE]
    }
    mc_mm = rand_vector(fulllist, model_entries_mul)
    if(any(x > 1) || any(x < -1)) {
      if(infval) {
        return(Inf)
      } else {
        return(10000)
      }
    }
    100 * (ncol(designmm)) / (nrow(designmm) * max(diag(mc_mm %*% solve(t(designmm) %*% designmm) %*% t(mc_mm))))
  }
  modelentries = names(calculate_level_vector(design,get_attribute(design,"model"), FALSE))
  factorvars = attr(design,"contrastslist")
  variables = all.vars(get_attribute(design,"model"))
  variables = variables[!variables %in% names(factorvars)]
  if(calculation_type == "random") {
    vals = list()
    lowest = 100
    for(i in 1:randsearches) {
      vals[[i]] = calculate_optimality_for_point(2*runif(length(variables))-1)
      if(vals[[i]] < lowest) {
        lowest = vals[[i]]
      }
    }
    return(min(unlist(vals)))
  }
  if(calculation_type == "optim") {
    vals = list()
    lowest = 100
    for(i in 1:randsearches) {
      vals[[i]] = optim(2*runif(length(variables))-1,calculate_optimality_for_point, method = "SANN")$value
      if(vals[[i]] < lowest) {
        lowest = vals[[i]]
      }
    }
    return(min(unlist(vals)))
  }
  if(calculation_type == "custom" && !is.null(design_space_mm)) {
    return(GEfficiency(designmm, design_space_mm))
  }
}
