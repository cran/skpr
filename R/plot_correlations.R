#'@title Plots design diagnostics
#'
#'@description Plots design diagnostics
#'
#'@param genoutput The output of either gen_design or eval_design/eval_design_mc
#'@param model Default `NULL`. Defaults to the model used in generating/evaluating
#'the design, augmented with 2-factor interactions. If specified, it will override the default
#'model used to generate/evaluate the design.
#'@param customcolors A vector of colors for customizing the appearance of the colormap
#'@param pow Default 2. The interaction level that the correlation map is showing.
#'@param custompar Default NULL. Custom parameters to pass to the `par` function for base R plotting.
#'@param standardize Default `TRUE`. Whether to standardize (scale to -1 and 1 and center) the continuous numeric columns. Not
#'standardizing the numeric columns can increase multi-collinearity (predictors that are correlated with other predictors in the model).
#'@param plot Default `TRUE`. If `FALSE`, this will return the correlation matrix.
#'@return Silently returns the correlation matrix with the proper row and column names.
#'@import graphics grDevices
#'@export
#'@examples
#'#We can pass either the output of gen_design or eval_design to plot_correlations
#'#in order to obtain the correlation map. Passing the output of eval_design is useful
#'#if you want to plot the correlation map from an externally generated design.
#'
#'#First generate the design:
#'
#'candidatelist = expand.grid(cost = c(15000, 20000), year = c("2001", "2002", "2003", "2004"),
#'                            type = c("SUV", "Sedan", "Hybrid"))
#'cardesign = gen_design(candidatelist, ~(cost+type+year)^2, 30)
#'plot_correlations(cardesign)
#'
#'#We can also increase the level of interactions that are shown by default.
#'
#'plot_correlations(cardesign, pow = 3)
#'
#'#You can also pass in a custom color map.
#'plot_correlations(cardesign, customcolors = c("blue", "grey", "red"))
#'plot_correlations(cardesign, customcolors = c("blue", "green", "yellow", "orange", "red"))
plot_correlations = function(
  genoutput,
  model = NULL,
  customcolors = NULL,
  pow = 2,
  custompar = NULL,
  standardize = TRUE,
  plot = TRUE
) {
  #Remove skpr-generated REML blocking indicators if present
  if (!is.null(attr(genoutput, "splitanalyzable"))) {
    if (attr(genoutput, "splitanalyzable")) {
      allattr = attributes(genoutput)
      remove_cols = which(colnames(genoutput) %in% allattr$splitcolumns)
      if (length(remove_cols) > 0) {
        genoutput = genoutput[, -remove_cols, drop = FALSE]
        allattr$names = allattr$names[-remove_cols]
      }
      attributes(genoutput) = allattr
    }
  }
  if (!is.null(attr(genoutput, "splitcolumns"))) {
    allattr = attributes(genoutput)
    genoutput = genoutput[,
      !(colnames(genoutput) %in% attr(genoutput, "splitcolumns")),
      drop = FALSE
    ]
    allattr$names = allattr$names[
      !allattr$names %in% attr(genoutput, "splitcolumns")
    ]
    attributes(genoutput) = allattr
  }
  if (!is.null(attr(genoutput, "augmented"))) {
    if (attr(genoutput, "augmented")) {
      allattr = attributes(genoutput)
      genoutput = genoutput[,
        !(colnames(genoutput) %in% "Block1"),
        drop = FALSE
      ]
      allattr$names = allattr$names[!allattr$names %in% "Block1"]
      attributes(genoutput) = allattr
    }
  }
  if (is.null(attr(genoutput, "variance.matrix"))) {
    genoutput = eval_design(genoutput, model, 0.2)
  }
  V = attr(genoutput, "variance.matrix")
  if (is.null(model)) {
    if (!is.null(attr(genoutput, "runmatrix"))) {
      variables = paste0("`", colnames(attr(genoutput, "runmatrix")), "`")
    } else {
      variables = paste0("`", colnames(genoutput), "`")
    }
    linearterms = paste(variables, collapse = " + ")
    linearmodel = paste0(c("~", linearterms), collapse = "")
    model1 = as.formula(paste(
      c(
        linearmodel,
        as.character(aliasmodel(as.formula(linearmodel), power = pow)[2])
      ),
      collapse = " + "
    ))
    if (!is.null(attr(genoutput, "generating.model"))) {
      modelfactors = colnames(attr(
        terms.formula(attr(genoutput, "generating.model")),
        "factors"
      ))
      quadmodelfactors = colnames(attr(
        terms.formula(model1, "factors"),
        "factors"
      ))
      otherterms = modelfactors[!modelfactors %in% quadmodelfactors]
      model = as.formula(paste(c(model1, otherterms), collapse = " + "))
    } else {
      model = model1
    }
  }
  if (!is.null(attr(genoutput, "runmatrix"))) {
    genoutput = attr(genoutput, "runmatrix")
  }
  factornames = colnames(genoutput)[
    unlist(lapply(genoutput, class)) %in% c("factor", "character")
  ]
  if (length(factornames) > 0) {
    contrastlist = list()
    for (name in 1:length(factornames)) {
      contrastlist[[factornames[name]]] = contr.simplex
    }
  } else {
    contrastlist = NULL
  }
  #------Normalize/Center numeric columns ------#
  if (standardize) {
    for (column in 1:ncol(genoutput)) {
      if (is.numeric(genoutput[, column])) {
        midvalue = mean(c(max(genoutput[, column]), min(genoutput[, column])))
        genoutput[, column] = (genoutput[, column] - midvalue) /
          (max(genoutput[, column]) - midvalue)
      }
    }
  }

  mm = model.matrix(model, genoutput, contrasts.arg = contrastlist)
  #Generate pseudo inverse as it's likely the model matrix will be singular with extra terms
  cormat = abs(cov2cor(getPseudoInverse(t(mm) %*% solve(V) %*% mm))[-1, -1])
  if (!plot) {
    return(cormat)
  }

  if (is.null(customcolors)) {
    imagecolors = viridis::viridis(101)
  } else {
    imagecolors = colorRampPalette(customcolors)(101)
  }
  if (is.null(custompar)) {
    par(mar = c(5, 3, 7, 0))
  } else {
    do.call(par, custompar)
  }
  image(
    t(cormat[ncol(cormat):1, , drop = FALSE]),
    x = 1:ncol(cormat),
    y = 1:ncol(cormat),
    zlim = c(0, 1),
    asp = 1,
    axes = F,
    col = imagecolors,
    xlab = "",
    ylab = ""
  )
  axis(
    3,
    at = 1:ncol(cormat),
    labels = colnames(mm)[-1],
    pos = ncol(cormat) + 1,
    las = 2,
    hadj = 0,
    cex.axis = 0.8
  )
  axis(
    2,
    at = ncol(cormat):1,
    labels = colnames(mm)[-1],
    pos = 0,
    las = 2,
    hadj = 1,
    cex.axis = 0.8
  )

  legend(
    length(colnames(mm)[-1]) + 1,
    length(colnames(mm)[-1]),
    c("0", "", "", "", "", "0.5", "", "", "", "", "1.0"),
    title = "|r|\n",
    fill = imagecolors[c(seq(1, 101, 10))],
    xpd = TRUE,
    bty = "n",
    border = NA,
    y.intersp = 0.3,
    x.intersp = 0.1,
    cex = 1
  )
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  results = t(cormat[ncol(cormat):1, , drop = FALSE])
  colnames(results) = rev(colnames(mm)[-1])
  rownames(results) = colnames(mm)[-1]
  invisible(results)
}
