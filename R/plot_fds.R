#'@title Fraction of Design Space Plot
#'
#'@description Creates a fraction of design space plot
#'
#'@param genoutput The design, or the output of the power evaluation functions. This can also be a list
#'of several designs, which will result in all of them being plotted in a row (for easy comparison).
#'@param model Default `NULL`. The model, if `NULL` it defaults to the model used in `eval_design` or `gen_design`.
#'@param continuouslength Default `11`. The precision of the continuous variables. Decrease for faster (but less precise) plotting.
#'@param plot Default `TRUE`. Whether to plot the FDS, or just calculate the cumulative distribution function.
#'@param sample_size Default `10000`. Number of samples to take of the design space.
#'@param yaxis_max Default `NULL`. Manually set the maximum value of the prediction variance.
#'@param description Default `Fraction of Design Space`. The description to add to the plot. If a vector and multiple designs
#'passed to genoutput, it will be the description for each plot.
#'@return Plots design diagnostics, and invisibly returns the vector of values representing the fraction of design space plot. If multiple
#'designs are passed, this will return a list of all FDS vectors.
#'@import graphics grDevices
#'@export
#'@examples
#'#We can pass either the output of gen_design or eval_design to plot_correlations
#'#in order to obtain the correlation map. Passing the output of eval_design is useful
#'#if you want to plot the correlation map from an externally generated design.
#'
#'#First generate the design:
#'
#'candidatelist = expand.grid(X1 = c(1, -1), X2 = c(1, -1))
#'
#'design = gen_design(candidatelist, ~(X1 + X2), 15)
#'
#'plot_fds(design)
plot_fds = function(
  genoutput,
  model = NULL,
  continuouslength = 1001,
  plot = TRUE,
  sample_size = 10000,
  yaxis_max = NULL,
  description = "Fraction of Design Space"
) {
  if (inherits(genoutput, "list") && length(genoutput) > 1) {
    old.par = par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(1, length(genoutput)))
    fds_values = list()
    if (!plot && !is.null(yaxis_max)) {
      warning(
        "`plot = FALSE` but `yaxis_max` non-NULL. Setting `yaxis_max` to NULL"
      )
      yaxis_max = NULL
    }
    if (is.null(yaxis_max)) {
      for (i in 1:length(genoutput)) {
        fds_values[[i]] = plot_fds(
          genoutput[[i]],
          model = model,
          continuouslength = continuouslength,
          plot = FALSE
        )
      }
      yaxis_max = max(unlist(fds_values)) + max(unlist(fds_values)) / 20
    }
    if (length(description) == 1) {
      description = rep(description, length(genoutput))
    }
    if (plot) {
      for (i in 1:length(genoutput)) {
        fds_values[[i]] = plot_fds(
          genoutput[[i]],
          model = model,
          continuouslength = continuouslength,
          plot = plot,
          yaxis_max = yaxis_max,
          description = description[i]
        )
      }
    }
    return(invisible(fds_values))
  }
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

  Iopt = attr(genoutput, "I")
  if(is.null(Iopt)) {
    stop("No I-optimality value found in design--was your design generated outside of skpr? If so, pass in a high resolution candidate set to `high_resolution_candidate_set` to ensure I-optimality is computed.")
  }
  V = attr(genoutput, "variance.matrix")

  if (is.null(model)) {
    if (!is.null(attr(genoutput, "generating.model"))) {
      model = attr(genoutput, "generating.model")
    } else {
      model = ~.
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
      contrastlist[[factornames[name]]] = "contr.sum"
    }
  } else {
    contrastlist = NULL
  }
  sample_list = list()

  for (col in 1:ncol(genoutput)) {
    if (inherits(genoutput[, col], c("factor", "character"))) {
      vals = unique(genoutput[, col])
    }
    if (is.numeric(genoutput[, col])) {
      vals = seq(-1, 1, length.out = continuouslength)
    }
    sample_list[[colnames(genoutput)[col]]] = vals[sample(
      seq_len(length(vals)),
      size = sample_size,
      replace = TRUE
    )]
  }
  samples = as.data.frame(sample_list)

  #------Normalize/Center numeric columns ------#
  for (column in 1:ncol(genoutput)) {
    if (is.numeric(genoutput[, column])) {
      midvalue = mean(c(max(genoutput[, column]), min(genoutput[, column])))
      genoutput[, column] = (genoutput[, column] - midvalue) /
        (max(genoutput[, column]) - midvalue)
    }
  }
  mm = model.matrix(model, genoutput, contrasts.arg = contrastlist)
  samplemm = model.matrix(model, samples, contrasts.arg = contrastlist)

  testcor = solve(t(mm) %*% solve(V) %*% mm)

  v = list()

  for (i in 1:nrow(samplemm)) {
    xi = samplemm[i, ]
    v[[i]] = t(xi) %*% testcor %*% xi
  }

  vars = do.call(rbind, v)
  varsordered = vars[order(vars)]
  meanindex = which(
    abs(mean(varsordered) - varsordered) ==
      min(abs(mean(varsordered) - varsordered))
  )

  scale = varsordered[meanindex]
  if (length(scale) > 1) {
    scale = scale[1]
  }
  varsorderedscaled = varsordered / scale * Iopt
  midval = varsorderedscaled[sample_size / 2]
  if (is.null(yaxis_max)) {
    maxyaxis = max(varsorderedscaled) + max(varsorderedscaled) / 20
  } else {
    maxyaxis = yaxis_max
  }
  if (plot) {
    plot(
      1:length(varsorderedscaled) / length(varsorderedscaled),
      varsorderedscaled,
      ylim = c(0, maxyaxis),
      type = "n",
      xlab = description,
      ylab = "Prediction Variance",
      xlim = c(0, 1),
      xaxs = "i",
      yaxs = "i"
    )
    abline(v = 0.5, untf = FALSE, lty = 2, col = "red", lwd = 2)
    abline(h = midval, untf = FALSE, lty = 2, col = "red", lwd = 2)
    lines(
      1:length(varsorderedscaled) / length(varsorderedscaled),
      varsorderedscaled,
      lwd = 2,
      col = "blue"
    )
    invisible(varsorderedscaled)
  } else {
    return(varsorderedscaled)
  }
}
