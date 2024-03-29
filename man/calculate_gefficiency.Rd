% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_gefficiency.R
\name{calculate_gefficiency}
\alias{calculate_gefficiency}
\title{Calculate G Efficiency}
\usage{
calculate_gefficiency(
  design,
  calculation_type = "random",
  randsearches = 10000,
  design_space_mm = NULL
)
}
\value{
Normalized run matrix
}
\description{
Either calculates G-Efficiency by Monte Carlo sampling from the design space (ignoring constraints),
searching for the maximum point (slower but higher quality), or by using a user-specified candidate set
for the design space (fastest).
}
\keyword{internal}
