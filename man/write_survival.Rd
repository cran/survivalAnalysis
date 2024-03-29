% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{write_survival}
\alias{write_survival}
\title{Print the essentials of a SurvivalAnalysisUnivariateResult.}
\usage{
write_survival(
  ...,
  file,
  label = NULL,
  p_precision = 3,
  hr_precision = 2,
  time_precision = 1,
  include_end_separator = FALSE,
  timespan_unit = c("days", "months", "years")
)
}
\arguments{
\item{...}{Results generated by \code{\link{analyse_survival}},
or \code{\link{analyse_multivariate}}, or lists of such objects}

\item{file}{A connection, or a character string naming the file to print to.
(see \code{\link{cat}})}

\item{label}{A label describing the result,
or a vector of the same size as results in ... (will then be mapped 1:1)}

\item{p_precision, hr_precision, time_precision}{Precision with which to print floating point values}

\item{include_end_separator}{Boolean:
Append "\\n---\\n" as separator?
Comes handy if printing multiple results following each other}

\item{timespan_unit}{Unit for time spans: "days", "months" or "years"}
}
\value{
None (invisible NULL).
}
\description{
Write complete textual information for one or multiple survival
analysis results in a text file.
}
\details{
As write_survival takes potentially multiple objects, it cannot
return its input in a cleanly defined way.
You can still elegantly combine \code{write_survival} in a pipe followed by
\code{\link{kaplan_meier_plot}} or \code{\link{kaplan_meier_grid}}
for a single input object if you apply the
tee pipe operator \code{\%T>\%} in front of \code{write_survival}.
}
