% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_functions.R
\name{plot_all}
\alias{plot_all}
\title{Visualize the results in different types of plots}
\usage{
plot_all(dep, plots = c("volcano", "heatmap", "single", "freq",
  "comparison"))
}
\arguments{
\item{dep}{SummarizedExperiment,
Data object which has been generated by \code{\link{analyze_dep}}
or the combination of \code{\link{test_diff}}
and \code{\link{add_rejections}}.}

\item{plots}{"volcano", "heatmap", "single", "freq" and/or "comparison",}
}
\value{
Pdfs containg the desired plots.
}
\description{
\code{plot_all} visualizes the results of the differential protein
expression analysis in different types of plots. These are (1)
volcano plots, (2) heatmaps, (3) single protein plots, (4) frequency plots
and/or (5) comparison plots.
}
\examples{
# Load datasets
data <- UbiLength
exp_design <- UbiLength_ExpDesign

# Import and process data
se <- import_MaxQuant(data, exp_design)
processed <- process(se)

# Differential protein expression analysis
dep <- analyze_dep(processed, "control", "Ctrl")

\dontrun{
# Plot all plots
plot_all(dep)
}
}
