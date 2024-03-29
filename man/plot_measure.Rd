% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_measure}
\alias{plot_measure}
\title{Box plot/Violin plot of gene expressions or meta measures}
\usage{
plot_measure(
  dataset,
  measures,
  group_by,
  split_by = NA,
  pal_setup = "Set2",
  plot_type = "combined"
)
}
\arguments{
\item{dataset}{A Seurat object.}

\item{measures}{A character vector - names of genes or meta measures to plot.}

\item{group_by}{A string - the grouping variable in the metadata.}

\item{split_by}{A string - the splitting variable in the metadata.}

\item{pal_setup}{A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used; 

Or a character vector of \code{RColorBrewer} palettes. Multiple palettes can be specified, in case of many levels to be colored;
 
Or a character vector of colors.}

\item{plot_type}{A string - type of the plot. Should be values in \code{c("box", "violin", "combined")}.}
}
\value{
A plot.
}
\description{
Box plot/Violin plot of gene expressions or meta measures
}
