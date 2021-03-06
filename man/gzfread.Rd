% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.R
\name{gzfread}
\alias{gzfread}
\title{gzfread}
\usage{
gzfread(path, sep = "\\t", out_dir = NULL)
}
\arguments{
\item{path}{a character}
}
\value{
data.table
}
\description{
Import text from gzipped or gunzipped file to data.table
}
\details{
Read gzipped file into data.table.
}
\author{
http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
}
\seealso{
Other utility functions: \code{\link{get_bed}},
  \code{\link{get_file_paths}},
  \code{\link{get_variant_table}}
}

