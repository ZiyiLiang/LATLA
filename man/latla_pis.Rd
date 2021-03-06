% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latla_funcs.R
\name{latla_pis}
\alias{latla_pis}
\title{Estimated sparsity level}
\usage{
latla_pis(x, d, pval, tau = 0.9, h = "auto", eps = 0.1)
}
\arguments{
\item{x}{Vector of the primary statistics.}

\item{d}{Distance matrix.}

\item{pval}{p-values of the primary statistics}

\item{tau}{A threshold to approximate the 'null set', roughly speaking,
we deem p-values greater than tau as null points. Default is
set to 0.9, see reference paper for more discussions on
parameter selection.}

\item{h}{Bandwidth for kernel estimation, default is set to 'auto',
where the bandwidth is automatically selected by build-in
method 'ste'(solve the equation).}

\item{eps}{Determines the size of the neighborhood, only m^(1-eps) data
points will be used for kernel estimation. eps should range
from 0 to 1.}
}
\value{
The estimated sparsity levels
}
\description{
The is a kernel-based function that estimates the local sparsity levels
which are subsequently used to construct oracle-assisted weights.
}
