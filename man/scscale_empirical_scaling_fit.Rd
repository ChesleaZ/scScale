\name{scscale_empirical_scaling_fit}
\alias{scscale_fit_empirical_cell_law}
\alias{scscale_fit_empirical_umi_law}
\alias{scscale_fit_empirical_joint_law}
\alias{scscale_empirical_global_ceiling}
\alias{scscale_empirical_cell_curve}
\alias{scscale_empirical_umi_curve}
\alias{scscale_empirical_joint_curve}
\title{Fit Direct Empirical Scaling Laws}
\usage{
scscale_fit_empirical_cell_law(data, n_col = "n", mi_col = "I_empirical",
  start = NULL, control = list(maxit = 10000))

scscale_fit_empirical_umi_law(data, U_col = "U", mi_col = "I_empirical",
  start = NULL, control = list(maxit = 10000))

scscale_fit_empirical_joint_law(data, n_col = "n", U_col = "U",
  mi_col = "I_empirical", cell_fit = NULL, umi_fit = NULL,
  mode = c("hybrid", "transfer", "free"), start = NULL,
  control = list(maxit = 20000))

scscale_empirical_global_ceiling(cell_fit, umi_fit, n0, U0)

scscale_empirical_cell_curve(n, I_infinity, C_n, A_n)

scscale_empirical_umi_curve(U, I_infinity, C_U, A_U, U_offset)

scscale_empirical_joint_curve(n, U, I_infinity, C_n, A_n, C_U, A_U,
  U_offset)
}
\arguments{
\item{data}{Data frame containing empirical MI observations.}
\item{n_col, U_col, mi_col}{Column names for cell number, UMI depth, and empirical MI.}
\item{start}{Optional named list of starting values.}
\item{control}{Control list passed to \code{\link[stats]{optim}}.}
\item{cell_fit, umi_fit}{Marginal empirical scaling fits used to initialize or transfer joint-law parameters.}
\item{mode}{Joint fitting mode. \code{"transfer"} fixes marginal gap parameters and fits only the global ceiling; \code{"hybrid"} fixes scale parameters and refits the ceiling and amplitudes; \code{"free"} refits all joint parameters.}
\item{n, U}{Cell number and UMI depth values for prediction.}
\item{n0, U0}{Reference cell number and UMI depth for extrapolating the global ceiling from marginal fits.}
\item{I_infinity, C_n, A_n, C_U, A_U, U_offset}{Scaling-law parameters.}
}
\description{
Fit direct nonlinear scaling laws to fixed-rank empirical subspace MI values,
without spike calling or Marchenko--Pastur bulk fitting.  The cell law is
\deqn{I(n)=I_{\infty|U_0}-C_n\log(1+A_n/n).}
The UMI law is
\deqn{I(U)=I_{\infty|n_0}-C_U\log(1+A_U/(U+U_0)).}
The joint law is additive in the two empirical information gaps:
\deqn{I(n,U)=I_\infty-C_n\log(1+A_n/n)-C_U\log(1+A_U/(U+U_0)).}
}
\value{
The fitting functions return an object of class
\code{scscale_empirical_scaling_fit}, with fitted coefficients, residuals,
fit metrics, and \code{predict()} support. \code{scscale_empirical_global_ceiling}
returns the two marginally implied global ceilings and their average.
}
\keyword{models}
