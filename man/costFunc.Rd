% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/costFuncR6.R
\docType{class}
\name{costFunc}
\alias{costFunc}
\title{\code{costFunc} class}
\description{
An \code{R6} class specifying a cost function
}
\details{
Creates an instance of \code{costFunc} \code{R6} class, used in initialisation of change-point detection modules. Currently
supports the following cost functions:
\itemize{
\item \strong{L1 cost function}:
\deqn{c_{L_1}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \tilde{y}_{(a+1)...b} \|_1}
where \eqn{\tilde{y}_{(a+1)...b}} is the coordinate-wise median of the segment. If \eqn{a \ge b - 1}, return 0.
\item \strong{L2 cost function}:
\deqn{c_{L_2}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \bar{y}_{(a+1)...b} \|_2^2}
where \eqn{\bar{y}_{(a+1)...b}} is the empirical mean of the segment. If \eqn{a \ge b - 1}, return 0.
\item \strong{SIGMA cost function}:
\deqn{c_{\sum}(y_{(a+1)...b}) := (b - a)\log \det \hat{\Sigma}_{(a+1)...b}} where \eqn{\hat{\Sigma}_{(a+1)...b}} is
the empirical covariance matrix of the segment without Bessel's correction. Here, if \code{addSmallDiag = TRUE}, a small
bias \code{epsilon} is added to the diagonal of estimated covariance matrices to improve numerical stability. \cr
\cr
By default, \code{addSmallDiag = TRUE} and \code{epsilon = 1e-6}. In case \code{addSmallDiag = TRUE}, if the computed determinant of covariance matrix is either 0 (singular)
or smaller than \code{p*log(epsilon)} - the lower bound, return \code{(b - a)*p*log(epsilon)}, otherwise, output an error message.
\item \strong{VAR(r) cost function}:
\deqn{c_{\mathrm{VAR}}(y_{(a+1)...b}) := \sum_{t = a+r+1}^{b} \left\| y_t - \sum_{j=1}^r \hat A_j y_{t-j} \right\|_2^2}
where \eqn{\hat A_j} are the estimated VAR coefficients, commonly estimated via the OLS criterion. If system is singular,
\eqn{a-b < p*r+1} (i.e., not enough observations), or \eqn{a \ge n-p} (where \code{n} is the time series length), return 0.
\item \strong{"LinearL2"} for piecewise linear regression process with \strong{constant noise variance}
\deqn{c_{\text{LinearL2}}(y_{(a+1):b}) := \sum_{t=a+1}^b \| y_t - X_t \hat{\beta} \|_2^2} where \eqn{\hat{\beta}} are OLS estimates on segment \eqn{(a+1):b}. If segment is shorter than the minimum number of
points needed for OLS, return 0.
}

If active binding \verb{$costFunc} Can be accessed or modified by assigning to \verb{$costFunc}, the default parameters will be used.
}
\section{Methods}{

\describe{
\item{\code{$new()}}{Initialises a \code{costFunc} object.}
\item{\code{$describe()}}{Describes the \code{costFunc} object.}
\item{\code{$clone()}}{Clones the \code{costFunc} object.}
}
}

\author{
Minh Long Nguyen \email{edelweiss611428@gmail.com}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{costFunc}}{Character. Cost function. Can be accessed or modified via \verb{$costFunc}. If \code{costFunc} is modified
and required parameters are missing, the default parameters are used.}

\item{\code{pVAR}}{Integer. Vector autoregressive order. Can be accessed or modified via \verb{$pVAR}.}

\item{\code{addSmallDiag}}{Logical. Whether to add a bias value to the diagonal of estimated covariance matrices to stabilise matrix operations. Can be accessed or modified via \verb{$addSmallDiag}.}

\item{\code{epsilon}}{Double. A bias value added to the diagonal of estimated covariance matrices to stabilise matrix operations. Can be accessed or modified via \verb{$epsilon}.}

\item{\code{intercept}}{Logical. Whether to include the intercept in regression problems. Can be accessed or modified via \verb{$intercept}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-costFunc-new}{\code{costFunc$new()}}
\item \href{#method-costFunc-pass}{\code{costFunc$pass()}}
\item \href{#method-costFunc-clone}{\code{costFunc$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-costFunc-new"></a>}}
\if{latex}{\out{\hypertarget{method-costFunc-new}{}}}
\subsection{Method \code{new()}}{
Initialises a \code{costFunc} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{costFunc$new(costFunc, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{costFunc}}{Character. Cost function. Supported values include \code{"L2"}, \code{"VAR"},
and \code{"SIGMA"}. Default: \code{L2}.}

\item{\code{...}}{Optional named parameters required by specific cost functions. \cr
If any required parameters are missing or null, default values will be used.

For \code{"L1"} and \code{"L2"}, there is no extra parameter.

For \code{"SIGMA"}, supported parameters are:
\describe{
\item{\code{addSmallDiag}}{Logical. If \code{TRUE}, add a small value to the diagonal of estimated covariance matrices
to stabilise matrix operations. Default: \code{TRUE}.}
\item{\code{epsilon}}{Double. If \code{addSmallDiag = TRUE}, a small positive value added to the diagonal of estimated covariance matrices to stabilise
matrix operations. Default: \code{1e-6}.}
}

For \code{"VAR"}, \code{pVAR} is required:
\describe{
\item{\code{pVAR}}{Integer. Vector autoregressive order. Must be a positive integer. Default: \code{1L}.}
}

For \code{"LinearL2"}, \code{intercept} is required:
\describe{
\item{\code{intercept}}{Logical. Whether to include the intercept in regression problems. Default: \code{TRUE}.}
}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-costFunc-pass"></a>}}
\if{latex}{\out{\hypertarget{method-costFunc-pass}{}}}
\subsection{Method \code{pass()}}{
Returns a list of configuration parameters to initialise \code{detection} modules.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{costFunc$pass()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-costFunc-clone"></a>}}
\if{latex}{\out{\hypertarget{method-costFunc-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{costFunc$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
