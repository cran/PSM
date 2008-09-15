\name{PSM.plot}
\alias{PSM.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Basic plots of data and output}
\description{
  Create basic plots of data and state estimates in PSM. 
}

\usage{
  PSM.plot(Data, Smooth = NULL, indiv = NULL, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
    Data list, see description in \code{\link{PSM.estimate}}.
  }
  \item{Smooth}{
    Output from \code{\link{PSM.smooth}}.
  }
  \item{indiv}{
    A vector of integers with which individuals to include.
  }
  \item{type}{
    A vector of strings listing the types of plots to create. The
    possibilities are: 'X', 'longX', 'Xp', 'Xf', 'Xs', 'Yp', 'Ys',
    'Yp.Y', 'Ys.Y', 'res', 'acf', 'eta'. If a
    string is preceeded by 'logx.', 'logy.' or 'logx.logy.' the
    corresponding axis is shown on log-scale. An example is: type=c('Xs','logy.Ys.Y','res','eta')
  }
}

\value{
  None (invisible \code{NULL}).
}

\note{
  For further details please also read the package vignette pdf-document
  by writing vignette("PSM") in R.
}

\references{Please visit \url{http://www.imm.dtu.dk/psm} or refer to the
  help page for \code{\link{PSM}}. } 

\author{ Stig B. Mortensen and Søren Klim }

\seealso{   \code{\link{PSM}}, \code{\link{PSM.estimate}}, \code{\link{PSM.smooth}},
  \code{\link{PSM.simulate}}, \code{\link{PSM.template}}
}
\examples{
cat("\nExamples are included in the package vignette.\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
