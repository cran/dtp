% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtp.R
\name{dtp}
\alias{dtp}
\title{Dynamic threshold panel function}
\usage{
dtp(f, data, index = c("state", "year"), maxlags, initnum, conf, conf1, conf2,
  graph = TRUE)
}
\arguments{
\item{f}{the formula}

\item{data}{the dataframe}

\item{index}{c=("id","time")}

\item{maxlags}{the instruments maximum lags number}

\item{initnum}{column number of the initial variable}

\item{conf}{the confidential interval}

\item{conf1}{the confidential interval}

\item{conf2}{the confidential interval}

\item{graph}{TRUE to show the threshold plot}
}
\description{
Dynamic threshold panel function
}
\examples{
############################################
# Fit the dynamic threshold panel data
############################################
# Load data
data(Mena)
############################################
# Estimate the model
############################################
reg<-dtp(GDPPC ~ FDI+OPEN|INF|INF,Mena,index=c("pays","ann"),4,2,0.95,0.8,1,graph = TRUE)
summary(reg)


}
