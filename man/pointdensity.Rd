% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pointdensity_0.3.5.R
\name{pointdensity}
\alias{pointdensity}
\title{Point density function for geospatial data}
\usage{
pointdensity(df, lat_col, lon_col, date_col = NULL, grid_size, radius)
}
\arguments{
\item{df}{Data frame minimally containing latitude and longitude of spatial point data}

\item{lat_col}{name of column in \code{df} that contains latitude or vertical dimension of data}

\item{lon_col}{name of column in \code{df} that contains longitude or horizontal dimension of data}

\item{date_col}{name of column in \code{df} that contains date associated with the event}

\item{grid_size}{distance in kilometers between the grid lines that will support discretization of data and density reference}

\item{radius}{distance in kilometers that represents the local neighborhood where an event adds density}
}
\description{
This function maps a dataset of geospatial points to a regular grid and calculates the density and temporal average of the points.
}
\details{
\code{pointdensity} returns a density count and the temporal average for every point in the original list.  The dataframe returned includes four columns:  lat, lon, count, and date_avg.  The "lat" column is the original latitude data; the "lon" column is the original longitude data; the "count" is the density count of the number of points within a defined radius (the neighborhood); and the date_avg column includes the average date of each point in the neighborhood.  Designed specifically for geospatial point processes and originally developed for military applications, this technique applies to any geospatial point process where there is a desire for an explainable measurement of density and maintaining fidelity of the original point locations.  Typical spatial density plotting algorithms, such as kernel density estimation, implement some type of smoothing function that often results in a density value that is difficult to interpret.  \code{pointdensity} was designed for ease of interpretation.  Potential applications include analysis of military events,  crime, and real estate transactions.  An example follows with the Arigon data using \pkg{ggmap} (recommended) for visualization: \cr\cr
\code{Arigon_density <- pointdensity(df = Arigon, lat_col = "latitude", lon_col = "longitude",} \cr
\code{date_col = "date", grid_size = 1, radius = 2)} \cr
\code{map_base <- qmap(location="44.12,-120.83", zoom = 7, darken=0.3)} \cr
\code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2,} \cr
\code{data = Arigon_density) + scale_colour_gradient(low = "green", high = "red")} \cr\cr

Here is another example using the crime dataset from \pkg{ggmap}:\cr\cr
\code{H_crime <- pointdensity(df = clean_crime, lat_col = "lat", lon_col = "lon",} \cr
\code{grid_size = 1, radius = 4)}\cr
\code{map_base <- qmap(location="29.76,-95.42", zoom = 11, darken=0.3)}\cr
\code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2,} \cr 
\code{data = H_crime) + scale_colour_gradient(low = "green", high = "red")}
}
\examples{
Arigon_test <- Arigon[1:1000,]
Arigon_density <- pointdensity(df = Arigon_test, lat_col = "latitude", 
lon_col = "longitude", date_col = "date", grid_size = 1, radius = 2)
}
\references{
Wand, M. P. (1994). Fast Computation of Multivariate Kernel Estimators. \emph{Journal of Computational and Graphical Statistics}, 3, 433-445.
}
\author{
Paul Evangelista \email{paul.evangelista@usma.edu}

David Beskow \email{david.beskow@usma.edu}
}
\keyword{density}
\keyword{spatial}
