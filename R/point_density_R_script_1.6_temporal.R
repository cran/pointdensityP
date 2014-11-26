#' Arigon dataset
#' 
#' A mock dataset containing meaningless events in a fictional state of Arigon (data overlays Oregon).
#' \itemize{
#'   \item latitude. Latitude of event.
#'   \item longitude. Longitude of event.
#'   \item date. Date of event.
#' }
#' @docType data
#' @keywords datasets
#' @name Arigon
#' @author LTC Steve Henderson and The Department of Systems Engineering at West Point
#' @usage data(Arigon)
#' @format A data frame with 80000 rows and 3 variables 
NULL
#' Houston crime dataset
#' 
#' Lightly cleaned Houston crime; no NA events included and all dates recognized by \code{pointdensity}; data from January 2010 to August 2010 geocoded with Google Maps and courtesy of \pkg{ggmap}
#' @docType data
#' @keywords datasets
#' @name clean_crime
#' @author Houston Police Department, City of Houston
#' @references http://www.houstontx.gov/police/cs/stats2.htm
NULL
#' Point density function for geospatial data  
#' 
#' This function maps a dataset of geospatial points to a a discrete grid and calculates the density and temporal average of the points.  
#'
#' \code{pointdensity} returns a density count and the temporal average for every point in the original list.  The dataframe returned includes four columns:  lat, lon, count, and date_avg.  The "lat" column is the original latitude data; the "lon" column is the original longitude data; the "count" is the density count of the number of points within a radius of radius*grid_size (the neighborhood); and the date_avg column includes the average date of each point in the neighborhood.  Designed specifically for geospatial point processes and originally developed for military applications, this technique applies to any geospatial point process where there is a desire for an explainable measurement of density and maintaining fidelity of the original point locations.  Typical spatial density plotting algorithms, such as kernel density estimation, implement some type of smoothing function that often results in a density value that is difficult to interpret.  Potential applications include analysis of military events,  crime, and real estate transactions.  An example follows with the Arigon data using \pkg{ggmap} (recommended) for visualization: \cr\cr
#' \code{Arigon_density <- pointdensity(df = Arigon, lat_col = "latitude", lon_col = "longitude", date_col = "date", grid_size = 0.01, radius = 2)} \cr
#' \code{map_base <- qmap(location="44.12,-120.83", zoom = 7, darken=0.3)} \cr
#' \code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2, data = Arigon_density) + scale_colour_gradient(low = "green", high = "red")} \cr\cr
#'
#' Here is another example using the crime dataset from \pkg{ggmap}:\cr\cr
#' \code{H_crime <- pointdensity(df = clean_crime, lat_col = "lat", lon_col = "lon", grid_size = 0.01, radius = 4)}\cr
#' \code{map_base <- qmap(location="29.76,-95.42", zoom = 11, darken=0.3)}\cr
#' \code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2, data = H_crime) + scale_colour_gradient(low = "green", high = "red")}
#'
#' @param df Data frame minimally containing latitude and longitude of spatial point data.
#' @param lat_col name of column in \code{df} that contains latitude or vertical dimension of data
#' @param lon_col name of column in \code{df} that contains longitude or horizontal dimension of data
#' @param date_col name of column in \code{df} that contains date associated with the event
#' @param grid_size distance in latitude of the grid that will support discretization of data and density reference
#' @param radius distance in grid steps that represents the local neighborhood where an event adds density 
#' @keywords spatial density
#' @export
#' @author Paul Evangelista \email{paul.evangelista@@usma.edu}
#' @author David Beskow \email{david.beskow@@usma.edu}
#' @references Wand, M. P. (1994). Fast Computation of Multivariate Kernel Estimators. \emph{Journal of Computational and Graphical Statistics}, 3, 433-445.
#' @examples 
#' Arigon_test <- Arigon[1:10000,]
#' Arigon_density <- pointdensity(df = Arigon_test, lat_col = "latitude", 
#' lon_col = "longitude", date_col = "date", grid_size = 0.01, radius = 2)

pointdensity <- function(df, lat_col, lon_col, date_col = NULL, grid_size, radius){
  h<-new.env()  ## hash that will store the density count
  avg_date<-new.env() ## hash that will store the average date
  
  lat_data <- df[,lat_col]
  lat<-lat_data*(1/grid_size)
  lat<-round(lat,0)
  lat<-lat*(grid_size)

  lon_data <- df[,lon_col]
  lon<-lon_data*(1/grid_size)
  lon<-round(lon,0)
  lon<-lon*(grid_size)

  if(is.null(date_col)){
    date <- rep(0,length(lon))
  }
  if(!is.null(date_col)){
    date <- as.Date(df[,date_col])
    date <- as.numeric(date)
  }

  key.vec<-paste(lat,lon,sep="-")
  counter<-0

  pb <- txtProgressBar(title="point density calculation progress", label="0% done", min=0, max=100, initial=0, style = 3)
  data_length <- length(lat)

  for(i in 1:data_length){
    counter <- counter + 1
    if(counter > 99){
      flush.console()
      counter <- 0
    }
    lat.vec<-seq(lat[i]-radius*grid_size,lat[i]+radius*grid_size,grid_size)
    for(lat.temp in lat.vec){
      t<-sqrt(round(((radius*grid_size)^2-(lat.temp-lat[i])^2),8))
      t<-t/cos(lat.temp*2*pi/360)
      t<-t/grid_size
      t<-round(t,0)
      t<-t*grid_size

      lon.vec<-seq(lon[i]-t,lon[i]+t,grid_size)
      for(lon.temp in lon.vec){
        key<-paste(lat.temp, lon.temp, sep="-")
	if(is.null(h[[key]])){
          h[[key]]=1
          avg_date[[key]]=date[i]
        }
        else{
          h[[key]]<-h[[key]]+1
          avg_date[[key]] = avg_date[[key]] + date[i]
	     }
      }
    }
   info <- sprintf("%d%% done", round((i/data_length)*100)) 
   #setWinProgressBar(pb, i/(data_length)*100, label=info)
   setTxtProgressBar(pb, i/(data_length)*100, label=info)
  }
  close(pb)

  count_val <- rep(0,length(key.vec))
  avg_date_val <- rep(0,length(key.vec))

  for(i in 1:length(key.vec)){
    count_val[i] <- h[[key.vec[i]]]	
    avg_date_val[i] <- avg_date[[key.vec[i]]]/count_val[i]
    #cat(i,avg_date_val[i],avg_date[[key.vec[i]]])

    #avg_date_val[i] <- avg_date[[key.vec[i]]]
  }
  final<-data.frame(lat=lat_data,lon=lon_data,count=count_val,dateavg = avg_date_val)
  final<-final[order(final$count),]
  return(final)
}

