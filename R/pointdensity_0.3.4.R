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
#' This function maps a dataset of geospatial points to a regular grid and calculates the density and temporal average of the points.  
#'
#' \code{pointdensity} returns a density count and the temporal average for every point in the original list.  The dataframe returned includes four columns:  lat, lon, count, and date_avg.  The "lat" column is the original latitude data; the "lon" column is the original longitude data; the "count" is the density count of the number of points within a defined radius (the neighborhood); and the date_avg column includes the average date of each point in the neighborhood.  Designed specifically for geospatial point processes and originally developed for military applications, this technique applies to any geospatial point process where there is a desire for an explainable measurement of density and maintaining fidelity of the original point locations.  Typical spatial density plotting algorithms, such as kernel density estimation, implement some type of smoothing function that often results in a density value that is difficult to interpret.  \code{pointdensity} was designed for ease of interpretation.  Potential applications include analysis of military events,  crime, and real estate transactions.  An example follows with the Arigon data using \pkg{ggmap} (recommended) for visualization: \cr\cr
#' \code{Arigon_density <- pointdensity(df = Arigon, lat_col = "latitude", lon_col = "longitude",} \cr
#' \code{date_col = "date", grid_size = 1, radius = 2)} \cr
#' \code{map_base <- qmap(location="44.12,-120.83", zoom = 7, darken=0.3)} \cr
#' \code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2,} \cr
#' \code{data = Arigon_density) + scale_colour_gradient(low = "green", high = "red")} \cr\cr
#'
#' Here is another example using the crime dataset from \pkg{ggmap}:\cr\cr
#' \code{H_crime <- pointdensity(df = clean_crime, lat_col = "lat", lon_col = "lon",} \cr
#' \code{grid_size = 1, radius = 4)}\cr
#' \code{map_base <- qmap(location="29.76,-95.42", zoom = 11, darken=0.3)}\cr
#' \code{map_base + geom_point(aes(x = lon, y = lat, colour = count), shape = 16, size = 2,} \cr 
#' \code{data = H_crime) + scale_colour_gradient(low = "green", high = "red")}
#'
#' @param df Data frame minimally containing latitude and longitude of spatial point data
#' @param lat_col name of column in \code{df} that contains latitude or vertical dimension of data
#' @param lon_col name of column in \code{df} that contains longitude or horizontal dimension of data
#' @param date_col name of column in \code{df} that contains date associated with the event
#' @param grid_size distance in kilometers between the grid lines that will support discretization of data and density reference
#' @param radius distance in kilometers that represents the local neighborhood where an event adds density 
#' @keywords spatial density
#' @import data.table
#' @export
#' @author Paul Evangelista \email{paul.evangelista@@usma.edu}
#' @author David Beskow \email{david.beskow@@usma.edu}
#' @references Wand, M. P. (1994). Fast Computation of Multivariate Kernel Estimators. \emph{Journal of Computational and Graphical Statistics}, 3, 433-445.
#' @examples 
#' Arigon_test <- Arigon[1:1000,]
#' Arigon_density <- pointdensity(df = Arigon_test, lat_col = "latitude", 
#' lon_col = "longitude", date_col = "date", grid_size = 1, radius = 2)

pointdensity <- function(df, lat_col, lon_col, date_col = NULL, grid_size, radius){
  
  #density calculation function - uses vector and matrix operations for fast calculation of local density
  calc_density <- function (lati, loni, grid_size, radius, count, sumDate){
    #vector that contains all latitude grids in neighborhood
    lat.vec <- seq(lati - radius * grid_size, lati + radius * grid_size, grid_size)
    
    #lat.vec.t is the tolerance vector, indicating whether or not a certain lat/lon is within radius
    lat.vec.t <- acos(cos(radius*grid_size)/cos(lat.vec - lati))
    lat.vec.t <- lat.vec.t/cos(lat.vec * 2 * pi/360)
    lat.vec.t <- round(lat.vec.t/grid_size,0)*grid_size
    
    #vector that contains all longitude grids in neighborhood
    lon.vec <- seq(loni - radius * grid_size, loni + radius * grid_size, grid_size)
    
    #matrix that contains lon position of every grid in neighborhood
    lon.mat <- matrix(lon.vec,nrow = length(lon.vec),ncol = length(lon.vec))
    tlon.mat <- abs(lon.mat - loni)
    tlon.mat <- t(tlon.mat)
    
    #apply latitude tolerance, zero-out all points not within neighborhood
    temp <- lat.vec.t - tlon.mat
    temp[temp < (grid_size-(1E-6))] <- 0
    temp[temp > 0] <- 1
    temp2 <- temp*t(lon.mat)
    
    #matrix containing latitude of all grids in neighborhood
    lat.mat <- matrix(rev(lat.vec),nrow = length(lat.vec),ncol = length(lat.vec))
    
    #convert matrices into vectors to support tidy data structure
    lat.vec <- c(lat.mat)
    lon.vec <- c(temp2)
    count.vec <- rep(count,length(lat.vec))
    sumDate.vec <- rep(sumDate,length(lat.vec))
    return.mat <- cbind(lat.vec,lon.vec,count.vec,sumDate.vec)
    
    #eliminate all rows not in neighborhood (0 value for lon)
    row_sub = apply(return.mat, 1, function(x) all(x[1]*x[2] !=0 ))
    return.mat <- return.mat[row_sub,]
    return(return.mat)
  }
  
  lat_c<-lon_c<-cat_r<-count<-sumDate<-ind<-NULL
  
  grid_size <- round(grid_size/111.2, digits = 3)
  rad_km <- radius 			## initial radius measurement in km
  rad_dg <- rad_km/111.2  		## radius as a latitudinal distance
  rad_steps <- round(rad_dg/grid_size)  ## number of steps within grid 
  rad_km <- rad_steps * grid_size * 111.2	## radius rounded to nearest grid step
  cat("\nThe radius was adjusted to ",rad_km,"km in order to accomodate the grid size\n\n") 
  cat("\nThe grid size is ",grid_size," measured in degrees\n\n")
  
  radius <- rad_steps  			## assign to original variable
  
  #round all latitude data to nearest grid
  lat_data <- df[,lat_col]
  lat <- lat_data * (1/grid_size)
  lat <- round(lat, 0)
  lat <- lat * (grid_size)
  lat <- round(lat,3)
  
  #round all longitude data to nearest grid
  lon_data <- df[,lon_col]
  lon <- lon_data * (1/grid_size)
  lon <- round(lon, 0)
  lon <- lon * (grid_size)
  lon <- round(lon,3)
  
  if (is.null(date_col)) {
    date <- rep(0, length(lon))
  }
  if (!is.null(date_col)) {
    date <- as.Date(df[, date_col])
    date <- as.numeric(date)
  }
  
  olat_olon <- paste(lat_data,lon_data,sep = "_")
  rlat_rlon <- paste(lat,lon,sep = "_")

  o_DT <- data.table(lat_c = lat, lon_c = lon, date = date, lat_o = lat_data, lon_o = lon_data, cat_o = olat_olon, cat_r = rlat_rlon, ind = rep(1,length(lon_data)))
  o_DT <- o_DT[order(o_DT$cat_r),]
  
  yy <- o_DT[,list(lat_c = max(lat_c), lon_c = max(lon_c), sumDate=sum(date), count=length(date)), by="cat_r"]
  yy <- yy[order(yy$cat_r),]
  yy_orig <- yy
  yy[,cat_r:=NULL]
  idxs <- seq(1,length(yy$lat_c),1)
  yy <- cbind(yy, ind = idxs)
  
  total_measures <- length(yy$lat_c)*(2*radius +1)*(2*radius +1)
  cat("There are ", length(yy$lat_c)," unique grids that require ", total_measures, " measurements...\n\n")
  
  pb <- txtProgressBar(title="point density calculation progress", label="0% done", min=0, max=100, initial=0, style = 3)
  
  inventory.mat <- matrix(nrow = 1E6, ncol = 4)	
  inventory = 0
  temp_inventory = 0
  inventory.dt <- data.table(lat_c = numeric(0), lon_c = numeric(0), count = numeric(0), sumDate = numeric(0))
  for(i in 1:length(yy$lat_c)){
    #for(i in 1:100000){
    temp3 <- calc_density(lati = yy$lat_c[i], loni = yy$lon_c[i], grid_size = grid_size, radius = radius, count = yy$count[i], sumDate = yy$sumDate[i])
    inventory <- inventory + length(temp3[,1])
    temp_inventory <- temp_inventory +1
    if(temp_inventory == 100){
      info <- sprintf("%d%% done", round((i/length(yy$lat_c))*100)) 
      setTxtProgressBar(pb, i/(length(yy$lat_c))*100, label=info)
      temp_inventory = 0
    }
    if(inventory > 1E6){
      inventory.mat <- inventory.mat[(1:(inventory - length(temp3[,1]))),]
      temp.dt <- data.table(lat_c = inventory.mat[,1], lon_c = inventory.mat[,2], count = inventory.mat[,3], sumDate = inventory.mat[,4])
      temp.dt <- temp.dt[,list(count = sum(count),sumDate = sum(sumDate)), by = "lat_c,lon_c"]
      inventory.dt <- rbind(inventory.dt, temp.dt)
      inventory.mat <- matrix(nrow = 1E6, ncol = 4)	
      inventory.mat[1:length(temp3[,1]),] <- temp3
      inventory = 0
    }		
    else{
      inventory.mat[(inventory - length(temp3[,1])+1):inventory,] <- temp3
    }
  }
  
  #consolidate all density counts
  inventory.mat <- inventory.mat[(1:inventory),]
  temp.dt <- data.table(lat_c = inventory.mat[,1], lon_c = inventory.mat[,2], count = inventory.mat[,3], sumDate = inventory.mat[,4])
  temp.dt <- temp.dt[,list(count = sum(count),sumDate = sum(sumDate)), by = "lat_c,lon_c"]
  inventory.dt <- rbind(inventory.dt, temp.dt)
  
  #append inventory count with original point locations, mark each original row location with ind=1
  yy_ind <- data.table(lat_c = yy$lat_c, lon_c = yy$lon_c, count = rep(0,length(yy$lat_c)), sumDate = rep(0,length(yy$lat_c)), ind = yy$ind)
  #create data table with all density information, calculation row location marked with ind=0
  inventory.dt_ind <- data.table(lat_c = inventory.dt$lat_c, lon_c = inventory.dt$lon_c, count = inventory.dt$count, sumDate = inventory.dt$sumDate, ind = rep(0,length(inventory.dt$lat_c)))
  
  #this is necessary to clean up minor round-off error and ensure all points are on the grid
  inventory2.dt_ind <- data.table(lat_c = round(inventory.dt_ind$lat_c,3), lon_c = round(inventory.dt_ind$lon_c,3), count = inventory.dt_ind$count, sumDate = inventory.dt_ind$sumDate, ind = inventory.dt_ind$ind)
  inventory.dt_ind <- inventory2.dt_ind
  
  inventory.dt <- rbind(inventory.dt_ind, yy_ind)
  inventory.dt<-inventory.dt[order(inventory.dt$lat_c, inventory.dt$lon_c),]
  
  #final consolidation of all density
  inventory3.dt <- inventory.dt[,list(sumDate = sum(sumDate), count = sum(count), ind = sum(ind)), by = "lat_c,lon_c"]
  
  #retain only grid points associated with original points
  final_inventory.dt <- inventory3.dt[inventory3.dt$ind > 0,]
  #calculate the average date of every density location
  date_avg = final_inventory.dt$sumDate / final_inventory.dt$count
  #the final df contains the density for every grid point, but does not contain original locations
  final<-data.table(lat=final_inventory.dt$lat_c,lon=final_inventory.dt$lon_c,count=final_inventory.dt$count,dateavg = date_avg, ind = final_inventory.dt$ind)
  final<-final[order(final$ind),]
  final.1 <- cbind(final,yy_count = yy_orig$count)
  nfinal <- final.1[rep(seq(1,nrow(final.1)), final.1$yy_count)]
  o_DT.check <- cbind(o_DT,nfinal)
  
  #below sorting is necessary only for troubleshooting / viewing grid locations
  
  #final df that contains the original locations, hence the "o"
  #o_final <- data.frame(lat = lat_data, lon = lon_data, count = rep(0,length(lat_data)), dateavg = rep(0,length(lat_data)))
  o_final <- data.table(lat = o_DT$lat_o, lon = o_DT$lon_o, count = nfinal$count, dateavg = nfinal$dateavg)
  
  o_final<-o_final[order(o_final$count),]
  
  cat("done...\n\n")	
  return(o_final)
  }
  