################################################################################
#' Extract data from the ncdf file 
#' @param .nc  Open netCDF file (results from after open it with `ncdf4::nc_open`)
#' @param .dname Name of the variable to be extracted  
#' @return Array with the variable data 

get_ncdf_data <- function(.nc, .dname, .start = NA, .count = NA){
  array <- ncdf4::ncvar_get(.nc, .dname, .start, .count)
  dlname <- ncdf4::ncatt_get(.nc, .dname,"standard_name")
  dunits <- ncdf4::ncatt_get(.nc, .dname,"units")
  fillvalue <- ncdf4::ncatt_get(.nc, .dname,"_FillValue")
  # replace netCDF fill values with NA's
  array[array == fillvalue$value] <- NA
  return(array)
}

################################################################################
#' Convert wind speed from "x_wind_" and "y_wind_" components, 
#' to wind magnitude (m/s) and direction (0 - 360 degrees) 
#' @param .r SpatRaster with "x_wind_" and "y_wind_" wind speed (e.g., NORA3 data) 
#' @param .height Height to get the data (i.e., 20, 50, 100, 250, 750)
#' @return SpatRaster with wind speed magnitude (_mag) and direction (_dir)

get_wind_height <- function(.r, .height = 100){
  
  # Velocity vector (y, y)
  ux = terra::subset(.r, paste0("x_wind_", .height, "m"))
  uy = terra::subset(.r, paste0("y_wind_", .height, "m"))
  
  # Calculate magnityd (mag) and direction (dir)
  u_mag = sqrt(ux^2 + uy^2)
  names(u_mag) <- "magnitude"
  u_dir = terra::atan2(y = uy, x = ux) * 180/pi
  names(u_dir) <- "direction"
  
  # Generate raster 
  u = c(u_mag, u_dir) 
  
  # Output as one raster
  return(u)
  
}

################################################################################
#' Function for downloading wind data from NORA3 (https://thredds.met.no/thredds/projects/nora3.html)
#' The function extract the data from the file "fc<YYYYMMDDHH>_<leadtime>_fp.nc"
#' @param .year Character indicating the year of the data (from "1959" to "2024")
#' @param .month Character indicating the month ("01", "02","03", ...., "12")
#' @param .day Character indicating the day ("01", "02", "03", ....)
#' @param .hour_group Character indicating the six hours period (four options: "00", "06", "12", "18")
#' @param .lead_time Character indicating lead time (six options: "004", "005", "006", "007", "008", "009")
#' @results Returns a single raster file (SpatRaster) with the magnitude ("_mag") and direction ("_dir") of the wind 
#' at different heights (z = 10, 20, 50, 100, 250, 750)  

get_wind_z <- function(.year,
                       .month,
                       .day,
                       .hour_group,
                       .lead_time){
  
  # URL of the data
  nora3_url <- paste0("https://thredds.met.no/thredds/dodsC/nora3/",
                      .year,
                      "/",
                      .month,
                      "/",
                      .day, 
                      "/",
                      .hour_group,
                      "/fc", 
                      .year,
                      .month, 
                      .day,
                      .hour_group,
                      "_",
                      .lead_time,
                      "_fp.nc")
  
  # Open the netCDF file
  ncin <- ncdf4::nc_open(nora3_url)
  
  # Get values of the variables
  # Coordinate  (x, y)
  # Time (time)
  # wind speed at 10 m above ground - height4 ("wind_speed")
  # Direction at 10 m above ground - height4 (wind_direction)
  # x_wind_z[x,y,height2,time] ("x_wind_z")
  # y_wind_z[x,y,height2,time] ("y_wind_z")
  ff <- function(.x) { get_ncdf_data(.nc = ncin, .dname = .x) } 
  data_list <- c("x", "y", "time", "wind_speed", "wind_direction", "x_wind_z", "y_wind_z") |> 
    purrr::map(ff)
  names(data_list) <- c("lon", "lat", "time", "wind_speed", "wind_direction", "x_wind_z", "y_wind_z")
  
  # Close link
  ncdf4::nc_close(ncin)
  
  # Create dataframe with values
  df <- expand.grid(data_list$lon ,data_list$lat) |> 
    tibble::as_tibble() |> 
    dplyr::rename_with(~ c("lon", "lat"), 1:2) |> 
    # Add wind speed and direction at 10 m
    dplyr::mutate(wind10_mag = as.vector(data_list$wind_speed),
                  wind10_dir = as.vector(data_list$wind_direction)) |> 
    # Add wind speed at h2 
    dplyr::bind_cols(
      matrix(data_list$x_wind_z, ncol = 6) |> 
        tibble::as_tibble(.name_repair = "unique") |> 
        dplyr::rename_with(~c("x_wind_20m",
                              "x_wind_50m",
                              "x_wind_100m", 
                              "x_wind_250m", 
                              "x_wind_500m",
                              "x_wind_750m"), 
                           .cols = everything())
    ) |> 
    dplyr::bind_cols(
      matrix(data_list$y_wind_z, ncol = 6) |> 
        tibble::as_tibble(.name_repair = "unique") |> 
        dplyr::rename_with(~c("y_wind_20m", 
                              "y_wind_50m", 
                              "y_wind_100m", 
                              "y_wind_250m", 
                              "y_wind_500m", 
                              "y_wind_750m"), 
                           .cols = everything())
    )
  
  # Create Raster with all data
  r_crs <- "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs"
  r <- tidyterra::as_spatraster(df, crs = r_crs, digits = 4)
  # Add time 
  terra::time(r) <- rep(lubridate::as_datetime(data_list$time, tz = "UTC"), times = length(names(r)))
  
  # get wind speed and direction at all heights
  ff <- function(.height) { get_wind_height(r, .height) }
  u_height_list <- c(20, 50, 100, 250, 750)
  u_height <- purrr::map(u_height_list, ff)
  names(u_height) <- c("wind20", "wind50", "wind100", "wind250", "wind750")
  u_height <- u_height |>
    terra::rast() |> 
    dplyr::rename_with( ~ gsub("_1", "_mag", .x, fixed = TRUE)) |> 
    dplyr::rename_with( ~ gsub("_2", "_dir", .x, fixed = TRUE))
  
  # Return SpatRaster object
  rr <- c(tidyterra::select(r, wind10_mag,  wind10_dir), u_height)
  
  return(rr)
  
} 


#' @example
ws <- get_wind_z(.year = "2018",
                 .month = "11",
                 .day = "08",
                 .hour_group = "00",
                 .lead_time = "004")
ws







################################################################################
#' Function for downloading the wind profile at a special location
#' @param .year Character indicating the year of the data (from "1959" to "2024")
#' @param .month Character indicating the month ("01", "02","03", ...., "12")
#' @param .day Character indicating the day ("01", "02", "03", ....)
#' @param .hour_group Character indicating the six hours period (four options: "00", "06", "12", "18")
#' @param .lead_time Character indicating lead time (six options: "004", "005", "006", "007", "008", "009")
#' @param .long Longitude coordinate (WGS84)
#' @param .lat Latitude coordinate (WGS84)
#' @results Data frame with wind profile (z = 10, 20, 50, 100, 250, 750)     

get_wind_profile <- function(.year,
                             .month,
                             .day,
                             .hour_group,
                             .lead_time,
                             .long,
                             .lat){
  
  # URL of the data
  nora3_url <- paste0("https://thredds.met.no/thredds/dodsC/nora3/",
                      .year,
                      "/",
                      .month,
                      "/",
                      .day, 
                      "/",
                      .hour_group,
                      "/fc", 
                      .year,
                      .month, 
                      .day,
                      .hour_group,
                      "_",
                      .lead_time,
                      "_fp.nc")
  
  # Open the netCDF file
  ncin <- ncdf4::nc_open(nora3_url)
  
  # Point as SpatVector
  point <- cbind(.long, .lat) |> 
    terra::vect(crs="+proj=longlat +datum=WGS84") |> 
    terra::project("+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs")
  
  # Get coordinate  variables (subset around point - 6000 m) 
  point_xy <- tibble::as_tibble(terra::geom(point))
  LonIdx <- which(ncin$dim$x$vals > (point_xy$x - 6000) &
                    ncin$dim$x$vals < (point_xy$x + 6000) 
  )
  LatIdx <- which(ncin$dim$y$vals > (point_xy$y - 6000) &
                    ncin$dim$y$vals < (point_xy$y + 6000)
  )
  lon <- ncin$dim$x$val[LonIdx]
  lat <- ncin$dim$y$val[LatIdx]
  
  # Get time
  time <- ncdf4::ncvar_get(ncin,"time")
  
  # Get values of the variables
  # Coordinate  (x, y)
  # Time (time)
  # wind speed at 10 m above ground - height4 ("wind_speed")
  # Direction at 10 m above ground - height4 (wind_direction)
  # x_wind_z[x,y,height2,time] height2: 20  50 100 250 500 750 m above ground
  # y_wind_z[x,y,height2,time] ("y_wind_z")
  ff <- function(.x) { 
    get_ncdf_data(.nc = ncin,
                  .dname = .x,
                  .start = c(LonIdx[1], LatIdx[1], 1, 1),
                  .count = c(length(LonIdx), length(LatIdx), 6, 1)
                  )
  } 
  data_list <- c("x_wind_z", "y_wind_z") |> 
    purrr::map(ff)
  names(data_list) <- c("x_wind_z", "y_wind_z")
  
  # Close link
  ncdf4::nc_close(ncin)
  
  # Create dataframe with values
  df <- expand.grid(lon, lat) |> 
    tibble::as_tibble() |> 
    dplyr::rename_with(~ c("lon", "lat"), 1:2) |> 
    # # Add wind speed and direction at 10 m
    # dplyr::mutate(wind10_mag = as.vector(data_list$wind_speed),
    #               wind10_dir = as.vector(data_list$wind_direction)) |> 
    # Add wind speed at h2 
    dplyr::bind_cols(
      matrix(data_list$x_wind_z, ncol = 6) |> 
        tibble::as_tibble(.name_repair = "unique") |> 
        dplyr::rename_with(~c("x_wind_20m",
                              "x_wind_50m",
                              "x_wind_100m", 
                              "x_wind_250m", 
                              "x_wind_500m",
                              "x_wind_750m"), 
                           .cols = everything())
    ) |> 
    dplyr::bind_cols(
      matrix(data_list$y_wind_z, ncol = 6) |> 
        tibble::as_tibble(.name_repair = "unique") |> 
        dplyr::rename_with(~c("y_wind_20m", 
                              "y_wind_50m", 
                              "y_wind_100m", 
                              "y_wind_250m", 
                              "y_wind_500m", 
                              "y_wind_750m"), 
                           .cols = everything())
    )
  
  # Create Raster with all data
  r_crs <- "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs"
  r <- tidyterra::as_spatraster(df, crs = r_crs, digits = 4)
  # Add time 
  terra::time(r) <- rep(lubridate::as_datetime(time, tz = "UTC"), times = length(names(r)))
  
  # get wind speed and direction at all heights
  ff <- function(.height) { get_wind_height(r, .height) }
  u_height_list <- c(20, 50, 100, 250, 750)
  u_height <- purrr::map(u_height_list, ff)
  names(u_height) <- c("wind20", "wind50", "wind100", "wind250", "wind750")
  u_height <- u_height |>
    terra::rast() |> 
    dplyr::rename_with( ~ gsub("_1", "_mag", .x, fixed = TRUE)) |> 
    dplyr::rename_with( ~ gsub("_2", "_dir", .x, fixed = TRUE))
  
  # Extract values
  u_profile <- terra::extract(u_height, point) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(
      # long = .long,
      # lat = .lat,
      time = lubridate::as_datetime(time, tz = "UTC")
    ) |> 
    relocate(c(
      # long, lat,
      time), .after = "ID")
  
  return(u_profile)
  
} 

#' @example
wp <- get_wind_profile(.year = "2018",
                       .month = "11",
                       .day = "08",
                       .hour_group = "00",
                       .lead_time = "004",
                       .long = 4.899842,
                       .lat = 59.34659)
wp




########################################################################
#' Interpolate wind speed profile
#' @param .x data frame with wind speed at different heights (20, 50, 100, 250, 750))

get_inter_wind_profile <- function(x){
  
  # Vector with alpha depending on height
  alpha <- rep(NA, 5)
  for(i in seq_along(alpha)) {
    alpha[i] = as.numeric( (log(x[i, 2]/x[i+1, 2]) / log(x[i, 1]/x[i+1, 1])) )
  }
  
  # Create data frame for interpolation (steps = 1) 
  wind_inter <- tibble(height_m = seq(10, 750, 1)) |> 
    mutate(wind_m_s = NA) |> 
    mutate(alpha = case_when(
      height_m <= 20 ~ alpha[1],
      height_m >  20 &  height_m <=  50 ~ alpha[2],
      height_m >  50 &  height_m <= 100 ~ alpha[3],
      height_m > 100 &  height_m <= 250 ~ alpha[4],
      height_m > 250 &  height_m <= 750 ~ alpha[5]
    ))  
  
  # Interpolation (based on u [m/s] at 10 m)
  wind_inter$wind_m_s[1] = x[x["height_m"] == 10, ]$wind_m_s
  for(i in 2:length(wind_inter$height_m)) {
    wind_inter$wind_m_s[i] = wind_inter$wind_m_s[i-1] * (wind_inter$height_m[i] / wind_inter$height_m[i-1])^wind_inter$alpha[i]
  }
  
  return(wind_inter)
  
}












################################################################################
#' Function for getting NORA3 data at 100 metres in a point
#' @param .year Year to get the data
#' @param .month
#' @param .day 
#' @param .hour_group
#' @param .lead_time
#' @param .long
#' @param .lat

get_wind_100m_point <- function(.year,
                                .month,
                                .day,
                                .hour_group,
                                .lead_time,
                                .long,
                                .lat){
  
  # URL of the data
  nora3_url <- paste0("https://thredds.met.no/thredds/dodsC/nora3/",
                      .year,
                      "/",
                      .month,
                      "/",
                      .day, 
                      "/",
                      .hour_group,
                      "/fc", 
                      .year,
                      .month, 
                      .day,
                      .hour_group,
                      "_",
                      .lead_time,
                      "_fp.nc")
  
  # Open the netCDF file
  ncin <- ncdf4::nc_open(nora3_url)
  
  # Point as SpatVector
  point <- cbind(.long, .lat) |> 
    vect(crs="+proj=longlat +datum=WGS84") |> 
    project("+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs")
  
  # Get coordinate  variables (subset around point - 10000 m) 
  point_xy <- as_tibble(geom(point))
  LonIdx <- which(ncin$dim$x$vals > (point_xy$x - 10000) &
                    ncin$dim$x$vals < (point_xy$x + 10000) 
  )
  LatIdx <- which(ncin$dim$y$vals > (point_xy$y - 10000) &
                    ncin$dim$y$vals < (point_xy$y + 10000)
  )
  lon <- ncin$dim$x$val[LonIdx]
  lat <- ncin$dim$y$val[LatIdx]
  
  # Get time
  time <- ncdf4::ncvar_get(ncin,"time")
  
  # x_wind_z[x,y,height2,time] -- height2: 20  50 [100] 250 500 750 m above ground
  dname <- "x_wind_z"
  xh2_array <- ncdf4::ncvar_get(ncin,
                                dname,
                                start = c(LonIdx[1], LatIdx[1], 3, 1),
                                count = c(length(LonIdx), length(LatIdx), 1, 1)
  )
  dlname <- ncdf4::ncatt_get(ncin,dname,"standard_name")
  dunits <- ncdf4::ncatt_get(ncin,dname,"units")
  fillvalue <- ncdf4::ncatt_get(ncin,dname,"_FillValue")
  # replace netCDF fill values with NA's
  xh2_array[xh2_array == fillvalue$value] <- NA
  
  # y_wind_z[x,y,height2,time] -- -- height2: 20  50 100 250 500 750 m above ground
  dname <- "y_wind_z"
  yh2_array <- ncdf4::ncvar_get(ncin,
                                dname,
                                start = c(LonIdx[1], LatIdx[1], 3, 1),
                                count = c(length(LonIdx), length(LatIdx), 1, 1)
  )
  dlname <- ncdf4::ncatt_get(ncin,dname,"standard_name")
  dunits <- ncdf4::ncatt_get(ncin,dname,"units")
  fillvalue <- ncdf4::ncatt_get(ncin,dname,"_FillValue")
  # replace netCDF fill values with NA's
  yh2_array[yh2_array == fillvalue$value] <- NA
  
  # create dataframe with values
  df <- expand.grid(lon,lat) |> 
    as_tibble() |> 
    dplyr::rename_with(~ c("long", "lat"), 1:2) |> 
    # Add wind speed at 100 m 
    dplyr::mutate(
      # x_wind_z
      x_wind_100m = as.vector(xh2_array),
      # y_wind_z 
      y_wind_100m = as.vector(yh2_array),
    )
  
  # Create Raster with all data
  r_crs <- "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs"
  r <- tidyterra::as_spatraster(df, crs = r_crs, digits = 4) |> 
    # Calculate magnitude and direction
    mutate(u_mag = sqrt(x_wind_100m^2 + y_wind_100m^2),
           u_dir = terra::atan2(y = y_wind_100m, x = x_wind_100m) * 180/pi) |> 
    select(u_mag, u_dir)
  
  # Extract values
  p_wind_100m <- terra::extract(r, point) |> 
    as_tibble() |> 
    mutate(time = as_datetime(time, tz = "UTC"))
  
  ncdf4::nc_close(ncin)
  
  return(p_wind_100m)
  
}



