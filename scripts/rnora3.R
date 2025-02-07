
# Aggregated data
nora3_agg_url <- "https://thredds.met.no/thredds/dodsC/nora3agg/nora3hindcastaggregated.ncml"
ncin <- ncdf4::nc_open(nora3_agg_url)


# function for getting NORA3 data at 100 m in a point ----
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
  
  # y_wind_z[x,y,height2,time] -- -- height2: 20  50 [100] 250 500 750 m above ground
  dname <- "y_wind_z"
  yh2_array <- ncdf4::ncvar_get(ncin,
                                dname,
                                start = c(LonIdx[1], LatIdx[1], 3, 1),
                                count = c(length(LonIdx), length(LatIdx), 1,1 )
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

# Use the function for getting for 30 days at the same point: 2018-01-01 to 2018-01-31
df <- tibble(
  time = seq.Date(from = as.Date("2018-01-01"),
                  to = as.Date("2018-01-31"),
                  by = "day"),
  data = list(
    tibble(hour_group = rep(c("00", "06", "12", "18"), each = 6),
           lead_time = rep(c("004", "005", "006", "007", "008", "009"), time = 4) 
    )
  )
) |> 
  separate(time, c("year", "month", "day")) |> 
  unnest(data) |> 
  mutate(long = 4.899842,
         lat = 59.346588)

wind_100m <- df |> 
  rowwise() |> 
  mutate(wing_100m = get_wind_100m_point(.year = year,
                                         .month = month,
                                         .day = day,
                                         .hour_group = hour_group,
                                         .lead_time = lead_time,
                                         .long = long,
                                         .lat = lat)) |>  
  unnest(wing_100m) 


saveRDS(wind_100m, "data/wind_100m.rds")


