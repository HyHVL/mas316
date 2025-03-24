
# Download data

source("scripts/rnora3.r")

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

wind_profile <- df |> 
  rowwise() |> 
  mutate(wind_100m = get_wind_profile(.year = year,
                                         .month = month,
                                         .day = day,
                                         .hour_group = hour_group,
                                         .lead_time = lead_time,
                                         .long = long,
                                         .lat = lat)) |>  
  unnest(wind_100m) 


saveRDS(wind_100m, "data/wind_100m.rds")


######################################################################
# wind profile in one point along time

df <- tibble(
  time = seq.Date(from = as.Date("2018-01-01"),
                  to = as.Date("2018-01-01"),
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

wind_profile <- df |> 
  rowwise() |> 
  mutate(wind_z = get_wind_profile(.year = year,
                                      .month = month,
                                      .day = day,
                                      .hour_group = hour_group,
                                      .lead_time = lead_time,
                                      .long = long,
                                      .lat = lat)) |>  
  unnest(wind_z, names_repair = "unique") |> 
  
  
 wind_profile |> unnest(wind_z)






# ##########################
# nora3_agg_url <- "https://thredds.met.no/thredds/dodsC/nora3agg/nora3hindcastaggregated.ncml"
# ncin <- ncdf4::nc_open(nora3_agg_url)
