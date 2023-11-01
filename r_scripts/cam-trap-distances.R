library(raster)
library(sf)
library(dplyr)

records <- read.csv("raw_data/raw_data_gams_pa.csv")

# new column for each grid
records$sess <- paste0(records$data_source, "_", records$block)
# remove zoi
records <- filter(records, !(data_source == "zoi"))

# get number of times grid was surveyed
surveys <- records %>%
  group_by(sess) %>%
  mutate(surveys = n_distinct(as.character(year)),
         min_year = min(year), 
         max_year = max(year)) %>%
  distinct(sess, .keep_all = TRUE) %>%
  ungroup() %>%
  dplyr::select(region, sess, treatment, surveys, min_year, max_year) %>%
  mutate()

# filter to just sites in each grid
mydata <- distinct(records, station, .keep_all = TRUE)

# make spatial 
mydata <- st_as_sf(mydata, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 3111)


mydata_glenelg_ark_annya   <- filter(mydata, sess == unique(mydata$sess)[1])
mydata_glenelg_ark_cobb    <- filter(mydata, sess == unique(mydata$sess)[2])
mydata_glenelg_ark_hotspur <- filter(mydata, sess == unique(mydata$sess)[3])
mydata_glenelg_ark_lgnp_n  <- filter(mydata, sess == unique(mydata$sess)[4])
mydata_glenelg_ark_lgnp_s  <- filter(mydata, sess == unique(mydata$sess)[5])
mydata_glenelg_ark_mtclay  <- filter(mydata, sess == unique(mydata$sess)[6])
mydata_matt_annya          <- filter(mydata, sess == unique(mydata$sess)[7])
mydata_matt_cobb           <- filter(mydata, sess == unique(mydata$sess)[8])
mydata_matt_hotspur        <- filter(mydata, sess == unique(mydata$sess)[9])
mydata_matt_mtclay         <- filter(mydata, sess == unique(mydata$sess)[10])
mydata_matt_north          <- filter(mydata, sess == unique(mydata$sess)[11])
mydata_matt_south          <- filter(mydata, sess == unique(mydata$sess)[12])
mydata_otway_ark_FA1       <- filter(mydata, sess == unique(mydata$sess)[13])
mydata_otway_ark_FA2       <- filter(mydata, sess == unique(mydata$sess)[14])
mydata_otway_ark_FA3       <- filter(mydata, sess == unique(mydata$sess)[15])
mydata_otway_ark_FA4       <- filter(mydata, sess == unique(mydata$sess)[16])


#distmat <- st_distance(mydata)
#dim(distmat)

mindist <- function(x) {
  mindis <- vector(mode = "numeric", length = nrow(x))
  for(i in 1:nrow(x)){
    mindis[i] <- min(x[i, -i])
  }
  return(mindis)
}


# function to add columns to dataframe (by session - group_by() wasn't working!)
dist_function <- function(data){
  data <- mutate(data, 
                 sites = n(),
                 mindist = mindist(st_distance(data)),
                 maxdist = apply(st_distance(data), 1, max),
                 min_mindist = min(mindist),
                 mean_mindist = mean(mindist),
                 max_mindist = max(mindist))
  data <- slice(data, 1)
  data$geometry <- NULL
  data <- dplyr::select(data, sess, data_source, block, sites, min_mindist, mean_mindist, max_mindist)
  return(data)
}


cam_dist <- bind_rows(dist_function(mydata_glenelg_ark_annya), dist_function(mydata_glenelg_ark_cobb)) %>%
  bind_rows(., dist_function(mydata_glenelg_ark_hotspur )) %>%
  bind_rows(., dist_function(mydata_glenelg_ark_lgnp_n  )) %>%
  bind_rows(., dist_function(mydata_glenelg_ark_lgnp_s  )) %>%
  bind_rows(., dist_function(mydata_glenelg_ark_mtclay  )) %>%
  bind_rows(., dist_function(mydata_matt_annya          )) %>%
  bind_rows(., dist_function(mydata_matt_cobb           )) %>%
  bind_rows(., dist_function(mydata_matt_hotspur        )) %>%
  bind_rows(., dist_function(mydata_matt_mtclay         )) %>%
  bind_rows(., dist_function(mydata_matt_north          )) %>%
  bind_rows(., dist_function(mydata_matt_south          )) %>%
  bind_rows(., dist_function(mydata_otway_ark_FA1       )) %>%
  bind_rows(., dist_function(mydata_otway_ark_FA2       )) %>%
  bind_rows(., dist_function(mydata_otway_ark_FA3       )) %>%
  bind_rows(., dist_function(mydata_otway_ark_FA4       ))


# add in survey session
cam_dist <- left_join(cam_dist, surveys)

# reorder dataframe
cam_dist <- relocate(cam_dist, region, sess, data_source, block, treatment, surveys, min_year, max_year, sites)

# make pretty
cam_dist$region <- if_else(cam_dist$region == "otways", "Otways", "Glenelg")
cam_dist$data_source <- if_else(cam_dist$data_source == "glenelg_ark", "Glenelg Ark", cam_dist$data_source)
cam_dist$data_source <- if_else(cam_dist$data_source == "otway_ark", "Otway Ark", cam_dist$data_source)
cam_dist$data_source <- if_else(cam_dist$data_source == "matt", "MWR PhD", cam_dist$data_source)

cam_dist$block <- if_else(cam_dist$block == "annya", "Annya", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "cobb", "Cobboboonee", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "hotspur", "Hotspur", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "mtclay", "Mt Clay", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "north", "North", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "south", "South", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "lgnp_n", "LGNP N", cam_dist$block)
cam_dist$block <- if_else(cam_dist$block == "lgnp_s", "LGNP S", cam_dist$block)

# save
write.csv(cam_dist, "derived_data/survey_summary_distances.csv")


# END