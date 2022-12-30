library(ctmm)
library(tidyverse)
library(amt)
library(sf)
library(lubridate)
library(ggplot2)
library(viridis)
library(hrbrthemes)

## read data and change column names according to movebank format --------------------------------
animal_582 <- read_csv("./data/GPS_ID44582_062722_145829.csv") %>% 
  rename (individual.local.identifier = 'Device ID',
          timestamp = 'Date & Time [Local]',
          location.long = Longitude,
          location.lat = Latitude) %>%
  select(individual.local.identifier, timestamp, location.long, location.lat) %>%
  filter(timestamp >  '2018-4-19' & timestamp < '2018-5-18') %>% # collection started April 19, 2018
  mutate(individual.local.identifier = "pig_582",
         timestamp = force_tz(timestamp, "Asia/Makassar")) ### <<<------------- Dave check time zone

animal_583 <- read_csv("./data/GPS_ID44583_062722_145140.csv") %>% 
  rename (individual.local.identifier = 'Device ID',
          timestamp = 'Date & Time [Local]',
          location.long = Longitude,
          location.lat = Latitude) %>%
  select(individual.local.identifier, timestamp, location.long, location.lat) %>%
  filter(timestamp >  '2022-4-8' & timestamp <  '2022-5-22') %>%
  mutate(individual.local.identifier = "pig_583",
         timestamp = force_tz(timestamp, "Asia/Makassar")) ### <<<------------- Dave check time zone

# # of days monitored ---------------------------------------------------------------
difftime(range(animal_582$timestamp)[2], range(animal_582$timestamp)[1]) # 28 days 

difftime(range(animal_583$timestamp)[2], range(animal_583$timestamp)[1]) # 39 days


## variograms - visualizing movement autocorrelation --------------------------------
# The variogram represents the average square distance traveled (vertical axis) within some time lag (horizontal axis).
# dt <- c(5, 60) %#% "minute" # account for sampling interval changes 
# SVF <- variogram(animal_582, dt = dt)
# level <- c(0.5,0.95) # 50% and 95% CIs
# xlim <- c(0, 30 %#% "day") # 0-12 hour window
# plot(SVF,xlim=xlim,level=level)
# variogram.fit(SVF)

############################################################################################
## HR estimate -----------------------------------------------------------------------------
############################################################################################

# ------------------------------------------------------------------------------------------------
# animal 582 -------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

animal = animal_582  %>% 
  as.telemetry(., timezone = 'Asia/Makassar') # convert to telementry object    ### <<<------------- Dave check time zone
# plot(animal_582)

GUESS <- ctmm.guess (animal, interactive = FALSE)
M.OUF <- ctmm.fit(animal, GUESS) #continuous-velocity OUF model
# weighted AKDE - account for sample schedule change,  https://doi.org/10.1002/eap.1704
wAKDE <- akde(animal,M.OUF,weights=TRUE) 
#EXT <- extent(wAKDE,level=0.95)
# plot(animal,UD=wAKDE,xlim=EXT$x,ylim=EXT$y)  ## <<---- this give a plot with locations + HR
# title(expression("animal 582 weighted OUF AKDE"["C"]))
wAKDE_582 <- wAKDE
summary(wAKDE_582, level.UD = 0.95, units = F)
summary(wAKDE_582, level.UD = 0.50, units = F)

# ------------------------------------------------------------------------------------------------
# animal 583 -------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

animal = animal_583%>% 
  as.telemetry(., timezone = 'Asia/Makassar') # convert to telementry object    <<< --------- Dave check time zone
# plot(animal_583)

GUESS <- ctmm.guess (animal, interactive = FALSE)
M.OUF <- ctmm.fit(animal, GUESS) #continuous-velocity OUF model
# weighted AKDE - account for sample schedule change,  https://doi.org/10.1002/eap.1704
wAKDE <- akde(animal,M.OUF,weights=TRUE) 
#EXT <- extent(wAKDE,level=0.95)
# plot(animal,UD=wAKDE,xlim=EXT$x,ylim=EXT$y)  ## <<---- this give a plot with locations + HR
# title(expression("animal 583 weighted OUF AKDE"["C"]))
wAKDE_583 <- wAKDE
summary(wAKDE_583, level.UD = 0.95, units = F)
summary(wAKDE_583, level.UD = 0.50, units = F)

# ------------------------------------------------------------------------------------------------
# summary ----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# <<<< For 582 (monitored for 28 days), 95% UD is 0.89 [0.69, 1.11] km2, 50% UD is 0.15 [0.11, 0.18] km2. 
# <<<< For 583 (monitored for 39 days), 95% UD is 10.79 [3.65, 21.73] km2, 50% UD is 2.92 [0.99, 5.88] km2. 

############################################################################################
## Movement patterns  ----------------------------------------------------------------------
############################################################################################

# ------------------------------------------------------------------------------------------------
# animal 582 -------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
animal = animal_582 %>% drop_na() # 35 NA rows dropped
animal.XY <- animal %>%   
  mutate(longitude = animal$location.long, latitude = animal$location.lat) %>%
  st_as_sf(., coords = c("location.long", "location.lat"), crs = 4326) %>%
  st_transform(., 32650) %>%
  cbind(., st_coordinates(.)) %>%
  st_set_geometry(NULL)  # turn into a df with both lat long and x y

traj <- make_track(animal.XY, .x = X, .y = Y, .t = timestamp, 
                   crs = sp::CRS("+init=epsg:32650"))

# summarize sampling rate -----------------------------------------------------------------------------
summarize_sampling_rate(traj) 
# min    q1 median  mean    q3   max    sd     n unit             # n is total relocations
# <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr>
#   1  2.18  4.93   5.07  36.2  59.4 2120.  147.  1135 min  

# now focus on the 5 min period ----------------------------------------------------------------------
# now filter data to only the 5 min ones 
traj_5min <- traj %>% 
  track_resample(rate = minutes(5), tolerance = minutes(5)) %>%
  filter_min_n_burst(.) 
# length(unique(traj_5min$burst_))  # total 47 5-min bursts

# turn tracks into steps to calculate turning angles and movement speeds
traj_5min_steps <- traj_5min %>% 
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = T) %>%#calculate if a location was taken during the day or night 
  mutate(tod_end_ = factor(tod_end_, levels = c("dawn", "day", "dusk", "night")))
# ta_ - turning angles, sl_ step lengths

# summarize day-night # of relocations ---------------------------------------------------------------
summary(traj_5min_steps$tod_end_) 
#   day  dusk night  dawn 
#   433    86   199    76 

# visualize step length (i.e.distance moved every 5 min) by time of day  -----------------------------
traj_5min_steps %>% ggplot( aes(x=tod_end_, y=sl_, fill=tod_end_, color=tod_end_)) +
  geom_violin(width=2.1, size=0.2, alpha=.8) +    
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()  # <----- seems like not too different at different time of day 

# visualize step length (i.e.distance moved every 5 min) by date  ----------------------------------
traj_5min_steps %>% ggplot ( aes(x = t1_, y = sl_)) +
  geom_point() +
  theme_ipsum()   # seems like the animal move more in the later part of the capture 

# visualize heading direction by time of day  -----------------------------------------------------
traj_5min_steps %>%
  ggplot(aes(x = direction_p)) +
  geom_histogram() +
  coord_polar() +
  facet_wrap(facets = traj_5min_steps$tod_end_) +
  theme_ipsum()


# now focus on the 1 h general movement ----------------------------------------------------------------------
traj_1h <- traj %>% 
  track_resample(rate = minutes(60), tolerance = minutes(15)) %>%
  filter_min_n_burst() 
nrow(traj_1h) # total 378 h 

# turn tracks into steps to calculate turning angles and movement speeds
traj_1h_steps <- traj_1h %>% 
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = F) #calculate if a location was taken during the day or night, did not include dawn and dusk this time due to temporal interval
# ta_ - turning angles, sl_ step lengths

# summarize day-night # of relocations ---------------------------------------------------------------
summary(traj_1h_steps$tod_end_) 
#   day  night  
#   170   166   

# visualize step length (i.e.distance moved every 5 min) by time of day  -----------------------------
traj_1h_steps %>% ggplot( aes(x=tod_end_, y=sl_, fill=tod_end_, color=tod_end_)) +
  geom_violin(width=2.1, size=0.2, alpha=.8) +    
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()  # <---- still similar 

# visualize step length (i.e.distance moved every 5 min) by date  ----------------------------------
traj_1h_steps %>% ggplot ( aes(x = t1_, y = sl_)) +
  geom_point() +
  theme_ipsum()   # still seems like the animal move more in the later part of the capture 

# summarize step length by day (daily step length)
traj_1h_steps %>% mutate(date = date(t1_)) %>%
  group_by(date) %>% 
  dplyr::summarise(daily_sl_ = sum(sl_)) %>%
  ggplot ( aes(x = date, y = daily_sl_)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_ipsum()  # daily movement distance confirm that animal 582 move more in the later part of the capture 

traj_1h_steps %>% mutate(date = date(t1_)) %>%
  group_by(date) %>% 
  dplyr::summarise(daily_sl_ = sum(sl_)) %>%
  ggplot ( aes(x = date, y = daily_sl_)) +
  geom_point() +
  geom_smooth(method = loess) +
  theme_ipsum()  # loess regression shows that there might not be clear movement "changing point"

# visualize NSD
traj_1h$nsd <- nsd(traj_1h)
# hourly NSD
traj_1h %>% ggplot ( aes(x = t_, y = nsd)) +
  geom_line(color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=1.5) +
  theme_ipsum()  # same pattern - animal started to move around more after around May 2nd. NSD also showed that the animal slowly disperse away from the original point.
  
# ------------------------------------------------------------------------------------------------
# animal 583 -------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
animal = animal_583 %>% drop_na() # 96 NA rows dropped
animal.XY <- animal %>%   
  mutate(longitude = animal$location.long, latitude = animal$location.lat) %>%
  st_as_sf(., coords = c("location.long", "location.lat"), crs = 4326) %>%
  st_transform(., 32650) %>%
  cbind(., st_coordinates(.)) %>%
  st_set_geometry(NULL)  # turn into a df with both lat long and x y

traj <- make_track(animal.XY, .x = X, .y = Y, .t = timestamp, 
                   crs = sp::CRS("+init=epsg:32650"))

# summarize sampling rate -----------------------------------------------------------------------------
summarize_sampling_rate(traj) 

# min    q1 median  mean    q3   max    sd     n unit 
# <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr>
#   1   2.7  4.92    5.1  69.9  59.4 6642.  404.   812 min  

# now focus on the 5 min period ----------------------------------------------------------------------
# now filter data to only the 5 min ones 
traj_5min <- traj %>% 
  track_resample(rate = minutes(5), tolerance = minutes(5)) %>%
  filter_min_n_burst(.) 
# length(unique(traj_5min$burst_))  # total 47 5-min bursts

# turn tracks into steps to calculate turning angles and movement speeds
traj_5min_steps <- traj_5min %>% 
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = T) %>%#calculate if a location was taken during the day or night 
  mutate(tod_end_ = factor(tod_end_, levels = c("dawn", "day", "dusk", "night"))) 
# ta_ - turning angles, sl_ step lengths

# summarize day-night # of relocations ---------------------------------------------------------------
summary(traj_5min_steps$tod_end_) 
# day  dusk night  dawn 
# 354    55   120    43 

# visualize step length (i.e.distance moved every 5 min) by time of day  -----------------------------
traj_5min_steps %>% ggplot( aes(x=tod_end_, y=sl_, fill=tod_end_, color=tod_end_)) +
  geom_violin(width=2.1, size=0.2, alpha=.8) +    
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()  # <----- seems like a bit more active during the day??

traj_5min_steps2 <- traj_5min %>% 
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = F) 

traj_5min_steps2 %>% ggplot( aes(x=tod_end_, y=sl_, fill=tod_end_, color=tod_end_)) +
  geom_boxplot(width=2.1, size=0.2, alpha=.8) +    
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()  # < ------ the day night differece is not statiscally significant 

# visualize step length (i.e.distance moved every 5 min) by date  ----------------------------------
traj_5min_steps %>% ggplot ( aes(x = t1_, y = sl_)) +
  geom_point() +
  theme_ipsum()   # does not seem to have much activity differences

# visualize heading direction by time of day  -----------------------------------------------------
traj_5min_steps %>%
  ggplot(aes(x = direction_p)) +
  geom_histogram() +
  coord_polar() +
  facet_wrap(facets = traj_5min_steps$tod_end_) +
  theme_ipsum()  # most data from day so hard to say 

# now focus on the 1 h general movement ----------------------------------------------------------------------
traj_1h <- traj %>% 
  track_resample(rate = minutes(60), tolerance = minutes(15)) %>%
  filter_min_n_burst() 
nrow(traj_1h) # total 250 h 

# turn tracks into steps to calculate turning angles and movement speeds
traj_1h_steps <- traj_1h %>% 
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = F) #calculate if a location was taken during the day or night, did not include dawn and dusk this time due to temporal interval
# ta_ - turning angles, sl_ step lengths

# summarize day-night # of relocations ---------------------------------------------------------------
summary(traj_1h_steps$tod_end_) 
#   day  night  
#   149   68   

# visualize step length (i.e.distance moved every 5 min) by time of day  -----------------------------
traj_1h_steps %>% ggplot( aes(x=tod_end_, y=sl_, fill=tod_end_, color=tod_end_)) +
  geom_violin(width=2.1, size=0.2, alpha=.8) +    
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()  # based on day movement it is more obvious that day night activity level differ 

day.sl <- (traj_1h_steps %>% filter (tod_end_ == "day"))$sl_
night.sl <- (traj_1h_steps %>% filter (tod_end_ == "night"))$sl_
t.test(day.sl, night.sl)
# Welch Two Sample t-test
# 
# data:  day.sl and night.sl
# t = 3.2713, df = 211.85, p-value = 0.00125    # and it is statistically significant. -----------------------
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   19.94261 80.41963
# sample estimates:
#   mean of x mean of y 
# 99.22307  49.04196  

# visualize step length (i.e.distance moved every 5 min) by date  ----------------------------------
traj_1h_steps %>% ggplot ( aes(x = t1_, y = sl_)) +
  geom_point() +
  theme_ipsum()   # still seems like the animal move more in the later part of the capture 

# summarize step length by day (daily step length)
traj_1h_steps %>% mutate(date = date(t1_)) %>%
  group_by(date) %>% 
  dplyr::summarise(daily_sl_ = sum(sl_)) %>%
  ggplot ( aes(x = date, y = daily_sl_)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_ipsum()  # 583 actually decrease in its daily movement distance, but just a little bit (the slope is not statistically significant)

traj_1h_steps %>% mutate(date = date(t1_)) %>%
  group_by(date) %>% 
  dplyr::summarise(daily_sl_ = sum(sl_)) %>%
  ggplot ( aes(x = date, y = daily_sl_)) +
  geom_point() +
  geom_smooth(method = loess) +
  theme_ipsum()  # loess regression shows that the daily movement distance actually goes up and down. may be relevant to habitat (will test later). also likely that the increased activity happens when the animal move from one staging area to another.

# visualize NSD
traj_1h$nsd <- nsd(traj_1h)
# hourly NSD
traj_1h %>% ggplot ( aes(x = t_, y = nsd)) +
  geom_line(color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=1.5) +
  theme_ipsum()  # seems like 583 has two general "staging" area. 

############################################################################################
## energy use   ---------------------------------------------------------------------------- 
############################################################################################ 
# <<< ------ accelerometer metadata 


############################################################################################
## HR use patterns  ----------------------------------------------------------------------- 
############################################################################################    
# this will require LC extraction  ### <<<< ------- Dave LC maps?


# for 583, test whether daily activity is dependent on the habitat it is in.

# interesting movement visualization (making maps)