library(ctmm)
library(tidyverse)
library(amt)
library(sf)

## read data and change column names according to movebank format --------------------------------
animal_582 <- read_csv("./data/GPS_ID44582_062722_145829.csv") %>% 
  rename (individual.local.identifier = 'Device ID',
          timestamp = 'Date & Time [GMT]',
          location.long = Longitude,
          location.lat = Latitude) %>%
  select(individual.local.identifier, timestamp, location.long, location.lat) %>%
  filter(timestamp >  '2018-4-19' & timestamp < '2018-5-18') %>% # collection started April 19, 2018
  mutate(individual.local.identifier = "pig_582")

animal_583 <- read_csv("./data/GPS_ID44583_062722_145140.csv") %>% 
  rename (individual.local.identifier = 'Device ID',
          timestamp = 'Date & Time [GMT]',
          location.long = Longitude,
          location.lat = Latitude) %>%
  select(individual.local.identifier, timestamp, location.long, location.lat) %>%
  filter(timestamp >  '2022-4-8' & timestamp <  '2022-5-22') %>%
  mutate(individual.local.identifier = "pig_583") 

# # of days monitored
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
# animal 582 ----------------
animal = animal_582  %>% 
  as.telemetry(., timezone = 'GMT') # convert to telementry object 
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

# animal 583 ----------------
animal = animal_583%>% 
  as.telemetry(., timezone = 'GMT') # convert to telementry object 
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

# summary: 
# <<<< For 582 (monitored for 28 days), 95% UD is 0.89 [0.69, 1.11] km2, 50% UD is 0.15 [0.11, 0.18] km2. 
# <<<< For 583 (monitored for 39 days), 95% UD is 10.79 [3.65, 21.73] km2, 50% UD is 2.92 [0.99, 5.88] km2. 

############################################################################################
## Movement patterns  ----------------------------------------------------------------------
############################################################################################

animal = animal_582 %>% drop_na() # 35 NA rows dropped
animal.XY <- animal %>%   
  mutate(longitude = animal$location.long, latitude = animal$location.lat) %>%
  st_as_sf(., coords = c("location.long", "location.lat"), crs = 4326) %>%
  st_transform(., 32650) %>%
  cbind(., st_coordinates(.)) %>%
  st_set_geometry(NULL)  # turn into a df with both lat long and x y

traj <- make_track(animal.XY, .x = X, .y = Y, .t = timestamp, 
                   crs = sp::CRS("+init=epsg:32650"))

# visualize speed 

# visualize NSD

# behavior change detection (is there a "typical" movement state?)



############################################################################################
## HR use patterns  -----------------------------------------------------------------------
############################################################################################
# this will require LC extraction 


# interesting movement visualization (making maps)