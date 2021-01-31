# prepare NYC taxi data for non-parametric Bayesian updating

setwd("~/Documents/Shop-window-NOC/non-parametric-updating/NYCtaxis")
datafolder<-"~/Documents/NYCtaxis/"

library(readr)
library(dplyr)
library(lubridate)


# split up 1st trip data file

currentfile<-read_csv(paste0(datafolder,"split1.csv"),
                      col_types=cols())

# drop variables
currentfile<-select(currentfile,
                    medallion,
                    pickup_datetime,
                    pickup_longitude,
                    pickup_latitude)

# generate datetime variables (wkday 1=Monday)
currentfile<-mutate(currentfile,
                    day=day(pickup_datetime),
                    month=month(pickup_datetime),
                    hour=hour(pickup_datetime),
                    wkday=wday(pickup_datetime,
                               week_start=1))

for(i in 1:31) {
  write_csv(filter(currentfile,day==i),
            paste0("trip_data_2013_1_",i,".csv"))
}

# then loop over the other files and allocate to the days
# you may see warnings about parsing the 
#    "store_and_fwd_flag" column but you can ignore these
pb<-txtProgressBar(min=0,max=296,style=3)
for(j in 2:296) {
  setTxtProgressBar(pb,j)
  currentfile<-read_csv(paste0(datafolder,"split",j,".csv"),
                        col_types=cols())
  # drop variables
  currentfile<-select(currentfile,
                      medallion,
                      pickup_datetime,
                      pickup_longitude,
                      pickup_latitude)
  # generate datetime variables (wkday 1=Monday)
  currentfile<-mutate(currentfile,
                      day=day(pickup_datetime),
                      month=month(pickup_datetime),
                      hour=hour(pickup_datetime),
                      wkday=wday(pickup_datetime,
                                 week_start=1))
  for(i in 1:31) {
    write_csv(filter(currentfile,day==i),
              paste0("trip_data_2013_1_",i,".csv"),
              append=TRUE)
  }
}
close(pb)

#####################################################

# split up 1st fare data file

currentfile<-read_csv(paste0(datafolder,"split1.csv"),
                      col_types=cols())

# drop variables
currentfile<-select(currentfile,
                    medallion,
                    pickup_datetime,
                    total_amount)

# generate datetime variables 
currentfile<-mutate(currentfile,
                    day=day(pickup_datetime),
                    logtotal=log(total_amount))

for(i in 1:31) {
  write_csv(filter(currentfile,day==i),
            paste0("fare_data_2013_1_",i,".csv"))
}

# then loop over the other files and allocate to the days
pb<-txtProgressBar(min=0,max=296,style=3)
for(j in 2:296) {
  setTxtProgressBar(pb,j)
  currentfile<-read_csv(paste0(datafolder,"split",j,".csv"),
                        col_types=cols())
  # drop variables
  currentfile<-select(currentfile,
                      medallion,
                      pickup_datetime,
                      total_amount)
  # generate datetime variables (wkday 1=Monday)
  currentfile<-mutate(currentfile,
                      day=day(pickup_datetime),
                      logtotal=log(total_amount))
  for(i in 1:31) {
    write_csv(filter(currentfile,day==i),
              paste0("fare_data_2013_1_",i,".csv"),
              append=TRUE)
  }
}
close(pb)
rm(currentfile)

#####################################################


# one day at a time, merge trip and fare data
for(i in 1:31) {
  tripfile<-read_csv(paste0("trip_data_2013_1_",i,".csv"),
                     col_types=cols())
  farefile<-read_csv(paste0("fare_data_2013_1_",i,".csv"),
                     col_types=cols())
  taxidata<-inner_join(farefile,tripfile)
  
  # drop variables
  taxidata<-select(taxidata,
                   logtotal,
                   pickup_longitude,
                   pickup_latitude,
                   month,
                   hour,
                   wkday)
  
  # write
  write_csv(taxidata,
            paste0("taxi_data_2013_1_",i,".csv"))
}
rm(farefile)
rm(tripfile)
