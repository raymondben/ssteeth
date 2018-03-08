library(raster)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(rworldmap)

data_dir <- "/rdsi/PUBLIC/raad/data"

cmap <- c("#7599D5FF","#6784CFFF","#5A6FC7FF","#535EC1FF","#504FB9FF","#5146B2FF","#5A46ABFF","#6C47A3FF","#844B9BFF","#9E5492FF","#B95E8AFF","#D26982FF","#E5747BFF","#F17F76FF","#F88B72FF","#FD996FFF","#FEA86CFF","#FFB96BFF","#FFC769FF","#FED367FF","#FEDB66FF")

roi <- c(94,180,-60,-30) ## region of interest E W S N

mapcoast <- getMap(resolution = "low") ## coastline data for plotting
## crop to our region of interest
mapcoast <- crop(mapcoast,extent(roi))

## read chronology data
x <- read_csv("VHamilton_whale chrons.csv")
## convert to long format
x <- x %>% gather(key="series",value="index_value",-year)

## so to extract e.g. albany series
## x %>% filter(series=="albany")

## ersst data
sst <- stack(file.path(data_dir, "ftp.cdc.noaa.gov/Datasets/noaa.ersst/sst.mnmean.v4.nc"))

sst <- crop(sst, extent(roi))

sstx <- as_tibble(values(sst)) ## put data into cat called mr tibbles
sstx$lon <- coordinates(sst)[,1]
sstx$lat <- coordinates(sst)[,2]
sstx <- sstx %>% gather(key="date",value="sst",-lon,-lat)
## is equivalent to
##sstx <- gather(sstx, key="date",value="sst",-lon,-lat)

sstx <- sstx %>% mutate(date=ymd(sub("^X","",date))) %>%
    mutate(year=year(date), month=month(date))

sstx_seasonal <- sstx %>% mutate(season=case_when(month %in% c(12,1,2)~"summer",
                                                  month %in% c(3,4,5)~"autumn",
                                                  month %in% c(6,7,8)~"winter",
                                                  month %in% c(9,10,11)~"spring")) %>%
    group_by(lon,lat,year,season) %>% summarize(sst=mean(sst)) %>%
    ungroup
## adjust the "year" of winter and spring to be the previous year
## so that year X of tooth time series matches to the previous year's winter or spring SST
sstx_seasonal <- sstx_seasonal %>% mutate(year=case_when(season %in% c("winter","spring")~year-1,
                                                         TRUE~year))

this_series <- "all_tas"
this_series <- "albany"

target_season <- "autumn"
target_season <- "spring"
target_season <- "winter"
target_season <- "summer"

thisx <- x %>% filter(series==this_series) %>% na.omit

rmap <- sst[[1]] ## use this as a template
values(rmap) <- NA_real_ ## but chuck away all the values in it
pmap <- rmap
ll_grid <- coordinates(rmap)
temp <- sstx_seasonal %>% filter(year %in% thisx$year & season==target_season)
for (i in seq_len(prod(dim(rmap)))) {
    this_sst <- temp %>% filter(lon==ll_grid[i,1] & lat==ll_grid[i,2]) %>% arrange(year)
    try({
        this_r <- cor.test(this_sst$sst, thisx$index_value, method="spearman")
        rmap[i] <- this_r$estimate
        pmap[i] <- this_r$p.value
        }, silent=TRUE)
}

plot(rmap, col=cmap, main=paste0(this_series," ",target_season))
plot(mapcoast, xlim = range(ll_grid[,1]), ylim = range(ll_grid[,2]), col = "grey50", add = TRUE)

plot(pmap<0.05, col=cmap, main=paste0(this_series," ",target_season))


## check one cell
this_sst <- temp %>% filter(lon==120 & lat==-50) %>% arrange(year)
plot(this_sst$sst,thisx$index_value)
