library(raster)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(rworldmap)

data_dir <- "/rdsi/PUBLIC/raad/data"

cmap <- c("#7599D5FF", "#6784CFFF", "#5A6FC7FF", "#535EC1FF", "#504FB9FF", "#5146B2FF", "#5A46ABFF", "#6C47A3FF", "#844B9BFF", "#9E5492FF", "#B95E8AFF", "#D26982FF", "#E5747BFF", "#F17F76FF", "#F88B72FF", "#FD996FFF", "#FEA86CFF", "#FFB96BFF", "#FFC769FF", "#FED367FF", "#FEDB66FF")

roi <- c(93, 181, -61, -29) ## region of interest E W S N ##*** expanded slightly

mapcoast <- getMap(resolution = "low") ## coastline data for plotting
## crop to our region of interest
mapcoast <- crop(mapcoast, extent(roi)+c(-1, 1, -1, 1)) ##*** add a bit around the edges

## read chronology data
x <- read_csv("VHamilton_whale chrons.csv")
## convert to long format
x <- x %>% gather(key="series", value="index_value", -year)

## so to extract e.g. albany series
## x %>% filter(series=="albany")

## choose which env data to use
env_data_to_use <- "ersstv5" ## can be "ersstv4" "ersstv5" "ncep_wind"

if (env_data_to_use=="ersstv4") {
    ## ersst v4 data
    envdat <- stack(file.path(data_dir, "ftp.cdc.noaa.gov/Datasets/noaa.ersst/sst.mnmean.v4.nc"))
} else if (env_data_to_use=="ersstv5") {
    ## ersst v5 data
    ersst5_files <- dir(file.path(data_dir, "ftp.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf"), pattern="nc$", full.names=TRUE)
    ersst5_dates <- ymd(paste0(sub("\\.nc$", "", sub(".*ersst\\.v5\\.", "", ersst5_files)), "01")) ## extract dates from file names
    ## subset to years of interest plus or minus a few
    idx <- year(ersst5_dates)>=(min(x$year)-2) & year(ersst5_dates)<=(max(x$year)+2)
    ersst5_files <- ersst5_files[idx]
    ersst5_dates <- ersst5_dates[idx]
    envdat <- stack(as.list(ersst5_files), varname="sst") ## also have "ssta" for anomaly available here
    ## set the names of the layers in this brick to be the corresponding dates
    names(envdat) <- paste0("X", ersst5_dates)
} else if (env_data_to_use=="ncep_wind") {
    envdat <- stack(file.path(data_dir, "ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/vwnd.mon.mean.nc")) ## u-component is zonal, v- is meridional
    ## keep just the YMD component of layer names
    names(envdat) <- substr(names(envdat), 1, 11) ## first 11 chars of each will be Xyyyy.mm.dd
}

envdat <- crop(envdat, extent(roi+c(-3, 3, -3, 3))) ##*** crop with extra space around the edges

envdatx <- as_tibble(values(envdat)) ## extract data from raster object and put it into a tibble (data.frame)
## add coords
envdatx$lon <- coordinates(envdat)[, 1]
envdatx$lat <- coordinates(envdat)[, 2]
## long format
envdatx <- envdatx %>% gather(key="date", value="env_value", -lon, -lat)
## is equivalent to
##envdatx <- gather(envdatx, key="date", value="env_value", -lon, -lat)

## reconvert the date back to an actual date object and extract the year and month
envdatx <- envdatx %>% mutate(date=ymd(sub("^X", "", date))) %>%
    mutate(year=year(date), month=month(date))

## construct a seasonal version of environmental data by averaging months
envdatx_seasonal <- envdatx %>% mutate(season=case_when(month %in% c(12, 1, 2)~"summer",
                                                  month %in% c(3, 4, 5)~"autumn",
                                                  month %in% c(6, 7, 8)~"winter",
                                                  month %in% c(9, 10, 11)~"spring")) %>%
    group_by(lon, lat, year, season) %>% summarize(env_value=mean(env_value)) %>%
    ungroup
## adjust the "year" of winter and spring to be the previous year
## so that year X of tooth time series matches to the previous year's winter or spring
envdatx_seasonal <- envdatx_seasonal %>% mutate(year=case_when(season %in% c("winter", "spring")~year-1,
                                                         TRUE~year))

this_series <- "all_tas"
this_series <- "albany"

target_season <- "autumn"
target_season <- "spring"
target_season <- "winter"
target_season <- "summer"

thisx <- x %>% filter(series==this_series) %>% na.omit

rmap <- envdat[[1]] ## use this as a template
values(rmap) <- NA_real_ ## but chuck away all the values in it
pmap <- rmap
ll_grid <- coordinates(rmap)
temp <- envdatx_seasonal %>% filter(year %in% thisx$year & season==target_season) ## subset the seasonal data once, outside of the loop
for (i in seq_len(prod(dim(rmap)))) {
    this_envdat <- temp %>% filter(lon==ll_grid[i, 1] & lat==ll_grid[i, 2]) %>% arrange(year)
    try({
        ## join this_envdat onto our time series
        corx <- thisx %>% left_join(this_envdat %>% dplyr::select(year, env_value), by="year")
        ## this will have all rows in x but NA env_values where we don't have environmental data (e.g. NCEP wind prior to 1947 or so)
        this_r <- cor.test(corx$env_value, corx$index_value, method="spearman", na.action="na.omit")
        rmap[i] <- this_r$estimate
        pmap[i] <- this_r$p.value
        }, silent=TRUE)
}

plot(rmap, col=cmap, main=paste0(this_series, " ", target_season))
plot(mapcoast, xlim = range(ll_grid[, 1]), ylim = range(ll_grid[, 2]), col = "grey50", add = TRUE)

plot(pmap<0.05, col=cmap, main=paste0(this_series, " ", target_season))


## check one cell
this_envdat <- temp %>% filter(lon==120 & lat==-50) %>% arrange(year)
plot(this_envdat$env_value, thisx$index_value)

## ggplot maps
library(ggplot2)

trim_tiles_to_extent <- function(data, x, y, width, height, extent) {
    ## assumes all tiles are equally sized
    ## extent should be c(xmin, xmax, ymin, ymax)
    xmin <- pmax(data[[x]]-width/2, extent[1])
    xmax <- pmin(data[[x]]+width/2, extent[2])
    idx <- xmin==xmax
    xmin[idx] <- NA
    xmax[idx] <- NA
    ymin <- pmax(data[[y]]-height/2, extent[3])
    ymax <- pmin(data[[y]]+height/2, extent[4])
    idx <- ymin==ymax
    ymin[idx] <- NA
    ymax[idx] <- NA
    data %>% mutate(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
}

ggcoast <- fortify(mapcoast)

ggx <- as.data.frame(cbind(coordinates(rmap), values(rmap)))
names(ggx) <- c("long", "lat", "r")

## find the largest (absolute) value so that we have positive and negative colour limits matching
colour_max <- max(abs(ggx$r), na.rm=TRUE)

## helper functions to make labels look nicer
num2lon <- function(z) paste0(z, " \u00B0E")
num2lat <- function(z) paste0(abs(z), " \u00B0S")

ggx <- trim_tiles_to_extent(ggx, x="long", y="lat", width=res(rmap)[1], height=res(rmap)[2], extent=roi)

ggplot(ggx, aes(long, lat)) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=r)) +
    ## use a red-blue anomaly map so that positive values are red and negative ones blue
    scale_fill_distiller(palette="RdBu", na.value="white", name="Correlation", limits=c(-1, 1)*colour_max) +
    ## add coast
    geom_polygon(data=ggcoast, aes(group=group), color="black", fill="grey50") +
    coord_fixed() + ## fixed aspect ratio. Probably better to use coord_map() here but it is appallingly slow
    ## get rid of decorations we don't want, make font bigger
    theme_bw() + theme(axis.line = element_blank(), text = element_text(size=18)) +
    ## set axis limits and format the labels
    scale_x_continuous(limits=roi[1:2], expand=c(0,0), labels=num2lon) +
    scale_y_continuous(limits=roi[3:4], expand=c(0,0), labels=num2lat) +
    labs(x="", y="")
