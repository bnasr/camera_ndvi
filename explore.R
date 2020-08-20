

#
dukehw_bijan <- fread('data/summ/dukehw_DB_1000_ndvi_3day.csv')
dukehw_bijan[, date:=as.Date(date)]

tmp <- tempfile()
download.file(url = 'https://phenocam.sr.unh.edu/data/archive/dukehw/ROI/dukehw_DB_1000_ndvi_3day.csv', destfile = tmp)
dukehw_tom <- fread(tmp)
dukehw_tom[, date:=as.Date(date)]

tmp <- tempfile()
download.file(url = 'https://phenocam.sr.unh.edu/data/archive/dukehw/ROI/dukehw_DB_1000_ndvi.csv', destfile = tmp)
dukehw_tom_stats <- fread(tmp)
dukehw_tom_stats[, date:=as.Date(date)]
dukehw_tom_stats[, datetime:=as.POSIXct(paste(date, local_std_time))]

dukehw_bijan_stats <- fread('data/stats/dukehw_DB_1000_ndvi_roistats.csv')
dukehw_bijan_stats[, datetime:=as.POSIXct(datetime)]


dukehw_tom[, plot(date, ndvi_50)]
dukehw_bijan[, plot(date,ndvi_50)]

dt <- merge(dukehw_bijan_stats, dukehw_tom_stats, by ='datetime')

dt[, plot(NDVI_c.x, NDVI_c.y)]
dt[, plot(exposure_rgb.x, exposure_rgb.y)]
dt[, plot(exposure_ir.x, exposure_ir.y)]
