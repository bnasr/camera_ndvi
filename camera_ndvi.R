library(RCurl)
library(plotly)
library(phenocamapi)
library(data.table)

source('funcs.R')

phenocam_url <- 'https://phenocam.sr.unh.edu/data/archive/'

modis_data_dir<- 'data/modis_time_series/'
stats_data_dir <- 'data/stats/'
figs_data_dir <- 'data/figs/'
summ_data_dir <- 'data/summ/'

dir.create(figs_data_dir, showWarnings = FALSE)
dir.create(stats_data_dir, showWarnings = FALSE)
dir.create(summ_data_dir, showWarnings = FALSE)




# phenos <- get_phenos()
# rois <- get_rois()
# ir_rois <- rois[site %in% phenos[infrared=='Y', site]]
# ir_rois[, ir_url := paste0(phenocam_url, site, '/ROI/', roi_name, '_IR_roistats.csv')]
# w = mapply(FUN = url.exists, ir_rois$ir_url)
# write.csv(ir_rois[which(as.vector(w))], row.names = FALSE, file = 'data/ir_rois.csv')

ir_rois <- fread('data/ir_rois.csv')
i=1
# i = 83
# roi_name <- 'dukehw_DB_1000'

n <- nrow(ir_rois)

# mapply(download.file, ir_rois$one_day_summary, paste0(summ_data_dir, basename(ir_rois$one_day_summary)))
# 
# mapply(download.file, ir_rois$three_day_summary, paste0(summ_data_dir, basename(ir_rois$three_day_summary)))
for(qtl in c('5', '10', '25', '50', '75', '90', '95', 'mean')){
  
  for(i in 1:n){
    site <- ir_rois[i, site]
    roi_name <- ir_rois[i, roi_name]
    
    cat(roi_name, qtl, '\n')
    
    ndvi_stats <- try(get_ndvi(roi_name, cache = TRUE, qtl = qtl))
    
    if(class(ndvi_stats)[1]=='try-error') next()
    if(nrow(ndvi_stats)==0) next()
    
    ndvi_stats[, year := year(datetime)]
    ndvi_stats[, d1 := yday(datetime)]
    ndvi_stats[, d3 := ceiling(d1/3)*3-1]
    
    ndvi_d1 <- ndvi_stats[, .(ndvi_mean = mean(NDVI_c, na.rm = T),
                              ndvi_std = sd(NDVI_c, na.rm = T),
                              ndvi_50 = quantile(NDVI_c, na.rm = T, probs = 0.50),
                              ndvi_75 = quantile(NDVI_c, na.rm = T, probs = 0.75),
                              ndvi_90 = quantile(NDVI_c, na.rm = T, probs = 0.90),
                              
                              evi2_mean = mean(EVI2_c, na.rm = T),
                              evi2_std = sd(EVI2_c, na.rm = T),
                              evi2_50 = quantile(EVI2_c, na.rm = T, probs = 0.50),
                              evi2_75 = quantile(EVI2_c, na.rm = T, probs = 0.75),
                              evi2_90 = quantile(EVI2_c, na.rm = T, probs = 0.90)),
                          .(year, doy = d1)]
    
    ndvi_d3 <- ndvi_stats[, .(ndvi_mean = mean(NDVI_c, na.rm = T),
                              ndvi_std = sd(NDVI_c, na.rm = T),
                              ndvi_50 = quantile(NDVI_c, na.rm = T, probs = 0.50),
                              ndvi_75 = quantile(NDVI_c, na.rm = T, probs = 0.75),
                              ndvi_90 = quantile(NDVI_c, na.rm = T, probs = 0.90),
                              
                              evi2_mean = mean(EVI2_c, na.rm = T),
                              evi2_std = sd(EVI2_c, na.rm = T),
                              evi2_50 = quantile(EVI2_c, na.rm = T, probs = 0.50),
                              evi2_75 = quantile(EVI2_c, na.rm = T, probs = 0.75),
                              evi2_90 = quantile(EVI2_c, na.rm = T, probs = 0.90)),
                          .(year, doy = d3)]
    
    x <- try({
      modis <- fread(paste0(modis_data_dir, site, '_gee_subset.csv'))
      modis[, date:=as.Date(date)]
      modis[, ndvi := (Nadir_Reflectance_Band2 - Nadir_Reflectance_Band1) /
              (Nadir_Reflectance_Band2 + Nadir_Reflectance_Band1)]
      
      modis[, evi := 2.5 * (Nadir_Reflectance_Band2 - Nadir_Reflectance_Band1) / 
              (Nadir_Reflectance_Band2 + 6 * Nadir_Reflectance_Band1 - 7.5 * Nadir_Reflectance_Band3 + 1)]
      
      modis[, evi2 := 2.5 * (Nadir_Reflectance_Band2 - Nadir_Reflectance_Band1) / 
              (Nadir_Reflectance_Band2 + 2.4 * Nadir_Reflectance_Band1 + 1)]
      
      # modis[abs(evi)>1, evi := NA]
      # modis[abs(evi2)>1, evi2 := NA]
    })
    
    ndvi_d1[, date := as.Date(doy, origin = paste0(year -1 , '-12-31'))]
    ndvi_d3[, date := as.Date(doy, origin = paste0(year -1 , '-12-31'))]
    
    write.csv(ndvi_stats, file = paste0(stats_data_dir, roi_name, '_ndvi_roistats.csv'), row.names = FALSE)
    write.csv(ndvi_d1, file = paste0(summ_data_dir, roi_name, '_ndvi_1day.csv'), row.names = FALSE)
    write.csv(ndvi_d3, file = paste0(summ_data_dir, roi_name, '_ndvi_3day.csv'), row.names = FALSE)
    
    gcc_d1 <- fread(paste0(summ_data_dir, roi_name, '_1day.csv'))
    gcc_d3 <- fread(paste0(summ_data_dir, roi_name, '_3day.csv'))
    
    gcc_d1[, date := as.Date(date)]
    gcc_d3[, date := as.Date(date)]
    
    
    
    # svg(filename = paste0(figs_data_dir, roi_name, '_', ifelse(qtl=='5','05',qtl),'.svg'), width = 8, height = 6)
    png(filename = paste0(figs_data_dir, roi_name, '_', ifelse(qtl=='5','05',qtl),'.png'), width = 12, height = 6, units = 'in', res = 600)
    layout(matrix(1:8, 4,2, byrow = T), widths = c(7,2))
    par(mar=c(0,0,0,6), oma = c(5,0,2,0), bty ='n')
    
    xlim <- range(gcc_d3$date)
    
    for(metric in c('mean', '50', '75', '90')){
      plot(gcc_d3$date, gcc_d3$gcc_90, type = 'l', col = 3, xlim = xlim, yaxt = 'n', xaxt = 'n', ylab ='', lwd =2)
      axis(4, line = -1.5, col = 3, col.ticks = 3, col.axis = 3)
      
      par(new = T)
      plot(ndvi_d3$date, ndvi_d3[[paste0('ndvi_',metric)]], type = 'l', col = 1, xlim = xlim, yaxt = 'n', xaxt = 'n', ylab ='', lwd = 2)
      axis(4, line = 0.5, col = 1, col.ticks = 1, col.axis = 1)
      
      par(new = T)
      plot(ndvi_d3$date, ndvi_d3[[paste0('evi2_',metric)]], type = 'l', col = 4, lty=1, xlim = xlim, yaxt = 'n', xaxt = 'n', ylab ='', lwd= 1)
      axis(4, line = 2.5, col = 4, col.ticks = 4, col.axis = 4)
      
      par(new = T)
      plot(modis$date, modis$ndvi, type = 'l', col = 2, xlim = xlim, yaxt = 'n', xaxt = 'n', ylab ='', lwd=1)
      axis(4, line = 4.5, col = 2, col.ticks = 2, col.axis = 2)
      
      par(new = T)
      plot(modis$date, modis$evi, type = 'l', col = 6, xlim = xlim, yaxt = 'n', xaxt = 'n', ylab ='', lwd=1, ylim = quantile(modis$evi, probs = c(0.01, 0.99), na.rm = T))
      # axis(4, line = 6.5, col = 6, col.ticks = 6, col.axis = 6)
      
      par(new = T)
      plot(modis$date, modis$evi2, type = 'l', col = 7, xlim = xlim, yaxt = 'n', xaxt = ifelse(metric =='90','s','n'), ylab ='', lwd=1, ylim = quantile(modis$evi2, probs = c(0.01, 0.99), na.rm = T))
      # axis(4, line = 8.5, col = 7, col.ticks = 7, col.axis = 7)
      
      plot(type='n', 1, 1, axes =F, xlab ='', ylab='')
      legend('center', legend = c('Gcc-90', paste0(c('NDVI-', 'EVI2-'), metric), 'MODIS-NDVI', 'MODIS-EVI', 'MODIS-EVI2'), bty= 'n', col = c(3,1,4,2,6,7), lty = c(1,1,1,1,1,1), lwd = 2)
    }
    
    mtext(outer = T, text = paste(roi_name, qtl))
    dev.off()
  }
  
}