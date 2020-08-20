get_ndvi <- function(roi_name, 
                     cache = TRUE, 
                     solar_elev_threshold = 5, 
                     qtl = 'mean',
                     lag_threshold_sec = 600){
  
  roi_name_parsed <- strsplit(roi_name, split = '_')[[1]]
  
  
  if(qtl%in%c('5', '10', '25', '50', '75', '90', '95')){
    qtl = paste0(qtl, '_qtl')
  }else if(qtl!='mean') {
    stop('the index argument is not valid!')
  }
  
  rgb_url <- sprintf('%s%s/ROI/%s_roistats.csv', phenocam_url, roi_name_parsed[1], roi_name)
  ir_url <- sprintf('%s%s/ROI/%s_IR_roistats.csv', phenocam_url, roi_name_parsed[1], roi_name)
  
  if(!url.exists(rgb_url)) stop(paste(rgb_url, 'does not exist!'))
  if(!url.exists(ir_url)) stop(paste(ir_url, 'does not exist!'))
  
  rgb_file <- paste0(stats_data_dir, basename(rgb_url))
  ir_file <- paste0(stats_data_dir, basename(ir_url))
  
  if(cache==FALSE){
    download.file(rgb_url, rgb_file)
    download.file(ir_url, ir_file)
  }
  
  if(!file.exists(rgb_file)) download.file(rgb_url, rgb_file)
  if(!file.exists(ir_file)) download.file(ir_url, ir_file)
  
  rgb_stats <- fread(rgb_file)
  ir_stats <- fread(ir_file)
  
  rgb_stats[, datetime := as.POSIXct(paste(date, local_std_time))]
  ir_stats[, datetime := as.POSIXct(paste(date, local_std_time))]
  
  rgb_data <- rgb_stats[,c('datetime', 'solar_elev', 'exposure', paste0('r_', qtl),  paste0('g_', qtl), paste0('b_', qtl)), with=FALSE]
  ir_data <- ir_stats[,c('datetime', 'exposure', paste0('ir_', qtl)), with=FALSE]
  
  rgb_data[, key_code := floor(as.integer(datetime) / lag_threshold_sec) * lag_threshold_sec]
  ir_data[, key_code := floor(as.integer(datetime) / lag_threshold_sec ) * lag_threshold_sec]
  
  dt <- merge(rgb_data, 
              ir_data[, c(paste0('ir_', qtl), 'key_code', 'exposure'), with = FALSE], 
              by = 'key_code')
  
  colnames(dt) <- c('key_code', 'datetime', 'solar_elev', 'exposure_rgb', 'r', 'g', 'b', 'ir', 'exposure_ir')
  
  dt$key_code <- NULL
  dt <- dt[solar_elev<5]
  
  dt[, Y := 0.3 * r + 0.59 * g + 0.11 * b]
  dt[, Z_prime := ir / sqrt(exposure_ir)]
  dt[, R_prime := r / sqrt(exposure_rgb)]
  dt[, Y_prime := Y / sqrt(exposure_rgb)]
  dt[, X_prime := Z_prime - Y_prime]
  dt[, NDVI_c := (X_prime - R_prime) / (X_prime + R_prime)]
  dt[, EVI2_c := 2.5 * (X_prime - R_prime) / (X_prime + 2.4 * R_prime + 1)]
  
}
