lunar_irrad_DNB <- function(date_time){
  print("Input Date and Time (Format = YYYYMMDDHHmm):")
  print("(YYYY= year (2010-2030), MM=month (01-12), DD=day-of-month, HH=UTC-hour (00-23), mm=minue (00-59)")
  print("NOTE VALID RANGE: From 201001010000 (Jan 1, 2010, 0000 UTC) to 203012312300 (Dec 31, 2030, 2300 UTC)")
  print(paste0("User-entered date_time = ",date_time))
  
  test_year = as.numeric(as.numeric(date_time)/1.0e+08)
  if((test_year < 2010) | (test_year > 2030)){
    print(paste0('Illegal Year (Valid for 2010-2030): ', test_year))
  }
  lun_irrad_scl(as.numeric(date_time))

}

lun_irrad_scl <- function(date_time){
  
  mean_earthsun_dist=149597870.700
  mean_earthmoon_dist=384402.0
  radius_earth=6378.140
  #mean_earthsun_dist=149598016
  #mean_earthmoon_dist=384400.0
  #radius_earth=6378.140
  
  
  #Shorten date/time to match table and separate minutes out
  
  egroup <- as.integer(trunc(date_time/100))
  
  mins <- date_time - egroup*100 
  
  if (mins == 0) {
    exact =  TRUE 
  }else{
    exact =  FALSE 
  }
  
  unit20=read.table("DIST_2010-2030.dat")
  colnames(unit20) <- c('group','phase_angle','earthsun_dist','earthmoon_dist')
  
  
  #Find revelant moon phase angle, earth-sun and earth moon distance from table
  
  group1 <- unit20[unit20$group == egroup,]$group
  phase_angle1 <- unit20[unit20$group == egroup,]$phase_angle
  earthsun_dist1 <- unit20[unit20$group == egroup,]$earthsun_dist
  earthmoon_dist1 = unit20[unit20$group == egroup,]$earthmoon_dist
  
  if (exact){
    group = group1
    current_phase_angle = phase_angle1
    current_earthsun_dist = earthsun_dist1
    current_earthmoon_dist = earthmoon_dist1
  }
  
  if (!exact){
    queryrow = unit20[which(unit20$group == group1)+1,]
    phase_angle2 <- unit20[which(unit20$group == group1)+1,]$phase_angle
    earthsun_dist2 <- unit20[which(unit20$group == group1)+1,]$earthsun_dist
    earthmoon_dist2 <- unit20[which(unit20$group == group1)+1,]$earthmoon_dist
    frac = mins/60.0
    delta = phase_angle2 - phase_angle1
    current_phase_angle = phase_angle1 + (delta * frac)
    delta = earthsun_dist2 - earthsun_dist1
    current_earthsun_dist = earthsun_dist1 + (delta * frac)
    delta = earthmoon_dist2 - earthmoon_dist1
    current_earthmoon_dist = earthmoon_dist1 + (delta * frac)
  }
    
    
  unit22=read.table("lunar_irrad_1AU_MeanME_CONVOLVED_DNB.dat",skip = 3)
  
  
  #unit22=pd.read_csv("lunar_irrad_1AU_MeanME_CONVOLVED_DNB.dat",sep="\t",skiprows = 3, names = 'x')
  colnames(unit22) = c("phase","lunar_irrad")
  
  
  #Determine convolved irradiance from table
  queryrow = unit22[unit22$phase == ceiling(current_phase_angle),]
  queryrow_prev = unit22[unit22$phase == ceiling(current_phase_angle)-1,]
  frac = (queryrow$phase - current_phase_angle) / (queryrow$phase - queryrow_prev$phase)
  delta = queryrow$lunar_irrad - queryrow_prev$lunar_irrad
  lunar_irrad_dnb_interp = queryrow$lunar_irrad - (delta * frac)
  
  
  #Compute scaling factor
  cos_phase_angle = cos(current_phase_angle * pi/180)
  T1 = mean_earthsun_dist**2.0 + mean_earthmoon_dist**2.0 + 2.0*mean_earthmoon_dist*mean_earthsun_dist*cos_phase_angle
  T2 = current_earthsun_dist**2.0 + current_earthmoon_dist**2.0 + 2.0*current_earthmoon_dist*current_earthsun_dist*cos_phase_angle
  T3 = ((mean_earthmoon_dist-radius_earth) / (current_earthmoon_dist - radius_earth))**2.0
  
  SCALE_FACTOR = (T1/T2)*T3
  
  #Compute lunar irradiance
  
  lun_irrad_scl = lunar_irrad_dnb_interp * SCALE_FACTOR
  
  print(paste0("Current Earth/Moon Distance (km) = ",current_earthmoon_dist))
  print(paste0(" % of mean = ", 100.0*current_earthmoon_dist/mean_earthmoon_dist))
  print(paste0("Current Earth/Sun Distance (km) = ",current_earthsun_dist))
  print(paste0(" % of mean = ", 100.0*current_earthsun_dist/mean_earthsun_dist))
  print(paste0("Current lunar phase angle (degrees) = ",current_phase_angle))
  print(paste0("SCALE_FACTOR = ",SCALE_FACTOR))
  print(paste0("Lunar Irradiance (mW/m^2-micron):", lun_irrad_scl))
}
