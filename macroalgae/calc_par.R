### Script to calculate theroteitcal par atthe sea surface for any latitude

get_theoretical_par<-function(rlat){
  #using the Gregg and Carder Scheme as implemented in the atmos package provided by Bernard Gentili <gentili@obs-vlfr.f
  # calculate the daily total par irradiance for the latitude given fora whole year 
  
  #from G+C - to calculate integrated total par
  
  
  day_total<-function(jday,rlat){
    day<-seq(from=0,to=23.9,by=0.1)
    totalpar<-function(jday,rlat,hr){
      h = 6.6256E-34
      c = 2.998E8
      hc = 1.0/(h*c)
      x<-GreggCarder.f(jday=jday,rlat=rlat,rlon=0,hr=hr)
      1.0E-9 * hc/6.023E17 * sum(x$Ed * x$lam)
    }
    y<-sapply(day,totalpar,jday=jday,rlat=rlat)
    sum(y,na.rm = TRUE)/length(day)
  }
  
  year<-1:365
  sapply(year,day_total,rlat=rlat)
}


# example code to generate csv file for latitude of bantry bay

# bantry_par<-get_theoretical_par(51.65)
# write.csv(file='bantry_par_theoretical.csv',data.frame(jday=1:365,par_GC=bantry_par),row.names = FALSE)