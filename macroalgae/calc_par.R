### Script to calculate theroteitcal par atthe sea surface for any latitude

library(atmos)
modelled_par<-function(rlat){
  #using the Gregg and Carder Scheme as implemented in the atmos package provided by Bernard Gentili <gentili@obs-vlfr.f
  # calculate the daily total par irradiance for the latitude given fora whole year 
  
  #from G+C - to calculate integrated total par
  
  
  day_total<-function(jday,rlat){
    day<-seq(from=0.4,to=23.9,by=0.5)
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

# bantry_par<-modelled_par(51.65)
# write.csv(file='bantry_par_theoretical.csv',data.frame(jday=1:365,par_GC=bantry_par),row.names = FALSE)


generate_par_table<-function(){
  #generate table of par data for each day for every latitude 0-90 and write it out as a csv
  latitudes<-0:90
  x<-sapply(latitudes,modelled_par)
  write.table(x=x,file = 'PARTable.csv',row.names=FALSE, col.names=FALSE,sep=',')
  return('done')
}


# uncomment line below to generate new partable (note this takes a good few minutes)
# generate_par_table()

substitute_par<-function(observed_par,modelled_par){
  #function to subsitute where PAR values are missing - modelled par is clear-sky so use scaling factor to give representative par value
  #this function expects inputs which are 365 (days) long. Anything else will make bad things happen so better handle that outside of this function..
  normalised_par<-observed_par/modelled_par
  
  indexer<-which(is.na(observed_par))
  subs<-rnorm(n=length(indexer),mean = mean(normalised_par,na.rm=TRUE),sd = sd(normalised_par,na.rm=TRUE))
  observed_par[indexer]<-pmin(modelled_par[indexer],modelled_par[indexer]*subs)
  observed_par
  }
