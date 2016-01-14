# --------------------------------------------------
# Utility function to create the input data set for the fitting of ET
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# September 2015
# -------------------------------

# merges the data and puts in et.period to allow 8 day aggregation
ETa.merge <- function(Flowdata,ETdata, fill=0) {
  #browser()
  # Flowdata is the standard hydromad input
  # ETdata is the zoo object with 8 day ET data
  # fill is the value to put in the missing dates for ETa
  # this can be NA or 0
  # calculate days.cor    
  days.cor <- c(diff(time(ETdata)),6)
  # Now we need to replace the -360 values at the end of the year
  # run leap.fun to account for leap years
  # calculate et.period
  et.period <- zoo(rep(time(ETdata), times=days.cor),
                   order.by= seq.Date(as.Date(time(ETdata)[1]),
                                      by=1,length=sum(days.cor)))
  # now merge all
  out.zoo <- merge(Flowdata ,aET=ETdata, et.period=et.period,all=T)
  # insert fill values for missing aET data
  out.zoo[is.na(out.zoo$aET),"aET"] <- fill
  return(out.zoo)
}
