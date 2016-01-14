# --------------------------------------------------
# hydromad objective functions
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# September 2015
# -------------------------------

# use Viney's objective function(includes Bias), to fit 
# see http://hydromad.catchment.org/#hydromad.stats
hydromad.stats("viney" = function(Q, X, ...) {
  hmadstat("r.squared")(Q, X, ...) -
    5*(abs(log(1+hmadstat("rel.bias")(Q,X)))^2.5)})

# This is one way to do this, but uses the NSE
# Define an objective function that aggregates ET based on specified periods
# use buildTsObjective, but this uses the NSE
hydromad.stats("ETfun" = function(...,DATA,U) {
  .(buildTsObjective(DATA$aET, groups=DATA$et.period, 
                     FUN=sum))(DATA$aET,U$ET,...)})

# Define an objective function that aggregates ET based on specified periods
#This can be any objective function, specify using objf
hydromad.stats("ETaggrViney" = function(..., DATA=DATA,U=U,objf=hmadstat("viney")) {
  #This should be mean if the observed aET is repeated for each point in the period
  # using sum because inserted 0 values in data (ETa.merge)
  # inserted "coredata" statement after discussion with J Guillaume
  # relates to how bias is calculated
  aET.fin <- aggregate(DATA$aET,list(date=coredata(DATA$et.period)),sum)
  ET.fin <- aggregate(U$ET,list(date=coredata(DATA$et.period)),sum)

  obj <- objf(coredata(aET.fin),coredata(ET.fin),...)
  return(obj)
})


# Build a master objective function that has w as a further parameter
hydromad.stats("JointQandET" = function(Q,X,...,w,DATA,U,objf=hmadstat("viney")) {
#browser()
  # prepare the ET data
  aET.fin <- aggregate(coredata(DATA$aET),list(date=coredata(DATA$et.period)),sum)[,2]
  ET.fin <- aggregate(U$ET,list(date=coredata(DATA$et.period)),sum)
  # calculate overall objective function
  obj <- w*objf(Q,X,...) + (1-w)*objf(coredata(aET.fin),coredata(ET.fin),...)
  return(obj)

})