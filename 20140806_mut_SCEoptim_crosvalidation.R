
## Muttama
## New version February 2014
## SCE optimisation

## different objective functions: using Viney's as well as rsq and rsqlog
## also include shape


library(hydromad)
library(lattice)
library(latticeExtra)
library(zoo)

rootdir <- "x:/vervoort/research/susie_Pub/my_project"

setwd(paste(rootdir,"/Data/Australia/Muttuma and Jugiong",sep=""))

#----------------------------
##        READ IN DATA       ####
#-----------------------------

## First need to have a data frame with P, Q, MODIS ET (named E)

#require(hydromad)

# flow data
Flow <- read.csv("Flow_Muttuma_01031970_01032012.csv")
Flow$Date <- as.Date(Flow[,1],"%d/%m/%Y")

#convert Flow from ML/day to mm/day
Flow$newFlow <- convertFlow(Flow$Flow, from="ML", area.km2=1159.8) # Area from ArcGIS

# rainfall data
Rain <- read.csv("Rainfall_Coota_Av_01011970_29022012.csv")
# colnames(Rain)[2] <- "Rain"
Rain$Date <- as.Date(Rain[,1],"%d/%m/%Y")

# set missing rainfall data to 0
Rain[,2] <- replace(Rain[,2], is.na(Rain[,2])==T,0)

# temperature
Temp <- read.csv("maxT_Coota_01011957_17032012.csv")
Temp$Date <- as.Date(Temp[,1],"%d/%m/%Y")

head(Temp)
nrow(Temp)
# Interpolate Temperature to deal with missing data
length(na.omit(Temp[,2]))
n <- 1:nrow(Temp)
Temp$MaxT2 <- loess.smooth(n, Temp$MaxT, span = 0.02, degree = 2, evaluation = length(Temp$MaxT))$y

foo <- rep(0,length(Temp$MaxT))
foo <- replace(foo,is.na(Temp$MaxT)==T,1)
foo1 <- replace(Temp$MaxT,is.na(Temp$MaxT)==T,0)
plot(foo)
foo.b <- foo*Temp$MaxT2+foo1
#plot(foo.b)


#read in MODIS ET
# Two different files to read in: logmean and normal mean. Only spline fit

# use trigger to decide
logged <- FALSE


ETdir <- paste(rootdir,"/Data/MODIS/Interpolated Data and plots",sep="")
# not logged
# MODIS ET
ET <- read.csv(paste(ETdir,"/Muttama_mean_smooth.csv",sep=""))
ET$zooDates <- as.Date(ET[,2])

# Raw
ETraw <- read.csv(paste(ETdir,"/Muttama_raw.csv",sep=""))
ETraw$zooDate <-as.Date(ETraw[,2])
ETraw

# linear
linear.et <- read.csv(paste(ETdir,"/Muttama_mean_linear.csv",sep=""))
linear.et$zooDate <- as.Date(linear.et[,2])
#head(linear.et)


# logged
#ETlog <- read.csv(paste(ETdir,"/logmean/Jugiong_logmean_smooth.csv",sep=""))
#ETlog$Dates <-as.Date(ETlog[,2])
# Compare both on one plot
#plot(as.Date(ET$zooDate),ET[,3],type="l",xlab="Date",ylab="ET (mm/day)",
#	main="Muttama catchment")
#lines(as.Date(ETlog$Dates),ETlog[,3],col="blue")
#legend("topright",c("reg. mean","log mean"),col=c(1,"blue"),lty=1)

#if (logged == TRUE) {
#	ET <- ETlog
#	ET$zooDate <- ETlog$Dates
#}


## use zoo to create timeseries (match dates of Q, P, T, aET)
## need to have four columns of Q, P, E and aET
require(zoo)
tsQ <- zoo(Flow$newFlow,Flow$Date,frequency=1)
tsP <- zoo(Rain[, 2],Rain$Date,frequency=1)
tsET <- zoo(ET[, 3],ET$zooDates,frequency=1)
#tsET <- zoo(linear.et$line.interp, linear.et$zooDate, frequency=1)
tsT <- zoo(foo.b, Temp$Date, frequency=1)
Mut <- merge(P=tsP,Q=tsQ,E=tsT, aET=tsET, all=FALSE) # aET is MODIS ET

#xyplot(Mut, main="Muttama")
#xyplot(tsQ, ps=18)




### DEFINE CALIBRATION AND VALIDATION PERIODS ####
data.cal <- window(Mut, start = "2000-01-01",end = "2006-12-31",)
# xyplot(data.cal, main="Muttama Calibration (2000-2006)")
# ok<-complete.cases(data.cal)
# summary(ok) #GOOD few missing values (9)

data.cal.mut <- data.cal


data.val <- window(Mut, start = "2007-01-01", end = "2011-12-31")

#**********************************************************************************************
#-----------------------------
#    SETUP HYDROMAD ####
#-----------------------------

# SETUP HYDROMAD
## reset f parameter
hydromad.options(cmd=list(f=c(0.01,3)))

## make sure return_state carries through when rfit is specified
# 1)
#trace(hydromad:::doRoutingFit,edit=T)
# This should open an editor with the function doRoutingFit in it.
# 2) Replace the line
# object$parlist <- as.list(coef(object, which = "sma", warn = FALSE))
# with
# object$parlist <- as.list(coef(object, which = "sma", warn = 
# FALSE,etc=TRUE))
# 3) Close and save the changes

## trace warnings - to check when SRIV fails
## If SRIV fails early in the iteration it should be ok. https://github.com/josephguillaume/hydromad/issues/4
hydromad.options(trace=TRUE)
options(warn=1)


# ## fix fitting of (1,0) model allowing v_q -> 0
# # insert
# # if(!"v_q" %in% names(model$coefficients)) model$coefficients["v_q"]<-0
# # to the second last line (i.e. above "model")
# trace(expuh.sriv.fit,edit=T)
# trace(expuh.ls.fit,edit=T)
# trace(expuh.inverse.fit, edit=T)
# 





##------------------------
##    MODELLING ####
##------------------------

## specify delay during calibration period
delay <- estimateDelay(data.cal, rises=T, plot=T)


## MODEL 1 - calibrate on Qobs ####
# try model orders indicates little difference between (1,1) and (1,0)
# use (1,0)
hydromad.options(order=c(1,1))
     mutmod.Q <- hydromad(DATA=data.val,
          sma="cmd", d=200, e=c(0,25),#shape=c(0,2),
            routing="expuh",v_s = c(0,1), tau_s = c(3,1000),delay=delay,
            return_state=TRUE)


# # rsq and rsqlog
#        mutfit.Q.r <- fitByOptim(mutmod.Q, samples=500, objective=hmadstat("r.squared"),
# 		 method="PORT")
#        summary(mutfit.Q.r)
#        coef(mutfit.Q.r)
# # agrees with Susie's results
# # rsqlog
#        mutfit.Q.rl <- fitByOptim(mutmod.Q, samples=500, objective=hmadstat("r.sq.log"),
# 		 method="PORT")
#        summary(mutfit.Q.rl)
#        coef(mutfit.Q.rl)
# Agrees with Susie's results
# use Viney's (includes Bias), see http://hydromad.catchment.org/#hydromad.stats
hydromad.stats("viney" = function(Q, X, ...) {
    hmadstat("r.squared")(Q, X, ...) -
      5*(abs(log(1+hmadstat("rel.bias")(Q,X)))^2.5)})

#        mutfit.Q.v <- fitByOptim(mutmod.Q, samples=500, objective=~hmadstat("viney")(Q, X),
# 		 method="PORT")
#        summary(mutfit.Q.v)
#        coef(mutfit.Q.v)
# do this 10 times and save results
Store <- data.frame(rel.bias = numeric(length=10),r.squared=numeric(length=10),
                    r.sq.sqrt=numeric(length=10),r.sq.log=numeric(length=10),
                    f=numeric(length=10),e=numeric(length=10),shape=numeric(length=10),
                    v_s=numeric(length=10),tau_s=numeric(length=10))
for (i in 1:10) {
  mutfit.Q.v <- fitBySCE(mutmod.Q,  objective=~hmadstat("viney")(Q, X),
                         control=list(ncomplex=20))
  
  s <- summary(mutfit.Q.v)
  Store[i,1:4] <- do.call(rbind,s[7:10])[,1]
  Store[i,5:9] <-coef(mutfit.Q.v)[c(1,2,4:6)] 
} 
# now what do we do with Store?
Store$model <- "mutfit.Q.v"
write.csv(Store,"20140806_ResultsMuttamaCrossVal.csv",row.names=F)
#save.image("x:/vervoort/research/misc/current.Rdata")
##-------------------------------------------
#### Do again, but now no shape #### 
# -------------------------------------
#    mutmod.Q <- hydromad(DATA=data.cal,
#                            sma="cmd", d=200, e=c(0,100), routing="expuh",rfit=list("sriv", order=c(1,0),delay=delay),
#                            return_state=TRUE)
# ------------------------------------------------

# # Test more complex model
# mutmod.Q2 <- hydromad(DATA=data.cal,
#                      sma="cmd", d=200, e=c(0,5),shape=c(0,2),
#                      routing="expuh",rfit=list("sriv", order=c(2,1),delay=delay),
#                      return_state=TRUE)
# mutfit.Q.v2 <- fitByOptim(mutmod.Q2, samples=500, objective=~hmadstat("viney")(Q, X),
#                          method="PORT")
# summary(mutfit.Q.v2)
# coef(mutfit.Q.v2)

# define model based on Viney's
mutfit.Q <- hydromad(DATA = data.val, d = 200, sma = "cmd",return_state = TRUE,
                      routing = "expuh",  f = mean(Store[,5]),
				e = mean(Store[,6]), shape=mean(Store[,7]), v_s = mean(Store[,8]),
                     tau_s =  mean(Store[,9]))
summary(mutfit.Q)
coef(mutfit.Q)
hmadstat("viney")(Q=data.val$Q,X=mutfit.Q$fitted.values)
#0.2506236
xyplot(mutfit.Q)
xyplot(mutfit.Q$U)

## MODEL 2 ### 
       mod.5<-hydromad(DATA=data.val,
                       sma="cmd",d=200, e=c(0,5), #shape=c(0,2),
                       routing="expuh", v_s = c(0,1), tau_s = c(3,1000),delay=delay,
                       return_state=TRUE)
w=0.5
for (i in 1:10) {
  mutfit.5.v <- fitBySCE(mod.5, objective=~w*hmadstat("viney")(DATA$aET,U$ET) + 
                           (1-w)*hmadstat("viney")(Q,X),control=list(ncomplex=20))
  
  
  s <- summary(mutfit.5.v)
  Store[i,1:4] <- do.call(rbind,s[7:10])[,1]
  Store[i,5:9] <-coef(mutfit.5.v)[c(1,2,4:6)] 
} 

# Store the results
Store$model <- "mutfit.5.v"
write.table(Store,"20140807_ResultsModel2MuttamaCrossval.csv",row.names=F)



# two solutions
# solution 1, 9 out of 10
mutfit.5a <- hydromad(DATA = data.val, d = 200, return_state = T, sma = "cmd", 
                     routing = "expuh", f =  mean(Store[,5]), 
                      e =  mean(Store[,6]), shape= mean(Store[,7]),
                      v_s = mean(Store[,8]), tau_s=mean(Store[,9]))
summary(mutfit.5a)
coef(mutfit.5a)
w*hmadstat("viney")(data.val$aET,mutfit.5a$U$ET) + 
  (1-w)*hmadstat("viney")(data.val$Q,mutfit.5a$fitted.values)
# -0.05834951
mutfit.5 <- mutfit.5a
xyplot(mutfit.5)
xyplot(mutfit.5a$U)
#Solution 2
mutfit.5b <- hydromad(DATA = data.cal, d = 200, return_state = T, sma = "cmd", 
                      routing = "expuh", f =  mean(Store[1,5]), 
                      e =  mean(Store[1,6]), shape= mean(Store[1,7]),
                      v_s = mean(Store[1,8]), tau_s=mean(Store[1,9]))
summary(mutfit.5b)
coef(mutfit.5b)
w*hmadstat("viney")(data.cal$aET,mutfit.5b$U$ET) + 
  (1-w)*hmadstat("viney")(data.cal$Q,mutfit.5b$fitted.values)
# -01.62292

save.image("x:/vervoort/research/misc/current.Rdata")

## MODEL 3 ####
        # # Calibrate SMA on ET, use neighbours (Jugiong) routing
        # #                 f        e   d shape    tau_s v_s v_q delay
        # # Model 1 0.4499936 0.9794826 200 0.7010715 5.321812   0.8973   0.1027     0
        # 
        # # define SMA
sma.ET.v <- hydromad(DATA=data.cal, sma="cmd", d=200,e=c(0,5), return_state=T)
hydromad.options(objective=~hmadstat("viney")(DATA$aET,U$ET))

#fitsma.ET.v <- fitByOptim(sma.ET.v,samples=500,method="PORT")



for (i in 1:10) {
  fitsma.ET.v <- fitBySCE(sma.ET.v, control=list(ncomplex=20))  
  mod <- update(fitsma.ET.v, routing="expuh", tau_s= 5.322, v_s= 0.8973,
                v_q= 0.1027, delay= 0.0000000)
  
  s <- summary(mod)
  Store[i,1:4] <- do.call(rbind,s[7:10])[,1]
  Store[i,5:9] <-coef(mod)[c(1,2,4:6)] 
} 
Store$model <- "mutfit.ET"
write.table(Store,"20140321_ResultsMuttama.csv",row.names=F,append=T,sep=",",col.names=F)

# #add routing from jugfit.Q
#define mod based on 3 replicates from 3
mutfit.ET <- hydromad(DATA = data.cal, d = 200, return_state = T, sma = "cmd", 
                      f = mean(Store[,5]), e = mean(Store[,6]), shape=mean(Store[,7]),
				tau_s = 5.322, v_s = 0.8973, v_q = 0.1027, 
                      delay = 0, routing = "expuh")
summary(mutfit.ET)
coef(mutfit.ET)
hmadstat("viney")(Q=data.cal$Q,X=mutfit.ET$fitted.values)
# -72.26694
save.image("x:/vervoort/research/misc/current.Rdata")


## MODEL 4 ####
# Jugiong
#f           e           d       shape       tau_s         v_s 
#0.4499608   0.9797882 200.0000000   0.0000000   5.3235623   0.8972692 # use neighbours SMA and routing (mutfit.Q)
pred.mut <- hydromad(DATA = data.cal, d = 200, return_state = TRUE, f = 0.4499608, 
               e = 0.9797882, shape = 0, sma = "cmd", 
			tau_s = 5.3235623, v_s = 0.8972692, v_q = 0, delay = 0, routing = "expuh")
summary(pred.mut)
coef(pred.mut)





##-------------------------
## COMPARE ALL MODS ####
##-------------------------
all.mods.cal <- runlist("Model 1"=mutfit.Q, "Model 2"=mutfit.5, "Model 3"=mutfit.ET, "Model 4" = pred.mut)
round(summary(all.mods.cal),4)

coef(all.mods.cal)

trellis.par.set(strip.background = list(col =
                                          "gray75"),add.text=list(font=2,cex=1.2))

scales=list(cex=1.3, relation="same")


streamflow=xyplot(all.mods.cal, main="Muttama Streamflow Calibration", 
                  scales=scales,
                  ylab=list(label="Streamflow (mm)", cex=1.3),
                  xlab=list(label="Date", cex=1.3),
                  col=c("black","Gray50"),lty=c(1,2),lwd=c(1,2))
scales$y$log <- FALSE
rainfall=xyplot(observed(all.mods.cal[[1]], select = "P", all = FALSE), 
                scales=scales,superpose = TRUE, type = "h",col="black",
                lty=1,lwd=1)

c(streamflow,rainfall=rainfall,layout=c(1,NA),x.same=NA,y.same=NA)


## EVALUATION ####
mutfit.Q.val <- update(mutfit.Q, newdata=data.cal)
summary(mutfit.Q.val)
mutfit.5.val <- update(mutfit.5, newdata=data.cal)
summary(mutfit.5.val)

mutfit.ET.val <- update(mutfit.ET, newdata=data.val)
pred.mut.val <- update(pred.mut, newdata=data.val)
all.mods.val <- runlist("Model 1"=mutfit.Q.val, "Model 2"=mutfit.5.val, "Model 3"=mutfit.ET.val, "Model 4" = pred.mut.val)
round(summary(all.mods.val),4)
coef(all.mods.val)

streamflow=xyplot(all.mods.val, main="Muttama Streamflow Evaluation", 
                  scales=scales,
                  ylab=list(label="Streamflow (mm)", cex=1.3),
                  xlab=list(label="Date", cex=1.3),
                  col=c("black","Gray50"),lty=c(1,2),lwd=c(1,2))
scales$y$log <- FALSE
rainfall=xyplot(observed(all.mods.val[[1]], select = "P", all = FALSE), 
                scales=scales,superpose = TRUE, type = "h",col="black",
                lty=1,lwd=1)

c(streamflow,rainfall=rainfall,layout=c(1,NA),x.same=NA,y.same=NA)


##---------------------
##     SOME PLOTS
##---------------------

##---------------
##      Q
##---------------

## check sum[X]=sum[U]
Q <- sum(mutfit.Q$data$Q, na.rm=T)
X <- sum(mutfit.Q$fitted.values, na.rm=T)
U <- sum(mutfit.Q$U$U, na.rm=T)
as.data.frame(c(Q, X, U), row.names=c("Q", "X", "U"))

Q2 <- sum(mutfit.5$data$Q, na.rm=T)
X2 <- sum(mutfit.5$fitted.values, na.rm=T)
U2 <- sum(mutfit.5$U$U, na.rm=T)
as.data.frame(c(Q2, X2, U2), row.names=c("Q", "X", "U"))

Q3 <- sum(mutfit.ET$data$Q, na.rm=T)
X3 <- sum(mutfit.ET$fitted.values, na.rm=T)
U3 <- sum(mutfit.ET$U$U, na.rm=T)
as.data.frame(c(Q3, X3, U3), row.names=c("Q", "X", "U"))

Q4 <- sum(pred.mut$data$Q, na.rm=T)
X4 <- sum(pred.mut$fitted.values, na.rm=T)
U4 <- sum(pred.mut$U$U, na.rm=T)
as.data.frame(c(Q4, X4, U4), row.names=c("Q", "X", "U"))





## look at U predictions
ee <- mutfit.Q$U$U
ff <- mutfit.5$U$U
gg <- mutfit.ET$U$U
hh <- pred.mut$U$U

plot(gg, main="Modelled U - Muttama (calibration) (rsquared)", col="magenta", 
     xlim=c(as.Date("2002-04-10"),as.Date("2005-12-31")), xlab="Time", ylab="Modelled U (mm)",
     lwd=1)
lines(ff, col="blue", lwd=1)
lines(ee, col="orange", lwd=1)
par(new = TRUE)
plot(data.cal$Q, axes=FALSE, bty = "n", xlim=c(as.Date("2002-04-10"),as.Date("2005-12-31")),
     xlab = "", ylab="", col="black", lwd=1)
axis(side=4, at = pretty(range(data.cal$Q)))
par(mar=c(5, 4, 4, 4) + 0.1)
mtext("Streamflow (mm)",side=4,line=3)
lgd.txt=c("Model 3", "Model 2", "Model 1", "Qobs")
legend("topleft",lgd.txt,lty=c(1,rep(1,2)),lwd=2,
       col=c("magenta","blue","orange", "black"))



## remove warmup period from data
## warmup period is 100 days which is 10-4-2000
data.cal.lim <- window(data.cal, start="2000-04-10", end="2006-12-31")
ETmod1 <- window(mutfit.Q$U$ET, start="2000-04-10", end="2006-12-31")
ETmod2 <- window(mutfit.5$U$ET, start="2000-04-10", end="2006-12-31")
ETmod3 <- window(mutfit.ET$U$ET, start="2000-04-10", end="2006-12-31")
ETmod4 <- window(pred.mut$U$ET, start="2000-04-10", end="2006-12-31")

## set the zero values of ETpred to NA
ETmod1 <- ifelse(ETmod1==0, NA, ETmod1)
ETmod2 <- ifelse(ETmod2==0, NA, ETmod2)
ETmod3 <- ifelse(ETmod3==0, NA, ETmod3)
ETmod4 <- ifelse(ETmod4==0, NA, ETmod4)

## Plot aET vs Etpred
par(mfrow=c(1,3))
plot(data.cal$aET, ETmod1, main="r.squared - Q calibration (Muttama)")
abline(0,1)
plot(data.cal$aET, ETmod2, main="r.squared - w=0.5, Q, ET (Muttama)")
abline(0,1)
plot(data.cal$aET, ETmod3, main="r.squared - ET calibration (Muttama)")
abline(0,1)

require(epiR)
## lins correlation concordance
a<-epi.ccc(data.cal.lim$aET, ETmod1)
a$rho.c #0.5100196
b<-epi.ccc(data.cal.lim$aET, ETmod2)
b$rho.c #0.5831868
c<-epi.ccc(data.cal.lim$aET, ETmod3)
c$rho.c #0.5316669
d<-epi.ccc(data.cal.lim$aET, ETmod4)
d$rho.c #0.5867814

cal <- data.frame(loc=rep("Muttama",4),per=rep("cal",4),
     rbind(as.numeric(a$rho.c),as.numeric(b$rho.c),
	as.numeric(c$rho.c),as.numeric(d$rho.c)))
write.table(cal,file="20140327_mut_rho_values.csv",row.names=F,col.names=F,sep=",")


## Plot aET vs Etpred
par(mfrow=c(1,4), mar=c(4, 4, 2, 0))
plot(data.cal$aET, ETmod1,xlab="ETobs (mm)", ylab="ETpred (mm)",main="Model 1",xlim=c(0,5),ylim=c(0,5))
abline(0,1)
par(mar=c(4, 3, 2, 1))
plot(data.cal$aET, ETmod2, xlab="ETobs (mm)", main="Model 2",xlim=c(0,5),ylim=c(0,5))
abline(0,1)
par(mar=c(4, 3, 2, 1))
plot(data.cal$aET, ETmod3, xlab="ETobs (mm)",main="Model 3",xlim=c(0,5),ylim=c(0,5))
abline(0,1)
par(mar=c(4, 3, 2, 1))
plot(data.cal$aET, ETmod4, xlab="ETobs (mm)",main="Model 4",xlim=c(0,5),ylim=c(0,5))
abline(0,1)

par(mfrow=c(1,3), mar=c(4, 4.5, 2, 0), ps=25)
plot(data.cal$aET, ETmod1,xlab="ETobs (mm)", ylab="ETpred (mm)",main="Model 1",xlim=c(0,5),ylim=c(0,5))
abline(0,1, col="red", lwd=1.5)
par(mar=c(4, 3, 2, 0))
plot(data.cal$aET, ETmod2, xlab="ETobs (mm)", main="Model 2",xlim=c(0,5),ylim=c(0,5))
abline(0,1, col="red", lwd=1.5)
par(mar=c(4, 3, 2, 1))
plot(data.cal$aET, ETmod3, xlab="ETobs (mm)",main="Model 3",xlim=c(0,5),ylim=c(0,5))
abline(0,1, col="red", lwd=1.5)


## evaluation


ETmod1.val <- mutfit.Q.val$U$ET
ETmod2.val <- mutfit.5.val$U$ET
ETmod3.val <- mutfit.ET.val$U$ET
ETmod4.val <- pred.mut.val$U$ET

## set zero in ETpred to NA
ETmod1.val <- ifelse(ETmod1.val==0, NA, ETmod1.val)
ETmod2.val <- ifelse(ETmod2.val==0, NA, ETmod2.val)
ETmod3.val <- ifelse(ETmod3.val==0, NA, ETmod3.val)
ETmod4.val <- ifelse(ETmod4.val==0, NA, ETmod4.val)

a<-epi.ccc(data.val$aET, ETmod1.val)
a$rho.c #0.1476603
b<-epi.ccc(data.val$aET, ETmod2.val)
b$rho.c #0.2364832
c<-epi.ccc(data.val$aET, ETmod3.val)
c$rho.c #0.4361659
d<-epi.ccc(data.val$aET, ETmod4.val)
d$rho.c #0.3161434

val <- data.frame(loc=rep("Muttama",4),per=rep("val",4),
     rbind(as.numeric(a$rho.c),as.numeric(b$rho.c),
	as.numeric(c$rho.c),as.numeric(d$rho.c)))
write.table(val,file="20140327_mut_rho_values.csv",append=T,row.names=F,col.names=F,sep=",")


#++++++++++++++++++++++++++++++++++++++++++++++++++
# Residual analysis
# --------------------------------------------


hmadstat("r.squared")(data.cal$aET,mutfit.Q$U$ET)
hmadstat("r.squared")(data.cal$aET,mutfit.5$U$ET)
hmadstat("r.squared")(data.cal$aET,mutfit.ET$U$ET)
hmadstat("rel.bias")(data.cal$aET,mutfit.Q$U$ET)
hmadstat("rel.bias")(data.cal$aET,mutfit.5$U$ET)
hmadstat("rel.bias")(data.cal$aET,mutfit.ET$U$ET)
hmadstat("r.sq.log")(data.cal$aET,mutfit.Q$U$ET)
hmadstat("r.sq.log")(data.cal$aET,mutfit.5$U$ET)
hmadstat("r.sq.log")(data.cal$aET,mutfit.ET$U$ET)



# first in calibration period, more important for model structure
Qres.mod1 <- data.cal$Q - mutfit.Q$fitted.values
Qres.mod2 <- data.cal$Q - mutfit.5$fitted.values
Qres.mod3 <- data.cal$Q - mutfit.ET$fitted.values
Qres.mod4 <- data.cal$Q - pred.mut$fitted.values

# make standard model error plots
# Function to make plots
erplot <- function(data,predicted,name) {
  par(mfrow=c(2,2))
  resid = data - predicted
  plot(density(na.omit(as.numeric(resid))), 
       main=paste(name,"density of residuals"))
  plot(data,resid,pch=16, col="red", 
       xlab="observed", ylab="Residual", main=name)
  qqplot(data,predicted, main=paste("Q-Q plot",name), 
         xlab="Observed", ylab="Predicted", pch=16, col="blue")
  qqnorm(resid,main= "Normal Q-Q Plot residuals")
  par(mfrow=c(1,1))
  
}

# standard error plots
erplot(data.cal$Q,mutfit.Q$fitted.values,"Model 1")
erplot(data.cal$Q,mutfit.5$fitted.values,"Model 2")
erplot(data.cal$Q,mutfit.ET$fitted.values,"Model 3")


# Now plot the ET residuals
# first in calibration period, more important for model structure
ETres.mod1 <- data.cal$aET - mutfit.Q$U$ET
ETres.mod2 <- data.cal$aET - mutfit.5$U$ET
ETres.mod3 <- data.cal$aET - mutfit.ET$U$ET
ETres.mod4 <- data.cal$aET - pred.mut$U$ET

# make standard model error plots
erplot(data.cal$aET,mutfit.Q$U$ET,"Model 1")
erplot(data.cal$aET,mutfit.5$U$ET,"Model 2")
erplot(data.cal$aET,mutfit.ET$U$ET,"Model 3")

# Combined pairs plot
com.df <- cbind(data.cal$Q,data.cal$aET,Qres.mod1,ETres.mod1,
                Qres.mod2,ETres.mod2,Qres.mod3,ETres.mod3)
pairs(com.df)

## Now do some plotting
plot(ETres.mod1, col="black", ylim=c(-4,6), main="residuals ET prediction",
     xlab="Date", ylab="Residuals and ET estimate")#, xlim=c(as.Date("2000-04-10"),as.Date("2004-12-31")))
lines(ETres.mod2, col="red")
lines(ETres.mod3, col="green")
#lines(ETres.mod4, col="blue")
lines(data.cal$aET, col="blue")
lines(0.15*data.cal$E, col="pink")

lgd.txt=c("Model 1", "Model 2", "Model 3",  "observed (MODIS)", "model input~0.15*MaxT")
legend("topleft",lgd.txt,lty=1,lwd=1,
       col=c("black","red","green", "blue","pink"))
