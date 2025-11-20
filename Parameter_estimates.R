
###### R Script for PNAS publication by Croll et. al. 2025 ######
# This R script fits the functions for temperature dependence and temperature 
# time series used for the model simulations.
#
# Output: Figure S1-S4

###### SETTINGS AND PACKAGES ######

# device specific settings

# load packages
library(rTPC)
library(ggplot2)
library(lubridate)
library(hms)
library(heatwaveR)
library(ggpubr)

Sys.setlocale("LC_ALL", "English")

# Arrhenius function
Arrfunc <- function(temp, sc, pA){
  k = 8.617*10^{-5}
  Tr = 8.588
  val <- sc*exp(-pA/k*(1/(temp+273.15)-1/(Tr+273.15)))
  val[which(is.infinite(val))] <- 10E10
  return(val)
}


###### PARAMETER TEMPERATURE DEPENDENCE #####


### RESOURCE PRODUCTIVITY ###
# Based on Copepod egg productivity fitted to graph data from Holste & Peck (2005)

holste2005 <- read.csv("Data/Holste2005.csv")
names(holste2005) <- c("temp","EP")

start_vals <- get_start_vals(holste2005$temp, holste2005$EP, model_name =
                               'sharpeschoolhigh_1981')

holste2005fit <- nls(EP~sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref=8.588 ), 
                     data = holste2005, start= start_vals, algorithm="port",
                     control=list(tol=1e-10,minFactor=1e-20), trace = TRUE)

holste2005$pred <- predict(holste2005fit)

holste2005_plot <- ggplot(data=holste2005,aes(x=temp))+
  geom_point(aes(y=EP))+
  geom_line(aes(y=pred), colour="red")+
  theme_classic()+
  labs(x="Temperature", y="Resource egg production")


### RESOURCE GENERATION TIME ###
## Zooplankton generation time approximation of the relationship fitted by Huntly & Lopez (1992)

huntly1992func <- function(temp){
  gentime <- 128.8*exp(-0.120*temp)
  return(gentime);
}

huntly1992diff <- function(x){
  sc <- x[1]
  pA <- x[2]
  
  val <- integrate(function(temp){ (huntly1992func(temp) - Arrfunc(temp, sc, pA) )^2  }, lower = 0, upper= 35)
  return(val$value)
}

huntly1992fit <- nlm(huntly1992diff,c(50,-0.8),print.level=2, steptol = 1e-10, iterlim=1000, stepmax = 10)

huntly1992predict <- data.frame(temp=seq(0,35,0.01), 
                                     huntly = huntly1992func(seq(0,35,0.01)), 
                                     Arr=Arrfunc(seq(0,35,0.01), huntly1992fit$estimate[1], huntly1992fit$estimate[2]))

huntly1992_plot <- ggplot(data=huntly1992predict)+
  geom_line(aes(x=temp, y=huntly))+
  geom_line(aes(x=temp, y=Arr), colour="red")+
  labs(x="Temperature", y="Resource generation time")+
  theme_classic()

### HERRING GASTRIC EVACUATION ###
# Gastric evacuation approximation of the relationship fitted by Bernreuter et al 2009

bernreuter2009func <- function(temp){
  GE <- 0.0177*exp(0.0775*temp)*(1-1/(1+exp(-0.659*(temp-23.989))))
  return(GE);
}


bernreuter2009diff <- function(x){
  sc <- x[1]
  pA <- x[2]
  pD <- x[3]
  pH <- x[4]
  Tr = 8.588
  
  val <- integrate(function(temp){ (bernreuter2009func(temp) - sharpeschoolhigh_1981(temp=temp, r_tref=sc, e=pA, eh=pD, th=pH, tref=Tr) )^2  }, lower = 5, upper= 25)
  
  return(val$value)
}

bernreuter2009fit <- nlm(bernreuter2009diff,c(1,1,1,22),print.level=2, steptol = 1e-10, iterlim=1000, stepmax = 10)

bernreuter2009predict <- data.frame(temp=seq(5,25,0.01), 
                                    bernreuter = bernreuter2009func(seq(5,25,0.01)), 
                                    SS=sharpeschoolhigh_1981(seq(5,25,0.01), bernreuter2009fit$estimate[1], bernreuter2009fit$estimate[2], bernreuter2009fit$estimate[3], bernreuter2009fit$estimate[4], tref=8.588)) 

bernreuter2009_plot<- ggplot(data=bernreuter2009predict)+
  geom_line(aes(x=temp, y=bernreuter))+
  geom_line(aes(x=temp, y=SS), colour="red")+
  labs(x="Temperature", y="Herring gastric evacuation")+
  theme_classic()


### HERRING HATCHING TIME ###
## Hatching time approximation of the relationship fitted by Peck et al 2012

peck2012HTfunc <- function(temp){
  HT <- 4461.9*temp^(-1.232)
  return(HT);
}


peck2012HTdiff <- function(x){
  sc <- x[1]
  pA <- x[2]
  HTmin <- x[3]
  Tr = 8.588
  
  val <- integrate(function(temp){ (peck2012HTfunc(temp) - HTmin - Arrfunc(temp, sc, pA) )^2  }, lower = 5, upper= 20)
  
  return(val$value)
}

peck2012HTfit <- nlm(peck2012HTdiff,c(218.152104,  -1.530587, 101.441589),print.level=2, steptol = 1e-10, iterlim=1000, stepmax = 100)

peck2012HTpredict <- data.frame(temp=seq(5,20,0.01), 
                                peck = peck2012HTfunc(seq(5,20,0.01)), 
                                Arr= peck2012HTfit$estimate[3]+ Arrfunc(seq(5,20,0.01), peck2012HTfit$estimate[1] , peck2012HTfit$estimate[2] ))

peck2012HT_plot <- ggplot(data=peck2012HTpredict)+
  geom_line(aes(x=temp, y=peck))+
  geom_line(aes(x=temp, y=Arr), colour="red")+
  labs(x="Temperature", y="Herring hatching time")+
  theme_classic()


### HERRING EGG SURVIVAL ###
# Egg survival fitted to graph data from Peck et al 2012

peck2012ES <- read.csv("Data/peck2012survival.csv")

names(peck2012ES) <- c("temp","ES")


start_vals <- get_start_vals(peck2012ES$temp, peck2012ES$ES, model_name =
                               'sharpeschoolhigh_1981')

peck2012ESfit <- nls(ES~sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref=8.588 ), 
                   data = peck2012ES, start= c(r_tref=10, e=0.0266, eh=2, th=2), algorithm="port",
                   control=list(tol=1e-10,minFactor=1e-20), trace = TRUE)

peck2012ESpred<-data.frame(temp = seq(3,22,0.5))
peck2012ESpred$pred <- predict(peck2012ESfit,newdata = peck2012ESpred )

peck2012ES_plot <- ggplot(data=peck2012ES, aes(x=temp))+
  geom_point(aes(y=ES))+
  geom_line(data = peck2012ESpred,aes(y=pred), colour="red")+
  labs(x="Temperature", y="Hering egg survival")+
  theme_classic()



### FIGURE ###

tempres_plot  <- ggarrange(holste2005_plot, huntly1992_plot, bernreuter2009_plot, peck2012HT_plot, peck2012ES_plot, ncol=2, nrow=3, align="v",legend ="none", 
                                       labels=c("a", "b", "c", "d", "e"))

ggsave("Figures/FigureS2.pdf", tempres_plot, width=15, height=23, units="cm")



###### TEMPERATURE TIMESERIES #####

### Load data
SMHI33002 <- read.table(file="Data/SMHI_tempdata_33002.txt", sep=";", skip= 6, header=TRUE, fill=TRUE,quote = "")

SMHI33002 <- SMHI33002[,c(1,2,4)] # select appropriate columns
names(SMHI33002) <- c("textdate","temp","depth")

SMHI33002 <- SMHI33002[which(SMHI33002$depth==20),]

# convert dates in various formats for computing
SMHI33002$ymdhms <- ymd_hms(SMHI33002$textdate)
SMHI33002$ymd <- date(SMHI33002$ymdhms)
SMHI33002$hms <- as_hms(SMHI33002$ymdhms)
SMHI33002$doy <- yday(SMHI33002$ymd)
SMHI33002$year <- year(SMHI33002$ymd)


### AGGREGATE DATA PER DAY ###

daydata33002 <- reshape(SMHI33002, v.names = "temp", idvar=c("ymd","depth","year"), drop=c("ymdhms","textdate"), timevar=c("hms"),direction="wide")

daydata33002$avgtemp <- rowMeans(daydata33002[,which(startsWith(names(daydata33002),"temp."))],na.rm=TRUE)
daydata33002 <- daydata33002[,-which(startsWith(names(daydata33002),"temp."))]

#### TRY VARIOUS PERIODIC MODELS ####

testmodel1 <- nls(avgtemp~avgT+amp*sin(2*pi/365*(doy-shift)), 
                  data=daydata33002,
                  start=list(avgT=8, amp=8, shift=160),
                  trace=TRUE)

testmodel2 <- nls(avgtemp~avgT+ampT*sin( 2*pi/365*(doy-shiftT+ampS*sin(2*pi/365*(doy-shiftS)))), 
                  data=daydata33002,
                  start=list(avgT=8, ampT=8, shiftT=160, ampS=45, shiftS=5),
                  trace=TRUE,
                  control=list(maxiter=200))

testmodel3 <- nls(avgtemp~avgT+(ampT+ampA*sin(2*pi/365*(doy-shiftA)))*sin( 2*pi/365*(doy-shiftT)), 
                  data=daydata33002,
                  start=list(avgT=8, ampT=8, shiftT=160, ampA=1, shiftA=160),
                  trace=TRUE,
                  control=list(maxiter=200))

testmodel4 <- nls(avgtemp~avgT+(ampT+ampA*sin(2*pi/365*(doy-shiftA)))*sin( 2*pi/365*(doy-(shiftT+ampS*sin(2*pi/365*(doy-shiftS))))), 
                  data=daydata33002,
                  start=list(avgT=8, ampT=8, shiftT=160, ampA=1, shiftA=160, ampS=45, shiftS=5),
                  trace=TRUE,
                  control=list(maxiter=200, minFactor=10^-5))

modelpredict1 <- data.frame(doy =1:365)
modelpredict1$predtemp <- predict(testmodel1, newdata = modelpredict1)

modelpredict2 <- data.frame(doy =1:365)
modelpredict2$predtemp <- predict(testmodel2, newdata = modelpredict2)

modelpredict3 <- data.frame(doy =1:365)
modelpredict3$predtemp <- predict(testmodel3, newdata = modelpredict3)

modelpredict4 <- data.frame(doy =1:365)
modelpredict4$predtemp <- predict(testmodel4, newdata = modelpredict4)

AIC(testmodel1,testmodel2, testmodel3,testmodel4)

temperture_plot <- ggplot(data=daydata33002,aes(x=doy,y=avgtemp))+
  geom_point(size=0.2)+
  geom_line(data=modelpredict1,aes(y=predtemp, colour="\nSine (AIC: 16112)\n"), size=1)+
  geom_line(data=modelpredict2,aes(y=predtemp, colour="\nSine\n+periodic shift (AIC: 16016)\n"), size=1)+
  geom_line(data=modelpredict3,aes(y=predtemp, colour="\nSine\n+periodic amplitude (AIC: 16027)\n"), size=1)+
  geom_line(data=modelpredict4,aes(y=predtemp, colour="\nSine\n+periodic shift\n+periodic amplitude (AIC: 15866)\n"), size=1)+
  theme_classic()+
  labs(x="Julian date", y="Daily temperature")+
  theme(legend.title = element_blank())

ggsave("Figures/FigureS1.pdf", plot=temperture_plot, width=15, height=7.5, units="cm")



### HEATWAVE PATTERNS

## create a data frame with all dates in the timeseries
alldays33002 <- data.frame(t=seq(from=min(daydata33002$ymd), to=max(daydata33002$ymd),by=1))

# add temperatures
alldays33002 <- merge(alldays33002, daydata33002, by.x="t", by.y="ymd", all.x=TRUE)
names(alldays33002) <- c("t","depth","doy","year","temp")

## Extract heatwaves

# Calculate climate

cl33002 <- ts2clm(data=alldays33002, climatologyPeriod = c(min(alldays33002$t),max(alldays33002$t)), pctile = 90)

# detect heatwaves
hw33002 <- detect_event(data=cl33002)

hwintensity_plot <- ggplot(data=hw33002$event, aes(x=yday(date_peak), y=intensity_mean))+
  geom_errorbar(aes(ymin=intensity_mean, ymax=intensity_max), colour="grey")+
  geom_errorbarh(aes(xmin=yday(date_start), xmax=yday(date_end)), colour="grey")+
  geom_point()+
  theme_classic()+
  labs(x="Julian date", y="Heatwave intensity (degrees)")+
  scale_y_continuous(limits=c(0,7))+
  scale_x_continuous(limits=c(0,365))

hwduration_plot <- ggplot(data=hw33002$event, aes(x=yday(date_peak), y=duration, size=intensity_mean))+
  geom_point()+
  theme_classic()+
  labs(x="Julian date", y="Heatwave duration (days)")+
  scale_y_continuous(limits=c(0,50))+
  scale_x_continuous(limits=c(0,365))+
  theme(legend.position="none")

ggsave("Figures/FigureS3.pdf", plot=hwintensity_plot, width=15, height=7.5, units="cm")
ggsave("Figures/FigureS4.pdf", plot=hwduration_plot, width=15, height=7.5, units="cm")










