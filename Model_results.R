
###### R Script for Ecology Letters publication by Croll et. al. 2026 ######
# This script contains the R-code to import and process the data form the model
# simulations and produce the figures corresponding figures
#
# Output: Figure 2-4 and S5-S9

###### SETTINGS AND PACKAGES ######

# device specific settings

#  Load needed packages
library(ggplot2)
library(PSPManalysis)
library(lubridate)
library(dplyr)
library(ggpubr)
library(ggnewscale)
library(scales)
library(patchwork)

source("csbtxtread.R")

Sys.setlocale("LC_ALL", "English")


###### MAIN BIFURCATION GRAPHS (FIG 2) ######

### LOAD SIMULATION DATA ###

# dataframe with varied variables for names of files
HWbif_vars <- expand.grid(DUR = c(5, 10, 20, 40), INC=1:6, KEEP.OUT.ATTRS = FALSE)

# load simulation data into list
HWbif_list <- apply(HWbif_vars, 1, function(x){
  name <- paste0("Model/HWbif/herringEBT_",x[1],"_",x[2],".out")
  bifdata <- read.delim(name,
             sep = "\t",
             header=FALSE,
             col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "HWSTART2"))
  return(as.data.frame(bifdata))
  } )

# combine dataframes
HWbif <- data.frame(do.call("rbind",HWbif_list))

# calculate actual time
HWbif$reltime <-  HWbif$time%%36500+125

# calculate the number of years after the heatwave
HWbif$hwyear <- (HWbif$reltime-HWbif$HWSTART )/365

# calculate total harvested biomass
HWbif$harvtot <- rowSums(HWbif[,17:56])

# Extract special points that need to be indicated in graph
HWSTART_pointdatasubset <- HWbif[ which( HWbif$HWDUR==10    & 
                                           (  (HWbif$HWINC==5 & match(HWbif$HWSTART, c(29340, 29363,29444,29232)))  | 
                                                (   (HWbif$HWINC==2 & HWbif$HWSTART==29363) ) 
                                           ) ),]

### LOAD SIMULATIONS FOR TEMPERATURE CURVE ###

# load baseline data without heatwaves
HWts_env_base <- read.delim(paste0("Model/HWts/herringEBT_base.out"),
                            sep = "\t",
                            header=FALSE,
                            col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harvmass",0:39)))

# load runs with heatwaves

# list with values needed for analysis per run
HWts_vars <- data.frame(RUNNAME = c("10_5_140","10_5_163", "10_5_244", "10_5_397", "10_2_163"),
                        HWSTART = c(29340, 29363, 29444, 29597, 29363 ),
                        HWINC = c(5, 5, 5, 5, 2),
                        HWDUR = c(10, 10, 10, 10, 10),
                        POSTHWTIME = c(29355-125, 29378-125, 29459-125, 29612-125, 29378-125 ),
                        PREREPTIME = 29689-125,
                        PERSTART = 29325,
                        PEREND = 29690)


# read environment output 
HWts_env_list <- apply(HWts_vars, 1, function(x){
  envfilename <- paste0("Model/HWts/herringEBT_",x[1],".out")
  envdata <- read.delim(envfilename,
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harvmass",0:39)))
  
  envdata$reltime <- envdata$time  + 125
  envdata <- cbind("base"=HWts_env_base,"hw"=envdata)
  envdata$HWDATE <- x[2]
  envdata$HWINC <- x[3]
  
  return(as.data.frame(envdata))
} )

HWts_env <- data.frame(do.call("rbind",HWts_env_list))


### FIGURE 2a ###

HWSSTART_temp_plot <-ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black", size=1.2)+
  geom_line(data=HWts_env, aes(x=as.Date(hw.reltime-29200), y=hw.temp, colour=as.factor(hw.HWSTART), linetype=as.factor(hw.HWINC)))+
  geom_line(data=HWts_env, aes(x=as.Date(hw.reltime-29200), y=base.temp))+
  theme_classic()+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_y_continuous(limits=c(0, 20))+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  labs(x="Date", y="Temperature (\u00B0C)")+
  geom_text(aes(x=as.Date(163.5), y=18.5), label="Spawning\nseason", size=3.5)+
  geom_text(aes(x=as.Date(125+365), y=18.5), label="Measuring\nmoment", size=3.5)+
  theme(legend.position = "none")


### FIGURE 2b ###

HWSTART_INCabs_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 9.070128E+08, linetype="dashed")+
  geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif[which(HWbif$HWDUR==10   & HWbif$hwyear >=0  &  HWbif$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=eggbuff, colour=as.factor(HWINC)), linewidth=0.5)+
  scale_x_date(date_labels = "%b %d", limits=as.Date(c(124.5,125.5+365)) )+
  scale_y_continuous( limits=c( 8.7E8, 9.5E8 ) )+
  scale_colour_manual(values=c("grey60","grey50","grey40","grey30","grey20","grey10"))+
  labs(x="Midpoint of the heatwave", y="Total egg production", colour="Temperature increase during heatwave")+
  geom_text(aes(x=as.Date(163.5), y=9.5E8), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  geom_point(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear >=0  &  HWSTART_pointdatasubset$hwyear <1),], 
             aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=eggbuff, colour=as.factor((HWSTART-125)%%365), shape=as.factor(HWINC)) )+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(unique(c(29340, 29363,29444,29232))),"%b %d"), guide="none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(legend.position = "top")

HWSTART_INCrel_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif[which(HWbif$HWDUR==10   & HWbif$hwyear >=0  &  HWbif$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y= (eggbuff-9.070128E+08)/9.070128E+08*100, colour=as.factor(HWINC)), linewidth=0.5)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_y_continuous( limits=c( -3, 5 ) )+
  scale_colour_manual(values=c("grey60","grey50","grey40","grey30","grey20","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Temperature increase during heatwave (\u00B0C)")+
  #geom_text(aes(x=as.Date(163.5), y=5), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  geom_point(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear >=0  &  HWSTART_pointdatasubset$hwyear <1),], 
             aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=(eggbuff-9.070128E+08)/9.070128E+08*100, colour=as.factor((HWSTART-125)%%365), shape=as.factor(HWINC)) )+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(unique(c(29340, 29363,29444,29232))),"%b %d"), guide="none")+
  scale_shape_manual(name="Heatwave increase (\u00B0C)", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5), axis.text.x = element_blank())


### FIGURE 2c ###

HWSTART_DURabs_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 9.070128E+08, linetype="dashed")+
  geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif[which(HWbif$HWINC==5   & HWbif$hwyear >=0  &  HWbif$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=eggbuff, colour=as.factor(HWDUR)), linewidth=0.5)+
  scale_y_continuous( limits=c( 8.2E8, 10.1E8 ) )+
  scale_x_date(date_labels = "%b %d", limits=as.Date(c(124.5,125.5+365)) )+
  scale_colour_manual(values=c("grey70","grey50","grey30","grey10"))+
  labs(x="Midpoint of the heatwave", y="Total egg production", colour="Duration of the heatwave (days)")+
  geom_text(aes(x=as.Date(163.5), y=10.1E8), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  geom_point(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear >=0  &  HWSTART_pointdatasubset$hwyear <1 & HWSTART_pointdatasubset$HWINC == 5),], aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=eggbuff, colour=as.factor((HWSTART-125)%%365)) )+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(unique(c(29340, 29363,29444,29232))),"%b %d"), guide="none")+
  theme_classic()+
  theme(legend.position = "top")

HWSTART_DURrel_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif[which(HWbif$HWINC==5   & HWbif$hwyear >=0  &  HWbif$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=(eggbuff-9.070128E+08)/9.070128E+08*100, colour=as.factor(HWDUR)), linewidth=0.5)+
  scale_y_continuous( limits=c( -10, 13) )+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_colour_manual(values=c("grey70","grey50","grey30","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Duration of heatwave (days)")+
  #geom_text(aes(x=as.Date(163.5), y=13), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  geom_point(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear >=0  &  HWSTART_pointdatasubset$hwyear <1 & HWSTART_pointdatasubset$HWINC == 5),], aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=(eggbuff-9.070128E+08)/9.070128E+08*100, colour=as.factor((HWSTART-125)%%365)) )+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(unique(c(29340, 29363,29444,29232))),"%b %d"), guide="none")+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5))


HWSTARTrel_plot <- wrap_plots(list(HWSSTART_temp_plot, plot_spacer(),
                                   HWSTART_INCrel_plot, plot_spacer() , 
                                   HWSTART_DURrel_plot), 
                              ncol = 1 , 
                              axis_titles = "collect_x", 
                              tag_level = 'new', 
                              heights=c(0.6, -0.1, 1,-0.1,1))+plot_annotation(tag_levels="a")

ggsave("Figures/Figure2.pdf", HWSTARTrel_plot, width=15, height=20, units="cm")



### FIGURE 3 ###




HWlt_egg_plot <- ggplot(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear>=0  ),] )+
  geom_line(aes(x=floor(((time+125)%%36500)/365-floor((HWSTART-125)/365))-1, y=(eggbuff-9.070128E+08)/9.070128E+08*100, colour=as.factor((HWSTART-125)%%365), linetype=as.factor(HWINC)))+
  theme_classic()+  
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(name="% deviation in\ntotal egg production")+
  scale_x_continuous(limits=c(0,12), name="Years after heatwave", n.breaks = 15)+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="top", 
        axis.title.x =element_blank(),
        axis.text.x=element_blank())

HWlt_zoop_plot<- ggplot(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear>=0  ),] )+
  geom_line(aes(x=floor(((time+125)%%36500)/365-floor((HWSTART-125)/365))-1, y=(zoop-0.1689956)/0.1689956*100, colour=as.factor((HWSTART-125)%%365), linetype=as.factor(HWINC)))+
  theme_classic()+    
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  geom_hline(yintercept =0, linetype="dashed")+
  scale_y_continuous(name="% deviation in\nzooplankton")+
  scale_x_continuous(limits=c(0,12), name="Years after heatwave", n.breaks = 15)+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="top",
        axis.title.x =element_blank(),
        axis.text.x=element_blank())


HWlt_harv_plot <- ggplot(data=HWSTART_pointdatasubset[which(HWSTART_pointdatasubset$hwyear>=0  ),] )+
  geom_line(aes(x=floor(((time+125)%%36500)/365-floor((HWSTART-125)/365))-1, y=(harvtot-1163241)/1163241*100, colour=as.factor((HWSTART-125)%%365), linetype=as.factor(HWINC)))+
  theme_classic()+    
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(name="% deviation in\nyield")+
  scale_x_continuous(limits=c(0,12), name="Years after heatwave", n.breaks = 15)+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="top")


HWlt_plot <- wrap_plots(list(guide_area(), HWlt_egg_plot, plot_spacer(), HWlt_zoop_plot, plot_spacer(),  HWlt_harv_plot), ncol = 1 , guides="collect", tag_level = 'new', heights = c(0.1, 1, -0.2, 1, -0.2, 1))+plot_annotation(tag_levels="a")

ggsave("Figures/Figure3.pdf", HWlt_plot, width=15, height=15, units="cm")



###### INDIVIDUAL HEATWAVE RUNS (FIG 4) ######


## LOAD SIMULATED DATA ##

# load baseline data without heatwaves
HWts_env_base <- read.delim(paste0("Model/HWts/herringEBT_base.out"),
                            sep = "\t",
                            header=FALSE,
                            col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harvmass",0:39)))

# load runs with heatwaves

# list with values needed for analysis per run
HWts_vars <- data.frame(RUNNAME = c("10_5_140","10_5_163", "10_5_244", "10_5_397", "10_2_163"),
                        HWSTART = c(29340, 29363, 29444, 29597, 29363 ),
                        HWINC = c(5, 5, 5, 5, 2),
                        HWDUR = c(10, 10, 10, 10, 10),
                        POSTHWTIME = c(29355-125, 29378-125, 29459-125, 29612-125, 29378-125 ),
                        PREREPTIME = 29689-125,
                        PERSTART = 29325,
                        PEREND = 29690)


# read environment output 
HWts_env_list <- apply(HWts_vars, 1, function(x){
  envfilename <- paste0("Model/HWts/herringEBT_",x[1],".out")
  envdata <- read.delim(envfilename,
                              sep = "\t",
                              header=FALSE,
                              col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harvmass",0:39)))
  
  envdata$reltime <- envdata$time  + 125
  envdata <- cbind("base"=HWts_env_base,"hw"=envdata)
  envdata$HWDATE <- x[2]
  envdata$HWINC <- x[3]
  
  return(as.data.frame(envdata))
  } )

HWts_env <- data.frame(do.call("rbind",HWts_env_list))

# read structure for posthw plot
posthw_struc_list <- apply(HWts_vars, 1, function(x){
  strucfilename <- paste0("Model/HWts/herringEBT_",x[1],".csb")
  
  hwstruc <- csbread(strucfilename, state = as.integer(x[5])+1)
  basestruc <- csbread("Model/HWts/herringEBT_base.csb", state = as.integer(x[5])+1) #unfortunately this package/function is very unstable and sometimes crashed R-studio
  
  hwstruc <- as.data.frame(hwstruc[[3]])
  basestruc <- as.data.frame(basestruc[[3]])
  
  names(hwstruc) <- c("dens","age","indx","indy","mass","len","cond","btime","bage","matage")
  names(basestruc) <- c("dens","age","indx","indy","mass","len","cond","btime","bage","matage")
  
  hwstruc$coho <- floor(hwstruc$btime/365)
  hwstruc$totage <- hwstruc$age + hwstruc$bage
  
  basestruc$coho <- floor(basestruc$btime/365)
  basestruc$totage <- basestruc$age + basestruc$bage
  
  hwstruc_summary <- summarise(.data = hwstruc, 
                                     totdens = sum(dens),
                                     maxcond=max(cond), 
                                     mincond=min(cond),
                                     meancond=weighted.mean(cond,dens),
                                     maxlen=max(len), 
                                     minlen=min(len),
                                     meanlen = weighted.mean(len,dens),
                                     maxage=max(totage), 
                                     minage=min(totage),
                                     meanage = weighted.mean(totage,dens),
                                     .by=c("coho"))
  
  basestruc_summary <- summarise(.data = basestruc, 
                               totdens = sum(dens),
                               maxcond=max(cond), 
                               mincond=min(cond),
                               meancond=weighted.mean(cond,dens),
                               maxlen=max(len), 
                               minlen=min(len),
                               meanlen = weighted.mean(len,dens),
                               maxage=max(totage), 
                               minage=min(totage),
                               meanage = weighted.mean(totage,dens),
                               .by=c("coho"))
  
  struc_summary <- cbind("base"=basestruc_summary ,"hw"=hwstruc_summary )
  
  struc_summary$HWDATE <- x[[2]]
  struc_summary$HWINC <- x[[3]]
  
  return(struc_summary)
  
} )

posthw_struc <- data.frame(do.call("rbind",posthw_struc_list))



# read structure for prereproduction plot
prerep_struc_list <- apply(HWts_vars, 1, function(x){
  strucfilename <- paste0("Model/HWts/herringEBT_",x[1],".csb")
  
  hwstruc <- csbread(strucfilename, state = as.integer(x[6])+1)
  basestruc <- csbread("Model/HWts/herringEBT_base.csb", state = as.integer(x[6])+1)
  
  hwstruc <- as.data.frame(hwstruc[[3]])
  basestruc <- as.data.frame(basestruc[[3]])
  
  names(hwstruc) <- c("dens","age","indx","indy","mass","len","cond","btime","bage","matage")
  names(basestruc) <- c("dens","age","indx","indy","mass","len","cond","btime","bage","matage")
  
  hwstruc$coho <- floor(hwstruc$btime/365)
  hwstruc$totage <- hwstruc$age + hwstruc$bage
  
  basestruc$coho <- floor(basestruc$btime/365)
  basestruc$totage <- basestruc$age + basestruc$bage
  
  hwstruc_summary <- summarise(.data = hwstruc, 
                               totdens = sum(dens),
                               maxcond=max(cond), 
                               mincond=min(cond),
                               meancond=weighted.mean(cond,dens),
                               maxlen=max(len), 
                               minlen=min(len),
                               meanlen = weighted.mean(len,dens),
                               maxage=max(totage), 
                               minage=min(totage),
                               meanage = weighted.mean(totage,dens),
                               .by=c("coho"))
  
  basestruc_summary <- summarise(.data = basestruc, 
                                 totdens = sum(dens),
                                 maxcond=max(cond), 
                                 mincond=min(cond),
                                 meancond=weighted.mean(cond,dens),
                                 maxlen=max(len), 
                                 minlen=min(len),
                                 meanlen = weighted.mean(len,dens),
                                 maxage=max(totage), 
                                 minage=min(totage),
                                 meanage = weighted.mean(totage,dens),
                                 .by=c("coho"))
  
  struc_summary <- cbind("base"=basestruc_summary ,"hw"=hwstruc_summary )
  
  struc_summary$HWDATE <- x[[2]]
  struc_summary$HWINC <- x[[3]]
  
  return(struc_summary)
  
} )

prerep_struc <- data.frame(do.call("rbind",prerep_struc_list))


## FIGURE 4 ##



zoop_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect(data=data.frame(HWDATE = as.numeric(unique(HWts_env$HWDATE))), aes(xmin=as.Date(HWDATE-1-29200), xmax=as.Date(HWDATE+10-1-29200), ymin=-Inf, ymax=Inf, fill=as.factor(HWDATE)), alpha = 0.2  )+
  geom_vline(data=data.frame(HWDATE = as.numeric(unique(HWts_env$HWDATE))), aes(xintercept=as.Date(HWDATE+15-29200), colour=as.factor(HWDATE)), linetype="dotted")+
  geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWts_env, aes(x=as.Date(hw.reltime-29200), y=(hw.zoop-base.zoop)/base.zoop*100, colour=as.factor(HWDATE), linetype = as.factor(HWINC)))+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_date( limits = as.Date(c(124.5,163+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-25, 400), "% deviation in\nzooplankton", breaks=c(-25,0,25,50,100))+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"), guide="none")+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"), guide="none")+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(163.5), y=130), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


dens_posthw_plot<- ggplot()+
  geom_line(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.totdens-base.totdens)/base.totdens*100, colour=as.factor(HWDATE),linetype = as.factor(HWINC)))+
  geom_point(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.totdens-base.totdens)/base.totdens*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-25, 10), "% deviation in\ncohort density")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(title="Directly after the heatwave")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "top",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_posthw_plot <- ggplot()+
  geom_line(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.meanlen-base.meanlen)/base.meanlen*100, colour=as.factor(HWDATE),linetype = as.factor(HWINC)))+
  geom_point(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.meanlen-base.meanlen)/base.meanlen*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-0.1, 6), "% deviation in\n individual length")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "top",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


cond_posthw_plot <- ggplot()+
  geom_line(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.meancond-base.meancond)/base.meancond*100, colour=as.factor(HWDATE),linetype = as.factor(HWINC)))+
  geom_point(data=posthw_struc, aes(x=hw.meanage/365, y=(hw.meancond-base.meancond)/base.meancond*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 5), limits=c(-3, 15), "% deviation in\nindividual condition")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "top",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))




dens_prerep_plot<- ggplot()+
  geom_line(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.totdens-base.totdens)/base.totdens*100, colour=as.factor(HWDATE), linetype = as.factor(HWINC)))+
  geom_point(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.totdens-base.totdens)/base.totdens*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-25, 10), "% deviation in\ncohort density")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(title="Before following spawning season")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))

length_prerep_plot <- ggplot()+
  geom_line(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.meanlen-base.meanlen)/base.meanlen*100, colour=as.factor(HWDATE),linetype = as.factor(HWINC)))+
  geom_point(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.meanlen-base.meanlen)/base.meanlen*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-0.1, 6), "% deviation in\n individual length")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )

cond_prerep_plot <- ggplot()+
  geom_line(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.meancond-base.meancond)/base.meancond*100, colour=as.factor(HWDATE), linetype = as.factor(HWINC)))+
  geom_point(data=prerep_struc, aes(x=hw.meanage/365, y=(hw.meancond-base.meancond)/base.meancond*100, colour=as.factor(HWDATE), shape= as.factor(HWINC)))+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 5), limits=c(-3, 15), "% deviation in\nindividual condition")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(        axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                legend.position = "top")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )


  
HWstruc_plot<- wrap_plots(A=zoop_plot, B= dens_posthw_plot, C=dens_prerep_plot, D=length_posthw_plot, E=length_prerep_plot, G=cond_posthw_plot, H=cond_prerep_plot,
           design = "AA
           ##
           BC
           ##
           DE
           ##
           GH",
           ncol = 2 , guides="collect", tag_level = 'new', heights = c(1, -0, 1, -0, 1, -0, 1))+
  plot_annotation(tag_levels="a")&
  theme(legend.position = "top")


ggsave(paste0("Figures/Figure4.pdf"), HWstruc_plot , width=15, height=25, units="cm")







###### SUPPLEMENT FIGURE S7, FGIURE 3 WITH ABOSOLUTE VALUES (REP 29340) #####

zoop_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_vline(data=data.frame(HWDATE = as.numeric(unique(HWts_env$HWDATE))), aes(xintercept=as.Date(HWDATE+15-29200), colour=as.factor(HWDATE)), linetype="dotted")+
  geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWts_env, aes(x=as.Date(hw.reltime-29200), y=(hw.zoop-base.zoop)/base.zoop*100, colour=as.factor(HWDATE), linetype = as.factor(HWINC)))+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_date( limits = as.Date(c(124.5,163+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(-25, 400), "% deviation in\nzooplankton", breaks=c(-25,0,25,50,100))+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"), guide="none")+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"), guide="none")+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(163.5), y=130), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))

zoop_plot_abs1 <- ggplot()+
  geom_rect(aes(xmin=as.Date(29325), xmax=as.Date(29325+77 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(aes(xmin=as.Date(29325+365), xmax=as.Date(29325+38+365 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(data=data.frame(HWDATE = 29340), aes(xmin=as.Date(HWDATE-1), xmax=as.Date(HWDATE+11), ymin=-Inf,ymax=Inf) ,fill="#F8766D",  alpha=0.1 )+
  geom_vline(data=data.frame(HWDATE = 29340), aes(xintercept=as.Date(HWDATE+15)), colour="#F8766D", linetype="dotted")+
  geom_vline(xintercept=as.Date(29689), linetype="dotted")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29340),], aes(x=as.Date(hw.reltime), y=hw.zoop), colour="#F8766D", linetype ="solid")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29340),], aes(x=as.Date(hw.reltime), y=base.zoop), colour="black", linetype = "dashed")+
  theme_classic()+
  scale_x_date( limits = as.Date(c(29325,29325+38+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(limits=c(0, 0.5), "Relative zooplankton\nbiomass (gram)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(29363.5), y=0.5), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


dens_posthw_plot_abs1 <- ggplot()+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.totdens), colour="black", linetype = "dashed")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.totdens), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.totdens), colour="#F8766D", linetype ="solid")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.totdens, ), colour="#F8766D")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Directly after the heatwave")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_posthw_plot_abs1 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#F8766D", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#F8766D")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#F8766D")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


cond_posthw_plot_abs1 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#F8766D", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meancond), colour="#F8766D")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meancond), colour="#F8766D")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))




dens_prerep_plot_abs1 <- ggplot()+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE),linetype = as.factor(HWINC)), colour="black", linetype = "dashed")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE), shape= as.factor(HWINC)), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.totdens), colour="#F8766D", linetype ="solid")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.totdens), , colour="#F8766D")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Before following spawning season")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_prerep_plot_abs1 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#F8766D", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#F8766D")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#F8766D")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )



cond_prerep_plot_abs1 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#F8766D", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meancond), colour="#F8766D")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29340),], aes(x=hw.meanage/365, y=hw.meancond), colour="#F8766D")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(        axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )




HWstruc_plot_abs1<- wrap_plots(A=zoop_plot_abs1, B= dens_posthw_plot_abs1, C=dens_prerep_plot_abs1, D=length_posthw_plot_abs1, E=length_prerep_plot_abs1, G=cond_posthw_plot_abs1, H=cond_prerep_plot_abs1,
                          design = "AA
           ##
           BC
           ##
           DE
           ##
           GH",
                          ncol = 2 , guides="collect", tag_level = 'new', heights = c(1, -0, 1, -0, 1, -0, 1))+
  plot_annotation(tag_levels="a")


ggsave(paste0("Figures/FigigureS7R.pdf"), HWstruc_plot_abs1 , width=15, height=22.5, units="cm")


###### SUPPLEMENT FIGURE S8, FGIURE 3 WITH ABOSOLUTE VALUES (REP 29340) #####

zoop_plot_abs2 <- ggplot()+
  geom_rect(aes(xmin=as.Date(29325), xmax=as.Date(29325+77 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(aes(xmin=as.Date(29325+365), xmax=as.Date(29325+38+365 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(data=data.frame(HWDATE = 29363), aes(xmin=as.Date(HWDATE-1), xmax=as.Date(HWDATE+11), ymin=-Inf,ymax=Inf), fill="#7CAE00", alpha=0.1 )+
  geom_vline(data=data.frame(HWDATE = 29363), aes(xintercept=as.Date(HWDATE+15)), colour="#7CAE00", linetype="dotted")+
  geom_vline(xintercept=as.Date(29689), linetype="dotted")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29363),], aes(x=as.Date(hw.reltime), y=hw.zoop, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29363),], aes(x=as.Date(hw.reltime), y=base.zoop), colour="black", linetype = "dashed")+
  theme_classic()+
  scale_x_date( limits = as.Date(c(29325,29325+38+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(limits=c(0, 0.6), "Relative zooplankton\nbiomass (gram)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(29363.5), y=0.5), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


dens_posthw_plot_abs2 <- ggplot()+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.totdens), colour="black", linetype = "dashed")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.totdens), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.totdens, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.totdens, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Directly after the heatwave")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_posthw_plot_abs2 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen,  group=as.factor(HWINC)), fill="#7CAE00", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


cond_posthw_plot_abs2 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, ymin=hw.mincond, ymax=hw.maxcond, group=as.factor(HWINC)), fill="#7CAE00", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))




dens_prerep_plot_abs2 <- ggplot()+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.totdens), colour="black", linetype = "dashed")+
  geom_point(data=prerep_struc[which(posthw_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.totdens), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.totdens, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.totdens, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Before following spawning season")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))



length_prerep_plot_abs2 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen,  group=as.factor(HWINC)), fill="#7CAE00", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meanlen, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )



cond_prerep_plot_abs2 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, ymin=hw.mincond, ymax=hw.maxcond, group=as.factor(HWINC)), fill="#7CAE00", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, linetype = as.factor(HWINC)), colour="#7CAE00")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29363),], aes(x=hw.meanage/365, y=hw.meancond, shape= as.factor(HWINC)), colour="#7CAE00")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1.5, 4.5), "Individual relative\ncondition (-)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  theme(        axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )




HWstruc_plot_abs2<- wrap_plots(A=zoop_plot_abs2, B= dens_posthw_plot_abs2, C=dens_prerep_plot_abs2, D=length_posthw_plot_abs2, E=length_prerep_plot_abs2, G=cond_posthw_plot_abs2, H=cond_prerep_plot_abs2,
                               design = "AA
           ##
           BC
           ##
           DE
           ##
           GH",
                               ncol = 2 , guides="collect", tag_level = 'new', heights = c(1, -0, 1, -0, 1, -0, 1))+
  plot_annotation(tag_levels="a")


ggsave(paste0("Figures/FigureS8R.pdf"), HWstruc_plot_abs2 , width=15, height=22.5, units="cm")



###### SUPPLEMENT FIGURE S9, FGIURE 3 WITH ABOSOLUTE VALUES (REP 29444) #####

zoop_plot_abs3 <- ggplot()+
  geom_rect(aes(xmin=as.Date(29325), xmax=as.Date(29325+77 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(aes(xmin=as.Date(29325+365), xmax=as.Date(29325+38+365 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(data=data.frame(HWDATE = 29444), aes(xmin=as.Date(HWDATE-1), xmax=as.Date(HWDATE+11), ymin=-Inf,ymax=Inf) ,fill="#00BFC4",  alpha=0.1 )+
  geom_vline(data=data.frame(HWDATE = 29444), aes(xintercept=as.Date(HWDATE+15)), colour="#00BFC4", linetype="dotted")+
  geom_vline(xintercept=as.Date(29689), linetype="dotted")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29444),], aes(x=as.Date(hw.reltime), y=hw.zoop), colour="#00BFC4", linetype ="solid")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29444),], aes(x=as.Date(hw.reltime), y=base.zoop), colour="black", linetype = "dashed")+
  theme_classic()+
  scale_x_date( limits = as.Date(c(29325,29325+38+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(limits=c(0, 0.5), "Relative zooplankton\nbiomass (gram)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(29363.5), y=0.5), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


dens_posthw_plot_abs3 <- ggplot()+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.totdens), colour="black", linetype = "dashed")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.totdens), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.totdens), colour="#00BFC4", linetype ="solid")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.totdens, ), colour="#00BFC4")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Directly after the heatwave")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_posthw_plot_abs3 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#00BFC4", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#00BFC4")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#00BFC4")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


cond_posthw_plot_abs3 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#00BFC4", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meancond), colour="#00BFC4")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meancond), colour="#00BFC4")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))




dens_prerep_plot_abs3 <- ggplot()+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE),linetype = as.factor(HWINC)), colour="black", linetype = "dashed")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE), shape= as.factor(HWINC)), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.totdens), colour="#00BFC4", linetype ="solid")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.totdens), , colour="#00BFC4")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Before following spawning season")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_prerep_plot_abs3 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#00BFC4", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#00BFC4")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#00BFC4")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )



cond_prerep_plot_abs3 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#00BFC4", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meancond), colour="#00BFC4")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29444),], aes(x=hw.meanage/365, y=hw.meancond), colour="#00BFC4")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(        axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )




HWstruc_plot_abs3<- wrap_plots(A=zoop_plot_abs3, B= dens_posthw_plot_abs3, C=dens_prerep_plot_abs3, D=length_posthw_plot_abs3, E=length_prerep_plot_abs3, G=cond_posthw_plot_abs3, H=cond_prerep_plot_abs3,
                               design = "AA
           ##
           BC
           ##
           DE
           ##
           GH",
                               ncol = 2 , guides="collect", tag_level = 'new', heights = c(1, -0, 1, -0, 1, -0, 1))+
  plot_annotation(tag_levels="a")


ggsave(paste0("Figures/FigureS9R.pdf"), HWstruc_plot_abs3 , width=15, height=22.5, units="cm")




###### SUPPLEMENT FIGURE S10, FGIURE 3 WITH ABOSOLUTE VALUES (REP 29597) #####

zoop_plot_abs4 <- ggplot()+
  geom_rect(aes(xmin=as.Date(29325), xmax=as.Date(29325+77 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(aes(xmin=as.Date(29325+365), xmax=as.Date(29325+38+365 ), ymin=-Inf,ymax=Inf),fill="gray90")+
  geom_rect(data=data.frame(HWDATE = 29597), aes(xmin=as.Date(HWDATE-1), xmax=as.Date(HWDATE+11), ymin=-Inf,ymax=Inf) ,fill="#C77CFF",  alpha=0.1 )+
  geom_vline(data=data.frame(HWDATE = 29597), aes(xintercept=as.Date(HWDATE+15)), colour="#C77CFF", linetype="dotted")+
  geom_vline(xintercept=as.Date(29689), linetype="dotted")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29597),], aes(x=as.Date(hw.reltime), y=hw.zoop), colour="#C77CFF", linetype ="solid")+
  geom_line(data=HWts_env[which(HWts_env$HWDATE==29597),], aes(x=as.Date(hw.reltime), y=base.zoop), colour="black", linetype = "dashed")+
  theme_classic()+
  scale_x_date( limits = as.Date(c(29325,29325+38+365)),date_labels = "%b", name=element_blank() )+
  scale_y_continuous(limits=c(0, 0.5), "Relative zooplankton\nbiomass (gram)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_fill_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%d %b"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  geom_text(aes(x=as.Date(29363.5), y=0.5), label="Spawning\nseason", size=3.5)+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


dens_posthw_plot_abs4 <- ggplot()+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.totdens), colour="black", linetype = "dashed")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.totdens), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.totdens), colour="#C77CFF", linetype ="solid")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.totdens, ), colour="#C77CFF")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Directly after the heatwave")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_posthw_plot_abs4 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#C77CFF", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#C77CFF")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#C77CFF")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


cond_posthw_plot_abs4 <- ggplot()+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#C77CFF", alpha=0.2)+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meancond), colour="#C77CFF")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=posthw_struc[which(posthw_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meancond), colour="#C77CFF")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(legend.position = "none",
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))




dens_prerep_plot_abs4 <- ggplot()+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE),linetype = as.factor(HWINC)), colour="black", linetype = "dashed")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.totdens, colour=as.factor(HWDATE), shape= as.factor(HWINC)), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.totdens), colour="#C77CFF", linetype ="solid")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.totdens), , colour="#C77CFF")+
  scale_x_continuous(limits=c(0,20), name="Individual age (years)")+
  scale_y_continuous(trans="log10", limits=c(1E-11, 1E9), "Relative cohort\ndensity (number)")+
  scale_colour_discrete(name="Heatwave start", labels=format(as.Date(as.numeric(unique(HWts_env$HWDATE))),"%b %d"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  scale_shape_manual(name="Heatwave increase", values=c(1,16), guide = "none")+
  theme_classic()+
  labs(title="Before following spawning season")+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=10),
        plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4))


length_prerep_plot_abs4 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen, ymin=base.minlen, ymax=base.maxlen), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.minlen, ymax=hw.maxlen), fill="#C77CFF", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#C77CFF")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meanlen), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen), colour="#C77CFF")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(0,30), "Individual length\ndistribution (cm)")+
  theme_classic()+
  theme(axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )



cond_prerep_plot_abs4 <- ggplot()+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond, ymin=base.mincond, ymax=base.maxcond), fill="black", alpha=0.2)+
  geom_ribbon(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meanlen, ymin=hw.mincond, ymax=hw.maxcond), fill="#C77CFF", alpha=0.2)+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_line(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meancond), colour="#C77CFF")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=base.meanage/365, y=base.meancond), colour="black")+
  geom_point(data=prerep_struc[which(prerep_struc$HWDATE==29597),], aes(x=hw.meanage/365, y=hw.meancond), colour="#C77CFF")+
  scale_x_continuous(limits=c(-0.1,20), name="Individual age (years)")+
  scale_y_continuous(limits=c(1, 4.5), "Individual relative\ncondition (-)")+
  theme_classic()+
  theme(        axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                legend.position = "none")+
  theme(plot.tag.location="panel", plot.tag.position = "topleft", plot.tag = element_text(hjust=1.2, vjust=4) )




HWstruc_plot_abs4<- wrap_plots(A=zoop_plot_abs4, B= dens_posthw_plot_abs4, C=dens_prerep_plot_abs4, D=length_posthw_plot_abs4, E=length_prerep_plot_abs4, G=cond_posthw_plot_abs4, H=cond_prerep_plot_abs4,
                               design = "AA
           ##
           BC
           ##
           DE
           ##
           GH",
                               ncol = 2 , guides="collect", tag_level = 'new', heights = c(1, -0, 1, -0, 1, -0, 1))+
  plot_annotation(tag_levels="a")


ggsave(paste0("Figures/FigureS10R.pdf"), HWstruc_plot_abs4 , width=15, height=22.5, units="cm")






###### BIFURCATION WITH HIGH MUS (SUPPL FIG 11) ######

### LOAD SIMULATION DATA ###

# dataframe with varied variables for names of files
HWbif_highMUS_vars <- expand.grid(DUR = c("5", "10", "20", "40"), INC=1:6, SUB = c("a","b"), KEEP.OUT.ATTRS = FALSE)

# load simulation data into list
HWbif_highMUS_list <- apply(HWbif_highMUS_vars, 1, function(x){
  name <- paste0("Model/HWbif_highMUS/herringEBT_",x[1],"_",x[2],x[3],".out")
  
  if(!file.exists(name)) return(NA)
  
  bifdata <- read.delim(name,
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "HWSTART2"))
  return(as.data.frame(bifdata))
} )

# combine dataframes
HWbif_highMUS <- na.omit(data.frame(do.call("rbind",HWbif_highMUS_list)))

# calculate actual time
HWbif_highMUS$reltime <-  HWbif_highMUS$time%%36500+125

# calculate the number of years after the heatwave
HWbif_highMUS$hwyear <- (HWbif_highMUS$reltime-HWbif_highMUS$HWSTART )/365

# calculate total harvested biomass
HWbif_highMUS$harvtot <- rowSums(HWbif_highMUS[,17:56])




### FIGURE S11 ### 


HWSTART_INCrel_highMUS_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif_highMUS[which(HWbif_highMUS$HWDUR==10   & HWbif_highMUS$hwyear >=0  &  HWbif_highMUS$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y= (eggbuff-2.158150E+09)/2.158150E+09*100, colour=as.factor(HWINC)), linewidth=0.5)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_y_continuous( limits=c( -3, 4 ) )+
  scale_colour_manual(values=c("grey60","grey50","grey40","grey30","grey20","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Temperature increase during heatwave (\u00B0C)")+
  geom_text(aes(x=as.Date(163.5), y=4), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5), axis.text.x = element_blank())



HWSTART_DURrel_highMUS_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif_highMUS[which(HWbif_highMUS$HWINC==5   & HWbif_highMUS$hwyear >=0  &  HWbif_highMUS$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y=(eggbuff-2.158150E+09)/2.158150E+09*100, colour=as.factor(HWDUR)), linewidth=0.5)+
  scale_y_continuous( limits=c( -7, 7) )+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_colour_manual(values=c("grey70","grey50","grey30","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Duration of heatwave (days)")+
  #geom_text(aes(x=as.Date(163.5), y=13), label="Spawning\nseason", size=2.5)+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5))


HWSTARTrel_highMUS_plot <- wrap_plots(list(HWSTART_INCrel_highMUS_plot, plot_spacer() , 
                                   HWSTART_DURrel_highMUS_plot), 
                              ncol = 1 , 
                              axis_titles = "collect_x", 
                              tag_level = 'new', 
                              heights=c( 1,-0.1,1))+plot_annotation(tag_levels="a")



HWSTARTrel_highMUS_plot <- wrap_plots(list(HWSTART_INCrel_highMUS_plot, plot_spacer() , HWSTART_DURrel_highMUS_plot), ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1,1))+plot_annotation(tag_levels="a")

ggsave("Figures/FigureS11R.pdf", HWSTARTrel_highMUS_plot, width=15, height=15, units="cm")








###### BIFURCATION WITH HIGER MEAN TEMPERATURES (SUPPL FIG 13) ######

### LOAD SIMULATION DATA ###

# dataframe with varied variables for names of files
HWbif_TCC_vars <- expand.grid(INC=1:6, SUB = c("a","b"), KEEP.OUT.ATTRS = FALSE)

# load simulation data into list
HWbif_TCC_list <- apply(HWbif_TCC_vars, 1, function(x){
  name <- paste0("Model/HWbif_TCC/herringEBT_",x[1],x[2],".out")
  
  if(!file.exists(name)) return(NA)
  
  bifdata <- read.delim(name,
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "TCC","HWSTART2"))
  
  bifdata$eggbuffbase <- bifdata$eggbuff[2]
  
  return(as.data.frame(bifdata))
} )

# combine dataframes
HWbif_TCC <- na.omit(data.frame(do.call("rbind",HWbif_TCC_list)))

# calculate actual time
HWbif_TCC$reltime <-  HWbif_TCC$time%%36500+125

# calculate the number of years after the heatwave
HWbif_TCC$hwyear <- (HWbif_TCC$reltime-HWbif_TCC$HWSTART )/365

# calculate total harvested biomass
HWbif_TCC$harvtot <- rowSums(HWbif_TCC[,17:56])

### FIGURE S13 ### 


HWSTART_INCrel_TCC_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif_TCC[which(HWbif_TCC$HWDUR==10   & HWbif_TCC$hwyear >=0  &  HWbif_TCC$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y= (eggbuff-eggbuffbase)/eggbuffbase*100, colour=as.factor(TCC)), linewidth=0.5)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_y_continuous( limits=c( -3, 9 ) )+
  scale_colour_manual(values=c("grey60","grey50","grey40","grey30","grey20","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Center of seasonal temperature pattern (\u00B0C)")+
  geom_text(aes(x=as.Date(163.5), y=4), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5))


ggsave("Figures/FigureS13R.pdf", HWSTART_INCrel_TCC_plot, width=15, height=7.5, units="cm")




##### BIFRUCATION WITH ALternatIVE EGG MORTALITY (FIG S12)

### LOAD SIMULATION DATA ###

# dataframe with varied variables for names of files
HWbif_altegg_vars <- expand.grid(INC=1:6, KEEP.OUT.ATTRS = FALSE)

# load simulation data into list
HWbif_altegg_list <- apply(HWbif_altegg_vars, 1, function(x){
  name <- paste0("Model/HWbif_altegg/herringEBT_10_",x[1],".out")
  
  if(!file.exists(name)) return(NA)
  
  bifdata <- read.delim(name,
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1))))
  
  bifdata$eggbuffbase <- bifdata$eggbuff[2]
  
  return(as.data.frame(bifdata))
} )

# combine dataframes
HWbif_altegg <- na.omit(data.frame(do.call("rbind",HWbif_altegg_list)))

# calculate actual time
HWbif_altegg$reltime <-  HWbif_altegg$time%%36500+125

# calculate the number of years after the heatwave
HWbif_altegg$hwyear <- (HWbif_altegg$reltime-HWbif_altegg$HWSTART )/365

# calculate total harvested biomass
HWbif_altegg$harvtot <- rowSums(HWbif_altegg[,17:56])

### FIGURE S12 ### 
HWSTARTrel_altegg_plot <- ggplot()+
  geom_rect( aes(xmin=as.Date(124.5), xmax=as.Date(202.5), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_rect( aes(xmin=as.Date(124.5+365), xmax=as.Date(163+365), ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #geom_vline(xintercept=as.Date(125+365), linetype="dotted", colour="black")+
  geom_line(data=HWbif_altegg[which(HWbif_altegg$HWDUR==10   & HWbif_altegg$hwyear >=0  &  HWbif_altegg$hwyear <1),], 
            aes(x=as.Date((HWSTART-125+365+0.5*HWDUR)%%365+125), y= (eggbuff-eggbuffbase)/eggbuffbase*100, colour=as.factor(HWINC)), linewidth=0.5)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(124.5,163+365)), breaks=date_breaks(width="2 month")  )+
  scale_y_continuous( limits=c( -5, 5 ) )+
  scale_colour_manual(values=c("grey60","grey50","grey40","grey30","grey20","grey10"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(x="Midpoint of the heatwave", y="% deviation in total egg production", colour="Temperature increase during heatwave (\u00B0C)")+
  geom_text(aes(x=as.Date(163.5), y=4), label="Spawning\nseason", size=2.5)+
  ggnewscale::new_scale_colour()+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.title = element_text(hjust=0.5))

HWSTARTrel_altegg_comb  <- ggarrange(HWSTARTrel_altegg_plot, peck2012ES_plot_alt, ncol=2, nrow=1, align="h", 
                           labels=c("a", "b"), widths = c(2,1))


ggsave("Figures/FigureS12R.pdf", HWSTARTrel_altegg_comb, width=22, height=7.5, units="cm")




##### COMPARISION OF THE EFFECTS OF JUMPS (FIGS S5 AND S6) #####


# dataframe with varied variables for names of files


# load simulation data
  data_norm <- read.delim("Model/jumptest/herringEBT.out",
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","cohono","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "TCC","HWSTART2"))

data_norm$totharv <- rowSums(data_norm[,17:56])

data_jmp <- read.delim("Model/jumptest/herringEBT_jmptst.out",
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","cohono","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "TCC","HWSTART2"))

data_jmp$totharv <- rowSums(data_jmp[,17:56])

data_jmptst <- merge(data_norm,data_jmp, by="time", suffixes = c("_norm","_jmp") )


avg_cohonum = mean(data_jmptst$cohono_norm)

reltol_zoop_plot <- ggplot(data = data_jmptst , aes(x=time/365, y=(zoop_norm-zoop_jmp)/zoop_jmp))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line(colour="red")+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin zooplankton density")

reltol_egg_plot <- ggplot(data = data_jmptst[which(data_jmptst$time%%365==0),] , aes(x=time/365, y=(eggbuff_norm-eggbuff_jmp)/eggbuff_jmp))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line(colour="red")+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin egg production")

reltol_harv_plot <- ggplot(data = data_jmptst[which(data_jmptst$time%%365==0),] , aes(x=time/365, y=(totharv_norm-totharv_jmp)/totharv_jmp))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line(colour="red")+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin harvested biomass")


reltol_plot <- wrap_plots(list(reltol_zoop_plot, reltol_egg_plot, reltol_harv_plot), ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,1,1))+plot_annotation(tag_levels="a")


ggsave("Figures/FigureS5.pdf", reltol_plot, width=15, height=20, units="cm")



# load runs with heatwaves

# list with values needed for analysis per run
HWjmptst_vars <- data.frame(RUNNAME = c("10_5_140","10_5_163", "10_5_244", "10_5_397", "10_2_163"))


# read environment output 
HWjmptst_list <- apply(HWjmptst_vars, 1, function(x){
  data_norm <- read.delim(paste0("Model/jumptest/herringEBT_",x,".out"),
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","cohono","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "TCC","HWSTART2"))
  
  data_norm$totharv <- rowSums(data_norm[,17:56])
  
  data_jmp <- read.delim(paste0("Model/jumptest/herringEBT_jmptst_",x,".out"),
                         sep = "\t",
                         header=FALSE,
                         col.names=c("time","cohono","zoop","bent","eggbuff","egglarv","juvnum","juvmass","adltnum","adltmass","matage","maxlen","temp","HWSTART","HWDUR","HWINC",paste0("harv",seq(0,39,1)), "TCC","HWSTART2"))
  
  data_jmp$totharv <- rowSums(data_jmp[,17:56])
  
  data_jmptst <- merge(data_norm,data_jmp, by="time", suffixes = c("_norm","_jmp") )
  
  return(as.data.frame(data_jmptst))
} )

HWjmptst <- data.frame(do.call("rbind",HWjmptst_list))



HWjmptst$HWyear = floor((HWjmptst$time - HWjmptst$HWSTART_jmp)/365 )

HWreltol_zoop_plot<- ggplot(data = HWjmptst[which(HWjmptst$HWyear>=0),] , aes(x=time/365, y=(zoop_norm-zoop_jmp)/zoop_jmp, colour=as.factor((HWSTART_jmp-125)%%365), linetype=as.factor(HWINC_jmp)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line()+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin zooplankton density")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="none")
    

HWreltol_egg_plot <- ggplot(data = HWjmptst[which(HWjmptst$time%%365==0),] , aes(x=time/365, y=(eggbuff_norm-eggbuff_jmp)/eggbuff_jmp, colour=as.factor((HWSTART_jmp-125)%%365), linetype=as.factor(HWINC_jmp)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line()+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin egg production")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="none")

HWreltol_harv_plot <- ggplot(data = HWjmptst[which(HWjmptst$time%%365==0),] , aes(x=time/365, y=(totharv_norm-totharv_jmp)/totharv_jmp,colour=as.factor((HWSTART_jmp-125)%%365), linetype=as.factor(HWINC_jmp)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept=2*(avg_cohonum*3+3)*1E-10)+
  geom_hline(yintercept=-2*(avg_cohonum*3+3)*1E-10)+
  geom_line()+
  theme_classic()+
  labs(x="integration time (years)", y="relative deviation\nin harvested biomass")+
  scale_colour_discrete(name=element_blank(), labels=c("Heatwave during spawning season", "Heatwave during spawning peak","Heatwave after spawning season","Heatwave before spawning season"))+
  scale_linetype_manual(name="Heatwave increase", values=c("dashed","solid"), guide = "none")+
  theme(legend.position="none")

HWreltol_plot <- wrap_plots(list(HWreltol_zoop_plot, HWreltol_egg_plot, HWreltol_harv_plot), ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,1,1))+plot_annotation(tag_levels="a")

ggsave("Figures/FigureS6.pdf", HWreltol_plot, width=15, height=20, units="cm")
