#plot_treatment_timeline.R

library(tidyverse)
library(lubridate)
library(cowplot)


#funciton to read in and format hobo data
read_hobo_csv = function(infileName, group){
  read.csv(infileName) %>% 
    rename('Time'=1,
           'Temp'=2,
           'Light'=3) %>% 
    mutate(Time=mdy_hm(Time)) %>%
    mutate(Group = group) %>% 
    as_tibble()
}



# SELECT DATASET ----------------------------------------------------------

#execute the chunk for whichever experiment you want to plot

#chilled
source('beww_experiment/chilled/loadData.R')

#heat1
source('beww_experiment/heat/loadData.R')

#heat2
source('beww_experiment/three_way/loadHeat2.R')

#three-way  switches only
source('beww_experiment/three_way/loadSwitchesOnly.R')

#all three-way at once
source('beww_experiment/three_way/loadAllAtOnce.R')



# SET UP DATA -------------------------------------------------------------


#MERGE AND FORMAT WITH HOURS STARTING AT ZERO
dat =  rbind(cdat, tdat) %>% 
  filter(Time>experimentStart) %>% 
  filter(Time<=experimentEnd) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60)
         )

# PLOT --------------------------------------------------------------------

# #by time/date
# dat %>% 
#   ggplot(aes(x=Time, y=Temp, color=Group)) +
#   geom_line(lwd=2) +
#   geom_vline(xintercept=pdat$Time, linetype="dotted", lwd=1)+
#   geom_vline(xintercept=mdy_hm('2018-12-13 00:12:00'))


#by experimental hour
dat %>% 
  ggplot(aes(x=expHour, y=Temp, color=Group)) +
  # geom_vline(xintercept=round(seq(min(dat$expHour), max(dat$expHour), by = 0.5),1), lwd=0.1, color='grey') +
  geom_line(lwd=2) +
  labs(x='Experiment time (hours)', y=expression(paste("Temp ",degree,"C"))) +
  geom_vline(xintercept=pdat$expHour[pdat$fixedSamples], linetype="dotted", lwd=1) +
  scale_x_continuous(breaks = round(seq(min(dat$expHour), max(dat$expHour), by = 3),0)) +
  theme(legend.title=element_blank())



# RUN ALL AT ONCE ---------------------------------------------------------

sourceFiles = c('beww_experiment/chilled/loadData.R',
                'beww_experiment/heat/loadData.R',
                'beww_experiment/three_way/loadHeat2.R',
                'beww_experiment/three_way/loadSwitchesOnly.R')

plotList = list()
for (fileName in sourceFiles){
  source(fileName)
  dat =  rbind(cdat, tdat) %>% 
    filter(Time>experimentStart) %>% 
    filter(Time<=experimentEnd) %>% 
    mutate(expDay=day(Time)-day(experimentStart),
           expHour=
             24*expDay + 
             hour(Time) + 
             minute(Time)/60 - 
             (hour(experimentStart) + minute(experimentStart)/60)
    )
  
  plt=dat %>% 
    ggplot(aes(x=expHour, y=Temp, color=Group)) +
    # geom_vline(xintercept=round(seq(min(dat$expHour), max(dat$expHour), by = 0.5),1), lwd=0.1, color='grey') +
    geom_line(lwd=2) +
    labs(x='Experiment time (hours)', y=expression(paste("Temp ",degree,"C"))) +
    geom_vline(xintercept=pdat$expHour[pdat$fixedSamples], linetype="dotted", lwd=1) +
    scale_x_continuous(breaks = round(seq(min(dat$expHour), max(dat$expHour), by = 3),0)) +
    theme(legend.title=element_blank())
  plotList[[fileName]]=plt
}

#run salinity
source('')
plotList['salinity']=

#build them
quartz()
plotList



