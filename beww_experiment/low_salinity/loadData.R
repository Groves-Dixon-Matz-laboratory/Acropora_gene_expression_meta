#low salinity
experimentStart = mdy_hm('12/7/18 20:00')
experimentEnd = mdy_hm('12/8/18 18:01')
pdat = read.csv('beww_experiment/low_salinity/photoTimes.csv') %>% 
  mutate(Time=mdy_hm(Time)) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60),
         fixedSamples = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  ) %>% 
  filter(Time<experimentEnd)
pdat

#revise to include a return to 35 ppt point, assuming 5L per hour
pdat$Timepoint = c(0,1,3,4,5)
pdat$Time<-NULL
returnSaltLine = c(2,1, 7, FALSE)
names(returnSaltLine) = colnames(pdat)
pdat = rbind(pdat, returnSaltLine) %>% 
  mutate(fixedSamples=as.logical(fixedSamples))
  arrange(Timepoint)
pdat

#set pdat hours and add salinity measures
dat = pdat %>% 
  mutate(control = 35,
         treatment = c(15, 15, 35, 35, 35, 35)
         
  ) %>% 
  gather(key='Group', value = 'salinity', c('control', 'treatment'))


plt=dat %>% 
  ggplot(aes(x=expHour, y=salinity, color=Group)) +
  # geom_vline(xintercept=round(seq(min(dat$expHour), max(dat$expHour), by = 0.5),1), lwd=0.1, color='grey') +
  # geom_line(lwd=2) +
  geom_smooth(se=F, lwd=2, span=0.4) +
  labs(x='Experiment time (hours)', y=expression("Approx. Salinity (ppt)")) +
  geom_vline(xintercept=pdat$expHour[pdat$fixedSamples], linetype="dotted", lwd=1) +
  scale_x_continuous(breaks = round(seq(min(dat$expHour), max(dat$expHour), by = 3),0)) +
  theme(legend.title=element_blank())
