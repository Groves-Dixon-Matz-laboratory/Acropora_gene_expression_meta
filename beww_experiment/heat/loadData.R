#heat
controlIn = 'beww_experiment/heat/official_heat_CONTROL_12-10-18.csv'
treatIn = 'beww_experiment/heat/official_heat_TREATMENT_12-10-18.csv'
experimentStart = mdy_hm('12/10/18 10:25')
experimentEnd = mdy_hm('12/11/18 14:25')
cdat = read_hobo_csv(controlIn, 'control')
tdat = read_hobo_csv(treatIn, 'treatment')
pdat = read.csv('beww_experiment/heat/photoTimes.csv') %>% 
  mutate(Time=mdy_hm(Time)) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60),
         fixedSamples = c(FALSE, TRUE, TRUE)
  ) %>% 
  filter(Time<experimentEnd)
pdat