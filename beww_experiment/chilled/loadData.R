#chilled
controlIn = 'beww_experiment/chilled/chilled_official_CONTROL_12-5-18.csv'
treatIn = 'beww_experiment/chilled/chilled_official_TREATMENT_12-5-18.csv'
experimentStart = mdy_hm('12/5/18 17:12')
experimentEnd = mdy_hm('12/6/18 14:49')
cdat = read_hobo_csv(controlIn, 'control')
tdat = read_hobo_csv(treatIn, 'treatment')
pdat = read.csv('beww_experiment/chilled/photoTimes.csv') %>% 
  mutate(Time=mdy_hm(Time)) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60),
         fixedSamples = c(FALSE, TRUE, FALSE, FALSE)
  ) %>% 
  filter(Time<experimentEnd)
pdat


