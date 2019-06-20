#heat2
controlIn = 'beww_experiment/three_way/official_All_three_CONTROL_12-12-18.csv'
treatIn = 'beww_experiment/three_way/official_All_three_heat2_alongside_allThree_12-12-18.csv'
experimentStart = mdy_hm('12/12/18 11:30')
experimentEnd = mdy_hm('12/13/18 9:00') #end after photo time 2 bc heats died
cdat = read_hobo_csv(controlIn, 'control') %>% 
  filter(Time<experimentEnd)
tdat = read_hobo_csv(treatIn, 'treatment') %>% 
  filter(Time<experimentEnd)
pdat = read.csv('beww_experiment/three_way/photoTimes.csv') %>% 
  mutate(Time=mdy_hm(Time)) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60),
         fixedSamples = c(FALSE, FALSE, TRUE, FALSE, FALSE)
  ) %>% 
  filter(Time<experimentEnd)
pdat
