#three-way  switches only
controlIn = 'beww_experiment/three_way/official_All_three_CONTROL_12-12-18.csv'
CtoH = 'beww_experiment/three_way/official_All_three_Cold_to_hot_12-12-18.csv'
HtoC = 'beww_experiment/three_way/official_All_three_Hot_to_Cold_12-12-18.csv'
experimentStart = mdy_hm('12/12/18 11:25')
experimentEnd = mdy_hm('12/13/18 18:40')
cdat = read_hobo_csv(controlIn, 'control')
ch = read_hobo_csv(CtoH, 'Cold then Hot')
hc = read_hobo_csv(HtoC, 'Hot then Cold') 
tdat = rbind(ch, hc)
pdat = read.csv('beww_experiment/three_way/photoTimes.csv') %>% 
  mutate(Time=mdy_hm(Time)) %>% 
  mutate(expDay=day(Time)-day(experimentStart),
         expHour=
           24*expDay + 
           hour(Time) + 
           minute(Time)/60 - 
           (hour(experimentStart) + minute(experimentStart)/60),
         fixedSamples = c(FALSE, FALSE, TRUE, TRUE, FALSE)
  ) %>% 
  filter(Time<experimentEnd)
pdat