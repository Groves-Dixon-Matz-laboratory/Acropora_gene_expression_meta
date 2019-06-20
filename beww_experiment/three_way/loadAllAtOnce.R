#all three-way at once
controlIn = 'beww_experiment/three_way/official_All_three_CONTROL_12-12-18.csv'
CtoH = 'beww_experiment/three_way/official_All_three_Cold_to_hot_12-12-18.csv'
HtoC = 'beww_experiment/three_way/official_All_three_Hot_to_Cold_12-12-18.csv'
heat2 = 'beww_experiment/three_way/official_All_three_heat2_alongside_allThree_12-12-18.csv'
experimentStart = mdy_hm('12/12/18 10:30')
experimentEnd = mdy_hm('12/13/18 18:40')
cdat = read_hobo_csv(controlIn, 'control')
h2 = read_hobo_csv(heat2, 'Heat2') %>% 
  filter(Time<'2018-12-13 00:52:00')
ch = read_hobo_csv(CtoH, 'Cold then Hot')
hc = read_hobo_csv(HtoC, 'Hot then Cold') %>% 
  filter(Time<'2018-12-13 00:52:00')
tdat = rbind(ch, hc, h2)
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
