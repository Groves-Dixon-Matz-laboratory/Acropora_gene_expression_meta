#setwd("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts")




data1=read.table("MWU_MF_ethanol_ForMWU.csv",header=T)
data2=read.table("MWU_MF_methanol_ForMWU.csv",header=T)
data2=read.table("MWU_MF_stress_ForMWU.csv",header=T)

goods=intersect(data1$term,data2$term)
#goods=unique(as.character(c(data1$term[data1$p.adj<=0.1],data2$term[data2$p.adj<=0.1])))
length(goods)

data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]

# all overlapping GO terms
ress=merge(data1,data2,by="term") %>% 
  mutate(quad = if_else(delta.rank.x < 0 & delta.rank.y > 0,
                        'XdownYup',
                        'none'),
         quad = if_else(delta.rank.x > 0 & delta.rank.y > 0,
                        'XupYup',
                        quad),
         quad = if_else(delta.rank.x > 0 & delta.rank.y < 0,
                        'XupYdown',
                        quad),
         quad = if_else(delta.rank.x < 0 & delta.rank.y < 0,
                        'XdownYdown',
                        quad))
table(ress$quad)

plot(delta.rank.x~delta.rank.y,ress,xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms highly signifcant in any of the two datasets
sigs=(ress$p.adj.x<=0.01 | ress$p.adj.y<=0.01)
sum(sigs) # 71
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms signifcant in both datasets
sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1)
sum(sigs) # 20
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)



# look at genes -----------------------------------------------------------
dat = ress %>% 
  filter(pval.x < 0.05 & pval.y < 0.05)
#up in both
dat %>% 
  filter(quad=='XupYup')

#down in both
dat %>% 
  filter(quad=='XdownYdown')

#up in X down in Y
dat %>% 
  filter(quad=='XupYdown')

#up in Y down in X
dat %>% 
  filter(quad=='XdownYup')

