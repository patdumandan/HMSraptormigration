gams=read.csv(file.choose(), h=T)
gams$years=(gams$year-mean(gams$year))/(2*sd(gams$year))
gams$obs_days=log(gams$obs.hours/24)
gams$agri.cov=gams$pasture.cov+gams$crop.cov
str(gams)

#standardized

gams$urban.covs=(gams$urban.cov-mean(gams$urban.cov))/(2*sd(gams$urban.cov))
gams$forest.covs=(gams$forest.cov-mean(gams$forest.cov))/(2*sd(gams$forest.cov))
gams$agri.covs=(gams$agri.cov-mean(gams$agri.cov))/(2*sd(gams$agri.cov))

require(rstanarm)
require(dplyr)
amke=gams%>%filter(species=="AMKE")
baea=gams%>%filter(species=="BAEA")
blvu=gams%>%filter(species=="BLVU")
bwha=gams%>%filter(species=="BWHA")
coha=gams%>%filter(species=="COHA")
goea=gams%>%filter(species=="GOEA")
merl=gams%>%filter(species=="MERL")
nogo=gams%>%filter(species=="NOGO")
noha=gams%>%filter(species=="NOHA")
ospr=gams%>%filter(species=="OSPR")
pefa=gams%>%filter(species=="PEFA")
rlha=gams%>%filter(species=="RLHA")
rsha=gams%>%filter(species=="RSHA")
rtha=gams%>%filter(species=="RTHA")
ssha=gams%>%filter(species=="SSHA")
tuvu=gams%>%filter(species=="TUVU")

amke1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=amke)
summary(amke2)                 
amke2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=amke)

baea1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=baea)
summary(baea2)                 
baea2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=baea)
loo(baea2, cores=1)

blvu1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=blvu)
summary(blvu2)                
blvu2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=blvu)
loo(blvu2, cores=1)

bwha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=bwha)
summary(bwha2)                
bwha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=bwha)
loo(bwha2, cores=1)

coha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=coha)
summary(coha2)                
coha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=coha)
loo(coha2, cores=1)

goea1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=goea)
summary(goea2)                
goea2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=goea)
loo(goea2, cores=1)

merl1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=merl)
summary(merl2)                
merl2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=merl)
loo(merl2, cores=1)

nogo1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=nogo)
summary(nogo2)                
nogo2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=nogo)
loo(nogo2, cores=1)

noha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=noha)
summary(noha2)                
noha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=noha)
loo(noha2, cores=1)

ospr1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=ospr)
summary(ospr2)                
ospr2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=ospr)
loo(ospr2, cores=1)

pefa1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=pefa)
summary(pefa2)                
pefa2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=pefa)
loo(pefa2, cores=1)

rlha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=rlha)
summary(rlha2)                
rlha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=rlha)
loo(rlha2, cores=1)

rsha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=rsha)
summary(rsha2)                
rsha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=rsha)
loo(rsha2, cores=1)

rtha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=rtha)
summary(rtha2)                
rtha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=rtha)
loo(rtha2, cores=1)

ssha1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=ssha)
summary(ssha2)                
ssha2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=ssha)
loo(ssha2, cores=1)

tuvu1=stan_glm.nb(count~years+weighted.pasture+weighted.urban+obs_days, data=tuvu)
summary(tuvu2)                
tuvu2=stan_glm.nb(count~years+weighed.forest+weighted.urban+obs_days, data=tuvu)
loo(tuvu2, cores=1)
