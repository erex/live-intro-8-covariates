## ----setup, include=FALSE----------------------------------------------------------
library(tint)
library(extraDistr)
library(dsims)
library(Distance)
# invalidate cache when the package version changes
knitr::opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tint'))
options(htmltools.dir.version = FALSE)


## ---- echo=TRUE, eval=FALSE--------------------------------------------------------
## a.covariate <- ds(my.data, transect="line", key="hn", formula=~size)


## ----fig-nocap-margin-first, fig.width=8, echo=FALSE-------------------------------
group.size <- rtpois(1000,lambda = 10)
hist(group.size, xlab="Group size")


## ----sim10, cache=TRUE, echo=FALSE, message=FALSE, warning=FALSE-------------------
set.seed(2400112)

simreps <- 350
trunc.distance <- 80
eg.region <- make.region()
my.density <- make.density(region=eg.region, constant=1)
covariate.list <- list()
covariate.list$size <- list(distribution="ztruncpois", mean = 10)
pop.desc <- make.population.description(region=eg.region, 
                                        density=my.density,
                                        covariates = covariate.list, N=200)
cov.params <- list(size = c(.1))
detect <- make.detectability(scale.param = 10, 
                             cov.param = cov.params, 
                             truncation = trunc.distance)
plot(detect, pop.desc)
my.pop <- generate.population(pop.desc, detect, eg.region)
transects <- make.design(region = eg.region, spacing=200, truncation = trunc.distance)
plot(eg.region, generate.transects(transects), covered.area=TRUE)
detfn.size <- make.ds.analysis(dfmodel=list(~size),
                               key="hn", truncation = trunc.distance)
size.cov <- make.simulation(reps=simreps, design=transects, 
                            population.description = pop.desc, 
                            detectability = detect,
                            ds.analysis = detfn.size)
size.cov.sim <- run.simulation(size.cov, run.parallel = TRUE, max.cores=11)


## ---- echo=FALSE-------------------------------------------------------------------
# survey.design <- create.survey.results(size.cov, dht.tables = TRUE)
# my.survey <- data.for.distance(object = survey.design)
# scatter.smooth(my.survey$distance, my.survey$size)


## ----small-with-covar, fig.width=8, echo=FALSE-------------------------------------
hist(size.cov.sim@results$expected.size[1,1,1:simreps], main="Computed average group size\nGroup size covariate", xlab="Group size")
abline(v=covariate.list$size$mean, lwd=2, lty=3)
indiv <- size.cov.sim@results$individuals$N[1,1:6,(simreps+1)]

## ----nocov-10, echo=FALSE, message=FALSE, fig.width=8------------------------------
detfn.constant <- make.ds.analysis(dfmodel=list(~1),
                               key="hn", truncation = trunc.distance)
no.size.cov <- make.simulation(reps=simreps, design=transects, 
                            population.description = pop.desc, 
                            detectability = detect,
                            ds.analysis = detfn.size)
no.size.sim <- run.simulation(no.size.cov, run.parallel = TRUE, max.cores=11)
hist(no.size.sim@results$expected.size[1,1,1:simreps], main="Computed average group size\nNo covariate", xlab="Group size")
abline(v=covariate.list$size$mean, lwd=2, lty=3)
nosize.indiv <- no.size.sim@results$individuals$N[1,1:6,(simreps+1)]


## ----lognor-with, fig.margin=TRUE, echo=FALSE--------------------------------------
ml <- log(12)
sl <- log(3)
group.size <- rlnorm(1000, ml, sl)
mybr <- c(0,5,10,15,20,25,30,35,40,50,60,100,300,500)
hist(group.size, breaks = mybr)
covariate.list$size <- list(distribution="lognormal", meanlog = ml, sdlog = sl)
pop.desc <- make.population.description(covariates = covariate.list, N=200)
cov.params <- list(size = c(0.01))
detect <- make.detectability(scale.param = 30, 
                             cov.param = cov.params, 
                             truncation = trunc.distance)
plot(detect, pop.desc)
detfn.size <- make.ds.analysis(dfmodel=list(~scale(size)),
                               key="hn", truncation = trunc.distance)
tail.size.cov <- make.simulation(reps=simreps, design=transects, 
                                 pop=pop.desc, det=detect,
                                 ds.analysis = detfn.size)
tail.size.cov.sim <- run.simulation(tail.size.cov, run.parallel = TRUE, max.cores=11)
tail.indiv <- tail.size.cov.sim@results$individuals$N[1,1:6,(simreps+1)]


## ---- echo=FALSE-------------------------------------------------------------------
# survey.design <- create.survey.results(size.cov, dht.tables = TRUE)
# my.survey <- data.for.distance(object = survey.design)
# scatter.smooth(my.survey$distance, my.survey$size)


## ---- echo=FALSE, fig.margin=TRUE--------------------------------------------------
hist(tail.size.cov.sim@results$expected.size[1,1,1:simreps], 
     main="Computed average group size\nGroup size covariate", xlab="Group size")
abline(v=exp(ml+sl^2/2), lty=3)
est.tail.size.cov <- tail.size.cov.sim@results$individuals$N[1,1,1:simreps]
hist(est.tail.size.cov, main="Abundance of individuals\n~scale(size)")

## ---- lognor-without, echo=FALSE, fig.margin=TRUE----------------------------------
detfn.constant <- make.ds.analysis(dfmodel=list(~1),
                               key="hn", truncation = trunc.distance)
tail.size.NOcov <- make.simulation(reps=simreps, design=transects, 
                                   pop=pop.desc, det=detect,
                                   ds.analysis = detfn.constant)
tail.size.NOcov.sim <- run.simulation(tail.size.NOcov, run.parallel = TRUE, max.cores=11)
tail.indiv.NO <- tail.size.NOcov.sim@results$individuals$N[1,1:6,(simreps+1)]
hist(tail.size.NOcov.sim@results$expected.size[1,1,1:simreps], main="Computed average group size\nNo covariate", xlab="Group size")
abline(v=exp(ml+sl^2/2), lty=3)

## ---- echo=FALSE-------------------------------------------------------------------
est.tail.nocov <- tail.size.NOcov.sim@results$individuals$N[1,1,1:simreps]
hist(est.tail.nocov,
     main="Estimates of abundance (~1)", xlab="Abundance estimate",
     sub="Dotted vertical line is true value")
abline(v=exp(ml+sl^2/2)*200, lty=3)

