## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(kableExtra)


## ----------------------------------------------------------------------------------------------------
library(Distance)
data(amakihi)


## ---- fig.height=4-----------------------------------------------------------------------------------
hist(amakihi$distance, xlab="Radial distances (m)", main="Amakihi point transect data.")


## ---- fig.height=8-----------------------------------------------------------------------------------
# Plots of covariates
par(mfrow=c(2,2))
# Boxplots by obs
boxplot(amakihi$distance~amakihi$OBs, xlab="Observer", ylab="Distance (m)")
# Boxplots by hour after sunrise
boxplot(amakihi$distance~amakihi$HAS, xlab="Hour", ylab="Distance (m)")
# Plot of MAS vs distance (using dots)
plot(x=amakihi$MAS, y=amakihi$distance, xlab="Minutes after sunrise",
     ylab="Distance (m)", pch=20)
# Plot of HAS vs MAS (using dots)
plot(x=amakihi$HAS, y=amakihi$MAS, xlab="Hours after sunrise",
     ylab="Minutes after sunrise", pch=20)


## ----------------------------------------------------------------------------------------------------
# Adjusting the raw data
# Convert HAS to a factor
amakihi$HAS <- factor(amakihi$HAS)
# Set the reference level 
amakihi$OBs <- relevel(amakihi$OBs, ref="TKP")
amakihi$HAS <- relevel(amakihi$HAS, ref="5")


## ----------------------------------------------------------------------------------------------------
# Rescale MAS 
amakihi$MAS <- scale(amakihi$MAS)
summary(amakihi$MAS)


## ---- fig.height=4, message=FALSE, warning=FALSE-----------------------------------------------------
# Fit model selected by Marques et al (2007)
conv <- convert_units("meter", NULL, "hectare")
amak.hr.obs.mas <- ds(amakihi, transect="point", key="hr", formula=~OBs+MAS, convert.units = conv,
                      truncation=82.5)

# Plot selected model
plot(amak.hr.obs.mas, main="Model with OBs and MAS", pch=".", pdf=TRUE)


## ----------------------------------------------------------------------------------------------------
data(ETP_Dolphin)
head(ETP_Dolphin, n=3)
# Check conversion units
ETP_Dolphin_units


## ---- fig.height=4-----------------------------------------------------------------------------------
# Histogram of distances with lots of intervals
hist(ETP_Dolphin$distance, nclass=50, xlab="Distance (nm)",
     main="Tropical Pacific dolphin survey perpendicular distances")


## ---- fig.height=6-----------------------------------------------------------------------------------
# Boxplots of distances against factor covariates
par(mfrow=c(2,2))
# Search method
boxplot(ETP_Dolphin$distance~ETP_Dolphin$Search.method, xlab="Search method", 
        ylab="Distance (nm)")
# Cue 
boxplot(ETP_Dolphin$distance~ETP_Dolphin$Cue.type, xlab="Cue", ylab="Distance (nm)")
# Beaufort 
boxplot(ETP_Dolphin$distance~ETP_Dolphin$Beauf.class, xlab="Beaufort class", 
        ylab="Distance (nm)")
# Month
boxplot(ETP_Dolphin$distance~ETP_Dolphin$Month, xlab="Month", ylab="Distance (nm)")


## ---- message=FALSE----------------------------------------------------------------------------------
# Fit basic detection functions 
# Half normal
etp.hn <- ds(ETP_Dolphin, key="hn", adjustment=NULL)
# Hazard rate
etp.hr <- ds(ETP_Dolphin, key="hr", adjustment=NULL)
# Compare these fits
knitr::kable(as.data.frame(AIC(etp.hn, etp.hr))) %>%
    kable_styling(bootstrap_options = "condensed", full_width = F)  


## ---- message=FALSE----------------------------------------------------------------------------------
# Add covariates to hazard rate detection function
# Search method (factor)
etp.hr.search <- ds(ETP_Dolphin, key="hr", formula=~factor(Search.method))
# Cue type (factor)
etp.hr.cue <- ds(ETP_Dolphin, key="hr", formula=~factor(Cue.type))
# Beaufort class (factor)
etp.hr.bf <- ds(ETP_Dolphin, key="hr", formula=~factor(Beauf.class))
# Month (factor)
etp.hr.month <- ds(ETP_Dolphin, key="hr", formula=~factor(Month))

# Compare models (using pretty printing)
knitr::kable(summarize_ds_models(etp.hr, etp.hr.search, etp.hr.cue, etp.hr.bf, etp.hr.month),
               caption="ETP dolphin model selection.", digits=3) %>%
       kable_styling(bootstrap_options = "condensed", full_width = F)  


## ----------------------------------------------------------------------------------------------------
# Look at detection function part of the model object
etp.hr.search$ddf


## ---- fig.height=4-----------------------------------------------------------------------------------
# Plot search method detection function
plot(etp.hr.search, pch=".")


## ---- fig.width=4, fig.height=4----------------------------------------------------------------------
# Savannah sparrow 1980
data(Savannah_sparrow_1980)
# Check data
head(Savannah_sparrow_1980, n=3)
# Histogram of distances with lots of bins
hist(Savannah_sparrow_1980$distance, nclass=20, xlab="Distance (m)",
     main="Savannah sparrow radial distances '80")
conversion.factor <- convert_units("meter", NULL, "hectare")


## ---- eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------
# Fit different detection functions, truncation at 55m
# Half-normal 
Savannah_sparrow_1980.hn <- ds(data=Savannah_sparrow_1980, key="hn", adjustment="cos", truncation=55,
                 transect="point", convert.units=conversion.factor)
# Hazard
Savannah_sparrow_1980.hr <- ds(data=Savannah_sparrow_1980, key="hr", adjustment="cos", truncation=55,
                 transect="point", convert.units=conversion.factor)

# Half-normal with pasture covariate
Savannah_sparrow_1980.hn.region <- ds(data=Savannah_sparrow_1980, key="hn", truncation=55,
                        transect="point", convert.units=conversion.factor,
                        formula=~Region.Label)
# Hazard with pasture covariate
Savannah_sparrow_1980.hr.region <- ds(data=Savannah_sparrow_1980, key="hr", truncation=55,
                        transect="point", convert.units=conversion.factor,
                        formula=~Region.Label)
# Select between these models
AIC(Savannah_sparrow_1980.hn, Savannah_sparrow_1980.hr, Savannah_sparrow_1980.hn.region, Savannah_sparrow_1980.hr.region)


## ---- fig.cap="Note different PDF shapes caused by the pasture covariate."---------------------------
# Plot results of selected model
plot(Savannah_sparrow_1980.hn.region, pch=".", pdf=TRUE)


## ----------------------------------------------------------------------------------------------------
# Summarise results for selected model
summary(Savannah_sparrow_1980.hn.region)


## ----------------------------------------------------------------------------------------------------
# Savannah sparrow 1981
data(Savannah_sparrow_1981)
conversion.factor <- convert_units("meter", NULL, "hectare")


## ---- eval=T, message=FALSE, warning=FALSE-----------------------------------------------------------
# Fit alternative models 
# Half-normal detection function, truncation 55m 
Savannah_sparrow_1981.hn <- ds(data=Savannah_sparrow_1981, key="hn", adjustment="cos", truncation=55,
                 transect="point", convert.units=conversion.factor)
# Hazard rate
Savannah_sparrow_1981.hr <- ds(data=Savannah_sparrow_1981, key="hr", adjustment="cos", truncation=55,
                 transect="point", convert.units=conversion.factor)
# Half normal with pasture
Savannah_sparrow_1981.hn.region <- ds(data=Savannah_sparrow_1981, key="hn", truncation=55,
                        transect="point", convert.units=conversion.factor,
                        formula=~Region.Label)
# Hazard rate with pasture
Savannah_sparrow_1981.hr.region <- ds(data=Savannah_sparrow_1981, key="hr", truncation=55,
                        transect="point", convert.units=conversion.factor,
                        formula=~Region.Label)
# Compare models
AIC(Savannah_sparrow_1981.hn, Savannah_sparrow_1981.hr, Savannah_sparrow_1981.hn.region, Savannah_sparrow_1981.hr.region)


## ---- fig.height=4, fig.cap="Stronger influence of pasture covariate seen here."---------------------
# Plot results of selected model
plot(Savannah_sparrow_1981.hn.region, pch=".", pdf=TRUE, main="1981 sparrows by pasture.")

# Summary of results for selected model
summary(Savannah_sparrow_1981.hn.region)

