## ----setup, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(knitr)
library(kableExtra)


## ---------------------------------------------------------------------------
library(Distance)
data(amakihi)


## ---- fig.height=4----------------------------------------------------------
hist(amakihi$distance, xlab="Radial distances (m)", main="Amakihi point transect data.")


## ---- fig.height=8----------------------------------------------------------
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


## ---------------------------------------------------------------------------
# Adjusting the raw data
# Convert HAS to a factor
amakihi$HAS <- factor(amakihi$HAS)
# Set the reference level 
amakihi$OBs <- relevel(amakihi$OBs, ref="TKP")
amakihi$HAS <- relevel(amakihi$HAS, ref="5")


## ---- fig.height=5, message=FALSE, warning=FALSE----------------------------
# Fit model selected by Marques et al (2007)
conv <- convert_units("meter", NULL, "hectare")
amak.hr.obs.mas <- ds(amakihi, transect="point", key="hr", formula=~OBs+MAS, convert.units = conv,
                      truncation=82.5)

# Plot selected model

plot(amak.hr.obs.mas, showpoints=FALSE, main="Amakihi Observer and Minutes", pdf=TRUE)
sfzero <- data.frame(OBs="SGF", MAS=0)
sf180 <- data.frame(OBs="SGF", MAS=180)
t1zero <- data.frame(OBs="TJS", MAS=0)
t1180 <- data.frame(OBs="TJS", MAS=180)
t2zero <- data.frame(OBs="TKP", MAS=0)
t2180 <- data.frame(OBs="TKP", MAS=180)

add_df_covar_line(amak.hr.obs.mas, data=sfzero, lty=1, lwd=2,col="blue", pdf=TRUE)
add_df_covar_line(amak.hr.obs.mas, data=sf180, lty=2, lwd=2,col="blue", pdf=TRUE)
add_df_covar_line(amak.hr.obs.mas, data=t1zero, lty=1,lwd=2,col="red", pdf=TRUE)
add_df_covar_line(amak.hr.obs.mas, data=t1180, lty=2, lwd=2,col="red", pdf=TRUE)
add_df_covar_line(amak.hr.obs.mas, data=t2zero, lty=1,lwd=2,col="green", pdf=TRUE)
add_df_covar_line(amak.hr.obs.mas, data=t2180, lty=2, lwd=2,col="green", pdf=TRUE)

legend("topright", legend=c("SF, minutes=0",
                            "SF, minutes=180",
                            "TS, minutes=0",
                            "TS, minutes=180",
                            "TP, minutes=0",
                            "TP, minutes=180"),
       title="Covariate combination: observer and minutes",
       lty=rep(c(1,2),times=3), lwd=2, col=rep(c("blue","red","green"), each=2))


## ----distp, eval=FALSE------------------------------------------------------
## p_dist_table(amak.hr.obs.mas, proportion = TRUE)


## ----distp2, echo=FALSE-----------------------------------------------------
kable(p_dist_table(amak.hr.obs.mas, bins=seq(0, 0.6, 0.1), proportion = TRUE),
      digits = 3,
      caption="Distribution of $P_a(z_i)$ from preferred model when w=82.5") %>%
  kable_styling(full_width=FALSE) %>%
  row_spec(2, bold=TRUE, color="white", background="blue")


## ----moretrunc, echo=FALSE, eval=TRUE---------------------------------------
amak.hr.obs.mas.70 <- ds(amakihi, transect="point", key="hr", formula=~OBs+MAS, convert.units = conv,
                      truncation=70)
kable(p_dist_table(amak.hr.obs.mas.70, bins=seq(0, 0.6, 0.1), proportion = TRUE),
      digits = 3,
      caption="Distribution of $P_a(z_i)$ from preferred model when w=70") %>%
  kable_styling(full_width=FALSE) %>%
  row_spec(2, bold=TRUE, color="white", background="blue")


## ----estfigs, fig.cap="Density estimates for amakihi data set under a variety of detection function models and truncation distances.", echo=FALSE, out.width="80%"----
include_graphics("estimate-comparisons.png")


## ---------------------------------------------------------------------------
data(ETP_Dolphin)
head(ETP_Dolphin, n=3)
# Check conversion units
ETP_Dolphin_units[,1:2]


## ---- fig.height=4----------------------------------------------------------
# Histogram of distances with lots of intervals
hist(ETP_Dolphin$distance, nclass=50, xlab="Distance (nm)",
     main="Tropical Pacific dolphin survey perpendicular distances")


## ---- fig.height=6----------------------------------------------------------
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


## ---- message=FALSE---------------------------------------------------------
# Fit basic detection functions 
# Half normal
etp.hn <- ds(ETP_Dolphin, key="hn", adjustment=NULL)
# Hazard rate
etp.hr <- ds(ETP_Dolphin, key="hr", adjustment=NULL)
# Compare these fits
knitr::kable(as.data.frame(AIC(etp.hn, etp.hr))) %>%
    kable_styling(bootstrap_options = "condensed", full_width = F)  


## ---- message=FALSE---------------------------------------------------------
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


## ---------------------------------------------------------------------------
# Look at detection function part of the model object
etp.hr.search$ddf


## ---- fig.height=4----------------------------------------------------------
# Plot search method detection function
plot(etp.hr.search, pch=".")


## ---- colourful, fig.height=4, fig.cap="Detection function with cue type as covariate."----
plot(etp.hr.cue, main="ETP dolphin survey", showpoints=FALSE)
add_df_covar_line(etp.hr.cue, data = data.frame(Cue.type=1), col='red', lwd=2, lty=1)
add_df_covar_line(etp.hr.cue, data = data.frame(Cue.type=2), col='blue', lwd=2, lty=1)
add_df_covar_line(etp.hr.cue, data = data.frame(Cue.type=3), col='green', lwd=2, lty=1)
add_df_covar_line(etp.hr.cue, data = data.frame(Cue.type=4), col='purple', lwd=2, lty=1)
legend("topright", legend=c("Birds","Splashes","Unspecified","Floating objects"),
       col=c("red", "blue", "green", "purple"), lwd=2, title = "Cue type")


## ---- fig.width=4, fig.height=4---------------------------------------------
# Savannah sparrow 1980
data(Savannah_sparrow_1980)
# Check data
head(Savannah_sparrow_1980, n=3)
# Histogram of distances with lots of bins
hist(Savannah_sparrow_1980$distance, nclass=20, xlab="Distance (m)",
     main="Savannah sparrow radial distances '80")
conversion.factor <- convert_units("meter", NULL, "hectare")


## ---- eval=TRUE, message=FALSE, warning=FALSE-------------------------------
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


## ---- fig.cap="Note different PDF shapes caused by the pasture covariate."----
# Plot results of selected model
plot(Savannah_sparrow_1980.hn.region, pch=".", pdf=TRUE)


## ---------------------------------------------------------------------------
# Summarise results for selected model
summary(Savannah_sparrow_1980.hn.region)


## ---------------------------------------------------------------------------
# Savannah sparrow 1981
data(Savannah_sparrow_1981)
conversion.factor <- convert_units("meter", NULL, "hectare")


## ---- eval=T, message=FALSE, warning=FALSE----------------------------------
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


## ---- fig.height=6, fig.width=8, fig.cap="Stronger influence of pasture covariate seen here."----
pastures <- unique(Savannah_sparrow_1981$Region.Label)
plot(Savannah_sparrow_1981.hn.region, showpoints=FALSE, 
     main="Savannah sparrows with pasture covariate", pdf=TRUE)
k <- 1
for (i in pastures) {
  k <- k+1
  add_df_covar_line(Savannah_sparrow_1981.hn.region, 
                    data=data.frame(Region.Label=as.character(i)),
                    lty=1, col=k, lwd=3, pdf=TRUE)
}
legend("topright", legend=tolower(as.character(pastures)), 
       col=2:k, lwd=2, title = "Pastures")
text(-2,0.038, cex=0.9, pos=4,
     expression(widehat(sigma[p])==plain(exp)(widehat(beta[0]) + widehat(beta[1]) %.% p[1] + widehat(beta[2]) %.% p[2] + widehat(beta[3]) %.% p[3])))
library(plotrix)
parms <- data.frame(est=c(2.944, 0.736, 0.166, 0.271),
                    se=c(0.111, 0.373, 0.153, 0.179))
rownames(parms) <- c("b0", "b1", "b2", "b3")
addtable2plot(2, 0.027, parms, bty="o",
              display.rownames=TRUE,hlines=TRUE, cex=0.8,
              xpad=0.4, vlines=FALSE,title="Parameter estimates")
# Summary of results for selected model
summary(Savannah_sparrow_1981.hn.region)

