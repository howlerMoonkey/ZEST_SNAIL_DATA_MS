######################################################
######################################################
######## ZANZIBAR ZEST PROJECT SNAIL DATA ANALYSIS
######################################################
setwd("C:\\Users\\murr\\Documents")
getwd()

library(DHARMa)
library(bbmle)
library(effects)
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(ggthemes) ## not avail for R version 4
library(DHARMa) # for model diagnostics
library(emmeans)  # for post hoc test
library(DataExplorer)   # for data exploration
library(gridExtra)
library(lme4)
library(car)
library(reshape2)
library(ggfortify)
library(boot)
library(aod)
library(mgcv)
library(pscl)
source("HighstatLibV10.R")
#source("C:\\Users\\murr\\Documents\\HighstatLibV10.R")

#############################################
#############################################

data<-read.csv("ZEST_snail_data_2021-01-01.csv")

dim(data)
names(data)
# mutate- transform to as.factor etc where relevant
data1 <- data %>%
  mutate(Date_ed = lubridate::dmy(Date_ed), tz = "Africa/Tanzania",
         Date_appl_ed = lubridate::dmy(Date_appl_ed), tz = "Africa/Tanzania",
         duration = as.numeric(as.character(duration)),
         N_Duration_min = as.numeric(as.character(N_Duration_min)),
         month = as.factor(month),
         Year = lubridate::year(Date_ed),
         EMU_site_irn = as.factor(EMU_site_irn), 
         Visit = as.factor(Visit),
         niclo_visit = as.factor(niclo_visit),
         Depth_ed = as.factor(Depth_ed), ## orig < 0.5, 0.5-1m- transformed to category- 1, 2, 3, etc
         Flow_ed = as.factor(Flow_ed), ## orig < still, 0.5, 0.5-1m- transformed to category- 1, 2, 3, etc
         Wlevel_ed = as.factor(Wlevel_ed),
         Bul_pres_ed = as.factor(as.character(Bul_pres_ed)),
         pres_non_bulinids = as.factor(as.character(pres_non_bulinids)),
         Infsnails = as.factor(as.character(Infsnails)),
         Sprayed_Steffi = as.factor(as.character(Sprayed_Steffi)),
         sprayed_when_Bul_pres = as.factor(as.character(sprayed_when_Bul_pres)),
         Cond = scale(Cond, center = 0),
         TDS = scale(TDS, center = 0),
         Temp = scale(Temp, center = 0),
         Sal = scale(Sal, center = 0),
         Oxy = scale(Oxy, center = 0),
         PH = scale(PH, center = 0))

################################
## select relevant columns
data_ed <- data1 %>%
  dplyr::select(Bul_pres_ed, Cond, Date_appl_ed, Date_ed, Depth_ed,
                duration, EMU_site_irn, Flow_ed, Island, Infsnails, Infsnails_no,
                Oxy, month, N_Duration_min, niclo_visit, no_collected, include,
                Numberofcollectors, PH, pres_non_bulinids, Ricecombined_Steffi, Sal,
                Seas, Season, Shehia_cleaned, site_dry_ed, Sprayed_Steffi, 
                sprayed_when_Bul_pres, TDS, Temp, Village, Visit, wbtype_ed, wbtype_fin,
                Wbtype2_ed, Wlevel_ed, Year, Ricecombined_Steffi, Wbname_fin,
                Bglob, Bfor, Lym, Mel, Cleo, Pila, Lan, Cera, Neri, Thiara,
                mud, sand, gravel, rock, concrete, peat, roots, include,
                Lilies, rushes, rice, waterhyacinth, macrophyteweed, palmfronds,
                sedge, grass, Polygonum, none, Otherveg, monthly_precip, niclo_order, 
                month_survey_order, INF.ALL, TOT.ALL, INF.2017,INF.2016,INF.2015,INF.2014,
                INF.2013,INF.2012, TOT.2017,TOT.2016, TOT.2015, TOT.2014, TOT.2013, TOT.2012)
dim(data_ed)   
dim(data1)

# Remove NA values in duration, as well as predictors
#data_short <- na.omit(data_ed) ## omit all na values = ## removed about 1000 rows- was 3587, now 2471
data_short<- data_ed %>% 
  filter(!is.na(duration))

data_2<- data_ed %>% 
  filter(include !="no")
data_ed <- data_2

############################################
##### missing data

plot_missing(data_ed) + theme_few()
## remove- only niclo details- barrel drip etc; and Oxy which excluded anyway as correlation with salinity
## model with data fields where tlo missing data
data_ed$dur_missing <- ifelse(is.na(data_ed$duration), 1, 0)

miss_model <- glm(dur_missing ~ Temp + no_collected + month + Numberofcollectors + 
                    Wbtype2_ed +  Bul_pres_ed,
                  family = binomial,
                  data = data_ed)
summary(miss_model)
## negative correlation of no of collectors- remove term. otherwise no strong effects
#############################################

#' Overview over the data structure
plot_str(data_ed) + theme_few()
test <- as.data.frame(plot_str(data_ed))
write.csv(test, "test.csv")
introduce(data_ed)

#' Overview over predictor variables
plot_histogram(data_ed)  # Outliers in no_collectors

#' These outliers will drive results. Look at them critically
# 
##########################################
#' Are there any variance inflation factors (multicollinearity)? 
#' Check using a function from Zuur et al. 2010

pairs(data_ed[,c("Cond", "TDS", "Sal", "Oxy", "PH",
                 "Flow_ed", "Wlevel_ed", "Depth_ed","Temp","Depth_ed")],
      lower.panel = panel.cor)

corvif(data.table::as.data.table(data_ed)[, 
   c("Cond", "TDS", "Sal", "Oxy", "PH","Flow_ed", 
     "Wlevel_ed", "Depth_ed","Temp","Depth_ed"), with=FALSE])

# high GVIF values (10). For values of higher than 4, 
#only one of the two variables should be used in models, 
# to avoid multicollinearity.

## TDS and Conductvity- GVIF value of 5- use conductivity only

#' Look at general correlation matrix
plot_correlation(na.omit(data_ed), maxcat = 5L)  

#######################################
#######################################
##############STUDY DESIGN
## We have Bulinus count data, and alo Bulinus presence/absence data 
## (also infected snail counts and presence/absence)
## monthly surveys- with gaps- over 5 years (late 2011- early 2017)
## also data on molluscidiing with Niclosamide- 
## sites sprayed (generally) when Bulinus present- 
## no control sites- how best to assess any effect
## other data- water chemistry, water activities etc

#######################################

#' ### Make a GLMM

#' Make a maximum model. glmmTMB is a new package by Ben Bolker, that fits models faster, 
#' and allows to include arguments
# to account for zero-inflation, if needed.
# Include the sampling duration as an offset. It needs to be specified as log(), since we are 
# using a family distribution
# with log link (nbinom2)

## compare family for model
## negative binomial likely approperiate- count data, overdispersed (Bolker et al, 2009)
## not converging unless single terms- run only with monthly precipitation

## neg binomial 2
nbinom2_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                        monthly_precip + 
                        offset(log(duration)),
                      data=data_ed,
                      family=nbinom2)
## neg binomial 1
nbinom1_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                        monthly_precip + 
                        offset(log(duration)),
                      data=data_ed,
                      family=nbinom1)

## poisson
poisson_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                        monthly_precip + 
                        offset(log(duration)),
                      data=data_ed,
                      family=poisson)

## compare models
Anova.glmmTMB(poisson_m1) ## cant find function
glmmTMB:::Anova.glmmTMB(nbinom2_m1)
glmmTMB:::Anova.glmmTMB(nbinom1_m1)
glmmTMB:::Anova.glmmTMB(poisson_m1)

AICtab(nbinom2, nbinom1, poisson_m1)
#Error in UseMethod("logLik") : 
#no applicable method for 'logLik' applied to an object of class "function"
AICc(nbinom2, nbinom1, poisson_m1, modnames = NULL,
     second.ord = TRUE, nobs = NULL, sort = TRUE)
summary(nbinom2)$AICtab ## not working

## model diganostics
## all- no deviation in QQ plot but deviation in residuals
summary(nbinom1_m1)                    
sim_residualsnbinom1_m1 <- DHARMa::simulateResiduals(nbinom1_m1, 1000)  
plot(sim_residualsnbinom1_m1)  
testZeroInflation(sim_residualsnbinom1_m1)

###############################################
## inbuilt tests- test generic summary stats etc in DHARMa
#https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#zero-inflation-k-inflation-or-deficits
testResiduals(sim_residualsNiclo_m1)
testTemporalAutocorrelation(sim_residualsNiclo_m1)
testSpatialAutocorrelation(sim_residualsNiclo_m1)
testDispersion(sim_residualsNiclo_m1)
plotResiduals(sim_residualsNiclo_m1, data_ed$Temp)
testGeneric(sim_residualsNiclo_m1)
## Checking for missing predictors or quadratic effects
sim_residualsNiclo_m1$scaledResiduals ## checking for missing predictors or quadratic effects
plot(sim_residualsNiclo_m1, quantreg = T) ## general plot
## difficult to see here- becomes clear if plot against the environment
par(mfrow = c(1,2))
plotResiduals(sim_residualsNiclo_m1, Niclo_m1$Shehia_cleaned)
plotResiduals(sim_residualsNiclo_m1, Niclo_m1$Date_ed)
plotResiduals(sim_residualsNiclo_m1, Niclo_m1$EMU_site_irn)
plotResiduals(sim_residualsNiclo_m1, Niclo_m1$Monthly.precip)
countOnes <- function(x) sum(x == 1) # testing for no of 1s
testGeneric(sim_residualsNiclo_m1, summary = countOnes, alternative = "greater") # 1-inflation
ziformula = ~1 ## zero inflation testing

######################################### 14 dec

## interactions- site type and flow- eg river vs pond and flow
## prevalence models- add data
## do as binomial models

### 29 DEC 2020
## reshape prevalence data



