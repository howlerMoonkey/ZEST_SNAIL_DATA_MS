######################################################
######################################################
######## ZANZIBAR ZEST PROJECT SNAIL DATA ANALYSIS
######################################################
setwd("C:\\Users\\murr\\Desktop\\CURRENT_MANUSCRIPTS_2020\\ZEST_SNAIL_DATA_MS_2020\\ZANZIBAR_SNAIL_DATA_MS")
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
## Metadata

#data<-read.csv("ZEST_snail_data_2021-12-31.csv")
#list <- lapply(data, class)
#as.data.frame(list)
#write.csv(list, "zest_list.csv")
#zest <- read.csv("ZEST_list.csv")  
#names(zest)
#zest_meta <-
#  zest %>%
#  gather(names, values, orig_row_rec:INF.2012)
#write.csv(zest_meta, "zest_meta.csv")
################################################
################################################
## Prepare dataset for analysis

data<-read.csv("ZEST_snail_data_2021-01-11.csv")
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
         tot_pres_other_snail = as.factor(tot_pres_other_snail),
         Infsnails = as.factor(as.character(Infsnails)),
         Sprayed_Steffi = as.factor(as.character(Sprayed_Steffi)),
         sprayed_when_Bul_pres = as.factor(as.character(sprayed_when_Bul_pres)),
         Cond = scale(Cond, center = 0),
         TDS = scale(TDS, center = 0),
         Temp = scale(Temp, center = 0),
         Sal = scale(Sal, center = 0),
         Oxy = scale(Oxy, center = 0),
         PH = scale(PH, center = 0))%>%
    group_by(EMU_site_irn) %>% 
    mutate(latitude = mean(latitude, na.rm = TRUE),
         longitude = mean(longitude, na.rm = TRUE))

################################
## select relevant columns
data_ed <- data1 %>%
  dplyr::select(Bul_pres_ed, Cond, Oxy, Date_appl_ed, Date_ed, Depth_ed, N_Duration_min,
                duration, EMU_site_irn, Flow_ed, Island, Infsnails, Infsnails_no,
                month, Year, N_Duration_min, niclo_visit, no_collected, include,
                Numberofcollectors, PH, pres_non_bulinids, Ricecombined_Steffi, Sal,
                Seas, Season, Shehia_cleaned, site_dry_ed, Sprayed_Steffi, niclo_order, 
                sprayed_when_Bul_pres, sprayed_but_no_bul, TDS, Temp, Village, Visit, 
                wbtype_ed, wbtype_fin,Wbtype2_ed, Wlevel_ed, Ricecombined_Steffi, Wbname_fin,
                Bglob, Bfor, Lym, Mel, Cleo, Pila, Lan, Cera, Neri, Thiara,
                mud, sand, gravel, rock, concrete, peat, roots, include, tot_pres_other_genera,
                Lilies, rushes, rice, waterhyacinth, macrophyteweed, palmfronds,
                sedge, grass, Polygonum, none, Otherveg, monthly_precip, latitude, longitude,
                Dishes, bathing, Clothes, Carbike, Collectingwater, fording, swimmingplaying, 
                fishing,Ricecult, otherfarming, sanitation, ablution, Otheract, 
                month_survey_order, INF.ALL, TOT.ALL, INF.2017,INF.2016,INF.2015,INF.2014,
                INF.2013,INF.2012, TOT.2017,TOT.2016, TOT.2015, TOT.2014, TOT.2013, TOT.2012)
dim(data_ed)   
dim(data1)
data_2<- data_ed %>% 
  filter(include !="no") ## remove duplicate site visits
data_ed <- data_2

# Remove NA values in duration, as well as predictors
#data_short <- na.omit(data_ed) ## omit all na values = ## removed about 1000 rows- was 3587, now 2471
data_short<- data_ed %>% 
  filter(!is.na(duration))

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

## randomised selection of 15 shehias per island of 45 tot for ZEST:
## aim- testing different treatment arms- MDA only, MDA + snail control, MDA + WASH
## monthly surveys in waterbodies- with gaps- over 5yrs (late 2011- early 2017)
## We have Bulinus count and presence/absence data (& for infected Bulinus
## other data- water chemistry, eology- substrate and water plants, water activities etc

## data on molluscidiing with Niclosamide- 
## sites sprayed (gen) when Bulinus pres, also when not (also occas not sprayed WHEN pres)
## no control sites- how best to assess any effect
## short term and long term impacts if any

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

## Anova
Anova.glmmTMB(poisson_m1) ## cant find function
glmmTMB:::Anova.glmmTMB(nbinom2_m1) 
##Error in Anova.glmmTMB(Niclo_m1) : could not find function "Anova.glmmTMB"
glmmTMB::Anova(nbinom2_m1)
#Error: 'Anova' is not an exported object from 'namespace:glmmTMB'
Anova.glmmTMB(nbinom2_m1)
#Error: 'Anova' is not an exported object from 'namespace:glmmTMB'

car::Anova(nbinom2_m1) ## working. why above not working

## model diganostics
## all- no deviation in QQ plot but deviation in residuals
summary(nbinom1_m1)                    
sim_residualsnbinom1_m1 <- DHARMa::simulateResiduals(nbinom1_m1, 1000)  
plot(sim_residualsnbinom1_m1)  
testZeroInflation(sim_residualsnbinom1_m1)

######################################### 
### 29 DEC 2020 interactions- site type and flow- eg river vs pond and flow
## prevalence models- add data; do as binomial models
## reshape prevalence dataset

data <- read.csv("Ung_2011_tot.csv")
data_spread <- data %>%
  spread(age_group, tot.2011)
write.csv(data_spread, "ung_tot-2011.csv")
## write function as for loop with purrr to iterate for each year

##and merge with main dataset
data<-read.csv("ZEST_snail_data_key.csv")
data2 <- read.csv("ZEST_Prevalence_Data_2020-12-31.CSV")
merga <- merge(data, data2, by = "Shehia_cleaned", all.x = T)
write.csv(merga, "merge_zest.csv")

######################################### 
################ MODELS BY DATA TYPE- do presence apbse & count models

################ FULL MODEL BINOMIAL

Niclo_m1 <- glmmTMB(Bul_pres_ed ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                      month + Cond  + pres_non_bulinids + monthly_precip + Temp + Depth_ed  + 
                      wbtype_ed + Flow_ed + Wlevel_ed + Infsnails + 
                      offset(log(duration)),
                    data=data_ed,
                    family=binomial)

## originally converging- not now
######################################

### 2 JAN prevalence models
prev_m1 <- glmmTMB(INF.ALL/TOT.ALL ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                     wbtype_fin + Bul_pres_ed + no_collected + month + sanitation +
                     offset(log(duration)),
                   weights = TOT.ALL,
                   data=data_ed,
                   family=binomial)

## remove duration here- add offset of niclosamide duration to niclosamide models
## results in non-convergence
## significant deviation ks test (much more so with random effects, log duration removed)

#####################################
## Niclosamide data model
## add offset of duration of niclosamiding also? offset(log(N_Duration_min))
## many other factors- method, volume, no of sprayers etc
Niclo_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                        niclo_visit + 
                        offset(log(duration)),
                      data=data_ed,
                      family=nbinom2)

#####################################
##model for substrate data
Niclo_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                      Lilies + rushes + rice + waterhyacinth +  macrophyteweed + palmfronds +
                      grass + Polygonum +  
                      offset(log(duration)),
                    data=data_ed,
                    family=nbinom2)

#####################################
## model for plant data
## add interactions eg waterbody type?
Niclo_m1 <- glmmTMB(Bul_pres_ed ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                      Lilies + rushes + rice + waterhyacinth +  macrophyteweed + palmfronds +
                      grass + Polygonum +  
                      offset(log(duration)),
                    data=data_ed,
                    family=binomial)

#####################################
## model for water activities (more relevant for prevalence model but may illuminate 
## where are transmission sites- where is more usage- therefore likely transmission)
Niclo_m1 <- glmmTMB(Bul_pres_ed ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                      Dishes +bathing + Clothes + Carbike + Collectingwater + fording +
                      swimmingplaying + fishing +Ricecult + otherfarming + sanitation + 
                      ablution + Otheract +
                    offset(log(duration)),
                    data=data_ed,
                    family=binomial)

## add interaction of waterbody tpye & water activity eg appears-
## swimming/collecting water in large ponds and ablution not in this wbtype- 
## poss more likely in streams w tree cover- privacy, flowing water, avoid wbs where water collected etc
## add rice-steffi?

Niclo_m1 <- glmmTMB(no_collected ~ (1|Island/Shehia_cleaned/EMU_site_irn) + (1|month/Date_ed) +
                      Collectingwater + swimmingplaying + sanitation + 
                      ablution + wbtype_fin*ablution + wbtype_fin*Collectingwater + 
                      wbtype_fin*swimmingplaying +
                      offset(log(duration)),
                    data=data_ed,
                    family=nbinom2)
## ablution and wbtype_fin = significant positive correlation

#####################################
## model for water chemistry data


#####################################
## model for physical water data

glmmTMB::Anova(Niclo_m1)
Anova.glmmTMB(Niclo_m1)
glmmTMB::Anova.glmmTMB(Niclo_m1) 

car::Anova(Niclo_m1)
summary(Niclo_m1)
sim_residualsNiclo_m1 <- DHARMa::simulateResiduals(Niclo_m1, 1000)  
plot(sim_residualsNiclo_m1)  
testZeroInflation(sim_residualsNiclo_m1) 

## add no of collectors as offset?
