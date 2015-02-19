microsome analyses
========================================================
Statistical analysis and results presentation and discussion for total concentrations 
of parent active 
ingredients (atrazine, imidacloprid, fipronil, triadimenon, pendimethalin).

```{r eval=TRUE, echo=FALSE}

if(Sys.info()[4]=="DC2626UTPURUCKE"){
  frogsoildir <- "C:\\stp_drop\\Dropbox\\amphib_dermalexposure\\DATA\\PLEvsOLS2014\\"
}

if(Sys.info()[4]=="stp-air.local"){
  frogsoildir <- path.expand("~/Dropbox/amphib_dermalexposure/DATA/PLEvsOLS2014/")
}

library(MASS)
##################
#the data sets
##################
#everything
frog.soil <- read.table(paste(frogsoildir,"RDATA.csv",sep=""), header = TRUE, sep = ",")
colnames(frog.soil)
class(frog.soil$Day)
frog.soil$Day <- as.factor(frog.soil$Day)
class(frog.soil$Row)
frog.soil$Row <- as.factor(frog.soil$Row)
class(frog.soil$Column)
class(frog.soil$Pesticide)
class(frog.soil$SoilType)
class(frog.soil$Soil)
class(frog.soil$BB)
class(frog.soil$Weight)
class(frog.soil$Total) <- as.factor(frog.soil$Total)
class(frog.soil$Formulation) <- as.factor(frog.soil$Formulation)
class(frog.soil$Parent) <- as.factor(frog.soil$Parent)
frog.soil$bowlbcf <- frog.soil$BB/frog.soil$Soil
frog.soil$surface_area_total <- 1.131 * frog.soil$Weight^0.579
frog.soil$surface_area_footprint <- 0.425 * frog.soil$Weight^0.85
```

```{r eval=TRUE, echo=FALSE}
dim(frog.soil)
colnames(frog.soil)

#total only but with formulations
frog.soil.total <- frog.soil[which(frog.soil$Total==T),]
#View(frog.soil.total)
dim(frog.soil.total)
colnames(frog.soil.total)

#Formulation is unbalanced so drop formulation, parent ais only
frog.soil.total.ai <- frog.soil.total[which(frog.soil.total$Formulation==F & 
                                              frog.soil.total$Parent==T),]
#View(frog.soil.total.ai)
dim(frog.soil.total.ai)
colnames(frog.soil.total.ai)

#drop imidacloprid from ai analysis
frog.soil.total.noimid <- frog.soil.total.ai[-which(frog.soil.total.ai$Pesticide=="Imid"),]
#View(frog.soil.total.noimid)
dim(frog.soil.total.noimid)
colnames(frog.soil.total.noimid)

##################
#analyses
##################
#randomized block design for bowlbcfs
bowlbcf.total.aov <- aov(bowlbcf ~ Pesticide + SoilType + Formulation + surface_area_total, 
                           data = frog.soil.total)
summary(bowlbcf.total.aov)
bowlbcf.total.ai.aov <- aov(bowlbcf ~ Pesticide + SoilType + surface_area_total, 
                              data = frog.soil.total.ai)
summary(bowlbcf.total.ai.aov)
bowlbcf.total.noimid.aov <- aov(bowlbcf ~ Pesticide + SoilType + surface_area_total, 
                                  data = frog.soil.total.noimid)
summary(bowlbcf.total.noimid.aov)

#randomized block design for body burdens
bodyburden.total.aov <- aov(BB ~ Pesticide + SoilType + Soil + Formulation + surface_area_total, 
                           data = frog.soil.total)
summary(bodyburden.total.aov)
bodyburden.total.ai.aov <- aov(BB ~ Pesticide + SoilType + Soil + Formulation + surface_area_total, 
                            data = frog.soil.total.ai)
summary(bodyburden.total.ai.aov)
bodyburden.total.noimid.aov <- aov(BB ~ Pesticide + SoilType + Soil + Formulation + surface_area_total, 
                               data = frog.soil.total.noimid)
summary(bodyburden.total.noimid.aov)

# 5 t-tests for total parents
pesticides <- unique(frog.soil.total.ai$Pesticide)
length(pesticides)
for(pesticide in pesticides){
  frog.soil.total.ai.tt <- frog.soil.total.ai[which(frog.soil.total.ai$Pesticide==pesticide),]
  ols.bcf <- frog.soil.total.ai.tt[which(frog.soil.total.ai.tt$SoilType=="OLS"),]
  ple.bcf <- frog.soil.total.ai.tt[which(frog.soil.total.ai.tt$SoilType=="PLE"),]
  soil.t.test <- t.test(ols.bcf$bowlbcf,ple.bcf$bowlbcf)
  print("ols")
  print(ols.bcf$bowlbcf)
  print("ple")
  print(ple.bcf$bowlbcf)
  print(pesticide)
  print(soil.t.test)
}

#switch to aic instead of significance testing for regressing on body burdens
#everything
total.lm <- lm(BB ~ Pesticide + Soil + Formulation + SoilType + surface_area_footprint, 
               data = frog.soil.total)
stepAIC(total.lm)

#without imidacloprid
colnames(frog.soil.total.noimid)
#View(frog.soil.total.noimid)
noimd.lm <- lm(BB ~ Pesticide + Soil + Formulation + SoilType + surface_area_footprint, 
               data = frog.soil.total.noimid)
stepAIC(noimd.lm)

colnames(frog.soil.total.noimid)
noimd.lm <- lm(BB ~ Pesticide + Soil + Formulation + SoilType + surface_area_footprint, 
               data = frog.soil.total.noimid)
stepAIC(noimd.lm)

#using a mean soil instead of the bowl ratios for bcfs


#lump ai and formulation



noimd.aov <- aov(BB ~ Pesticide + Soil + Formulation + Soil.1 + Weight + surface_area_footprint, data = frog.soil.total.noimid)
summary(noimd.aov)
TukeyHSD(noimd.aov)

#formulation factor does not make sense unless we hange formulation name to equal active ingredient nems=e

