
## Clearing the global environment
rm(list=ls(all=TRUE))
rm(list=ls())

## Setting up the working directory
getwd()
setwd("C:/Users/yismail/Desktop/RStudio Work Directory/Plant Density")

## Loading packages
install.packages("dfoptim")
install.packages("lmerTest")
install.packages("read")

library(dfoptim)
library(lmerTest)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(devtools)
library(ggbiplot)
library(xlsx)
library(readxl)
library(lattice)
library(car)
library(agricolae)
library(FactoMineR)
library(psycho)
library(rcompanion)
library(multcompView)
library(emmeans)
library(lme4)
library(psych)


## Loading in the data
data <- read.csv("phenotypes2020A.csv", header=T, stringsAsFactors=F)

## Setting variables as factors for Genotype, block, rep, and env

data$Ent=as.factor(data$Genotype)
data$Female=as.factor(data$Female)
data$Male=as.factor(data$Male)
data$PED=as.factor(data$PED)
data$Env=as.factor(data$Env)
data$Rep=as.factor(data$Rep)
data$Block=as.factor(data$Block)

str(data)
traits <- colnames(data)[8:ncol(data)]

for(trait in traits){
  pdf(paste0("plots_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(data[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across environments
  stripchart(data[,trait] ~data$Env,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Environment",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_", trait, ".txt"))}
  
# Summary statistics of the trait
summary <- summary(data[,trait], )
out <- capture.output(summary)
cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
# Test for normality
normality <- shapiro.test(data[,trait])
out <- capture.output(normality)
cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)

# Find Outliers
1- Using Histogram, Scatter plot, and boxplot


# Run an ANOVA (switched : in Env:Rep and Rep:Block for nesting)
model <- lm(get(trait) ~ Genotype + Env + Env/Rep + Env/Rep/Block + Genotype:Env, data=data)
anova <- anova(model)
out <- capture.output(anova)
cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
# Run a random effects model
model.1 <- lmer(get(trait) ~ (1|Genotype) + (1|Env) + (1|Env/Rep) + (1|Env/Rep/Block) + (1|Genotype:Env), data = data, REML = TRUE)
# Decreasing stopping tolerances
strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
 if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
summary(model, correlation=FALSE)
random_effects <- ranef(model)
# Write out BLUPs for Genotypes
write.table(random_effects$Genotype, paste0("blups_", trait, ".csv"), col.names=F, row.names=F, sep=",")
# Summary of random effects
summary <- summary(model, correlation=FALSE)
out <- capture.output(summary)
cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
# Write out residuals from ANOVA
write.table(resid(model), paste0("resids_", trait, ".csv"), col.names=F, row.names=F, sep=",")
# Calculate hertiability 
model_variances <- as.data.frame(VarCorr(model))
h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/5)+(model_variances$vcov[9]/10))
out <- capture.output(h2)
cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
pdf(paste0("assumptions_", trait, ".pdf"), width = 15, height = 5)
par(mfrow=c(1,3))
  
# Model Fit with REML
plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
abline(h=0,col="red", lwd=1, lty=1)
# histogram of residuals
hist(residuals(model),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
x=seq(-5e-15,9e-15,5e-15)
curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
# qq plot
qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
  
dev.off()
  
}

#########################################################################
  
  
# I removed these traits from the data set because they caused issues during the loop run. 
# pheno = dplyr::select(pheno, -SDL, -TWA) 


# Have a look at the data!
pheno1 = pheno %>% drop_na(PHT)

(pheno1 %>% dplyr::group_by(Environ) %>% dplyr::summarise( avg=mean(PHT)) %>% arrange(desc(avg)))


(pheno.gy.env = pheno %>% filter(Environ == c("MF2012", "MF2013", "SF2012","SF2013", "Mon2013" ))  %>% 
               ggplot(aes(TGY, fill=Environ)) + 
                geom_histogram())


(ggplot(pheno, aes(x = Environ, TGY)) + 
    stat_summary(fun.y = mean, geom = "point", shape=18, color="red") +
    geom_boxplot() +
    theme_stata() )




# Fitting Mixed Linear Models ==========================================


# fit1 <-     asreml(fixed = GY  ~ 1 + Environ + Rep:Environ,
#                    random =    ~  Block:Rep:Environ + Entries + Entries:Environ,
#                    residual =  ~ idv(units),
#                    na.action = na.method(x="include"),
#                    data=pheno)
# plot(fit1)
# fit1$loglik
# summary(fit1)$varcomp
# wald(fit1)
# 
# summary(fit1, coef=TRUE)$coef.random
# summary(fit1, coef=TRUE)$coef.fixed


## fit2 is similar to the model used by Kadam et al. 2016

fit2 <-     asreml(fixed = PHT  ~ 1 + Environ + Rep:Environ,
                   random =    ~  Block:Rep:Environ + idv(Ent) + Ent:Environ,
                   residual =  ~ dsum( ~ idv(units) | Environ),
                   na.action = na.method(x="include", y ="include"),
                   data=pheno1)


# fit2b <-     asreml(fixed = GY  ~ 1 + Environ + Rep:Environ,
#                    random =    ~  Block:Rep:Environ + idv(Entries),
#                    residual =  ~ dsum( ~ idv(units) | Environ),
#                    na.action = na.method(x="include", y ="include"),
#                    data=pheno1)


plot(fit2)
fit2$loglik
summary(fit2)$varcomp
wald(fit2)

summary(fit2, coef=TRUE)$coef.random
summary(fit2, coef=TRUE)$coef.fixed

# lrt(fit2, fit2b)

pred.fit2g = predict(fit2, "Entries")
blup.g = pred.fit2g$pvals
ggplot(blup.g, aes(predicted.value)) + geom_histogram()

pred.fit2ge = predict(fit2, "Entries:Environ")
blup.ge = pred.fit2ge$pvals
ggplot(blup.ge, aes(x=Environ, predicted.value)) + geom_boxplot()
ggplot(blup.ge, aes(fill =(Environ !="SF2012"), predicted.value)) + geom_histogram()

# Export of BLUPs
# write.csv(pred.fit2, "NIFA2015Density.csv")


# Calculation of Heritability
vg.gy = fit2$vparameters[2]
vge.gy = fit2$vparameters[3]
ve = c(fit2$vparameters[5],fit2$vparameters[7],fit2$vparameters[9],fit2$vparameters[11],fit2$vparameters[13])
ve.gy = mean(ve)  

ne = 5

rep = data.frame(table(pheno1$Ent))
nr = harmonic.mean(rep$Freq, zero = FALSE)


(h2.gy=vg.gy/(vg.gy+vge.gy/ne+ve.gy/nr))


##########################################################################################
## Now let's loop through all the traits to estimate BLUPs, Varcomps, and heritabilities.


colnames(pheno)

df = pheno

l = nlevels(df$Ent)

#Create empty data frame for BLUP output
DataOutput <- data.frame(matrix(vector(),l,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-l so that the cbind will work later on
DataOutput$Entry <- unique(df[,6]) #fill in Entry numbers
DataOutput$Row <- c(1:l)

#this empty dataframe is for variance components
DataVarComp <- data.frame()
DataVarCompOutput <- data.frame()
HeritabilityData <- data.frame()

#this empty dataframe is for dropped variance components
drops <- c("var1","var2","sdcor") 

colnames(df)
str(df[,11:ncol(df)])
#take only the columns with numerical data
colnum=c(11:ncol(df))

for (i in 1:32){  #this second loop runs through each TRAIT, one at a time
  x=colnum[i]  #set the current [i] column as x
  trait=colnames(df)[x] #sets the current column header as the trait name
 
  pheno1 = df %>% drop_na(trait)
  
  df1 <- pheno1 #make a copy that we can use for our analysis
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  #our random effects mixed model
  model = asreml(fixed = y  ~ 1 + Environ + Rep:Environ,
         random =    ~  Block:Rep:Environ + idv(Ent) + Ent:Environ,
         residual =  ~ dsum( ~ idv(units) | Environ),
         na.action = na.method(x="include", y ="include"),
         data=df1)

#BLUPs   
pred.fit2g = predict(model, "Ent")
blup = data.frame(pred.fit2g$pvals$predicted.value)
colnames(blup) = trait
DataOutput <- cbind(DataOutput,blup)


# Variance Components
vg = model$vparameters[2]
vge = model$vparameters[3]
vex = c(model$vparameters[5],model$vparameters[7],model$vparameters[9],model$vparameters[11],model$vparameters[13])
ve = mean(vex)  

varcomp = data.frame(trait, vg, vge, ve)
DataVarComp = rbind(DataVarComp, varcomp)



}

write.csv(DataOutput, "NIFA_DensityBLUPs.csv")
  
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  blup = coef(model)$Entry #coef extracts model coefficients from lmer objects returned by modeling functions
  hist(blup[,1]) #plot it out
  colnames(blup) <- trait #rename the BLUP column by the trait in our loop
  #add columns to existing dataframe   
  DataOutput <- cbind(DataOutput,blup) #ammends our dataframe with the new BLUP column for the new trait.
  
  #Modify variance component df by
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  DataVarComp <- rbind(DataVarComp,varComp) 
}



#reshape our variance components dataframe so that we can run the heritability script
DataVarCompOutput <- reshape(DataVarComp, idvar = "Trait", timevar = "grp", direction = "wide")

#the broad sense heritability script
#Taken from Gioia et al. 2017
nloc = 14
HeritabilityData <- ((DataVarCompOutput[,3])) / (((DataVarCompOutput[,3])) + (((DataVarCompOutput[,5])) / (nloc)))

#create function to identify maximum value in a column
colMax <- function(data) sapply(data, max, na.rm = TRUE)
#create function to identify minimum value in a column
colMin <- function(data) sapply(data, min, na.rm = TRUE)
#create function to identify mean value in a column
colMean <- function(data) sapply(data, mean, na.rm = TRUE)
#create function to identify median value in a column
colMedian <- function(data) sapply(data, median, na.rm = TRUE)
#create function to identify standard dev value in a column
colStdev <- function(data) sapply(data, sd, na.rm = TRUE)

#summary statistics
DataColMax <- colMax(df[,colnum])
DataColMin <- colMin(df[,colnum])
DataColMean <- colMean(df[,colnum])
DataColMedian <- colMean(df[,colnum])
DataColStdev <- colStdev(df[,colnum])

#CVg <- ((DataVarCompOutput[,3])) / ((DataColMean)^2)
CVg <- (sqrt(DataVarCompOutput[,3])) / ((DataColMean))
#out <- LSD.test(model,"Entry", p.adj="bonferroni")

#bind the heritability to the variance component data frame

DataVarCompOutput <- cbind(DataVarCompOutput,HeritabilityData,CVg,DataColMin,DataColMax,DataColMean,DataColMedian,DataColStdev)
DataVarCompOutput[1:10,]

DataVarCompOutput %>% filter(Trait == "TRL") %>% select(HeritabilityData) 

#output that beast
write.csv(Day6DataOutput,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/BLUPsDay6_thinned_Nov14.csv")
write.csv(DataVarCompOutput,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/VarianceComponents_thinned_Nov14.csv")
? 2020 GitHub, Inc.
