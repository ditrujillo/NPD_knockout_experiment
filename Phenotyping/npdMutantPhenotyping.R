detach(NodArea)
detach(NodAreaAvg)
rm(list=ls())

###Load Packages and Data ------------------------------------------------------------------
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) }
if (!require("nlme")) { install.packages("nlme", dependencies = TRUE) }
if (!require("lme4")) { install.packages("lme4", dependencies = TRUE) }
if (!require("lmerTest")) { install.packages("lmerTest", dependencies = TRUE) }
if (!require("foreign")) { install.packages("foreign", dependencies = TRUE) }
if (!require("agricolae")) { install.packages("agricolae", dependencies = TRUE) }
if (!require("multcompView")) { install.packages("multcompView", dependencies = TRUE) }

library(ggplot2)
library(nlme)
library(lme4)
library(lmerTest)
library(foreign)
library(agricolae)
library(multcompView)

###Analyze Raw data -----------------------------------------------------------------------
NodArea <- read.delim("~/Data/GeneFams/ESPs/PhenotypingFinal/npdMutantPhenotyping.txt")
attach(NodArea)
View(NodArea)
Plant <- as.factor(Plant)

#Preliminary
str(NodArea)
formula(NodArea)
interaction.plot(Deletions, Plant, Nod)     		#Is the difference in Area btwn Genotypes significant?

#Test for equal variances then do t.test
var.test(PinkNod[Deletions=='3gene1'],PinkNod[Deletions=='5gene'], paired = FALSE)
t.test(PinkNod[Genotype=='Gm8.5'],PinkNod[Genotype=='Gm26.4'], paired = FALSE, var.equal = FALSE) #p-value = 9.416e-11

#Boxplots
boxplot(Nod~Deletions, lwd = 2, ylab = 'Area (mm2)', title='Total Nodule Area')
stripchart(Nod~Deletions, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = Plant)

boxplot(PinkNod~Deletions, lwd = 2, ylab = 'Pink Area (mm2)', title='Nodule Pink Area')
stripchart(PinkNod~Deletions, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = Plant)

ggplot(data=NodArea, aes(x=Deletions, y=Nod), na.rm = TRUE) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_dodge(0.5), aes(colour = as.factor(Plant)))+
  coord_cartesian(ylim = c(0, 4)) +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=28,face="bold"))

ggplot(data=NodArea, aes(x=Deletions, y=PinkNod), na.rm = TRUE) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_dodge(0.5), aes(colour = as.factor(Plant)))+
  coord_cartesian(ylim = c(0, 3)) +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=28,face="bold"))


##Test differences in the means: Mixed effects model using 'lme' ---------------------------
model1 <- aov(Nod ~ Deletions+Plant, data=NodArea) 
anova(model1)

#Nested with lme
model2 <- lme(Nod ~ Deletions, random =~1|Plant, method="REML", data=NodArea, na.action = na.exclude) 
anova(model2)
anova.lme(model2)
plot(model2)
hist(residuals(model2), col="darkgray", breaks=16)
plot(fitted(model2), residuals(model2)) 

Pinkmodel2 <- lme(PinkNod ~ Deletions, random =~1|Plant, method="REML", data=NodArea, na.action = na.exclude) 
anova(Pinkmodel2)
anova.lme(Pinkmodel2)
plot(Pinkmodel2)
hist(residuals(Pinkmodel2), col="darkgray", breaks=16)
plot(fitted(Pinkmodel2), residuals(Pinkmodel2)) 

#Nested with lmer
model3 = lmer(Nod ~ Deletions + (1|Plant), data=NodArea, REML=TRUE)
anova(model3)
rand(model3)
difflsmeans(model3, test.effs="Deletions")

Pinkmodel3 = lmer(PinkNod ~ Deletions + (1|Plant), data=NodArea, REML=TRUE)
anova(Pinkmodel3)
rand(Pinkmodel3)
difflsmeans(Pinkmodel3, test.effs="Deletions")


#Aggregate the raw nodule observations on a per plant basis --------------------------------
#Used for the final analysis in the paper:

aggregate(LeafLength~Plant+Genotype,data=NodArea,length) #LeafNum = Number of leaves
aggregate(LeafLength~Plant+Genotype,data=NodArea,sum) #Leaftot = Cummulative petiole lengths
aggregate(LeafLength~Plant+Genotype,data=NodArea,mean) #LeafLengthAvg = Average petiole length
aggregate(LeafLength~Plant+Genotype,data=NodArea,max) #Height = Longest petiole length ~ Plant height
aggregate(Nod~Plant+Genotype,data=NodArea,length) #NodNum = Number of nodules
aggregate(Nod~Plant+Genotype,data=NodArea,sum) #NodTot = Cumulative nodule area
aggregate(Nod~Plant+Genotype,data=NodArea,mean) #NodAvg = Average nodule area/size
aggregate(Nod~Plant+Genotype,data=NodArea,max) #NodMax = Largest nodule area/size
aggregate(PinkNod[PinkNod>0]~Plant[PinkNod>0]+Genotype[PinkNod>0],data=NodArea,length) #PinkNum = Number of pink nodules
aggregate(PinkNod~Plant+Genotype,data=NodArea,sum) #PinkTot = Cumulative pink nodule area
aggregate(PinkNod~Plant+Genotype,data=NodArea,mean) #PinkAvg = Average pink nodule area
aggregate(PinkNod~Plant+Genotype,data=NodArea,max) #PinkMax = Largest pink nodule area
aggregate(GreenNod[GreenNod>0]~Plant[GreenNod>0]+Genotype[GreenNod>0],data=NodArea,length) #GreenNum = Number of green nodules
aggregate(GreenNod~Plant+Genotype,data=NodArea,sum) #GreenTot = Cumulative green nodule area
aggregate(GreenNod~Plant+Genotype,data=NodArea,mean) #GreenAvg = Average green nodule area
aggregate(GreenNod~Plant+Genotype,data=NodArea,max) #GreenMax = Largest green nodule area

#All of the above calculations were compiled into a text file called:
#npdMutantPhenotypingAvg.txt
detach(NodArea)
rm(list=ls())

###Analyze compiled data on a per plant basis ----------------------------------------------
NodAreaAvg <- read.delim("~/Data/GeneFams/ESPs/PhenotypingFinal/npdMutantPhenotypingAvg.txt")
attach(NodAreaAvg)
View(NodAreaAvg)
Plant <- as.factor(Plant)
#Check for heteroscedasticity
aggregate(PinkTot~Deletions,data=NodAreaAvg,mean)
aggregate(PinkTot~Deletions,data=NodAreaAvg,var)

#Summarize data for the paper
#Leaftot = Cummulative petiole lengths
aggregate(.~Deletions,data=NodAreaAvg,mean)
aggregate(.~Deletions, NodAreaAvg, function(x) c(mean = mean(x), sd = sd(x)))

#Find p-values between plant lines
#Plant traits
TukeyHSD(aov(LeafNum~Deletions, data=NodAreaAvg)) #LeafNum = Number of leaves
TukeyHSD(aov(Leaftot~Deletions, data=NodAreaAvg)) #Leaftot = Cummulative petiole lengths
TukeyHSD(aov(LeafLengthAvg~Deletions, data=NodAreaAvg)) #LeafLengthAvg = Average petiole length
TukeyHSD(aov(Height~Deletions, data=NodAreaAvg)) #Height = Longest petiole length ~ Plant height
#Nodule traits
TukeyHSD(aov(NodNum~Deletions, data=NodAreaAvg)) #NodNum = Number of nodules
TukeyHSD(aov(NodTot~Deletions, data=NodAreaAvg)) #NodTot = Cumulative nodule area
TukeyHSD(aov(NodAvg~Deletions, data=NodAreaAvg)) #NodAvg = Average nodule area/size
TukeyHSD(aov(NodMax~Deletions, data=NodAreaAvg)) #NodMax = Largest nodule area/size
#Nodule pink measurements
TukeyHSD(aov(PinkNum~Deletions, data=NodAreaAvg))  #PinkNum = Number of pink nodules
TukeyHSD(aov(PinkTot~Deletions, data=NodAreaAvg)) #PinkTot = Cumulative pink nodule area
TukeyHSD(aov(PinkAvg~Deletions, data=NodAreaAvg)) #PinkAvg = Average pink nodule area
TukeyHSD(aov(PinkMax~Deletions, data=NodAreaAvg)) #PinkMax = Largest pink nodule area
#Nodule green measurements
TukeyHSD(aov(GreenNum~Deletions, data=NodAreaAvg))  #GreenNum = Number of green nodules
TukeyHSD(aov(GreenTot~Deletions, data=NodAreaAvg)) #GreenTot = Cumulative green nodule area
TukeyHSD(aov(GreenAvg~Deletions, data=NodAreaAvg))  #GreenAvg = Average green nodule area
TukeyHSD(aov(GreenMax~Deletions, data=NodAreaAvg))  #GreenMax = Largest green nodule area


#Draw plots

#Number of Leaves and Plant Height
boxplot(LeafNum~Deletions, lwd = 2, ylab = 'Number of Leaves')
stripchart(LeafNum~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(LeafNum~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(LeafNum~Deletions, data=NodAreaAvg))$Deletions))

boxplot(Height~Deletions, lwd = 2, ylab = 'Plant Height (mm)', ylim=c(10,70)) #Figure 4f
stripchart(Height~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(Height~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(Height~Deletions, data=NodAreaAvg))$Deletions))

#Cummulative leaf length and average leaf length
boxplot(Leaftot~Deletions, lwd = 2, ylab = 'Cummulative Leaf Height (mm)', ylim=c(10,470))
stripchart(Leaftot~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(Leaftot~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(Leaftot~Deletions, data=NodAreaAvg))$Deletions))

boxplot(LeafLengthAvg~Deletions, lwd = 2, ylab = 'Average Leaf Height (mm)')
stripchart(LeafLengthAvg~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(LeafLengthAvg~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(LeafLengthAvg~Deletions, data=NodAreaAvg))$Deletions))

#Number of Total, Pink or Senescing Nodules
boxplot(NodNum~Deletions, lwd = 2, ylab = 'Number of Nodules')
stripchart(NodNum~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(NodNum~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(NodNum~Deletions, data=NodAreaAvg))$Deletions))

boxplot(PinkNum~Deletions, lwd = 2, ylab = 'Number of Pink Nodules', ylim=c(0,30)) #Figure 4e
stripchart(PinkNum~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(PinkNum~Deletions, data=NodAreaAvg))
t.test(PinkNum[Deletions=='3gene1'],PinkNum[Deletions=='5gene'], paired = FALSE)
multcompLetters(extract_p(TukeyHSD(aov(PinkNum~Deletions, data=NodAreaAvg))$Deletions))

boxplot(GreenNum~Deletions, lwd = 2, ylab = 'Number of Senescing Nodules', ylim=c(0,11))
stripchart(GreenNum~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(GreenNum~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(GreenNum~Deletions, data=NodAreaAvg))$Deletions))

#Average Nodule Area, Pink Area, or Green Area
boxplot(NodAvg~Deletions, lwd = 2, ylab = 'Avg Nod Size (mm2)', ylim=c(0,2.5))
stripchart(NodAvg~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(NodAvg~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(NodAvg~Deletions, data=NodAreaAvg))$Deletions))

boxplot(PinkAvg~Deletions, lwd = 2, ylab = 'Avg Pink Area (mm2)', ylim=c(0,1.6))
stripchart(PinkAvg~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(PinkAvg~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(PinkAvg~Deletions, data=NodAreaAvg))$Deletions))

boxplot(GreenAvg~Deletions, lwd = 2, ylab = 'Avg Senescing Area (mm2)')
stripchart(GreenAvg~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(GreenAvg~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(GreenAvg~Deletions, data=NodAreaAvg))$Deletions))

#ggPlot
ggplot(data=NodAreaAvg, aes(x=Deletions, y=PinkAvg), na.rm = TRUE) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_dodge(0.5), aes(colour = as.factor(Plant)))+
  coord_cartesian(ylim = c(0,1.75)) + ylab("Average N-fixing area (mm2)")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"))

#Total Nodule Area, Pink Area, or Green Area
boxplot(NodTot~Deletions, lwd = 2, ylab = 'Total Nodule Area (mm2)', ylim=c(0,50))
stripchart(NodTot~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(NodTot~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(NodTot~Deletions, data=NodAreaAvg))$Deletions))

boxplot(PinkTot~Deletions, lwd = 2, ylab = 'Total Pink Area (mm2)', ylim=c(0,30))
stripchart(PinkTot~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(PinkTot~Deletions, data=NodAreaAvg))
plot(aov(PinkTot~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(PinkTot~Deletions, data=NodAreaAvg))$Deletions))
#TukeyHSD(aov(log(PinkTot+0.1)~Deletions, data=NodAreaAvg)) #Transformed data
#plot(aov(sqrt(PinkTot)~Deletions, data=NodAreaAvg))

boxplot(GreenTot~Deletions, lwd = 2, ylab = 'Total Senescing Area (mm2)', ylim=c(0,6))
stripchart(GreenTot~Deletions, vertical = TRUE, method = "jitter", jitter = 0.25, add = TRUE, pch = 20, col = Plant)
TukeyHSD(aov(GreenTot~Deletions, data=NodAreaAvg))
multcompLetters(extract_p(TukeyHSD(aov(GreenTot~Deletions, data=NodAreaAvg))$Deletions))


#Find correlations between data
plot(PinkTot~Height, ylab='Cumulative N-fixation zone size (mm2)', xlab='Plant height (mm)', col = Deletions, pch = 16, cex = 1.1)
legend("topleft",c(levels(Deletions), 'R2=0.63, p<2e-16'),col=c(1,2,3,4,5,6,7,8),pch=c(16,16,16,16,16,16,16,NA),ncol=2,cex=0.9,pt.cex=1.1)
abline(lm(PinkTot~Height), col='red')
summary(lm(PinkTot~Height))

plot(PinkTot~Leaftot, ylab='Cumulative N-fixation zone size (mm2)', xlab='Cummulative petiole lengths (mm)', col = Deletions, pch = 16, cex = 1.1)
legend("topleft",c(levels(Deletions), 'R2=0.82, p<2e-16'),col=c(1,2,3,4,5,6,7,8),pch=c(16,16,16,16,16,16,16,NA),ncol=2,cex=0.9,pt.cex=1.1)
abline(lm(PinkTot~Leaftot), col='red')
summary(lm(PinkTot~Leaftot))

plot(PinkTot~LeafNum, ylab='Cumulative N-fixation zone size (mm2)', xlab='Number of Leaves', col = Deletions, pch = 16, cex = 1.1)
legend("topleft",c(levels(Deletions), 'R2=0.71, p<2e-16'),col=c(1,2,3,4,5,6,7,8),pch=c(16,16,16,16,16,16,16,NA),ncol=2,cex=0.9,pt.cex=1.1)
abline(lm(PinkTot~LeafNum), col='red')
summary(lm(PinkTot~LeafNum))





