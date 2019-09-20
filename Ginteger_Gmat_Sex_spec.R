rm(list=ls())

library(lattice);library(tidyverse); library(lubridate);library(MCMCglmm);library(pedantics)
library(GGally); library(ggthemes); library(patchwork)


#### 0. Data import and cleanup ####

#Data<- read.csv("G.integer.Aim1.Full.csv",header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
Data<- read.csv("Royaute(PopComp)_Data.csv",header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
str(Data)


#### 0.1 Subset data by generation and F1 generation by Sex ####
Data_F0=subset(Data, Generation=="F0")
Data_F1=subset(Data, Generation=="F1")

Data_F1_F=subset(Data_F1, Sex=="F")
Data_F1_M=subset(Data_F1, Sex=="M")

Data_F=rbind(Data_F0,Data_F1_F)
Data_M=rbind(Data_F0,Data_F1_M)


#### 0.2 Pedigree generation ####
PooledPedigree2<-Data[,c(1,3,2)]
PooledPedigree.analysis<-fixPedigree(PooledPedigree2)


#### 1. Define MCMC chain iterations ####
Mult=1;NITT=120000;THIN=100;BURN=20000
#Mult=.5
Mult=40

#### 2. Gmat by sex ####
#### 2.1 PRIOR BASED ON LL estimates ####
Va=c(4*2.858,4*3.912,4*0.1912,4*0.04876,4*5.483,4*0.4556,4*0.03946)
Vr=c(31.434,94.845,1,1.47720,83.950,1,0.69877)

prior.7v<-
  list(R=list(V=diag(7)*Vr,nu=7),
       G=list(G1=list(V=diag(7)*Va,nu=7)))

#### 2.2 Run model with cross-sex G matrix ####
# remove non-adults
Data_SexSpec=subset(Data, Sex == "M" | Sex == "F")
str(Data_SexSpec)

Data_SexSpec$Latency.cen_M=with(Data_SexSpec, 
                                 ifelse(Sex=="M",Latency.cen,NA))
Data_SexSpec$Latency.cen_F=with(Data_SexSpec, 
                                 ifelse(Sex=="F",Latency.cen,NA))
Data_SexSpec$Latency.up_M=with(Data_SexSpec, 
                                ifelse(Sex=="M",Latency.up,NA))
Data_SexSpec$Latency.up_F=with(Data_SexSpec, 
                                ifelse(Sex=="F",Latency.up,NA))
Data_SexSpec$OF.Dist_M=with(Data_SexSpec, 
                                 ifelse(Sex=="M",OF.Dist,NA))
Data_SexSpec$OF.Dist_F=with(Data_SexSpec, 
                                 ifelse(Sex=="F",OF.Dist,NA))
Data_SexSpec$AP.Dist_M=with(Data_SexSpec, 
                                 ifelse(Sex=="M",AP.Dist,NA))
Data_SexSpec$AP.Dist_F=with(Data_SexSpec, 
                                 ifelse(Sex=="F",AP.Dist,NA))

Va_fm=c(4*2.858,4*3.912,4*5.483,
        4*2.858,4*3.912,4*5.483)
      
Vr_fm=c(31.434,94.845,1.47720,
        31.434,94.845,1.47720)

prior.6v<-
  list(R=list(V=diag(6)*Vr_fm,nu=6),
       G=list(G1=list(V=diag(6)*Va_fm,nu=6)))

GBmat.MCMC.5<-
  MCMCglmm(cbind(Latency.cen_F,Latency.up_F,sqrt(OF.Dist_F),sqrt(AP.Dist_F),
                 Latency.cen_M,Latency.up_M,sqrt(OF.Dist_M),sqrt(AP.Dist_M))~(trait-1)+
             trait:(Mass+Generation+Pop+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,3):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP),
           random = ~us(trait):animal,rcov=~idh(trait):units,
           pedigree = PooledPedigree.analysis,data=Data_SexSpec, 
           prior=prior.6v, verbose=F,
           family=c("cengaussian","gaussian","gaussian",
                    "cengaussian","gaussian","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr=T)#,
summary(GBmat.MCMC.5)
GBmat<-matrix(colMeans(GBmat.MCMC.5$VCV[,1:36]),nrow=6,
               dimnames = list(c("Latency_F","OF.Dist_F","AP.Dist_F","Latency_M","OF.Dist_M","AP.Dist_M"),
                               c("Latency_F","OF.Dist_F","AP.Dist_F","Latency_M","OF.Dist_M","AP.Dist_M")))
GBmat
HPDinterval(GBmat.MCMC.5$VCV[,1:36])
HPDinterval(GBmat.MCMC.5$VCV[,1:36],prob=.9)

GBmat.2<-matrix(colMeans(posterior.cor(GBmat.MCMC.5$VCV[,1:36])),nrow=6,
              dimnames = list(c("Latency_F","OF.Dist_F","AP.Dist_F","Latency_M","OF.Dist_M","AP.Dist_M"),
                              c("Latency_F","OF.Dist_F","AP.Dist_F","Latency_M","OF.Dist_M","AP.Dist_M")))
GBmat.2
HPDinterval(GBmat.MCMC.5$VCV[,1:36])
HPDinterval(GBmat.MCMC.5$VCV[,1:36],prob=.9)



#### 2.3 run Models ####
#Gmat.MCMC.F<-
#  MCMCglmm(cbind(Latency.cen,Latency.up,
#                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
#                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
#             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
#             at.level(trait,5):(Understimate.AP)+
#             at.level(trait,6):(Understimate.AP)+
#             at.level(trait,7):(Understimate.AP),
#           random = ~us(trait):animal,rcov=~us(trait):units,
#           pedigree = PooledPedigree.analysis,data=Data_F, 
#           prior=prior.7v, verbose=F,
#           family=c("cengaussian","gaussian","poisson","gaussian",
#                    "gaussian","poisson","gaussian"),
#           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr=T)#,
#summary(Gmat.MCMC.F)
#Gmat.F<-matrix(colMeans(Gmat.MCMC.F$VCV[,1:49]),nrow=7,
#               dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
#                               c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))
#Gmat.F
#HPDinterval(Gmat.MCMC.F$VCV[,1:49])
#HPDinterval(Gmat.MCMC.F$VCV[,1:49],prob=.9)
#plot(Gmat.MCMC.F$VCV)


#Gmat.MCMC.M<-
#  MCMCglmm(cbind(Latency.cen,Latency.up,
#                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
#                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
#             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
#             at.level(trait,5):(Understimate.AP)+
#             at.level(trait,6):(Understimate.AP)+
#             at.level(trait,7):(Understimate.AP),
#           random = ~us(trait):animal,rcov=~us(trait):units,
#           pedigree = PooledPedigree.analysis,data=Data_M, 
#           prior=prior.7v, verbose=F,
#           family=c("cengaussian","gaussian","poisson","gaussian",
#                    "gaussian","poisson","gaussian"),
#           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr=T)#,
#summary(Gmat.MCMC.M)
#Gmat.M<-matrix(colMeans(Gmat.MCMC.M$VCV[,1:49]),nrow=7,
#               dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
#                               c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))
#Gmat.M
#HPDinterval(Gmat.MCMC.M$VCV[,1:49])
#HPDinterval(Gmat.MCMC.M$VCV[,1:49],prob=.9)
#plot(Gmat.MCMC.M$VCV)

#### 3. Extract Analyses files ####
save.image(file = "Gmat_Sex_Spec.RData")

# Extract each model
#Gmat
#saveRDS(Gmat.MCMC.F, file = "Gmat.MCMC.F.rds")
#saveRDS(Gmat.MCMC.M, file = "Gmat.MCMC.M.rds")
saveRDS(GBmat.MCMC.5, file = "GBmat.MCMC.5.rds")


#### 4. Results: Plot  Data ####
library(tidyverse);library(ggthemes);library(patchwork);library(cowplot);library(ggbeeswarm);library(GGally)

Data_Plot=subset(Data, Sex == "M" | Sex == "F")
str(Data_Plot)

##### 4.1 Heritabilities by sex ####
GBmat

# Females
# Lat
posterior.mode(GBmat.MCMC.5$VCV[,1]/(GBmat.MCMC.5$VCV[,1]+GBmat.MCMC.5$VCV[,37]))
HPDinterval(GBmat.MCMC.5$VCV[,1]/(GBmat.MCMC.5$VCV[,1]+GBmat.MCMC.5$VCV[,37]))
# OF.Dist
posterior.mode(GBmat.MCMC.5$VCV[,8]/(GBmat.MCMC.5$VCV[,8]+GBmat.MCMC.5$VCV[,37]))
HPDinterval(GBmat.MCMC.5$VCV[,8]/(GBmat.MCMC.5$VCV[,8]+GBmat.MCMC.5$VCV[,37]))
# AP.Dist
posterior.mode(GBmat.MCMC.5$VCV[,15]/(GBmat.MCMC.5$VCV[,15]+GBmat.MCMC.5$VCV[,38]))
HPDinterval(GBmat.MCMC.5$VCV[,15]/(GBmat.MCMC.5$VCV[,15]+GBmat.MCMC.5$VCV[,38]))

# Males
# Lat
posterior.mode(GBmat.MCMC.5$VCV[,22]/(GBmat.MCMC.5$VCV[,22]+GBmat.MCMC.5$VCV[,39]))
HPDinterval(GBmat.MCMC.5$VCV[,22]/(GBmat.MCMC.5$VCV[,22]+GBmat.MCMC.5$VCV[,39]))
# OF.Dist
posterior.mode(GBmat.MCMC.5$VCV[,29]/(GBmat.MCMC.5$VCV[,29]+GBmat.MCMC.5$VCV[,40]))
HPDinterval(GBmat.MCMC.5$VCV[,29]/(GBmat.MCMC.5$VCV[,29]+GBmat.MCMC.5$VCV[,40]))
# AP.Dist
posterior.mode(GBmat.MCMC.5$VCV[,36]/(GBmat.MCMC.5$VCV[,36]+GBmat.MCMC.5$VCV[,41]))
HPDinterval(GBmat.MCMC.5$VCV[,15]/(GBmat.MCMC.5$VCV[,15]+GBmat.MCMC.5$VCV[,38]))


##### 4.2 Boxplot Sex-difference in behavior ####
bw.plot.Lat=ggplot(Data_Plot, aes(x= Sex, y=sqrt(Latency.Cens), fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Latency (s)") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
bw.plot.Lat

bw.plot.OF=ggplot(Data_Plot, aes(x= Sex, y=sqrt(OF.Dist), fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Open-field distance (cm)") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
bw.plot.OF

bw.plot.AP=ggplot(Data_Plot, aes(x= Sex, y=sqrt(AP.Dist), fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Antipredator response (cm)") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
bw.plot.AP

bw.plots= bw.plot.Lat + bw.plot.OF + bw.plot.AP + plot_layout(nrow = 1)
bw.plots
ggsave(bw.plots,filename = "bw.plots.jpeg")

##### 4.3 Boxplot breeding values differences by sex #####
# Extract BLUPS
library(tidybayes);library(emmeans);library(plyr)

Data_BLUPS=data.frame(posterior.mode(GBmat.MCMC.5$Sol))
str(Data_BLUPS)
Data_BLUPS$animal=factor(rownames(Data_BLUPS))
names(Data_BLUPS)=c("BLUPS","animal")
Data_BLUPS=Data_BLUPS[-c(1:80),]  # remove rows 1-80 containing all fixed effects coefficients
# Find row limits for each Trait
#Data_BLUPS[1038,] # Latency_F row 1:1037
#Data_BLUPS[2074,] # OF.Dist_F row 1037:2074
#Data_BLUPS[3111,] # AP.Dist_F row 2075:3111
#Data_BLUPS[4148,] # Latency_M row 3112:4148
#Data_BLUPS[5185,] # OF.Dist_M row 4149:5185
#Data_BLUPS[6222,] # AP.Dr_M ow 5186:6222

# Code each for its corresponding behavioral trait type
Data_BLUPS$Trait=factor(c(rep("Latency_F",1037),
                          rep("OF.Dist_F",1037),
                          rep("AP.Dist_F",1037),
                          rep("Latency_M",1037),
                          rep("OF.Dist_M",1037),
                          rep("AP.Dist_M",1037)),
                        levels = c("Latency_F","OF.Dist_F","AP.Dist_F",
                                   "Latency_M","OF.Dist_M","AP.Dist_M"))

#Data_BLUPS$Trait2=factor(c(rep("Latency",1037),
#                          rep("OF.Dist",1037),
#                          rep("AP.Dist",1037),
#                          rep("Latency",1037),
#                          rep("OF.Dist",1037),
#                          rep("AP.Dist",1037)),
#                        levels = c("Latency","OF.Dist","AP.Dist"))

# Recover only the individul unique ID code from the MCMCglmm BLUPS object
animal_split=str_split_fixed(Data_BLUPS$animal,"animal.", n=2) %>% data.frame()
names(animal_split)=c("Trait","animal")

animal_split$Trait=str_sub(animal_split$Trait,6, str_length(animal_split$Trait)-1) %>% as.factor()
animal_split$Trait=revalue(animal_split$Trait, c("Latency.cen_F"="Latency_F",
                                                 "Latency.cen_M"="Latency_M"))
animal_split$Trait=factor(animal_split$Trait, 
                          levels=c("Latency_F","OF.Dist_F","AP.Dist_F",
                                   "Latency_M","OF.Dist_M","AP.Dist_M"))
str(animal_split)
str(Data_BLUPS)


Data_BLUPS$animal=animal_split$animal
#Data_BLUPS$Trait=animal_split$Trait

# Dataframe containing the BLUPS of all individuals in the dataset for each behavioral trait
Data_BLUPS.2=spread(Data_BLUPS, Trait, BLUPS)
str(Data_BLUPS.2)
Data_BLUPS.2=merge(Data_BLUPS.2,Data[,c(1,6)])
Data_BLUPS.2=subset(Data_BLUPS.2, Sex == "M" | Sex == "F")

# Assign female breeding value if female and male breeding value if male
Data_BLUPS.2$Latency=with(Data_BLUPS.2, ifelse(Sex=="M", Latency_M, Latency_F))
Data_BLUPS.2$OF.Dist=with(Data_BLUPS.2, ifelse(Sex=="M", OF.Dist_M, OF.Dist_F))
Data_BLUPS.2$AP.Dist=with(Data_BLUPS.2, ifelse(Sex=="M", AP.Dist_M, AP.Dist_F))

# Add model intercept
Data_BLUPS.2$Latency=with(Data_BLUPS.2, ifelse(Sex=="M", Latency_M+mean(GBmat.MCMC.5$Sol[,4]), 
                                               Latency_F+mean(GBmat.MCMC.5$Sol[,1])))
Data_BLUPS.2$OF.Dist=with(Data_BLUPS.2, ifelse(Sex=="M", OF.Dist_M+mean(GBmat.MCMC.5$Sol[,5]),
                                               OF.Dist_F+mean(GBmat.MCMC.5$Sol[,2])))
Data_BLUPS.2$AP.Dist=with(Data_BLUPS.2, ifelse(Sex=="M", AP.Dist_M+mean(GBmat.MCMC.5$Sol[,6]),
                                               AP.Dist_F+mean(GBmat.MCMC.5$Sol[,3])))


blup.plot.Lat=ggplot(Data_BLUPS.2, aes(x= Sex, y=Latency, fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Latency") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
blup.plot.Lat

blup.plot.OF=ggplot(Data_BLUPS.2, aes(x= Sex, y=OF.Dist, fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Open-field distance") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
blup.plot.OF

blup.plot.AP=ggplot(Data_BLUPS.2, aes(x= Sex, y=AP.Dist, fill=Sex, color=Sex)) +
  geom_beeswarm(size=3, alpha = .5) + geom_violin(alpha = .25) +
  ylab("Antipredator response") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
blup.plot.AP

blup.plots= blup.plot.Lat + blup.plot.OF + blup.plot.AP + plot_layout(nrow = 1)
blup.plots

# facet plot version
Data_BLUPS.3=Data_BLUPS.2[,c(1,8:11)]
Data_BLUPS.3=Data_BLUPS.3 %>% gather(Trait, Value, Latency:AP.Dist, factor_key=TRUE)
names(Data_BLUPS.3)
# relabel Trait factor
library(plyr)
Data_BLUPS.3$Trait=revalue(Data_BLUPS.3$Trait,
                           c("OF.Dist"="Open-field activity",
                             "AP.Dist"="Antipredator response"))
Data_BLUPS.3$Sex=revalue(Data_BLUPS.3$Sex,
                           c("F"="Females",
                             "M"="Males"))

blup.plots.2=ggplot(Data_BLUPS.3, aes(x= Sex, y=Value, fill=Sex, color=Sex)) +
  geom_beeswarm(size=1, alpha = .25) + 
  #geom_violin(alpha = .25) +
  facet_wrap(~Trait) +
  xlab("Sex") + ylab("Breeding values") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() 

blup.plots.2=blup.plots.2+ theme(aspect.ratio=1,
                                 legend.position="none",
                                 axis.title = element_text(size=12),
                                 axis.text = element_text(size=12),
                                 strip.text.x = element_text(size=14))
  

# Plot all correlations
Gmf=ggpairs(Data_BLUPS.2[,c(2:7)],
            lower=list(continuous="smooth"),
            diag = list("continuous")) + theme_base()
Gmf


##### 4.3  Va compared by Sex ####
Va_bSex=data.frame(Trait=factor(rep(c("Latency","Open-field activity",
                                  "Antipredator response"),2),
                                levels = c("Latency","Open-field activity",
                                           "Antipredator response")),
                   Sex=factor(c("Females","Females","Females",
                                    "Males","Males","Males")),
                   mean=c(mean(GBmat.MCMC.5$VCV[,1]),
                          mean(GBmat.MCMC.5$VCV[,8]),
                          mean(GBmat.MCMC.5$VCV[,15]),
                          mean(GBmat.MCMC.5$VCV[,22]),
                          mean(GBmat.MCMC.5$VCV[,29]),
                          mean(GBmat.MCMC.5$VCV[,36])),
                   low.ci=c(HPDinterval(GBmat.MCMC.5$VCV[,1])[1],
                            HPDinterval(GBmat.MCMC.5$VCV[,8])[1],
                            HPDinterval(GBmat.MCMC.5$VCV[,15])[1],
                            HPDinterval(GBmat.MCMC.5$VCV[,22])[1],
                            HPDinterval(GBmat.MCMC.5$VCV[,29])[1],
                            HPDinterval(GBmat.MCMC.5$VCV[,36])[1]),
                   up.ci=c(HPDinterval(GBmat.MCMC.5$VCV[,1])[2],
                            HPDinterval(GBmat.MCMC.5$VCV[,8])[2],
                            HPDinterval(GBmat.MCMC.5$VCV[,15])[2],
                            HPDinterval(GBmat.MCMC.5$VCV[,22])[2],
                            HPDinterval(GBmat.MCMC.5$VCV[,29])[2],
                            HPDinterval(GBmat.MCMC.5$VCV[,36])[2]))
 
pd <- position_dodge(0.2)
Va_bSex.plot=ggplot(Va_bSex, aes(x= Sex, y=mean, fill=Sex, color=Sex)) +
  geom_point(size=5) + 
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=.75, position=pd) +
  facet_wrap(~Trait) +
  ylab("Genetic variance (Va)") + xlab("Sex") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() 
Va_bSex.plot= Va_bSex.plot + theme(aspect.ratio=1,legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size=14))

deltaVa=data.frame(Trait=factor(c(rep("Latency",1000),
                                  rep("Open-field activity",1000),
                                  rep("Antipredator response",1000)),
                                levels = c("Latency","Open-field activity",
                                           "Antipredator response")),
                   delta.dist=c(GBmat.MCMC.5$VCV[,1]-GBmat.MCMC.5$VCV[,22],
                                GBmat.MCMC.5$VCV[,8]-GBmat.MCMC.5$VCV[,29],
                                GBmat.MCMC.5$VCV[,15]-GBmat.MCMC.5$VCV[,36]))
str(deltaVa)
 

deltaVa.plot=ggplot(deltaVa, aes(x=delta.dist)) +
  facet_wrap(~Trait) +
  geom_histogram(aes(y=..density..), color="black",fill="white",) + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  xlab("deltaVa") + ylab("") +
  theme_bw() 

deltaVa.plot=deltaVa.plot + theme(aspect.ratio=1,legend.position="none",
                       axis.title = element_text(size=12),
                       axis.text = element_text(size=12),
                     strip.text.x = element_text(size=14))

##### 4.3.1  Figure 1: Breeding values and Va by sex  ####
mean.Va.plots=blup.plots.2/Va_bSex.plot/deltaVa.plot
mean.Va.plots 
ggsave(mean.Va.plots,filename = "mean.Va.plots.jpeg", scale = 1.5)

#### 4.3.2 Pmcmc for Figure 1 ####
# Phenotypic means
summary(GBmat.MCMC.5)

# Latency (Females - Males)
hist(GBmat.MCMC.5$Sol[,"traitLatency.cen_F"]-GBmat.MCMC.5$Sol[,"traitLatency.cen_M"])
posterior.mode(GBmat.MCMC.5$Sol[,"traitLatency.cen_F"]-GBmat.MCMC.5$Sol[,"traitLatency.cen_M"])
HPDinterval(GBmat.MCMC.5$Sol[,"traitLatency.cen_F"]-GBmat.MCMC.5$Sol[,"traitLatency.cen_M"])
table(GBmat.MCMC.5$Sol[,"traitLatency.cen_F"]-GBmat.MCMC.5$Sol[,"traitLatency.cen_M"]>0) # 0.95

# OF.Dist (Females - Males)
hist(GBmat.MCMC.5$Sol[,"traitOF.Dist_F"]-GBmat.MCMC.5$Sol[,"traitOF.Dist_M"])
posterior.mode(GBmat.MCMC.5$Sol[,"traitOF.Dist_F"]-GBmat.MCMC.5$Sol[,"traitOF.Dist_M"])
HPDinterval(GBmat.MCMC.5$Sol[,"traitOF.Dist_F"]-GBmat.MCMC.5$Sol[,"traitOF.Dist_M"])
table(GBmat.MCMC.5$Sol[,"traitOF.Dist_F"]-GBmat.MCMC.5$Sol[,"traitOF.Dist_M"]>0) # 0.95

# OF.Dist (Females - Males)
hist(GBmat.MCMC.5$Sol[,"traitAP.Dist_F"]-GBmat.MCMC.5$Sol[,"traitAP.Dist_M"])
posterior.mode(GBmat.MCMC.5$Sol[,"traitAP.Dist_F"]-GBmat.MCMC.5$Sol[,"traitAP.Dist_M"])
HPDinterval(GBmat.MCMC.5$Sol[,"traitAP.Dist_F"]-GBmat.MCMC.5$Sol[,"traitAP.Dist_M"])
table(GBmat.MCMC.5$Sol[,"traitAP.Dist_F"]-GBmat.MCMC.5$Sol[,"traitAP.Dist_M"]>0) # 0.95

# Breding values NOTE: Randomization of breeding values is probably a better solution
with(Data_BLUPS.2, t.test(Latency~Sex))
with(Data_BLUPS.2, t.test(OF.Dist~Sex))
with(Data_BLUPS.2, t.test(AP.Dist~Sex))

# Va
# Latency (Females - Males)
posterior.mode(GBmat.MCMC.5$VCV[,1]-GBmat.MCMC.5$VCV[,22])
HPDinterval(GBmat.MCMC.5$VCV[,1]-GBmat.MCMC.5$VCV[,22])
table(GBmat.MCMC.5$VCV[,1]-GBmat.MCMC.5$VCV[,22]>0) # 0.68

# OF.Dist (Females - Males)
posterior.mode(GBmat.MCMC.5$VCV[,8]-GBmat.MCMC.5$VCV[,29])
HPDinterval(GBmat.MCMC.5$VCV[,8]-GBmat.MCMC.5$VCV[,29])
table(GBmat.MCMC.5$VCV[,8]-GBmat.MCMC.5$VCV[,29]<0)# 0.95

# AP.Dist (Females - Males)
posterior.mode(GBmat.MCMC.5$VCV[,15]-GBmat.MCMC.5$VCV[,36])
HPDinterval(GBmat.MCMC.5$VCV[,15]-GBmat.MCMC.5$VCV[,36])
table(GBmat.MCMC.5$VCV[,15]-GBmat.MCMC.5$VCV[,36]<0)# 0.92

#### 4.4 Figure 2: Full Gmf matrix ####
library(ellipse)
GBmat.2
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")   
colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(length(GBmat.2[1,]))
colors = colramp(100)

tiff(filename="Corrplot_Gmf.tiff", units="in",width=4,height=4, res=1200)
plotcorr(GBmat.2, col=colors[((GBmat.2 + 1)/2) * 100], type= "lower", mar=c(0,0,0,0))
dev.off()



tiff(filename="Corrplot_Gmf.tiff", units="in",width=4,height=4, res=1200)
#plotcorr(xc, col=colors[5*xc + 6], type= "lower", mar=c(0,0,0,0))
plotcorr(GBmat.2, col=colors[5*GBmat.2 + 6], type= "lower", mar=c(0,0,0,0))
#plotcorr(GBmat, col=colors[5*GBmat + 6], type= "lower", mar=c(0,0,0,0))
dev.off()

tiff(filename="Corrplot_Gmf_upper.tiff", units="in",width=4,height=4, res=1200)
#plotcorr(xc, col=colors[5*xc + 6], type= "lower", mar=c(0,0,0,0))
plotcorr(GBmat.2, col=colors[5*GBmat.2 + 6], type= "upper", numbers = T, mar=c(0,0,0,0))
#plotcorr(GBmat, col=colors[5*GBmat + 6], type= "lower", mar=c(0,0,0,0))
dev.off()

# Corrplot version
library(corrplot)
library(RColorBrewer)

tiff(filename="Corrplot_Gmf_2.tiff", units="in",width=8,height=8, res=1200)
corrplot.mixed(GBmat.2, lower = "ellipse", upper = "number",number.cex = 2)
dev.off()

corrplot.mixed(GBmat.2, col = brewer.pal(n = 8, name = "RdBu"))

colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")   
cols=colorRampPalette(colors)
cols=brewer.pal(n = 8, name = "RdBu")

corrplot.mixed(GBmat.2, lower = "ellipse", upper = "number",number.cex = 2,  
               col = cols)


#### 4.4.1 Pmcmc for Figure 2 ####
# Sex-specific
# LatxOF
hist(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,2]-
       posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,23])
posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,2]-
                 posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,23])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,2]-
              posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,23])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,2]-
        posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,23]>0)

#LatxAP
hist(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,3]-
       posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,24])
posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,3]-
                 posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,24])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,3]-
              posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,24])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,3]-
        posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,24]>0)

#OFxAP
hist(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,9]-
       posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,30])
posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,9]-
                 posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,30])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,9]-
              posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,30])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,9]-
        posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,30]>0)

# B matrix
# Diagonal
posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,4])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,4])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,4]>0)

posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,11])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,11])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,11]>0)

posterior.mode(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,18])
HPDinterval(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,18])
table(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,18]>0)


#### 5.Additional analyses ####
##### 5.1 Compare delta r for sex-specific syndrome ####
delta.r=data.frame(Corr=factor(c(rep("Latency x Open-field activity",1000),
                                 rep("Latency x Antipredator response",1000),
                                 rep("Open-field activity x Antipredator response",1000)),
                                levels = c("Latency x Open-field activity",
                                           "Latency x Antipredator response",
                                           "Open-field activity x Antipredator response")),
                   delta.dist=c(posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,2]-
                                  posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,23],
                                posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,3]-
                                  posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,24],
                                posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,9]-
                                  posterior.cor(GBmat.MCMC.5$VCV[,1:36])[,30]))
str(delta.r)

delta.r.plot=ggplot(delta.r, aes(x=delta.dist)) +
  facet_wrap(~Corr) +
  geom_histogram(aes(y=..density..), color="black",fill="white",) + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  xlab("delta r") + ylab("") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     strip.text = element_text(size=12))
delta.r.plot

##### 5.2 Vector correlation #####
vec.corr<-function(z1=z1,z2=z2){
  (sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}

(FxM.vcorr.gmax<-abs(vec.corr(eigen(GBmat.5[1:3,1:3])$vectors[,1],
                              eigen(GBmat.5[4:6,4:6])$vectors[,1])))
(FxM.vcorr.g2<-abs(vec.corr(eigen(GBmat.5[1:3,1:3])$vectors[,2],
                              eigen(GBmat.5[4:6,4:6])$vectors[,2])))
(FxM.vcorr.g3<-abs(vec.corr(eigen(GBmat.5[1:3,1:3])$vectors[,3],
                            eigen(GBmat.5[4:6,4:6])$vectors[,3])))

FxM.vcorr.gmax.mcmc <- vector("numeric",1000)
for(i in 1:1000){
  Gmf=array(,dim=c(6,6,1000))
  Gf=array(,dim=c(3,3,1000))
  Gm=array(,dim=c(3,3,1000))
  Gmf[,,i]=matrix(GBmat.MCMC.5$VCV[,1:36][i,], nrow=6)
  Gf[,,i]=Gmf[1:3,1:3,i]
  Gm[,,i]=Gmf[4:6,4:6,i]
  
  FxM.vcorr.gmax.mcmc[i]=abs(vec.corr(eigen(Gf[,,i])$vectors[,1],
                       eigen(Gm[,,i])$vectors[,1]))
  }
mean(as.mcmc(FxM.vcorr.gmax.mcmc))
posterior.mode(as.mcmc(FxM.vcorr.gmax.mcmc))
HPDinterval(as.mcmc(FxM.vcorr.gmax.mcmc))
plot(as.mcmc(FxM.vcorr.gmax.mcmc))

mean(as.mcmc(acos(FxM.vcorr.gmax.mcmc)*(180/pi)))
posterior.mode(as.mcmc(acos(FxM.vcorr.gmax.mcmc)*(180/pi)))
HPDinterval(as.mcmc(acos(FxM.vcorr.gmax.mcmc)*(180/pi)))
plot(as.mcmc(acos(FxM.vcorr.gmax.mcmc)*(180/pi)))



# angle
acos(abs(FxM.vcorr.gmax))*(180/pi)
# Mantel test
mantel.test(GBmat.5[1:3,1:3],GBmat.5[4:6,4:6], graph = T)




##### 5.3 gmax alignment #####
# Project BLUPS on the g1 and g2 by sex
# F
Data_BLUPS.2_F=subset(Data_BLUPS.2, Sex=="F")
Gmax_F=as.data.frame(as.matrix(Data_BLUPS.2_F[,c(9:11)]) %*% eigen(GBmat.5[1:3,1:3])$vectors[,1:2])
# rotate female data by 35 degrees
library(spdep)
alpha=-mean(as.mcmc(acos(FxM.vcorr.gmax.mcmc)*(180/pi)))
alpha=alpha*90*pi/180
Gmax_F=Rotation(Gmax_F,alpha)
colnames(Gmax_F)=c("gmax","g2")

# M
Data_BLUPS.2_M=subset(Data_BLUPS.2, Sex=="M")
Gmax_M=as.data.frame(as.matrix(Data_BLUPS.2_M[,c(9:11)]) %*% eigen(GBmat.5[4:6,4:6])$vectors[,1:2])
names(Gmax_M)=c("gmax","g2")

Gmax_BLUPS=rbind(Gmax_F,Gmax_M)
Gmax_BLUPS$Sex=factor(c(rep("Females",length(Gmax_F[,1])),
                        rep("Males",length(Gmax_M[,1]))))

gmax.plot=ggplot(Gmax_BLUPS, aes(x= gmax, y=g2, group=Sex, fill=Sex, color=Sex)) +
  stat_ellipse(geom = "polygon", alpha = .45, size=1, color ="black") + 
  geom_point(size=3, alpha = .5) + 
  xlab(expression(bold(g[max]))) + 
  ylab(expression(bold(g[2]))) + 
  scale_color_wsj() + scale_fill_wsj() + coord_equal(ratio=1) +
  theme_test() + theme(aspect.ratio=1,
                       legend.position=c(.1, .9), legend.title = element_blank(),
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16),
                       strip.text = element_text(size=12))
gmax.plot



##### 5.4 Short presentation figure: Sex-specific syndrome and Shelter emergence cross sex ####
delta_r.gmax.plot=delta.r.plot/gmax.plot
ggsave(delta_r.gmax.plot, filename = "delta_r.gmax.plot.jpeg", width = 12, height = 12)

delta.r.plot.LxOF=ggplot(subset(delta.r,Corr=="Latency x Open-field activity"), aes(x=delta.dist)) +
  #facet_wrap(~Corr) +
  geom_histogram(aes(y=..density..), color="black",fill="white",) + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  xlab("delta r") + ylab("") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)

delta.r.plot.LxAP=ggplot(subset(delta.r,Corr=="Latency x Antipredator response"), aes(x=delta.dist)) +
  #facet_wrap(~Corr) +
  geom_histogram(aes(y=..density..), color="black",fill="white",) + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  xlab("delta r") + ylab("") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)

delta.r.plot.OFxAP=ggplot(subset(delta.r,Corr=="Open-field activity x Antipredator response"), aes(x=delta.dist)) +
  #facet_wrap(~Corr) +
  geom_histogram(aes(y=..density..), color="black",fill="white",) + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  xlab("delta r") + ylab("") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)


LatxOF=ggplot(Data_BLUPS.2, aes(x=OF.Dist, y=Latency, color=Sex)) +
  geom_point(size=3, alpha=.5) + scale_color_wsj() +
  xlab("Open-field activity") + ylab("Shelter emergence") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)
  
LatxAP=ggplot(Data_BLUPS.2, aes(x=AP.Dist, y=Latency, color=Sex)) +
  geom_point(size=3, alpha=.5) + scale_color_wsj() +
  xlab("Antipredator activity") + ylab("Shelter emergence") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)

OFxAP=ggplot(Data_BLUPS.2, aes(x=OF.Dist, y=AP.Dist, color=Sex)) +
  geom_point(size=3, alpha=.5) + scale_color_wsj() +
  ylab("Antipredator activity") + xlab("Open-field activity") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)
Fig2=(LatxOF+LatxAP+OFxAP)/(delta.r.plot.LxOF+delta.r.plot.LxAP+delta.r.plot.OFxAP)
Fig2

ggsave(Fig2, filename = "Fig2.jpeg", width = 12, height = 12)


Lat_FxLat_M=ggplot(Data_BLUPS.2, aes(y=Latency_F, x=Latency_M)) +
  geom_point(size=3, alpha=.5) + scale_color_wsj() +
  ylab("Shelter emergence (F)") + xlab("Shelter emergence (M)") +
  theme_bw() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)
Lat_FxLat_M
ggsave(Lat_FxLat_M, filename = "Lat_FxLat_M.jpeg", width = 6, height = 6)

Fig3=(LatxOF+LatxAP+OFxAP+Lat_FxLat_M)+plot_layout(ncol = 4)
Fig3

ggsave(Fig3, filename = "Fig3.jpeg", width = 14, height = 7)

#Intro slide
Fig_intro=ggplot(Data_BLUPS.2, aes(x=OF.Dist, y=Latency, color=Sex, group=Sex, fill=Sex)) +
  stat_ellipse(geom = "polygon", alpha = .7, size=1, color ="black") +
  geom_point(size=3, alpha=.5) +
  scale_color_wsj() + scale_fill_wsj() +
  xlab("Open-field activity") + ylab("Shelter emergence") +
  theme_base() + theme(legend.position="none",
                     axis.title = element_text(size=20),
                     axis.text = element_text(size=16),
                     aspect.ratio=1)
ggsave(Fig_intro, filename = "Fig_intro.jpeg", width = 4, height = 4)



##### Extra #####
#Pmat
Data_P_M=subset(Data,Sex=="M")
Data_P_F=subset(Data,Sex=="F")


Pmat_M=ggpairs(Data_P_M[,c("Latency.Cens","OF.Dist","AP.Dist")],
               lower=list(continuous="smooth"),
               diag = list("continuous")) + theme_base()
Pmat_M

Pmat_F=ggpairs(Data_P_F[,c("Latency.Cens","OF.Dist","AP.Dist")],
               lower=list(continuous="smooth"),
               diag = list("continuous")) + theme_base()
Pmat_F


Pmat_M=ggpairs(Data_P_M[,c("Latency.Cens","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")],
              lower=list(continuous="smooth"),
              diag = list("continuous")) + theme_base()
Pmat_M

Pmat_F=ggpairs(Data_P_F[,c("Latency.Cens","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")],
               lower=list(continuous="smooth"),
               diag = list("continuous")) + theme_base()
Pmat_F



