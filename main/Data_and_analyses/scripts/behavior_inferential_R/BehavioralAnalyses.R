#this is the script to reproduce results in the manuscript entitled Computational Underpinnings of Partisan Information Processing Biases and Associations with Depth of Cognitive Reasoning 

#specify the location of this script in your path, for me it's:
here::i_am("Data_and_analyses/scripts/behavior_inferential_R/BehavioralAnalyses.R")

#Load helper functions
source(here::here("Data_and_analyses", "scripts","behavior_inferential_R","helper_functions.R"))

#Import and load libraries
library(groundhog) # we are using groundhon to make this script reproducible: https://groundhogr.com/using/
pkgs <-  c("tidyverse","sjstats","ggplot2","lme4","lmerTest","Hmisc","car","lmtest","ROCR", "ggeffects", "sjPlot", "lavaan", "Rmisc", "here")
groundhog.day <- '2022-05-22'
groundhog.library(pkgs, groundhog.day)
#please follow instructions in console if packages do not load successfully. 

#load data  
#####
DF <- read.csv(here::here("Data_and_analyses","data","cleaned","behavioral_data.csv"), header=T, stringsAsFactors = FALSE, na.strings=c("","NA"))
DF <- DF[,-c(1)] #that one extra column that always shows up

#data for RT analysis (new) -- same as input for ddm models
DF_rt <- read.csv(here::here("Data_and_analyses","data","cleaned","rt_behavioral_data_sec.csv"), header=T, stringsAsFactors = FALSE, na.strings=c("","NA"))

#load DDM trace for visuals
DF_post = read.csv(here::here("Data_and_analyses","saved_hddm_models_and_parms" ,"model_outputs","trace_processed","FullModel_trace_all.csv"),header = F) 
colnames(DF_post) = c("a","t",
                 "z_int","v_int",
                 "v_cond","v_stim_n","v_stim_weight")
DF_post$v_bias <- DF_post$v_int + DF_post$v_cond #difference between ingroup and outgroup drift rate
#####

#demographics
#####
m_age <- mean(DF$age, na.rm = T)
sd_age <- sd(DF$age, na.rm = T)
DF_agg <- DF[!duplicated(DF$Participant),]
prop.table(table(DF_agg$gender))
gender <- prop.table(table(DF_agg$gender))
ethnicity<- prop.table(table(DF_agg$ethnicity))
#####

#center variables for evidence seen 
#####
DF$stim_weight  <- scale(DF$stim_weight , center=TRUE, scale=FALSE)  
DF$stim_n  <- scale(DF$stim_n , center=TRUE, scale=FALSE)  
DF$diff_in_min_out <- (DF$meanIngroup - DF$meanOutgroup) #negative values denote more evidence for outgroup
DF$diff_in_min_out_centered  <- scale(DF$diff_in_min_out , center=TRUE, scale=FALSE)  
DF$numDiff0 <-(DF$numIngroup - DF$numOutgroup) #negative values mean more outgroup
DF$numDiff_centered0  <- scale(DF$numDiff0 , center=TRUE, scale=FALSE) 
DF$varDiff <-(DF$SDIngroup - DF$SDOutgroup) 
DF$varDiff_centered  <- scale(DF$varDiff , center=TRUE, scale=FALSE) 
#####


##Behavioral Analyses
#####

### Partisan Motivations Bias Information Processing
#gchoice: ingroup coded 1, outgroup coded 0
M1 <- glmer(gChoice~numDiff_centered0+diff_in_min_out_centered+(numDiff_centered0+diff_in_min_out_centered|Participant), data = DF, family = "binomial")
M1.coef <- summary(M1)
M1.CI <- confint(M1)
M1.OR <- exp(fixef(M1))

### Partisans More Accurate When Ingroup is Better
#pre-registered analyses show similar results
accByCond <-DF %>% group_by(whosBetter, Participant) %>% dplyr::summarize(NumIn=n(),accuracy=sum(correct)/n())
pre_reg_Model <- t.test(accByCond$accuracy[which(accByCond$whosBetter=="ingroup")], accByCond$accuracy[which(accByCond$whosBetter=="outgroup")])

#analysis in manuscript
#performance based on condition 
M2 <- glmer(correct~as.factor(whosBetter)+(1|Participant), data = DF, family = "binomial")
M2.coef <-summary(M2)
M2.CI<- confint(M2)
M2.OR <- exp(fixef(M2))

#affiliation moderation
#polAffil2 = binary "democrat" vs. "republican" 
M3 <- glmer(correct~as.factor(whosBetter)*as.factor(polAffil2)+(1|Participant), data = DF, family = "binomial")
M3.coef <-summary(M3)
M3.CI<- confint(M3)
M3.OR <- exp(fixef(M3))

### Partisans Faster to Categorize Ingroup as More Honest
M4.a <- lmer(rt_log~as.factor(response) + (1|subj_idx), DF_rt)
M4.a.coef <-summary(M4.a) #significant main effect of response type, no effect of affiliation
M4.a.CI<- confint(M4.a)
#adding affiliation as covariate
M4.b <- lmer(rt_log~as.factor(response)+as.factor(Affil) + (1|subj_idx), DF_rt)
M4.b.coef <- summary(M4.b) #interaction term is not significant
#adding affiliation as interaction
M4.c <- lmer(rt_log~as.factor(response)*as.factor(Affil) + (1|subj_idx), DF_rt)
M4.c.coef <- summary(M4.c) #interaction term is not significant
#model comparison
lrt_M4<- anova(M4.a,M4.b, M4.c) #non significant LRT
#####

##Computational Results 
#####

#Starting point
options(scipen=999)
ci_hdi_z <- ci(DF_post$z_int, method = "HDI", ci = 0.95)
ci_hdi_z_low <- ci_hdi_z$CI_low
ci_hdi_z_high <- ci_hdi_z$CI_high
z_bias_mean <-mean(DF_post$z_int)
z_bias_prob <- pnorm(.5, mean(DF_post$z_int), sd(DF_post$z_int), lower.tail = F)

#Drift rate
ci_hdi_v <- ci(DF_post$v_bias , method = "HDI", ci = 0.95)
ci_hdi_v_low <- ci_hdi_v$CI_low
ci_hdi_v_high <- ci_hdi_v$CI_high
v_bias_mean <- mean(DF_post$v_bias)
v_bias_prob <- pnorm(0, mean(DF_post$v_bias), sd(DF_post$v_bias), lower.tail = F) #rounds up
#####

##Individual Differences
#####

##Confirmatory factor analysis
DFAgg <- DF[!duplicated(DF$Participant),] #this needs to be run on level 2, so remove duplicate ids

#run CFA for three measures related to reasoning
composite.model <- 'reasoning=~NFC+b1*CRT+b1*ANS
             reasoning~~reasoning '

composite.model <-'reasoning=~NFC+CRT+ANS
            reasoning~~1*reasoning'

composite.model <-'reasoning=~NFC+CRT+ANS
            reasoning~~reasoning
            CRT~~b1*CRT
            ANS~~b1*ANS'

composite.fit <- cfa(composite.model, data = DFAgg)
composite.fit<- summary(composite.fit, fit.measures = TRUE)

##Create composite score
toZ <- c("NFC", "CRT", "ANS")
DFAgg[toZ] <- lapply(toZ, function(x) scale(DFAgg[[x]])) #z-scoring prior to entering because each measure is on different scale

DFAgg$reasoning_comp <- rowMeans(DFAgg[toZ])
DF$reasoning_comp <- rowMeans(DFAgg[toZ])

DFAgg$reasoning_comp <- composite(DFAgg[,c("NFC", "CRT", "ANS")], Zitems = T) #z-scoring prior to entering because each measure is on different scale
DF$reasoning_comp <- composite(DF[,c("NFC", "CRT", "ANS")], Zitems = T)

#Main Effect of Cognitive Reasoning
M4 <- glmer(correct~scale(reasoning_comp)+(1|Participant), data = DF, family = "binomial")
M4.coef <-summary(M4)
M4.CI<- confint(M4)
M4.OR <- exp(fixef(M4))

#Cognitive Reasoning Moderates Behavioral Task Performance
M5 <- glmer(correct~as.numeric(affilStrength)*as.factor(whosBetter)*scale(reasoning_comp)+(1|Participant), data = DF, family = "binomial")
M5.coef <-summary(M5)
M5.CI<- confint(M5)
M5.OR <- exp(fixef(M5))

#all individual differences measures estimated in separate models
#Cognitive reflection
CRT.M <- glmer(correct~as.numeric(affilStrength)*as.factor(whosBetter)*scale(CRT)+(1|Participant), data = DF, family = "binomial")
CRT.M.coef <-summary(CRT.M)
CRT.M.CI<- confint(CRT.M)
CRT.M.OR <- exp(fixef(CRT.M))

#Numeracy
ANS.M <- glmer(correct~as.numeric(affilStrength)*as.factor(whosBetter)*scale(ANS)+(1|Participant), data = DF, family = "binomial")
ANS.M.coef <-summary(ANS.M)
ANS.M.CI<- confint(ANS.M)
ANS.M.OR <- exp(fixef(ANS.M))

#Need for cognition
NFC.M <- glmer(correct~as.numeric(affilStrength)*as.factor(whosBetter)*scale(NFC)+(1|Participant), data = DF, family = "binomial")
NFC.M.coef <-summary(NFC.M)
NFC.M.CI<- confint(NFC.M)
NFC.M.OR <- exp(fixef(NFC.M))
#####

##Cognitive Reasoning has Divergent Correlations with Drift Rate and Starting Point
colnames(DFAgg)
Cor_table <- apa.cor.table(DFAgg[c("NFC", "CRT", "ANS", "reasoning_comp", "V_int", "V_cond","v_bias","Z")], filename="CorTable_2022.doc", table.number=1)
cor_v_res <- cor.test(DFAgg$v_bias, DFAgg$reasoning_comp)
cor_z_res <- cor.test(DFAgg$Z, DFAgg$reasoning_comp)
cor_z__v <- cor.test(DFAgg$Z, DFAgg$v_bias)
#####

##Create Figures
#####

#Figure 3
#for visualization only (i.e., not for statistical interpretations)
DF$numDiff_centered0 <- as.numeric(DF$numDiff_centered0 )
DF$diff_in_min_out_centered <- as.numeric(DF$diff_in_min_out_centered )
#had to make numeric because I was running into issues (see: https://stackoverflow.com/questions/43218003/error-variables-were-specified-with-different-types-from-the-fit)
sub_avg = group_by(DF,Participant,numDiff_centered0,diff_in_min_out_centered, whosBetter) %>% 
  dplyr::summarise(Avg = mean(gChoice, na.rm=T)) #get averages for each variable
sub_avg <- sub_avg[!is.na(sub_avg$diff_in_min_out_centered),] #remove missing data as it generates error
res.glm = glm(gChoice ~ numDiff_centered0 + diff_in_min_out_centered, DF, family = binomial(link="probit"))#run probit model
pred.data = data.frame(numDiff_centered0 = rep(seq(-3, 3, len = 149),4), #generate data using predict, choices at all levels
                       whosBetter = c(rep("ingroup",298),rep("outgroup",298)),
                       diff_in_min_out_centered = rep(seq(-30, 30, len = 149),4))

pred.data$numDiff_centered0 <- as.numeric(pred.data$numDiff_centered0)#make numeric due to issue outlined above
pred.data$diff_in_min_out_centered <- as.numeric(pred.data$diff_in_min_out_centered)
gChoice = predict.glm(res.glm,pred.data,type = "response")

#run separately for ingroup and outgroup
res.glmIn = glm(gChoice ~ numDiff_centered0 + diff_in_min_out_centered+whosBetter, DF, family = binomial(link="probit"))
pred.dataIn = data.frame(numDiff_centered0 = rep(seq(-3, 3, len = 149),4),
                         whosBetter = c(rep("ingroup",596)),
                         diff_in_min_out_centered = rep(seq(-45, 45, len = 149),4))

pred.dataIn$numDiff_centered0 <- as.numeric(pred.dataIn$numDiff_centered0 )
pred.dataIn$diff_in_min_out_centered <- as.numeric(pred.dataIn$diff_in_min_out_centered )
gChoice_in = predict.glm(res.glmIn,pred.dataIn,type = "response")

res.glmOut = glm(gChoice ~ numDiff_centered0 + diff_in_min_out_centered+whosBetter, DF, family = binomial(link="probit"))
pred.dataOut = data.frame(numDiff_centered0 = rep(seq(-3, 3, len = 149),4),
                          whosBetter = c(rep("outgroup",596)),
                          diff_in_min_out_centered = rep(seq(-45, 45, len = 149),4))
pred.dataOut$numDiff_centered0 <- as.numeric(pred.dataOut$numDiff_centered0 )
pred.dataOut$diff_in_min_out_centered <- as.numeric(pred.dataOut$diff_in_min_out_centered )
gChoice_out = predict.glm(res.glmOut,pred.dataOut,type = "response")

#create figure
Figure3 <- ggplot(sub_avg, aes(color = as.factor(whosBetter))) +
  stat_summary(data = sub_avg, aes(x=diff_in_min_out_centered,y=Avg),
               fun=mean,geom="point",size=.3, alpha = 0.3) +
  geom_line(color="darkgoldenrod3",data = pred.dataIn, aes(x=(diff_in_min_out_centered),y=gChoice_in),size=1.5, alpha = 0.9) +
  geom_line(color="darkorchid4",data = pred.dataOut, aes(x=(diff_in_min_out_centered),y=gChoice_out),size=1.5, alpha = 0.9) +
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        legend.position = c(0.77, 0.19),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x  = element_text(size=13,color="black"),
        axis.text.y  = element_text(size=13,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 15))  +
  xlab("Evidence favoring more honest candidate") + ylab("P(Choose more honest candidate)") +  
  scale_x_continuous(breaks = c(-45,-30,-15,0, 15, 30, 45),labels=c("-45","-30", "-15", "0", "15","30", "45"))+ 
  coord_cartesian(xlim=c(-45, 45)) +
  scale_color_manual(name="Candidate Honesty",  labels=c("Outgroup more honest","Ingroup more honest"), values=c("darkorchid4", "darkgoldenrod3"))+
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

#RT figure
samplingDescr_Affil <- summarySEwithin(DF_rt, idvar='subj_idx', measurevar = 'rt', betweenvars = "Affil", withinvars = "response", na.rm = TRUE)
samplingDescr_Fav  <- summarySE(DF_rt, measurevar="rt", groupvars=c("response"))

#plot with SE within
rt_fav<- ggplot(samplingDescr_Fav, aes(x = response, y = rt))+
  geom_point(size = 4)+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=rt-se, ymax=rt+se), width=.1)+
  coord_cartesian(ylim = c(0, 3))+
  scale_fill_manual(values=c("lightblue","red"))+
  theme_classic() +
  coord_cartesian(ylim = c(2,3))+
  xlab("Response favorability") + ylab("Reaction time (seconds)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))+
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x  = element_text(size=13,color="black"),
        axis.text.y  = element_text(size=13,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 15))
#ggsave("rt_fav_2022.jpeg",  dpi = 700)

rt_fav_and_affil <- ggplot(samplingDescr_Affil, aes(x = response, y = rt, fill = Affil))+
  geom_point(size = 4)+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black') +
  geom_errorbar(aes(ymin=rt-se, ymax=rt+se), width=.1)+
  coord_cartesian(ylim = c(0, 3))+
  facet_wrap(~Affil)+
  scale_fill_manual(values=c("lightblue","red"))+
  theme_classic() +
  coord_cartesian(ylim = c(2,3))+
  xlab("Response favorability & political affiliation") + ylab("Reaction time (seconds)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))+
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x  = element_text(size=13,color="black"),
        axis.text.y  = element_text(size=13,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 15))
#ggsave("rt_fav_affil_2022.jpeg",  dpi = 700)

#Figure 4
model.DIC = NULL
model.DIC$Model = c('z & v','v', 'z')
model.DIC$DIC = c(-5589.12,-5399.21, -6.56) # See HDDM Analysis for DIC output
model.DIC = as.data.frame(model.DIC)
model.DIC$Model = factor(model.DIC$Model, levels = c('z & v','v','z'))
DIC.fig = ggplot(model.DIC, aes(x = Model, y = DIC)) + 
  geom_bar(stat = "identity",width = 0.3) +
  theme_classic()+
  theme(axis.text.x  = element_text(size=13,color="black"),
        axis.text.y  = element_text(size=13,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(face = "bold", size = 15))+
  scale_y_reverse()+ylab(paste("DIC", "\u394", "null"))+
  coord_cartesian(ylim = c(-50,-6000))

#Posterior parameters
z_int.fig = ggplot(data = DF_post, aes(x = z_int)) + 
  geom_histogram(aes(y = stat(ncount)), binwidth = .0015, fill=I("grey"), col=I("black"))+
  theme_classic()+ 
  ylab("Density") +  xlab("Starting point")+ theme(text=element_text(size=15, color = "black"))+
  theme(axis.text.x = element_text(size=15, color="black"), 
        axis.text.y = element_text(size=15,color="black"))+
  coord_cartesian(xlim = c(.5, .535), ylim = c(.03,1))+
  geom_vline(xintercept=c(.5), linetype="longdash", size = .7)

v_bias.fig = ggplot(data = DF_post, aes(x = v_bias)) + 
  geom_histogram(aes(y = stat(ncount)), binwidth = .0095, fill=I("grey"), col=I("black"))+
  theme_classic()+ 
  ylab("Density") +  xlab("Drift Rate \u394")+ theme(text=element_text(size=15, color = "black"))+
  theme(axis.text.x = element_text(size=15, color="black"), 
        axis.text.y = element_text(size=15,color="black"))+
  coord_cartesian(xlim = c(.1, .4), ylim = c(.03,1))
#add together into one figure
Figure4 <- ggarrange(DIC.fig, z_int.fig, v_bias.fig,labels = c("A", "B","C"), nrow = 1, font.label = list(size = 16, color = "black"))

#Figure 5
M4.model <- ggpredict(M4, terms = c("affilStrength", "reasoning_comp", "whosBetter"))
Figure5 <- ggplot(M4.model, aes(x = x, y = predicted, colour = group)) +
  geom_line() +
  geom_point()+
  facet_wrap(~facet)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.4,
                position=position_dodge(0.001)) +theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  xlab("Political Conviction") + 
  ylab("Accuracy (%)") +
  theme(legend.position = c(0.16, 0.13))+ labs(color='composite cognitive reasoning') + 
  scale_x_continuous(breaks = c(1,2,3),labels=c("Slightly", "Somewhat", "Very"))+ 
  scale_y_continuous(breaks = c(.45,.60, .75), labels = c("45%", "60%","75%"))+
  theme(axis.text.x = element_text(color="black", size=12), axis.text.y = element_text(color="black", size=12,  family="Times New Roman"),
        axis.title.x = element_text(color="black", size=16), axis.title.y = element_text(color="black", size=16,  family="Times New Roman"))+
  theme(text = element_text(size = 10,  family="Times New Roman"))+scale_color_hue(labels = c("-1 SD", "Mean", "+1 SD"))

#correlation figure between cog reasoning and DDM parameters

drift_cog <- ggplot(DFAgg, aes(x = reasoning_comp, y = v_bias))+geom_point()+
  geom_smooth(method =lm, colour = "grey3")+theme_classic() 
ggsave("drift_cog.jpeg", width = 4, height = 4, dpi = 700)

z_cog <-ggplot(DFAgg, aes(x = reasoning_comp, y = Z))+geom_point()+
  geom_smooth(method =lm, colour = "grey3")+theme_classic()
ggsave("Z_cog.jpeg", width = 4, height = 4, dpi = 700)

#####

##Print Results in the same order they appear in manuscript
print("Demographics")
print("Age")
print(m_age)
print("Gender")
print(gender)
print("Ethnicity")
print(ethnicity)
print("Political orientation")
print(polOrient)

print("Behavioral Analyses")
print("Partisan Motivations Bias Information Processing")
print(M1.coef)
print(M1.OR)
print("Figure 3")
print(Figure3)

print("Partisans More Accurate When Ingroup is Better")
print(M2.coef)
print(M2.OR)

print("Political afiliation moderates performance based on who is better")
print(M3.coef)
print(M3.OR)

print("Partisans Faster to Categorize Ingroup as More Honest")
print(M4.a.coef)
print(M4.a.CI)

print("Computational Analyses")
print("Starting point")
print("mean")
print(z_bias_mean)
print("lowerCI")
print(ci_hdi_z_low)
print("upperCI")
print(ci_hdi_z_high)
print("prob_test")
print(z_bias_prob)

print("Drift Rate")
print("mean")
print(v_bias_mean)
print("lowerCI")
print(ci_hdi_v_low)
print("upperCI")
print(ci_hdi_v_high)
print("prob_test")
print(v_bias_prob)
print("Figure 4")
print(v_bias_prob)

print("Individual Differences")
print("Main Effect of Cognitive Reasoning")
print(M4.coef)
print(M4.OR)
print("Cognitive Reasoning Moderates Behavioral Task Performance")
print(M5.coef)
print(M5.OR)
print(Figure5)

print("Correlation Matrix")
print(Cor_table)
print("Correlations between starting point and cognitive reasoning")
print(cor_v_res)
print("Correlations between drift rate and cognitive reasoning")
print(cor_z_res)
print("Correlations between drift rate and starting point")
print(cor_z__v)