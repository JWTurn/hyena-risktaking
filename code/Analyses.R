### Analyses ====
# Julie Turner
# prepped for publication: Sept. 11 2019


### Packages ----

libs <- c('data.table', 'lme4', 'bbmle', 'arm', 'car', 'survival', 'lmtest', 'foreign', 'MASS', 'Hmisc', 'reshape2')
lapply(libs, require, character.only = TRUE)


### Input data ----
path <- 'data/'


data <- readRDS(paste0(path, 'risktakn_pub.Rds'))
data.long <- readRDS(paste0(path, 'risktakn_long_pub.Rds'))

# we call the mock intruder 'target' because it is an archery target

data.target <- data[!(is.na(bold2))]


### risk-taking disturbance differences ----

#### Categorical analyses OLR ####

### exploration scores
data.target$exploration <- factor(data.target$exploration)
#subsetting to those with exploration score (i.e. approached within 30m)
data.target.exp <- data.target[exploration!=0]

olr.exp <- polr(exploration~ disturbance * (sex + rank_level + age )+ offset(log(dist_notice)), data.target.exp, Hess = T)
summary(olr.exp)

#calculating p-values
ctable.exp<- coef(summary(olr.exp))
p.exp <- pnorm(abs(ctable.exp[, "t value"]), lower.tail = FALSE) * 2
ctable.exp <- cbind(ctable.exp, "p value" = p.exp)


### boldness scores
data.target$bold2 <- factor(data.target$bold2, levels=c('fearful', 'neutral', 'bold', 'attack'))

olr.bold2 <- polr(bold2~ disturbance * (sex + rank_level + age )+ offset(log(dist_notice)), data.target, Hess = T)
summary(olr.bold2)
ctable.bold2<- coef(summary(olr.bold2))
p.bold2 <- pnorm(abs(ctable.bold2[, "t value"]), lower.tail = FALSE) * 2
ctable.bold2 <- cbind(ctable.bold2, "p value" = p.bold2)


#### ANOVA Continuous analyses ####
# closest distance
close.lm <- lm(log(closest_dist+1)~disturbance * sex + rank_level + disturbance:rank_level + age + disturbance:age + sex:rank_level
               + offset(dist_notice), data.target) 
close.aov <- aov(close.lm)
summary(close.aov)
TukeyHSD(close.aov)

# time to closest distance
tclose.lm <- lm(log(t_closest)~disturbance * sex + rank_level + disturbance:rank_level + age + disturbance:age + sex:rank_level
                + offset(dist_notice), data.target) 
tclose.aov <- aov(tclose.lm)
summary(tclose.aov)
TukeyHSD(tclose.aov)

# proportion of time within 30m
proptwin30.lm <- lm(proptwin30m~disturbance * sex + rank_level + disturbance:rank_level + age + disturbance:age + sex:rank_level
                    + offset(dist_notice), data.target) 
proptwin30.aov <- aov(proptwin30.lm)
summary(proptwin30.aov)
TukeyHSD(proptwin30.aov)



### Risk-taking survival analyses ----

data.surv <- data.target[sex == 'f' & !(is.na(longevity))]

data.surv[, uniqueN(id), by = c('disturbance', 'surv')] # n = 54, 15 dead

data.surv[,'ageYears'] <- (data.surv$date - data.surv$birthdate)/365
surv.obj<- with(data.surv, Surv(longevity, surv, type = 'right'))

# closest distance
surv.close <- coxph(surv.obj~  closest_dist * disturbance + closest_dist * rank_level + closest_dist:ageYears, data.surv)
sum.surv.close <- summary(surv.close) 
View(sum.surv.close$coefficients)

# time to closest distance
surv.tclose <- coxph(surv.obj~ t_closest * disturbance + t_closest * rank_level + t_closest:ageYears, data.surv)
sum.surv.tclose <- summary(surv.tclose) 
View(sum.surv.tclose$coefficients)

# proportion of time within 30m
surv.twin30m <- coxph(surv.obj~  proptwin30m * disturbance + proptwin30m * rank_level + proptwin30m:ageYears, data.surv)
sum.surv.twin30m <- summary(surv.twin30m)  
View(sum.surv.twin30m$coefficients)



### Risk-taking consistency ----

#### Target - baited box ####
# subsetting data to only model intruder and baited box tests
data.tarbox <- data.long[test == 'target' | test == 'box']
data.tarbox <- data.tarbox[id %chin% data.tarbox$id[duplicated(data.tarbox$id)]]
data.tarbox$id <- as.factor(data.tarbox$id)

#1 closest distance
# model with ID
box.close.ID<-lmer(closest_dist~ disturbance + sex + rank_level + test + (1|id), data.tarbox[complete.cases(data.tarbox[,'t_closest']),], REML = T)
# model without ID
box.close.noID<-lm(closest_dist~ disturbance + sex + rank_level + test, data.tarbox[complete.cases(data.tarbox[,'t_closest']),])

# testing which is better
lrtest(box.close.ID, box.close.noID) 
AICctab(box.close.ID, box.close.noID)


#2 time to closest distance
box.lat.ID<-lmer(t_closest~ disturbance + sex + rank_level + test + (1|id), data.tarbox[complete.cases(data.tarbox[,'t_closest']),], REML = F)
box.lat.noID<-lm(t_closest~ disturbance + sex + rank_level + test, data.tarbox[complete.cases(data.tarbox[,'t_closest']),])


lrtest(box.lat.ID, box.lat.noID) 
AICctab(box.lat.ID, box.lat.noID)




#### Target - Lion ####
# subsetting data to only model intruder and approaching lion tests
data.tarlion <- data.long[test == 'target' | test == 'lion']
data.tarlion <- data.tarlion[id %chin% data.tarlion$id[duplicated(data.tarlion$id)]]
data.tarlion$id <- as.factor(data.tarlion$id)


#1 closest distance
lion.close.ID <-lmer(log(closest_dist+1)~ disturbance + sex + rank_level + test + (1|id), data.tarlion, REML = T)
lion.close.noID <-lm(log(closest_dist+1)~ disturbance + sex + rank_level + test, data.tarlion)

lrtest(lion.close.ID, lion.close.noID) 
AICctab(lion.close.ID, lion.close.noID)


#2 Proportion of time w/in 30m of model & standardized distance to lions
lion.stdclose.ID <-lmer(log(stdClosest+1)~ disturbance + sex + rank_level + test + (1|id), data.tarlion, REML = T)
lion.stdclose.noID <-lm(log(stdClosest+1)~ disturbance + sex + rank_level + test, data.tarlion)

lrtest(lion.stdclose.ID, lion.stdclose.noID) 
AICctab(lion.stdclose.ID, lion.stdclose.noID)



#### Lion - baited box ####
# subsetting data to only approaching lions and baited box tests
data.lionbox <- data.long[test == 'lion' | test == 'box']
data.lionbox <- data.lionbox[id %chin% data.lionbox$id[duplicated(data.lionbox$id)]]
data.lionbox$id <- as.factor(data.lionbox$id)

#1 closest distance
libox.close.ID <-lmer(log(closest_dist+1)~ disturbance + sex + rank_level + test + (1|id), data.lionbox, REML = F)
libox.close.noID <-lm(log(closest_dist+1)~ disturbance + sex + rank_level + test, data.lionbox)

lrtest(libox.close.ID, libox.close.noID)  
AICctab(libox.close.ID, libox.close.noID)  




