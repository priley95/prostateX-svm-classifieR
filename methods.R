# SVM for prostate lesion classification: ProstateX #
## Libraries and data ----
### Libraries
library(tidyverse)
library(e1071) #svm, tune, tune.svm
library(kernlab)
library(mlr)
library(pROC)

### Data by finding 
# 1. T2-tra and ADC-tra ----
tra <- cbind(t2tra.df, tra.track) 
adc <- cbind(adc_transeq.df, adc_trans.track)

# Join data
tra.by.finding <- tra %>%
  inner_join(adc, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ProxID", "fid")) %>%  #, "AStra", "PZtra", "SVtra", "TZtra")) %>%#remove after being joined by
  rename(ClinSig = ClinSig.y) %>%
  mutate(ClinSig = as.factor(ClinSig)) %>%
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         adc1 = "1.y",  adc2 = "2.y", adc3 = "3.y", adc4 = "4.y", adc5 = "5.y", adc6 = "6.y", adc7 = "7.y", adc8 = "8.y", 
         adc9 = "9.y", adc10 = "10.y", adc11 = "11.y", adc12 = "12.y", adc13 = "13.y", adc14 = "14.y", adc15 = "15.y", adc16 = "16.y", 
         adc17 = "17.y", adc18 = "18.y", adc19 = "19.y", adc20 = "20.y", adc21 = "21.y", adc22 = "22.y", adc23 = "23.y", adc24 = "24.y", adc25 = "25.y")#,
         #AS = "ASadc", PZ = "PZadc", SV = "SVadc", TZ = "TZadc")
  
scaled.tra.by.finding <- scale(tra.by.finding[,1:50]) 

scaled.tra.by.finding <- cbind(scaled.tra.by.finding, as.factor(tra.by.finding$ClinSig))
colnames(scaled.tra.by.finding)[51] <- "ClinSig"
scaled.tra.by.finding$ClinSig <- scaled.tra.by.finding$ClinSig == 2

### *X-Validation* Data by finding - T2-tra and ADC 
# 10 fold X-Validation
xv <- makeResampleDesc("CV", iters = 10)

class.task1 <- makeClassifTask(id = "tra.ProstateX", data = scaled.tra.by.finding, target = "ClinSig", positive = "TRUE")

t.weight <- nrow(scaled.tra.by.finding)/(length(unique(scaled.tra.by.finding$ClinSig)) * nrow(filter(scaled.tra.by.finding, ClinSig == TRUE)))
f.weight <- nrow(scaled.tra.by.finding)/(length(unique(scaled.tra.by.finding$ClinSig)) * nrow(filter(scaled.tra.by.finding, ClinSig == FALSE)))
  
c.weight <- c(t.weight, f.weight)
names(c.weight) <- c("TRUE", "FALSE")

class.lrn <- makeLearner("classif.ksvm", class.weights = c.weight, predict.type = 'prob')

### Benchmarking - SVM parameter testing. Make the SVM great again.

discrete_ps = makeParamSet(
  makeDiscreteParam("C", values = c(0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 500, 5000, 7500)),
  makeDiscreteParam("sigma", values = c(10^-1, 10^-2, 0.015, 0.02, 10^-3, 10^-4, 10^-5))
)
ctrl = makeTuneControlGrid()

res1 = tuneParams(class.lrn, task = class.task1, resampling = xv, par.set = discrete_ps,
                 control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


class.lrn.lr <- makeLearner("classif.logreg", predict.type = "prob")
res1.lr <- crossval(class.lrn.lr, task = class.task1, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)



### 2. T2tra, ADCtra, zone ----
tra.zone <- cbind(t2tra.df, tra.track.zone)
#adc.zone <- cbind(adc_transeq.df, adc_trans.track)

tra.by.finding.zone <- tra.zone %>%
  inner_join(adc, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ProxID", "fid")) %>%  #, "AStra", "PZtra", "SVtra", "TZtra")) %>%#remove after being joined by
  rename(ClinSig = ClinSig.y) %>%
  mutate(ClinSig = as.factor(ClinSig)) %>%
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         adc1 = "1.y",  adc2 = "2.y", adc3 = "3.y", adc4 = "4.y", adc5 = "5.y", adc6 = "6.y", adc7 = "7.y", adc8 = "8.y", 
         adc9 = "9.y", adc10 = "10.y", adc11 = "11.y", adc12 = "12.y", adc13 = "13.y", adc14 = "14.y", adc15 = "15.y", adc16 = "16.y", 
         adc17 = "17.y", adc18 = "18.y", adc19 = "19.y", adc20 = "20.y", adc21 = "21.y", adc22 = "22.y", adc23 = "23.y", adc24 = "24.y", adc25 = "25.y",
         AS = "AStra", PZ = "PZtra", SV = "SVtra", TZ = "TZtra") %>%
  select(-AS, -PZ, -SV, -TZ, everything())

scaled.tra.by.finding.zone <- scale(tra.by.finding.zone[,1:50])

scaled.tra.by.finding.zone <- as.data.frame(cbind(scaled.tra.by.finding.zone, as.factor(tra.by.finding.zone[,51]), tra.by.finding.zone[,52:55]))
colnames(scaled.tra.by.finding.zone)[51] <- "ClinSig"
#scaled.tra.by.finding$ClinSig <- scaled.tra.by.finding$ClinSig == 2


# k-fold x-validation
class.task2 <- makeClassifTask(id = "tra.zone.ProstateX", data = scaled.tra.by.finding.zone, target = "ClinSig", positive = "TRUE")
#class weights are the same. Therefore, use same class weights and class.lrn (same learner)

res2 <- tuneParams(class.lrn, task = class.task2, resampling = xv, par.set = discrete_ps,
                       control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res2.lr <- crossval(class.lrn.lr, task = class.task2, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

# Generate a learner based on this optimal model
lrn2 <- setHyperPars(makeLearner("classif.ksvm", par.vals = res2$x, predict.type = "prob"))
m2 <- train(lrn2, class.task2)
m_p2 <- predict(m2, task = class.task2)

#table(m_p2$data$response, scaled.tra.by.finding.zone[,51])

# ROC curve
lrn2.roc <- generateThreshVsPerfData(m_p2, measures = list(fpr, tpr, mmce))
plotROCCurves(lrn2.roc)

# save this model
saveRDS(lrn2, "model2.RDS")
saveRDS(m2, "m2.RDS")

### 3. T2tra, ADCtra, zone, age, weight ----
tra.zone.ageweight <- cbind(t2tra.df, tra.track.age.weight)

tra.by.finding.zone.ageweight <- tra.zone.ageweight %>%
  inner_join(adc, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ProxID", "fid")) %>%  #, "AStra", "PZtra", "SVtra", "TZtra")) %>%#remove after being joined by
  rename(ClinSig = ClinSig.y) %>%
  mutate(ClinSig = as.factor(ClinSig)) %>%
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         adc1 = "1.y",  adc2 = "2.y", adc3 = "3.y", adc4 = "4.y", adc5 = "5.y", adc6 = "6.y", adc7 = "7.y", adc8 = "8.y", 
         adc9 = "9.y", adc10 = "10.y", adc11 = "11.y", adc12 = "12.y", adc13 = "13.y", adc14 = "14.y", adc15 = "15.y", adc16 = "16.y", 
         adc17 = "17.y", adc18 = "18.y", adc19 = "19.y", adc20 = "20.y", adc21 = "21.y", adc22 = "22.y", adc23 = "23.y", adc24 = "24.y", adc25 = "25.y",
         AS = "AStra", PZ = "PZtra", SV = "SVtra", TZ = "TZtra") %>%
  select(-AS, -PZ, -SV, -TZ, everything())

scaled.tra.by.finding.zone.ageweight <- scale(tra.by.finding.zone.ageweight[,1:52])

scaled.tra.by.finding.zone.ageweight <- as.data.frame(cbind(scaled.tra.by.finding.zone.ageweight, as.factor(tra.by.finding.zone.ageweight[,53]), tra.by.finding.zone.ageweight[,54:57]))
colnames(scaled.tra.by.finding.zone.ageweight)[53] <- "ClinSig"
#scaled.tra.by.finding$ClinSig <- scaled.tra.by.finding$ClinSig == 2

# k-fold x-validation
class.task3 <- makeClassifTask(id = "tra.zone.ageweight.ProstateX", data = scaled.tra.by.finding.zone.ageweight, target = "ClinSig", positive = "TRUE")
#class weights are the same. Therefore, use same class weights and class.lrn (same learner)

res3 <- tuneParams(class.lrn, task = class.task3, resampling = xv, par.set = discrete_ps,
                       control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res3.lr <- crossval(class.lrn.lr, task = class.task3, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


# Generate a learner based on this optimal model
lrn3 <- setHyperPars(makeLearner("classif.ksvm", par.vals = res3$x, predict.type = "prob"))
m3 <- train(lrn3, class.task3)
m_p3 <- predict(m3, task = class.task3)

# ROC curve
lrn3.roc <- generateThreshVsPerfData(m_p3, measures = list(fpr, tpr, mmce))
plotROCCurves(lrn3.roc)



# 4. T2sag, T2cor, T2tra ----
# T2sag, T2cor, T2tra 
### Data
sag <- cbind(t2sag.eq, sag.track)
cor <- cbind(t2cor.df, cor.track)
tra <- cbind(t2tra.df, tra.track)
data.by.finding <- sag %>%
  inner_join(cor, by = c("ProxID", "fid")) %>%
  inner_join(tra, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ClinSig.y", "ProxID", "fid")) %>% #remove after being joined by
  rename(sag1 = "1.x", sag2 = "2.x",sag3 = "3.x",sag4 = "4.x", sag5 = "5.x", sag6 = "6.x", sag7 = "7.x", sag8 = "8.x",
         sag9 = "9.x", sag10 = "10.x", sag11 = "11.x", sag12 = "12.x", sag13 = "13.x", sag14 = "14.x", sag15 = "15.x", sag16 = "16.x", 
         sag17 = "17.x", sag18 = "18.x", sag19 = "19.x", sag20 = "20.x", sag21 = "21.x", sag22 = "22.x", sag23 = "23.x", sag24 = "24.x", sag25 = "25.x",
         cor1 = "1.y",  cor2 = "2.y", cor3 = "3.y", cor4 = "4.y", cor5 = "5.y", cor6 = "6.y", cor7 = "7.y", cor8 = "8.y", 
         cor9 = "9.y", cor10 = "10.y", cor11 = "11.y", cor12 = "12.y", cor13 = "13.y", cor14 = "14.y", cor15 = "15.y", cor16 = "16.y", 
         cor17 = "17.y", cor18 = "18.y", cor19 = "19.y", cor20 = "20.y", cor21 = "21.y", cor22 = "22.y", cor23 = "23.y", cor24 = "24.y", cor25 = "25.y",
         tra1 = "1", tra2 = "2",tra3 = "3", tra4 = "4", tra5 = "5", tra6 = "6", tra7 = "7", tra8 = "8",
         tra9 = "9", tra10 = "10", tra11 = "11", tra12 = "12", tra13 = "13", tra14 = "14", tra15 = "15", tra16 = "16", 
         tra17 = "17", tra18 = "18", tra19 = "19", tra20 = "20", tra21 = "21", tra22 = "22", tra23 = "23", tra24 = "24", tra25 = "25")

#### Normalise the data - each input dimension: mean subtracted, divide by std dev.
scaled.data.by.finding <- scale(data.by.finding[,1:75])

scaled.data.by.finding <- as.data.frame(cbind(scaled.data.by.finding, data.by.finding[,76]))
colnames(scaled.data.by.finding)[76] <- "ClinSig"
scaled.data.by.finding$ClinSig <- scaled.data.by.finding$ClinSig == 1

# k-fold x validation
class.task4 <- makeClassifTask(id = "tra.zone.ageweight.ProstateX", data = scaled.data.by.finding, target = "ClinSig", positive = "TRUE")

## new class weights
t.weight.t2 <- nrow(scaled.data.by.finding)/(length(unique(scaled.data.by.finding$ClinSig)) * nrow(filter(scaled.data.by.finding, ClinSig == TRUE)))
f.weight.t2 <- nrow(scaled.data.by.finding)/(length(unique(scaled.data.by.finding$ClinSig)) * nrow(filter(scaled.data.by.finding, ClinSig == FALSE)))

c.weight.t2 <- c(t.weight.t2, f.weight.t2)
names(c.weight.t2) <- c("TRUE", "FALSE")

class.lrn.t2 <- makeLearner("classif.ksvm", class.weights = c.weight.t2, predict.type = 'prob')

res4 <- tuneParams(class.lrn.t2, task = class.task4, resampling = xv, par.set = discrete_ps,
                       control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res4.lr <- crossval(class.lrn.lr, task = class.task4, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


# 5. T2sag, T2cor, T2tra, zone ----
# use tra.track.zone

t2.zone.by.finding <- tra.zone %>%
  inner_join(sag, by = c("ProxID", "fid")) %>%
  inner_join(cor, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ClinSig.y", "ProxID", "fid")) %>% #remove after being joined by
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         sag1 = "1.y",  sag2 = "2.y", sag3 = "3.y", sag4 = "4.y", sag5 = "5.y", sag6 = "6.y", sag7 = "7.y", sag8 = "8.y", 
         sag9 = "9.y", sag10 = "10.y", sag11 = "11.y", sag12 = "12.y", sag13 = "13.y", sag14 = "14.y", sag15 = "15.y", sag16 = "16.y", 
         sag17 = "17.y", sag18 = "18.y", sag19 = "19.y", sag20 = "20.y", sag21 = "21.y", sag22 = "22.y", sag23 = "23.y", sag24 = "24.y", sag25 = "25.y",
         cor1 = "1", cor2 = "2",cor3 = "3", cor4 = "4", cor5 = "5", cor6 = "6", cor7 = "7", cor8 = "8",
         cor9 = "9", cor10 = "10", cor11 = "11", cor12 = "12", cor13 = "13", cor14 = "14", cor15 = "15", cor16 = "16", 
         cor17 = "17", cor18 = "18", cor19 = "19", cor20 = "20", cor21 = "21", cor22 = "22", cor23 = "23", cor24 = "24", cor25 = "25",
         AS = "AStra", PZ = "PZtra", SV = "SVtra", TZ = "TZtra") %>%
  select(-AS, -PZ, -SV, -TZ, everything())

scaled.t2.zone.by.finding <- scale(t2.zone.by.finding[,1:75])

scaled.t2.zone.by.finding <- as.data.frame(cbind(scaled.t2.zone.by.finding, t2.zone.by.finding[,76:80]))
colnames(scaled.data.by.finding)[76] <- "ClinSig"
#scaled.data.by.finding$ClinSig <- scaled.data.by.finding$ClinSig == 1

# k-fold x-validation
class.task5 <- makeClassifTask(id = "t2.zone.ProstateX", data = scaled.t2.zone.by.finding, target = "ClinSig", positive = "TRUE")
# same class weights as t2sag, t2cor, t2tra. Use class.lrn.t2
res5 <- tuneParams(class.lrn.t2, task = class.task5, resampling = xv, par.set = discrete_ps,
                     control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res5.lr <- crossval(class.lrn.lr, task = class.task5, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


# 6. T2sag, T2cor, T2tra, zone, age, weight ----
t2.zone.ageweight.by.finding <- tra.zone.ageweight %>%
  inner_join(sag, by = c("ProxID", "fid")) %>%
  inner_join(cor, by = c("ProxID", "fid")) %>%
  select(-one_of("ClinSig.x", "ClinSig.y", "ProxID", "fid")) %>% #remove after being joined by
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         sag1 = "1.y",  sag2 = "2.y", sag3 = "3.y", sag4 = "4.y", sag5 = "5.y", sag6 = "6.y", sag7 = "7.y", sag8 = "8.y", 
         sag9 = "9.y", sag10 = "10.y", sag11 = "11.y", sag12 = "12.y", sag13 = "13.y", sag14 = "14.y", sag15 = "15.y", sag16 = "16.y", 
         sag17 = "17.y", sag18 = "18.y", sag19 = "19.y", sag20 = "20.y", sag21 = "21.y", sag22 = "22.y", sag23 = "23.y", sag24 = "24.y", sag25 = "25.y",
         cor1 = "1", cor2 = "2",cor3 = "3", cor4 = "4", cor5 = "5", cor6 = "6", cor7 = "7", cor8 = "8",
         cor9 = "9", cor10 = "10", cor11 = "11", cor12 = "12", cor13 = "13", cor14 = "14", cor15 = "15", cor16 = "16", 
         cor17 = "17", cor18 = "18", cor19 = "19", cor20 = "20", cor21 = "21", cor22 = "22", cor23 = "23", cor24 = "24", cor25 = "25",
         AS = "AStra", PZ = "PZtra", SV = "SVtra", TZ = "TZtra") %>%
  select(-AS, -PZ, -SV, -TZ, everything())

scaled.t2.zone.ageweight.by.finding <- scale(t2.zone.ageweight.by.finding[,1:77])

scaled.t2.zone.ageweight.by.finding <- as.data.frame(cbind(scaled.t2.zone.ageweight.by.finding, t2.zone.ageweight.by.finding[,78:82]))
colnames(scaled.t2.zone.ageweight.by.finding)[78] <- "ClinSig"
#scaled.data.by.finding$ClinSig <- scaled.data.by.finding$ClinSig == 1

# k-fold x-validation
class.task6 <- makeClassifTask(id = "t2.zone.ageweight.ProstateX", data = scaled.t2.zone.ageweight.by.finding, target = "ClinSig", positive = "TRUE")
# same class weights as t2sag, t2cor, t2tra. Use class.lrn.t2
res6 <- tuneParams(class.lrn.t2, task = class.task6, resampling = xv, par.set = discrete_ps,
                          control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res6.lr <- crossval(class.lrn.lr, task = class.task6, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


# 7. T2sag, T2cor, T2tra, ADCtra, zone, age, weight ----
all.zone.ageweight.by.finding <- tra.zone.ageweight %>%
  inner_join(sag, by = c("ProxID", "fid")) %>%
  inner_join(cor, by = c("ProxID", "fid")) %>%
  inner_join(adc, by = c("ProxID", "fid")) %>%
  rename(tra1 = "1.x", tra2 = "2.x",tra3 = "3.x", tra4 = "4.x", tra5 = "5.x", tra6 = "6.x", tra7 = "7.x", tra8 = "8.x",
         tra9 = "9.x", tra10 = "10.x", tra11 = "11.x", tra12 = "12.x", tra13 = "13.x", tra14 = "14.x", tra15 = "15.x", tra16 = "16.x", 
         tra17 = "17.x", tra18 = "18.x", tra19 = "19.x", tra20 = "20.x", tra21 = "21.x", tra22 = "22.x", tra23 = "23.x", tra24 = "24.x", tra25 = "25.x",
         sag1 = "1.y",  sag2 = "2.y", sag3 = "3.y", sag4 = "4.y", sag5 = "5.y", sag6 = "6.y", sag7 = "7.y", sag8 = "8.y", 
         sag9 = "9.y", sag10 = "10.y", sag11 = "11.y", sag12 = "12.y", sag13 = "13.y", sag14 = "14.y", sag15 = "15.y", sag16 = "16.y", 
         sag17 = "17.y", sag18 = "18.y", sag19 = "19.y", sag20 = "20.y", sag21 = "21.y", sag22 = "22.y", sag23 = "23.y", sag24 = "24.y", sag25 = "25.y",
         cor1 = "1.x.x", cor2 = "2.x.x",cor3 = "3.x.x", cor4 = "4.x.x", cor5 = "5.x.x", cor6 = "6.x.x", cor7 = "7.x.x", cor8 = "8.x.x",
         cor9 = "9.x.x", cor10 = "10.x.x", cor11 = "11.x.x", cor12 = "12.x.x", cor13 = "13.x.x", cor14 = "14.x.x", cor15 = "15.x.x", cor16 = "16.x.x", 
         cor17 = "17.x.x", cor18 = "18.x.x", cor19 = "19.x.x", cor20 = "20.x.x", cor21 = "21.x.x", cor22 = "22.x.x", cor23 = "23.x.x", cor24 = "24.x.x", cor25 = "25.x.x",
         adc1 = "1.y.y",  adc2 = "2.y.y", adc3 = "3.y.y", adc4 = "4.y.y", adc5 = "5.y.y", adc6 = "6.y.y", adc7 = "7.y.y", adc8 = "8.y.y", 
         adc9 = "9.y.y", adc10 = "10.y.y", adc11 = "11.y.y", adc12 = "12.y.y", adc13 = "13.y.y", adc14 = "14.y.y", adc15 = "15.y.y", adc16 = "16.y.y", 
         adc17 = "17.y.y", adc18 = "18.y.y", adc19 = "19.y.y", adc20 = "20.y.y", adc21 = "21.y.y", adc22 = "22.y.y", adc23 = "23.y.y", adc24 = "24.y.y", adc25 = "25.y.y",
         ClinSig = "ClinSig.y.y", 
         AS = "AStra", PZ = "PZtra", SV = "SVtra", TZ = "TZtra") %>%
  select(-one_of("ClinSig.x", "ClinSig.y", "ClinSig.x.x", "ProxID", "fid")) %>% #remove after being joined by
  select(-AS, -PZ, -SV, -TZ, -ClinSig, everything())
  
scaled.all <- scale(all.zone.ageweight.by.finding[,1:102])
scaled.all <- as.data.frame(cbind(scaled.all, all.zone.ageweight.by.finding[,103:107]))

# k-fold x-validation
class.task7 <- makeClassifTask(id = "all.ProstateX", data = scaled.all, target = "ClinSig", positive = "TRUE")
#same class weights as class.lrn

res7 <- tuneParams(class.lrn, task = class.task7, resampling = xv, par.set = discrete_ps,
                                    control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

res7.lr <- crossval(class.lrn.lr, task = class.task7, iters = 10, stratify = TRUE, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)


# CHECK - zone, age, weight ----
ageweight.zone.check <- tra.track.age.weight[,3:9]
#ageweight.zone.check$ClinSig <- as.factor(ageweight.zone.check$ClinSig)
check <- as.data.frame(ageweight.zone.check)
class.task.check <- makeClassifTask(id = 'check', data = check, target = "ClinSig", positive = "TRUE")

## new class weights
t.weight.check <- nrow(check)/(length(unique(check$ClinSig)) * nrow(filter(check, ClinSig == TRUE)))
f.weight.check <- nrow(check)/(length(unique(check$ClinSig)) * nrow(filter(check, ClinSig == FALSE)))

c.weight.check <- c(t.weight.check, f.weight.check)
names(c.weight.check) <- c("TRUE", "FALSE")

## Learner
class.lrn.check <- makeLearner("classif.ksvm", class.weights = c.weight.check, predict.type = 'prob')

res.check <- tuneParams(class.lrn.check, task = class.task.check, resampling = xv, par.set = discrete_ps,
                     control = ctrl, measures = list(auc, setAggregation(auc, test.sd)), show.info = TRUE)

# Generate a learner based on this optimal model
check.lrn <- setHyperPars(makeLearner("classif.ksvm", par.vals = res.check$x))
m <- train(check.lrn, class.task.check)
m_p <- predict(m, task = class.task.check)
table(m_p$data$response, check$ClinSig)
