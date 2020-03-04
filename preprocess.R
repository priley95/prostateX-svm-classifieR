# Preprocessing script
# Using metadata, take the correct reference to the regions of interest for lesions and obtaining healthy tissue, to be used for source extraction.

# Load in packages
library(tidyverse)
library(oro.dicom)
library(IM) #histeq
library(mlr) #createDummyFeatures

# Metadata: .csv files for leison information. Join relevant metadata ----
findings <- read_csv('ProstateX-Findings-Train.csv')
images <- read_csv('ProstateX-Images-Train.csv') # Image metadata

# Join metadata from findings and images
joined <- images %>%
  select(ProxID:ijk, Dim) %>% 
  left_join(findings, by = "pos") %>%
  select(-one_of('WorldMatrix', 'ProxID.y', 'fid.y')) %>%
  rename(ProxID = ProxID.x, fid = fid.x) %>%
  separate(ijk, c('i_column', 'j_row', 'k_slice')) %>%
  mutate(i_column = as.numeric(i_column), j_row = as.numeric(j_row), k_slice = as.numeric(k_slice))

joined.by.modality <- split(joined, joined$Name) # 140 modalities?

# Seperate tibbles for the three t2- modalities to be studied, for joined data
t2.tse.sag0.data <- joined.by.modality$t2_tse_sag0 #358
t2.tse.sag0.data <- distinct(t2.tse.sag0.data) #343
t2.tse.sag0.data <- t2.tse.sag0.data[!duplicated(t2.tse.sag0.data[c("ProxID", "fid")]),] #330

t2.tse.cor0.data <- joined.by.modality$t2_tse_cor0 #length 339 findings with this modality.
t2.tse.cor0.data <- distinct(t2.tse.cor0.data) #338
t2.tse.cor0.data <- t2.tse.cor0.data[!duplicated(t2.tse.cor0.data[c("ProxID", "fid")]),] #330

t2.tse.tra0.data <- joined.by.modality$t2_tse_tra0 #412 -> Px-0191 has a t2_tse_tra under a different name. Ignore.
t2.tse.tra0.data <- distinct(t2.tse.tra0.data) #359
t2.tse.tra0.data <- t2.tse.tra0.data[!duplicated(t2.tse.tra0.data[c("ProxID", "fid")]),] #329
t2.tse.tra0.data <- t2.tse.tra0.data[!grepl("640x640", t2.tse.tra0.data$Dim), ] # caused problems - some have been removed (about 40)
#329

adc_trans.data <- joined.by.modality$ep2d_diff_tra_DYNDIST_ADC0 #265
adc_trans.data <- distinct(adc_trans.data) #265
#t2.tse.tra0.data <- t2.tse.tra0.data[!duplicated(t2.tse.tra0.data[c("ProxID", "fid")]),] #329

diff.data <- joined.by.modality$ep2d_diff_tra_DYNDIST0 #265
diff.data <- distinct(diff.data)

# Obtain DICOM data - header and image ---------
patientdir <- 'xyz/ProstateX data/PROSTATEx'

patients <- list.files(path = patientdir, recursive = T, pattern = ".dcm", full.names = T)

readDICOMfiles <- function(f){
  dcm.file <- readDICOMFile(f)
}

# DICOM data: header and image info
dicom.files <- map(patients, readDICOMfiles) 
saveRDS(dicom.files, "dicom_files.rds") # Save to file
dicom.files <- readRDS("dicom_files.rds") # Load back in

# Equalize the images. -----------
files.histeq <- list()

for (b in 1:length(dicom.files)){
  files.histeq[[b]] <- histeq(dicom.files[[b]]$img)
}
saveRDS(files.histeq, "files_histeq.RDS")
files.histeq <- readRDS("files_histeq.RDS")

# Create 'full' - containing header data, img data and equalized img data. ----
full <- list()
for (i in 1:length(dicom.files)){
  full[[i]] <- list(dicom.files[[i]]$hdr, dicom.files[[i]]$img, files.histeq[[i]])
  names(full[[i]]) <- c('hdr', 'img', 'eq')
}

saveRDS(full, "full.RDS")
full <- readRDS("full.RDS")

# Seperate dicom.files by different modalities --------
t2sag_index <- matrix(); t2cor_index <- matrix(); t2tra_index <- matrix(); adc_trans_index <- matrix(); diff_index <- matrix()
for (i in 1:length(full)){
t2sag_index[i] <- full[[i]]$hdr$value[26] == 't2_tse_sag' #indexed into list, says which are the value asked for!!
t2cor_index[i] <- full[[i]]$hdr$value[26] == 't2_tse_cor'
t2tra_index[i] <- full[[i]]$hdr$value[26] == 't2_tse_tra'
adc_trans_index[i] <- full[[i]]$hdr$value[26] == 'ep2d_diff_tra_DYNDIST_ADC' 
diff_index[i] <- full[[i]]$hdr$value[26] == 'ep2d_diff_tra_DYNDIST'
}
t2sag <- full[t2sag_index]; t2cor <- full[t2cor_index]; t2tra <- full[t2tra_index]; adc_trans <- full[adc_trans_index]; diff <- full[diff_index]

### note: InstanceNumber - 1 = k_slice ###

# Attain ROI and equalized ROI ----

# t2sag
a <- list()

for (m in (1:length(t2sag))){ # img data here
  for (n in (1:nrow(t2.tse.sag0.data))){ # reference to lesion here
    
    if ((t2.tse.sag0.data[n,]$ProxID == extractHeader(t2sag[[m]]$hdr, "PatientsName", numeric = FALSE)) 
        & (t2.tse.sag0.data[n,]$k_slice == (extractHeader(t2sag[[m]]$hdr, "InstanceNumber") - 1))
        & (substr(t2.tse.sag0.data[n,]$Name,1,10) == (extractHeader(t2sag[[m]]$hdr, "ProtocolName", numeric = FALSE)))){
      a[[n]] <- t2sag[[m]]
    }
  }
}

saveRDS(a, "a_sag_ref.RDS")
a <- readRDS("a_sag_ref.RDS")

sag.roi <- list(); sageq.roi <- list()#; age.sag <- matrix()

for (k in (1:length(a))){
    sag.roi[[k]] <- a[[k]]$img[((t2.tse.sag0.data[k,]$i_column) - 2) : ((t2.tse.sag0.data[k,]$i_column) + 2), ((t2.tse.sag0.data[k,]$j_row) - 2) : ((t2.tse.sag0.data[k,]$j_row) + 2)] 
    sageq.roi[[k]] <- a[[k]]$eq[((t2.tse.sag0.data[k,]$i_column) - 2) : ((t2.tse.sag0.data[k,]$i_column) + 2), ((t2.tse.sag0.data[k,]$j_row) - 2) : ((t2.tse.sag0.data[k,]$j_row) + 2)] 
    #age.sag[[k]] <- extractHeader(a[[k]]$hdr, "PatientsAge", numeric = FALSE)
} 


# t2sag track information
sag.track <- t2.tse.sag0.data[,c("ProxID", "fid", "ClinSig")]
zone.sag.dummy <- createDummyFeatures(t2.tse.sag0.data$zone)
sag.track.zone <- cbind(sag.track, zone.sag.dummy)

# Extract patient age and weight ----
age_sag <- NULL; weight_sag <- NULL; height_sag <- NULL; ProxID_sag <- NULL

for (k in 1:length(a)){
  age_sag[k] <- as.numeric(substr(extractHeader(a[[k]]$hdr, "PatientsAge", numeric= FALSE), 2, 3))
  weight_sag[k] <- extractHeader(a[[k]]$hdr, "PatientsWeight", numeric = TRUE)
  height_sag[k] <- extractHeader(a[[k]]$hdr, "PatientsSize", numeric = TRUE)
  ProxID_sag[k] <- extractHeader(a[[k]]$hdr, "PatientsName", numeric = FALSE)
}

patient_sag_ <- data.frame(ProxID_sag, age_sag, weight_sag, height_sag)
patient_sag <- subset(patient_sag_, !duplicated(patient_sag_$ProxID))

patient_sag <- patient_sag %>%
  rename(ProxID = ProxID_sag)
# add to a track data
sag.track.patient <- sag.track %>%
  left_join(patient_sag, by = "ProxID") %>%
  left_join(sag.track.zone, by = c("ProxID", "fid", "ClinSig"))



# Create sag - containing image data, ProxID, fid, age, weight, height, ClinSig, zone information
t2sag.eq <- NULL
for (i in 1:length(sageq.roi)){
  t2sag.eq[[i]] <- as.vector(sageq.roi[[i]])
}
t2sag.eq <- as.matrix(do.call(rbind, t2sag.eq))

sag <- cbind(t2sag.eq, sag.track.patient)

saveRDS(sag, "sag_data.RDS")

write.csv(sag, "t2sag_lesion.csv")


# cor
b <- list()
for (ma in (1:length(t2cor))){ #img data here
  for (na in (1:nrow(t2.tse.cor0.data))){ #reference to lesion here
    
    if ((t2.tse.cor0.data[na,]$ProxID == extractHeader(t2cor[[ma]]$hdr, "PatientsName", numeric = FALSE)) 
        & (t2.tse.cor0.data[na,]$k_slice == (extractHeader(t2cor[[ma]]$hdr, "InstanceNumber") - 1))
        & (substr(t2.tse.cor0.data[na,]$Name,1,10) == (extractHeader(t2cor[[ma]]$hdr, "ProtocolName", numeric = FALSE)))){
      b[[na]] <- t2cor[[ma]]
    }
  }
}
saveRDS(b, "b_cor_ref.RDS")
b <- readRDS("b_cor_ref.RDS")

cor.roi <- list(); coreq.roi <- list()
for (ka in (1:length(b))){
    cor.roi[[ka]] <- a[[ka]]$img[((t2.tse.cor0.data[ka,]$i_column) - 2) : ((t2.tse.cor0.data[ka,]$i_column) + 2), ((t2.tse.cor0.data[ka,]$j_row) - 2) : ((t2.tse.cor0.data[ka,]$j_row) + 2)] 
    coreq.roi[[ka]] <- a[[ka]]$eq[((t2.tse.cor0.data[ka,]$i_column) - 2) : ((t2.tse.cor0.data[ka,]$i_column) + 2), ((t2.tse.cor0.data[ka,]$j_row) - 2) : ((t2.tse.cor0.data[ka,]$j_row) + 2)] 
}

# t2cor track information
cor.track <- t2.tse.cor0.data[,c("ProxID", "fid", "ClinSig")]
zone.cor.dummy <- createDummyFeatures(t2.tse.cor0.data$zone)
cor.track.zone <- cbind(cor.track, zone.cor.dummy)

# Add to a track data
cor.track.patient <- cor.track %>%
  left_join(patient_sag, by = "ProxID") %>%
  left_join(cor.track.zone, by = c("ProxID", "fid", "ClinSig"))

# Create cor - containing image data, ProxID, fid, age, weight, height, ClinSig, zone information
t2cor.df <- NULL
for (i in 1:length(coreq.roi)){
  t2cor.df[[i]] <- as.vector(coreq.roi[[i]])
}
t2cor.df <- as.matrix(do.call(rbind, t2cor.df))
cor <- cbind(t2cor.df, cor.track.patient)

saveRDS(cor, "cor_data.RDS")
write.csv(cor, "t2cor_lesion.csv")

# tra
c <- list()
for (mb in (1:length(t2tra))){ #img data here
  for (nb in (1:nrow(t2.tse.tra0.data))){ #reference to lesion here
    
    if ((t2.tse.tra0.data[nb,]$ProxID == extractHeader(t2tra[[mb]]$hdr, "PatientsName", numeric = FALSE)) 
        & (t2.tse.tra0.data[nb,]$k_slice == (extractHeader(t2tra[[mb]]$hdr, "InstanceNumber") - 1))
        & (substr(t2.tse.tra0.data[nb,]$Name,1,10) == (extractHeader(t2tra[[mb]]$hdr, "ProtocolName", numeric = FALSE)))){
      c[[nb]] <- t2tra[[mb]]
    }
  }
}

# Some slices missing - remove. k_slice is higher than slice number available.
c[36] <- NULL #36
c[229] <- NULL #230 (-1 when removing 36)
c[229] <- NULL #231 (-2 after removing other two)

# Remove from metadata
t2.tse.tra0.data <- t2.tse.tra0.data[-c(36, 230, 231),]

saveRDS(c, "c_tra_ref.RDS")
c <- readRDS("c_tra_ref.RDS")


tra.roi <- list(); traeq.roi <- list()
for (kb in (1:length(c))){
    tra.roi[[kb]] <- c[[kb]]$img[((t2.tse.tra0.data[kb,]$i_column) - 2) : ((t2.tse.tra0.data[kb,]$i_column) + 2), ((t2.tse.tra0.data[kb,]$j_row) - 2) : ((t2.tse.tra0.data[kb,]$j_row) + 2)] 
    traeq.roi[[kb]] <- c[[kb]]$eq[((t2.tse.tra0.data[kb,]$i_column) - 2) : ((t2.tse.tra0.data[kb,]$i_column) + 2), ((t2.tse.tra0.data[kb,]$j_row) - 2) : ((t2.tse.tra0.data[kb,]$j_row) + 2)] 
}

tra.track <- t2.tse.tra0.data[,c("ProxID", "fid", "ClinSig")]
zone.tra.dummy <- createDummyFeatures(t2.tse.tra0.data$zone)
colnames(zone.tra.dummy) <- c("AStra", "PZtra", "SVtra", "TZtra")
tra.track.zone <- cbind(tra.track, zone.tra.dummy)

# add to a track data
tra.track.patient <- tra.track %>%
  left_join(patient_sag, by = "ProxID") %>%
  left_join(tra.track.zone, by = c("ProxID", "fid", "ClinSig"))

t2tra.df <- NULL
for (i in 1:length(traeq.roi)){
  t2tra.df[[i]] <- as.vector(traeq.roi[[i]])
}
t2tra.df <- as.matrix(do.call(rbind, t2tra.df))

tra <- cbind(t2tra.df, tra.track.patient)

saveRDS(tra, "tra_data.RDS")
write.csv(tra, "t2tra_lesion.csv")

# ADC trans files
adc_index <- list()
for (f in (1:length(adc_trans))){ #img data here
  for (g in (1:nrow(adc_trans.data))){ #reference to lesion here
    
    if ((adc_trans.data[g,]$ProxID == extractHeader(adc_trans[[f]]$hdr, "PatientsName", numeric = FALSE)) 
        & (adc_trans.data[g,]$k_slice == (extractHeader(adc_trans[[f]]$hdr, "InstanceNumber") - 1))
        & (substr(adc_trans.data[g,]$Name,1,13) == (extractHeader(adc_trans[[f]]$hdr, "ProtocolName", numeric = FALSE)))){
      adc_index[[g]] <- adc_trans[[f]]
    }
  }
}
saveRDS(adc_index, "adc_ref.RDS")
adc_index <- readRDS("adc_ref.RDS")

# remove out of bounds
adc_index[c(27,200)] <- NULL
adc_index[32] <- NULL
adc_index[c(25:32)] <- NULL #Px-0025
# remove from metadata
adc_trans.data <- adc_trans.data[-c(27,200),]
adc_trans.data <- adc_trans.data[-c(32),]
adc_trans.data <- adc_trans.data[-c(25:32),]

adc_trans.roi <- list(); adc_transeq.roi <- list()
for (ga in (1:length(adc_index))){
  adc_trans.roi[[ga]] <- adc_index[[ga]]$img[((adc_trans.data[ga,]$i_column) - 2) : ((adc_trans.data[ga,]$i_column) + 2), ((adc_trans.data[ga,]$j_row) - 2) : ((adc_trans.data[ga,]$j_row) + 2)] 
  adc_transeq.roi[[ga]] <- adc_index[[ga]]$eq[((adc_trans.data[ga,]$i_column) - 2) : ((adc_trans.data[ga,]$i_column) + 2), ((adc_trans.data[ga,]$j_row) - 2) : ((adc_trans.data[ga,]$j_row) + 2)] 
} 

adc_trans.track <- adc_trans.data[,c("ProxID", "fid", "ClinSig")]
zone.adc.dummy <- createDummyFeatures(adc_trans.data$zone)
colnames(zone.adc.dummy) <- c("ASadc", "PZadc", "SVadc", "TZadc")
adc_trans.track.zone <- cbind(adc_trans.track, zone.adc.dummy)
adc_trans.track.zone <- distinct(adc_trans.track.zone)

# add to a track data
adc.track.patient <- adc_trans.track %>%
  left_join(patient_sag, by = "ProxID") %>%
  left_join(adc_trans.track.zone, by = c("ProxID", "fid", "ClinSig"))

adc_trans.df <- NULL
for (i in 1:length(adc_transeq.roi)){
  adc_trans.df[[i]] <- as.vector(adc_transeq.roi[[i]])
}
adc_transeq.df <- as.matrix(do.call(rbind, adc_trans.df))
adc <- cbind(adc_transeq.df, adc.track.patient)

saveRDS(adc, "adc_data.RDS")
write.csv(adc, "adc_lesion.csv")


# diff
diff_index <- list()
for (h in (1:length(diff))){ #img data
  for(j in (1:nrow(diff.data))){ # ref to lesion
    
    if (((diff.data[j,]$ProxID == extractHeader(diff[[h]]$hdr, "PatientsName", numeric = FALSE))
         & (diff.data[j,]$k_slice == extractHeader(diff[[h]]$hdr, "InstanceNumber") - 1))
        & (substr(diff.data[j,]$Name, 1, 21) == (extractHeader(diff[[h]]$hdr, "SeriesDescription", numeric = FALSE)))){
      diff_index[[j]] <- diff[[h]]
    }
  }
}

# no data
#diff_index[33] <- NULL
#diff.data <- diff.data[-c(33),]

saveRDS(diff_index, "diff_ref.RDS")
diff_index <- readRDS("diff_ref.RDS")

# out of bounds
#diff_index[33] <- NULL
#diff.data <- diff.data[-c(33),]

diff.roi <- list(); diffeq.roi <- list()
for (ha in (1:length(diff_index))){
  diff.roi[[ha]] <- diff_index[[ha]]$img[((diff.data[ha,]$i_column) - 2) : ((diff.data[ha,]$i_column) + 2), ((diff.data[ha,]$j_row) - 2) : ((diff.data[ha,]$j_row) + 2)] 
  diffeq.roi[[ha]] <- diff_index[[ha]]$eq[((diff.data[ha,]$i_column) - 2) : ((diff.data[ha,]$i_column) + 2), ((diff.data[ha,]$j_row) - 2) : ((diff.data[ha,]$j_row) + 2)] 
}

# diff track information
diff.track <- diff.data[,c("ProxID", "fid", "ClinSig")]
zone.diff.dummy <- createDummyFeatures(diff.data$zone)
diff.track.zone <- cbind(diff.track, zone.diff.dummy)
diff.track.zone <- distinct(diff.track.zone)

# diff track data
diff.track.patient <- diff.track %>%
  left_join(patient_sag, by = "ProxID") %>%
  left_join(diff.track.zone, by = c("ProxID", "fid", "ClinSig"))

# Create sag - containing image data, ProxID, fid, age, weight, height, ClinSig, zone information
diff.eq <- NULL
for (i in 1:length(diffeq.roi)){
  diff.eq[[i]] <- as.vector(diffeq.roi[[i]])
}
diff.eq <- as.matrix(do.call(rbind, diff.eq))
diff <- cbind(diff.eq, diff.track.patient)

saveRDS(diff, "diff_data.RDS")
write.csv(diff, "diff_lesion.csv")