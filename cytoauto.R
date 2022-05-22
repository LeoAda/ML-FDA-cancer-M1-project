#############
##### 1 #####
#############
wd = "C:/Users/leo/Desktop/projet m1/WD"
setwd(wd)
#############
#############
#############


library(lattice)
library(fda.usc)
library(ggpubr)
library(data.table)
library(ddalpha)
library(FactoMineR)
library(factoextra)
library(e1071)
library(rstatix)
library(corrplot)
library(tidyverse)
library(gridExtra)
library(modelr)
library(stringr)
library(refund)
library(face)
library(fdapace)
rm(list = ls())


#############
##### 2 #####
#############
wd = "C:/Users/leo/Desktop/projet m1/WD"
source(paste(wd,"/cytofunctions.R",sep=""))
#############
#############
#############


T1<-Sys.time() 
options(digits=5)
use.quartz=FALSE




#############
##### 3 ##### Nombre de score dans score_all
############# 
scores.number=10
p.spline.e=10
#############
############# Nombre fonction de base
#############

#############
##### 4 ##### Score minimal (scores.number <= p.spline.e <= remove )
############# 
remove = 10
#############
#############
#############


#############
##### 5 ##### Liste des datasets
#############
#exp.list = c('20211123_cyto_E3MS_US_data_set')
exp.list = c('20210507_cyto_E2_MJ_data_set','20210507_E2_MUS_data_set','20210723_cyto_E2_BSMU_data_set','20210719_cyto_E2_BSMU_data_set')
#############
#############
#############


#############
##### 6 ##### Si on combine tous les scores_all à la fin dans un fichier
#############
combine_data = TRUE
comb_df_score <- data.frame()
comb_df_vpeak <- data.frame()
comb_df_maxmin <- data.frame()
#############
#############
#############









for (exp.name in exp.list){
print(exp.name)
  
  
exp.name.list = str_split(exp.name,"_")[[1]]
exp.name.prop = exp.name.list[str_starts(exp.name.list,"\\E")]
combined = str_detect(exp.name.prop ,"MS")
nbfreq = as.numeric(str_extract(exp.name.prop,"[[:digit:]]"))
work.directory = paste(wd,"/data_set/",exp.name,sep="")
setwd(work.directory)
cytoread(combined)
total.cell.number=dim(cellem.data)[1] ; total.cell.number # total number of cells
cellem.data$VLF.ampl=sqrt(cellem.data$VLF.real^2+cellem.data$VLF.imag^2)
cellem.data$VMF.ampl=sqrt(cellem.data$VMF.real^2+cellem.data$VMF.imag^2)
cellem.data$VHF.ampl=sqrt(cellem.data$VHF.real^2+cellem.data$VHF.imag^2)
cellem.data$opacity=cellem.data$VHF.ampl / cellem.data$VLF.ampl
VLF.phase = atan(cellem.data$VLF.imag/cellem.data$VLF.real) 
VHF.phase = atan(cellem.data$VHF.imag/cellem.data$VHF.real)
cellem.data$phase.ratio = VHF.phase / VLF.phase
cellem.data$contrast = cellem.data$phase.ratio / cellem.data$VLF.ampl 
cellem.data$cub.VLF.real = (cellem.data$VLF.real)^(1/3)
cellem.data$cub.VHF.imag = (cellem.data$VHF.imag)^(1/3)
files_wave=cellem.data$cell.name # all the cell names
if(!combined) { files_wavee=paste(files_wave,"_wavee.txt",sep="")
  wavee.data = lapply(files_wavee, read.table, header = FALSE,skip=1, sep=",") } # all waves
if(combined) { files_waveeme=paste(files_wave,"_waveeme.txt",sep="") 
               files_waveemm=paste(files_wave,"_waveemm.txt",sep="")
  wavee.data = lapply(files_waveeme, read.table, header = FALSE,skip=1, sep=",") # all waves
  wavem.data = lapply(files_waveemm, read.table, header = FALSE,skip=1, sep=",") } # all waves




#
# Correction erreur
#
error_index = which(sapply(wavee.data, nrow) <= remove)
if(length(error_index)>0 ){
  wavee.data = wavee.data[-error_index]
  
  print(nrow(cellem.data))
  
  cellem.data <- cellem.data[-c(error_index),]
  print(nrow(cellem.data))
  
  
  print("removed index :")
  print(error_index)
  total.cell.number = total.cell.number-length(error_index)
  error_index <- list()
}

cytowavee() # organize the electrical wave data for functional fit 
if(combined) {cytowavem() } # organize the mechanical wave data for functional fit 


tmin=0 ; tmax=3; 
b = create.bspline.basis(c(tmin,tmax),p.spline.e) 
coefmat.VLF.real = cytospline(wavee.data.time,wavee.data.VLF.real,tmin,tmax,p.spline.e)
coefmat.VLF.imag = cytospline(wavee.data.time,wavee.data.VLF.imag,tmin,tmax,p.spline.e)
coefmat.VMF.real = cytospline(wavee.data.time,wavee.data.VMF.real,tmin,tmax,p.spline.e)
coefmat.VMF.imag = cytospline(wavee.data.time,wavee.data.VMF.imag,tmin,tmax,p.spline.e)
coefmat.VHF.real = cytospline(wavee.data.time,wavee.data.VHF.real,tmin,tmax,p.spline.e)
coefmat.VHF.imag = cytospline(wavee.data.time,wavee.data.VHF.imag,tmin,tmax,p.spline.e)
Xfd.VLF.real = fd(coefmat.VLF.real,b)
Xfd.VLF.imag = fd(coefmat.VLF.imag,b)
Xfd.VMF.real = fd(coefmat.VMF.real,b)
Xfd.VMF.imag = fd(coefmat.VMF.imag,b)
Xfd.VHF.real = fd(coefmat.VHF.real,b)
Xfd.VHF.imag = fd(coefmat.VHF.imag,b)
if(combined) {   tmin=0 ; tmax=5; p.spline.m=10 
  c = create.bspline.basis(c(tmin,tmax),p.spline.m) 
  coefmat.MEC.ampl = cytospline(wavem.data.time,wavem.data.MEC.ampl,tmin,tmax,p.spline.m)
  Xfd.MEC.ampl = fd(coefmat.MEC.ampl,c) }



# ---------- Maxmin ----------
maxmin.VLF.real=c() ; maxmin.VLF.imag=c() ; maxmin.VMF.real=c() ; 
maxmin.VMF.imag=c() ; maxmin.VHF.real=c() ; maxmin.VHF.imag=c()
for (i in 1:total.cell.number) {
  maxmin.VLF.real[i]=max(Xfd.VLF.real[i]$coefs)-min(Xfd.VLF.real[i]$coefs)
  maxmin.VLF.imag[i]=max(Xfd.VLF.imag[i]$coefs)-min(Xfd.VLF.imag[i]$coefs) 
  maxmin.VMF.real[i]=max(Xfd.VMF.real[i]$coefs)-min(Xfd.VMF.real[i]$coefs) 
  maxmin.VMF.imag[i]=max(Xfd.VMF.imag[i]$coefs)-min(Xfd.VMF.imag[i]$coefs)
  maxmin.VHF.real[i]=max(Xfd.VHF.real[i]$coefs)-min(Xfd.VHF.real[i]$coefs) 
  maxmin.VHF.imag[i]=max(Xfd.VHF.imag[i]$coefs)-min(Xfd.VHF.imag[i]$coefs) }

if (combined) { maxmin.MEC.ampl=c() 
for (i in 1:total.cell.number) { 
  maxmin.MEC.ampl[i]=max(Xfd.MEC.ampl[i]$coefs)-min(Xfd.MEC.ampl[i]$coefs)}}
if(!combined) {maxmin_all=cytomaxmin.e(maxmin.VLF.real,maxmin.VLF.imag,
                                       maxmin.VMF.real,maxmin.VMF.imag,
                                       maxmin.VHF.real,maxmin.VHF.imag)  }
if(combined) {maxmin_all=cytomaxmin.m(maxmin.VLF.real,maxmin.VLF.imag,
                                      maxmin.VMF.real,maxmin.VMF.imag,
                                      maxmin.VHF.real,maxmin.VHF.imag, 
                                      maxmin.MEC.ampl) }
# Vpeak
if(!combined) { vpeaks_all=cellem.data[,c(8:13,2:4)] }
if(combined)  { vpeaks_all=cellem.data[,c(9:15,2:4)] }
# ---------- scores_all ----------

if(!combined) { scores_all=cytoscores.e(coefmat.VLF.real,coefmat.VLF.imag,
                      coefmat.VMF.real,coefmat.VMF.imag,
                      coefmat.VHF.real,coefmat.VHF.imag, 
                      cellem.data , scores.number) }
if(combined)  { scores_all=cytoscores.m(coefmat.VLF.real,coefmat.VLF.imag,
                                      coefmat.VMF.real,coefmat.VMF.imag,
                                      coefmat.VHF.real,coefmat.VHF.imag, 
                                      coefmat.MEC.ampl, 
                                      cellem.data , scores.number) }



scores_all <- transform(scores_all, type.num = as.numeric(type.num))
scores_all <- transform(scores_all, exp.number = as.numeric(exp.number))
if(combine_data){
  comb_df_score <- rbind(comb_df_score,scores_all)
  comb_df_vpeak <- rbind(comb_df_vpeak,vpeaks_all)
  comb_df_maxmin <- rbind(comb_df_maxmin,maxmin_all)
}

}
if(combine_data){
  write.csv(comb_df_score,paste(wd,"/comb_df_score.csv",sep=""))
  write.csv(comb_df_vpeak,paste(wd,"/comb_df_vpeak.csv",sep=""))
  write.csv(comb_df_maxmin,paste(wd,"/comb_df_maxmin.csv",sep=""))
}

print("fini")
