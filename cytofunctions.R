#-----------------------------------------------------------------------------
cytoread <- function(combined) {
#
# --- read the YYYYMMDD_HHMMSS_cytot (ee eme emm)  files
#  
  if(!combined) {    # electrical only
    files_cytotee = list.files(pattern="cytotee.txt",all.files = TRUE)
    cytotee.number=length(files_cytotee)
    cat("number of record file cytotee : ",cytotee.number )
#
# --- initialize the data.frame : cellem.data (elecal)
#
# - define the column (number and name) of cell.data
    cellem.data=data.frame(Characters=character(), Characters=character(),
                           Characters=character(), Characters=character() ,
                           Characters=character(), 
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           stringsAsFactors=FALSE)
    
    colnames(cellem.data) = c( 'cell.name' , 'cell.type' , 'type.num' ,
                               'exp.number' , 'exptype',
                               'passing.time', 'elec.duration',
                               'VLF.real', 'VLF.imag','VMF.real', 'VMF.imag' ,
                               'VHF.real', 'VHF.imag')
#
# - loop on the files to fill cellem.data
#
  for(ifile in 1:cytotee.number) {
      file.name=files_cytotee[ifile] # 
      file.data=read.table(file.name,sep=",") # read the record irecord
      cell.number=file.data[1,3];cell.number
      cell.type=file.data[1,1];cell.type
      exp.number=file.data[1,2] ; exp.number=substring(exp.number,5)
      type.num=0 # to be filled knowing all the cell types
      exptype=paste(exp.number,cell.type,sep="-")
# - loop on the cell.number for each record to fill the cell.data dataframe
    for(i in 1:cell.number) {
      cell.name=paste(substr(file.name,start=1,stop=16),as.character(i),sep="")
      cellem.data[nrow(cellem.data)+1,]=
        list(cell.name , cell.type , type.num , exp.number , exptype ,
             as.double(file.data[i+2,1]) , as.double(file.data[i+2,2]),
             as.double(file.data[i+2,3]) , as.double(file.data[i+2,4]),
             as.double(file.data[i+2,5]) , as.double(file.data[i+2,6]),
             as.double(file.data[i+2,7]),  as.double(file.data[i+2,8]))  }}
  total.cell.number=dim(cellem.data)[1] ; # total number of cells
  cell.type.list=unique(cellem.data$cell.type) ; # all type of cells
# - allocate type.num : index of type of cell for the classification 
  for(i in 1:total.cell.number) {
    cellem.data$type.num[i]=as.character(
      indexvalue(cellem.data$cell.type[i],cell.type.list) ) } }
  
  if(combined) { # electromechanical data
    files_cytotem = list.files(pattern="cytotem.txt",all.files = TRUE)
    cytotem.number=length(files_cytotem)
    cat("number of record file cytotem : ",cytotem.number )
#
# --- initialize the data.frame : cellem.data (electrical)
#
# - define the column (number and name) of cell.data
    cellem.data=data.frame(Characters=character(), Characters=character(),
                           Characters=character(), Characters=character() ,
                           Characters=character(), Doubles=double() , 
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),Doubles=double(),
                           Doubles=double(),stringsAsFactors=FALSE)
    colnames(cellem.data) = c( 'cell.name' , 'cell.type' , 'type.num' ,
                               'exp.number' , 'exptype',
                               'passing.time', 'elec.duration',
                               'mech.duration',
                               'VLF.real', 'VLF.imag','VMF.real', 'VMF.imag' ,
                               'VHF.real', 'VHF.imag', 'MEC.ampl')
    #
    # - loop on the files to fill cellem.data
    #
    for(ifile in 1:cytotem.number) {
      file.name=files_cytotem[ifile] # 
      file.data=read.table(file.name,sep=",") # read the record irecord
      cell.number=file.data[1,3];cell.number
      cell.type=file.data[1,1];cell.type
      exp.number=file.data[1,2] ; exp.number=substring(exp.number,5)
      type.num=0 # to be filled knowing all the cell types
      exptype=paste(exp.number,cell.type,sep="-")
# - loop on the cell.number for each record to fill the cell.data dataframe
      for(i in 1:cell.number) {
        cell.name=paste(substr(file.name,start=1,stop=16),as.character(i),sep="")
        cellem.data[nrow(cellem.data)+1,]=
          list(cell.name , cell.type , type.num , exp.number , exptype ,
               as.double(file.data[i+2,1]) , as.double(file.data[i+2,2]) ,
               as.double(file.data[i+2,9]) , 
               as.double(file.data[i+2,3]) , as.double(file.data[i+2,4]) ,
               as.double(file.data[i+2,5]) , as.double(file.data[i+2,6]) ,
               as.double(file.data[i+2,7]) , as.double(file.data[i+2,8]) ,
               as.double(file.data[i+2,10])   )  } }
    #
    total.cell.number=dim(cellem.data)[1] ; # total number of cells
    cell.type.list=unique(cellem.data$cell.type) ; # all type of cells
    # - allocate type.num : index of type of cell for the classification 
    for(i in 1:total.cell.number) {
      cellem.data$type.num[i]=as.character(
        indexvalue(cellem.data$cell.type[i],cell.type.list) ) }   }
assign("cellem.data", cellem.data, envir=.GlobalEnv)   }

#-----------------------------------------------------------------------------

indexvalue <- function(value,unique){
  # return the position of value in vector unique
  for( i in 1 : length(unique)) {
    if(value==unique[i]) {index=i} 
  }
  index   }

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# calculate the hertz model extraction and plot
cyto_hertz <- function(gap){
  cbind(MEC.ampl.hertz = 0, cellem.data) # create the hertz fit calculation
  
  maxtype=max(cellem.data$type.num) # for all the cell type
  #
  for (it in 1:maxtype) {
    cellem.data.typenum=cellem.data[cellem.data$type.num==it,] # isolate the cell type
    fit <- nls(MEC.ampl ~ hertz(size, a , b=gap), data=cellem.data.typenum, start=list(a=1))
    Eyoung=coef(fit)[[1]]  ;   Gap = gap 
# Calculate the results for the Hertz model for type cell it  for final plot   
    cellem.data$MEC.ampl.hertz[cellem.data$type.num==it] = 
          hertz(cellem.data$size[cellem.data$type.num==it], a=Eyoung , b=Gap)
    
    message(" R^2 regression = ",rsquare(fit, cellem.data.typenum)) 
    message(" E Young =  ",Eyoung)  
    message(" GAP =  ",Gap)
    
    p1 <- ggplot(data=cellem.data.typenum)
    p1=p1 + geom_point(aes(x=size, y=MEC.ampl),color="red", size=0.5)
    p1=p1+geom_function(fun = function(x) hertz(x, a=Eyoung , b=Gap)     )
    grid.arrange(p1,ncol=1) 
  }
  
  theme_set(theme_bw())
  p1<-ggplot(cellem.data) +
    geom_point(aes(x=size, y=MEC.ampl , color=cell.type) , size=0.1) +
    geom_line (aes(x=size, y=MEC.ampl.hertz , color=cell.type), size=0.5) +
    labs(x= "cubic.VLF.real")
  grid.arrange(p1,ncol=1)

  
  assign("cellem.data", cellem.data, envir=.GlobalEnv) }

# ---------------------------------------------------------------------------------------------------------------------------------------------------

hertz <- function(x, a, b) { 
  # x:cell diameter   a:young modulus   b:width of the constriction
  ifelse( x <= b, 0, a * 4/3 * (x/2)^(1/2) * ((x-b)/2)^(3/2))  }

# ---------------------------------------------------------------------------------------------------------------------------------------------------


cytowavee <- function() {
  ncell=length(wavee.data) # total number of cells
wavee.data.time    =lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,1]))
wavee.data.VLF.real=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,2]))
wavee.data.VLF.imag=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,3]))
wavee.data.VMF.real=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,4]))
wavee.data.VMF.imag=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,5]))
wavee.data.VHF.real=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,6]))
wavee.data.VHF.imag=lapply(1:ncell, function(i) as.numeric(wavee.data[[i]][,7]))

assign("wavee.data", wavee.data, envir=.GlobalEnv) 
assign("wavee.data.time", wavee.data.time, envir=.GlobalEnv) 
assign("wavee.data.VLF.real", wavee.data.VLF.real, envir=.GlobalEnv)
assign("wavee.data.VLF.imag", wavee.data.VLF.imag, envir=.GlobalEnv)
assign("wavee.data.VMF.real", wavee.data.VMF.real, envir=.GlobalEnv)
assign("wavee.data.VMF.imag", wavee.data.VMF.imag, envir=.GlobalEnv)
assign("wavee.data.VHF.real", wavee.data.VHF.real, envir=.GlobalEnv)
assign("wavee.data.VHF.imag", wavee.data.VHF.imag, envir=.GlobalEnv) }

# ------------------------------------------------------------

cytowavem <- function() {
  ncell=length(wavem.data) # total number of cells
  
    wavem.data.time     =lapply(1:ncell, function(i) as.numeric(wavem.data[[i]][,1]))
    wavem.data.MEC.ampl =lapply(1:ncell, function(i) as.numeric(wavem.data[[i]][,2]))
    
assign("wavem.data", wavem.data, envir=.GlobalEnv) 
assign("wavem.data.time", wavem.data.time , envir=.GlobalEnv) 
assign("wavem.data.MEC.ampl", wavem.data.MEC.ampl, envir=.GlobalEnv)  }

# ---------------------------------------------------

cytospline <- function(varnum1,varnum2,varnum1.min,varnum1.max,p){
  coefmat = NULL
  n=length(varnum1)
  for (j in 1:n){
    courbe = as.matrix(varnum2[[j]])
    t = as.vector(varnum1[[j]])
    b = create.bspline.basis(c(varnum1.min,varnum1.max),p)
    lambda=0
    fdParObj = fdPar(b,int2Lfd(2),lambda)
    fdSOpt = smooth.basis(t,courbe,fdParObj)
    coefmat  = cbind(coefmat,coef(fdSOpt))    }
  return(coefmat)   } 

#-------------------------------------------------------

cytoscores.e<- function(coefmat.VLF.real,coefmat.VLF.imag,
                      coefmat.VMF.real,coefmat.VMF.imag ,
                      coefmat.VHF.real,coefmat.VHF.imag,
                      cellem.data , scornum) {

scores.VLF.real=as.data.frame(t(coefmat.VLF.real)[,1:scornum])
scores.VLF.imag=as.data.frame(t(coefmat.VLF.imag)[,1:scornum])
scores.VMF.real=as.data.frame(t(coefmat.VMF.real)[,1:scornum])
scores.VMF.imag=as.data.frame(t(coefmat.VMF.imag)[,1:scornum])
scores.VHF.real=as.data.frame(t(coefmat.VHF.real)[,1:scornum])
scores.VHF.imag=as.data.frame(t(coefmat.VHF.imag)[,1:scornum])
scores_all=cbind(scores.VLF.real[,1:scornum],scores.VLF.imag[,1:scornum],
                 scores.VMF.real[,1:scornum],scores.VMF.imag[,1:scornum],
                 scores.VHF.real[,1:scornum],scores.VHF.imag[,1:scornum] )
scores_all$type.cell=cellem.data$cell.type
scores_all$type.num=cellem.data$type.num
scores_all$exp.number=cellem.data$exp.number
#
nomscol=c(1:scornum)
nomscol1=paste0("VLF.real.",nomscol)
nomscol2=paste0("VLF.imag.",nomscol)
nomscol3=paste0("VMF.real.",nomscol)
nomscol4=paste0("VMF.imag.",nomscol)
nomscol5=paste0("VHF.real.",nomscol)
nomscol6=paste0("VHF.imag.",nomscol)
colnames(scores_all)=c(nomscol1,nomscol2,nomscol3,nomscol4,nomscol5,nomscol6,
                       "cell.type","type.num", "exp.number") 
return(scores_all) }
  
# --------------------------------------

#-------------------------------------------------------

cytoscores.m<- function(coefmat.VLF.real,coefmat.VLF.imag,
                        coefmat.VMF.real,coefmat.VMF.imag ,
                        coefmat.VHF.real,coefmat.VHF.imag,
                        coefmat.MEC.ampl,cellem.data , scornum) {
  
  scores.VLF.real=as.data.frame(t(coefmat.VLF.real)[,1:scornum])
  scores.VLF.imag=as.data.frame(t(coefmat.VLF.imag)[,1:scornum])
  scores.VMF.real=as.data.frame(t(coefmat.VMF.real)[,1:scornum])
  scores.VMF.imag=as.data.frame(t(coefmat.VMF.imag)[,1:scornum])
  scores.VHF.real=as.data.frame(t(coefmat.VHF.real)[,1:scornum])
  scores.VHF.imag=as.data.frame(t(coefmat.VHF.imag)[,1:scornum])
  scores.MEC.ampl=as.data.frame(t(coefmat.MEC.ampl)[,1:scornum])
  scores_all=cbind(scores.VLF.real[,1:scornum],scores.VLF.imag[,1:scornum],
                   scores.VMF.real[,1:scornum],scores.VMF.imag[,1:scornum],
                   scores.VHF.real[,1:scornum],scores.VHF.imag[,1:scornum],
                   scores.MEC.ampl[,1:scornum] )              
  scores_all$type.cell=cellem.data$cell.type
  scores_all$type.num=cellem.data$type.num
  scores_all$exp.number=cellem.data$exp.number
  #
  nomscol=c(1:scornum)
  nomscol1=paste0("VLF.real.",nomscol)
  nomscol2=paste0("VLF.imag.",nomscol)
  nomscol3=paste0("VMF.real.",nomscol)
  nomscol4=paste0("VMF.imag.",nomscol)
  nomscol5=paste0("VHF.real.",nomscol)
  nomscol6=paste0("VHF.imag.",nomscol)
  nomscol7=paste0("MEC.ampl.",nomscol)
  colnames(scores_all)=c(nomscol1,nomscol2,nomscol3,nomscol4,nomscol5,nomscol6,
                         nomscol7,"cell.type","type.num", "exp.number") 
  return(scores_all) }

# --------------------------------------
  
cytomaxmin.e<- function(maxmin.VLF.real,maxmin.VLF.imag,
                      maxmin.VMF.real,maxmin.VMF.imag,
                      maxmin.VHF.real,maxmin.VHF.imag) {
  maxmin_all=as.data.frame(maxmin.VLF.real)
  colnames(maxmin_all)=c("VLF.real")
  maxmin_all$VLF.imag=maxmin.VLF.imag
  maxmin_all$VMF.real=maxmin.VMF.real
  maxmin_all$VMF.imag=maxmin.VMF.imag
  maxmin_all$VHF.real=maxmin.VHF.real
  maxmin_all$VHF.imag=maxmin.VHF.imag
  maxmin_all$cell.type=cellem.data$cell.type
  maxmin_all$type.num=cellem.data$type.num
  maxmin_all$exp.number=cellem.data$exp.number
  return(maxmin_all) }

# -----------------------------------------------

# --------------------------------------

cytomaxmin.m<- function(maxmin.VLF.real,maxmin.VLF.imag,
                        maxmin.VMF.real,maxmin.VMF.imag,
                        maxmin.VHF.real,maxmin.VHF.imag,maxmin.MEC.ampl) {
  maxmin_all=as.data.frame(maxmin.VLF.real)
  colnames(maxmin_all)=c("VLF.real")
  maxmin_all$VLF.imag=maxmin.VLF.imag
  maxmin_all$VMF.real=maxmin.VMF.real
  maxmin_all$VMF.imag=maxmin.VMF.imag
  maxmin_all$VHF.real=maxmin.VHF.real
  maxmin_all$VHF.imag=maxmin.VHF.imag
  maxmin_all$MEC.ampl=maxmin.MEC.ampl
  maxmin_all$cell.type=cellem.data$cell.type
  maxmin_all$type.num=cellem.data$type.num
  maxmin_all$exp.number=cellem.data$exp.number
  return(maxmin_all) }

# -----------------------------------------------

cytotraintest<- function (cellem.data,train.exp,train.percent,test.exp) {
  total.cell.number=dim(cellem.data)[1];total.cell.number  
  if(is.null(test.exp)) {
    train.sample=c(1:total.cell.number)[cellem.data$exp.number%in%train.exp]
    train.sample.number=length(train.sample)
    train.select.number=round(train.sample.number*train.percent/100)
    train.select=sample(train.sample,train.select.number)## random selection
    test.select=setdiff(train.sample,train.select)     }
  else{
    train.select=c(1:total.cell.number)[cellem.data$exp.number%in%train.exp]
    test.select= c(1:total.cell.number)[cellem.data$exp.number%in%test.exp]   }
  
  assign("train.select", train.select, envir=.GlobalEnv)
  assign("test.select", test.select, envir=.GlobalEnv)   }

# -----------------------------------------------------------------------

cytodf.cell.class.result <- function(train.select,test.select) {
  
  # - define and initialize the data frame to collect the cell classification results
  # - initiale with the number of cell (per cell.type) in the test
  
  total.cell.type=unique(cellem.data$cell.type) # the names of the different cell type in the exp
  tct=(length(total.cell.type)) # the number of different cell type in the exp
  cell.type.number=(length(total.cell.type)) # the number of different cell type in the exp
  
  all=cellem.data
  train=cellem.data[(train.select),]
  test=cellem.data[(test.select),]
  
  #- the number of cells per cell.type (for all, train and test ) 
  
  all.cell=lapply(1:tct, function(i) dim(all[all$cell.type==total.cell.type[i],])[1])
  train.cell=lapply(1:tct, function(i) dim(train[train$cell.type==total.cell.type[i],])[1])
  test.cell=lapply(1:tct, function(i) dim(test[test$cell.type==total.cell.type[i],])[1])

  
#  tcell=as.integer(test.cell[]) # number of cell per cell type in the test
  
  class.results=data.frame(Characters=character(), 
                           logical=logical(), logical=logical() ,
                           logical=logical(), logical=logical(),
                           Characters=character(),
                           Doubles=double(),Doubles=double(),
                           stringsAsFactors=FALSE)
  
  
  
  cell.class.result = data.frame(matrix(ncol = 9+cell.type.number, nrow = 0))
  colnames(cell.class.result) = c('SAMPLE',total.cell.type,
                                  'data.type', 'method' ,
                                  'class.prob','prediction',
                                  'LF' , 'MF', 'HF' , 'MEC')
  cell.class.result[nrow(cell.class.result)+1,]=
    c('ALL' , all.cell, ''  , '' , '' , '' , '' , '' , '' , '')
  cell.class.result[nrow(cell.class.result)+1,]=
    c('TRAIN' , train.cell, ''  , '' , '' , '' , '' , '' , '' , '')
  cell.class.result[nrow(cell.class.result)+1,]=
    c('TEST' , test.cell, ''  , '' , '' , '' , '' , '' , '' , '')
  
  assign("cellem.data", cellem.data, envir=.GlobalEnv)  
  assign("cell.class.result", cell.class.result, envir=.GlobalEnv)  }

# ----------------------------------------------- 

cytodf.class.result <- function(){
  # - define the data frame to collect the classification results

    class.results=data.frame(Characters=character(), 
                             logical=logical(), logical=logical() ,
                             logical=logical(), logical=logical(),
                             Characters=character(),
                             Doubles=double(),Doubles=double(),
                             stringsAsFactors=FALSE)
    
    colnames(class.results) = c( 'data.type' , 'LF' , 'MF', 'HF' , 'MEC' ,
                                 'method', 'class.prob', 'prediction') 
                             
  assign("class.results", class.results, envir=.GlobalEnv)  }

# -----------------------------------------------



cytoclassif<- function(sample_data,data.name,train.select,test.select,nbfreq,combined){
# Make the classification for all measurement and the 3 methods
  train_data=sample_data[(train.select),] # chose the training
  test_data =sample_data[(test.select),] # chose the test
  
  if(!combined) {

if(nbfreq==2) { 
  # cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE, MF=TRUE, HF=FALSE,MEC=FALSE)
  # cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE, MF=FALSE,HF=TRUE ,MEC=FALSE)
  cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE ,MF=TRUE, HF=FALSE,MEC=FALSE) }
  
if(nbfreq==3) { 
    
    cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE,HF=TRUE,MEC=FALSE)
#   cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=FALSE,HF=FALSE,MEC=FALSE)
#   cytoclassif.meas(train_data,test_data,data.name,combined,LF=FALSE,MF=TRUE,HF=FALSE,MEC=FALSE)
#   cytoclassif.meas(train_data,test_data,data.name,combined,LF=FALSE,MF=FALSE,HF=TRUE,MEC=FALSE)
    cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE,HF=FALSE,MEC=FALSE)
    cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=FALSE,HF=TRUE,MEC=FALSE)
    cytoclassif.meas(train_data,test_data,data.name,combined,LF=FALSE,MF=TRUE,HF=TRUE,MEC=FALSE) }}
 
  if(combined) {
    if(nbfreq==2) { 
      cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE, HF=FALSE,MEC=TRUE)
      cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE, HF=FALSE,MEC=FALSE) }
    
    if(nbfreq==3) { 
      
      cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE,HF=TRUE,MEC=FALSE)
      cytoclassif.meas(train_data,test_data,data.name,combined,LF=TRUE,MF=TRUE,HF=TRUE,MEC=TRUE ) }} 
  
  assign("class.results", class.results, envir=.GlobalEnv) 
  assign("trainX", trainX, envir=.GlobalEnv)  
  assign("trainY", trainY, envir=.GlobalEnv)
  assign("train_data", train_data, envir=.GlobalEnv)  
  assign("test_data", test_data, envir=.GlobalEnv)
  assign("testX", testX, envir=.GlobalEnv)  
  assign("testY", testY, envir=.GlobalEnv) }


# -----------------------------------------------

cytoclassif.meas<- function(train_data,test_data,data.name,combined,LF,MF,HF,MEC){
# compute the classification
  trainX = data.frame(matrix(ncol = 0, nrow = nrow(train_data)))
  testX =  data.frame(matrix(ncol = 0, nrow = nrow(test_data )))
#  LF=TRUE ; MF=TRUE ; HF=TRUE 
  print(data.name)
  if(data.name=="scores") {
    if (LF)  {  trainX = cbind(trainX,train_data[,c(1:10)])
                testX  = cbind(testX,test_data  [,c(1:10)])   }
    if (MF)  {  trainX = cbind(trainX,train_data[,c(11:20)])
                testX  = cbind(testX,test_data  [,c(11:20)])  }
    if (HF)  {  trainX = cbind(trainX,train_data[,c(21:30)])
                testX  = cbind(testX,test_data  [,c(21:30)])  }
    if (MEC) {  trainX = cbind(trainX,train_data[,c(31:35)])
                testX  = cbind(testX,test_data  [,c(31:35)])  }
    
    if(!combined) {  trainY=as.factor(train_data[,c(32)])  
                     testY =as.factor(test_data [,c(32)])     }
    if(combined) {   trainY=as.factor(train_data[,c(37)])  
                     testY =as.factor(test_data [,c(37)])     }
  }
  else {
    if (LF)  {  trainX = cbind(trainX,train_data[,c(1:2)])
                testX  = cbind(testX,test_data  [,c(1:2)])   }
    if (MF)  {  trainX = cbind(trainX,train_data[,c(3:4)])
                testX  = cbind(testX,test_data  [,c(3:4)])   }
    if (HF)  {  trainX = cbind(trainX,train_data[,c(5:6)])
                testX  = cbind(testX,test_data  [,c(5:6)])   } 
    if (MEC) {  trainX = cbind(trainX,train_data[,c(7)])
                testX  = cbind(testX,test_data  [,c(7)])     } 
    
    if(!combined) {  trainY=as.factor(train_data[,c(8)])  
                     testY =as.factor(test_data [,c(8)])     }
    if(combined) {   trainY=as.factor(train_data[,c(9)])  
                     testY =as.factor(test_data [,c(9)])     }
 }
  
  cytoclassif.meas.meth("lda",trainX,trainY,testX,testY,data.name,LF,MF,HF,MEC)
  cytoclassif.meas.meth("glm",trainX,trainY,testX,testY,data.name,LF,MF,HF,MEC)
  cytoclassif.meas.meth("gam",trainX,trainY,testX,testY,data.name,LF,MF,HF,MEC)
  
  assign("class.results", class.results, envir=.GlobalEnv) 
  assign("trainX", trainX, envir=.GlobalEnv)  
  assign("trainY", trainY, envir=.GlobalEnv)
  assign("testX", testX, envir=.GlobalEnv)  
  assign("testY", testY, envir=.GlobalEnv)
  # assign("pred1", pred1, envir=.GlobalEnv)  
  # assign("pred2", pred2, envir=.GlobalEnv) 
  # assign("pred3", pred3, envir=.GlobalEnv) 
  assign("cellem.data", cellem.data, envir=.GlobalEnv) } 



# -----------------------------------------------

cytoclassif.meas.meth<- function(meth,trainX,trainY,testX,testY,data.name,LF,MF,HF,MEC){

  test.number=length(testY) 
  ctrl = list(verbose = FALSE, draw = TRUE, alpha = 0.5, main="", xlab="", ylab="", xlim=c(0,0.2),xlim=c(0,0.1))
  
  out=classif.DD(trainY,trainX,depth="MhD",classif=meth, control=ctrl)
  summary(out)
  predY=predict(out,testX)
  pred.meth=sum(diag(table(predY,testY)))/test.number
  prob.class.meth=as.numeric(1.0-as.numeric(out[2]))

  prLF="X" ; if(!LF)  {prLF="  "}  
  prMF="X" ; if(!MF)  {prMF="  "}  
  prHF="X" ; if(!HF)  {prHF="  "}  
  prMEC="X" ; if(!MEC) {prMEC="  "}  
  # ------------------------------------------------ 
  print(unique(cellem.data$cell.type))
  print(table(testY))
  print(table(predY))
  # ----------------------------------

  class.results[nrow(class.results)+1,]=
    list(data.name , prLF , prMF, prHF, prMEC ,meth, prob.class.meth , pred.meth )

  total.cell.type.num=unique(cellem.data$type.num) # the names of the different cell type in the exp
  tct=(length(total.cell.type.num)) # the number of different cell type in the exp
  predict.cell=lapply(1:tct, function(i) length(predY[predY==total.cell.type.num[i]])) 

    cell.class.result[nrow(cell.class.result)+1,]=
      c('PREDICT' , predict.cell, data.name, meth , 
        round(prob.class.meth,digits=5) , round(pred.meth,digits=5) , 
        prLF , prMF, prHF, prMEC)
    
  assign("class.results", class.results, envir=.GlobalEnv)
  assign("trainX", trainX, envir=.GlobalEnv)
  assign("trainY", trainY, envir=.GlobalEnv)
  assign("testX", testX, envir=.GlobalEnv)
  assign("testY", testY, envir=.GlobalEnv)
  assign("cell.class.result", cell.class.result, envir=.GlobalEnv)
  assign("cellem.data", cellem.data, envir=.GlobalEnv) }



cyplot.lin <- function(dataset,combined,varnum1,varname1,type){
  #
  # - linear plot of 6 freq peaks vs varnum1 according to type
  #
  theme_set(theme_bw())
  p1<-ggplot(dataset, aes(x=varnum1, y=VLF.real , color=type)) +
    geom_point(size=0.1) +
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  p2<-ggplot(dataset, aes(x=varnum1, y=VLF.imag , color=type)) +
    geom_point(size=0.1) +
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  p3<-ggplot(dataset, aes(x=varnum1, y=VMF.real , color=type)) +
    geom_point(size=0.1) + 
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  p4<-ggplot(dataset, aes(x=varnum1, y=VMF.imag , color=type)) +
    geom_point(size=0.1) + 
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  p5<-ggplot(dataset, aes(x=varnum1, y=VHF.real , color=type)) +
    geom_point(size=0.1) + 
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  p6<-ggplot(dataset, aes(x=varnum1, y=VHF.imag , color=type)) +
    geom_point(size=0.1) + 
    labs(x =varname1) +
    geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2) 
  if(combined){ 
    p7<-ggplot(dataset, aes(x=varnum1, y=MEC.ampl , color=type)) +
      geom_point(size=0.1) + 
      labs(x =varname1) +
      geom_smooth(method=lm, se=FALSE , fullrange=FALSE)
    grid.arrange(p7,ncol=1) } }

# -------------------------------------------------------------

cyplot.box <- function(dataset,combined,type,typename){
  theme_set(theme_bw())
  p1<-ggplot(dataset, aes(type, VLF.real)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p2<-ggplot(dataset, aes(type, VLF.imag)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p3<-ggplot(dataset, aes(type, VMF.real)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p4<-ggplot(dataset, aes(type, VMF.imag)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p5<-ggplot(dataset, aes(type, VHF.real)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p6<-ggplot(dataset, aes(type, VHF.imag)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
  if(combined){ 
    p7<-ggplot(dataset, aes(type, MEC.ampl)) +
      geom_boxplot(aes(fill=type)) + 
      labs(x =NULL) +
      theme(axis.text.x = element_text(angle=65, vjust=0.6))
    grid.arrange(p7,ncol=1) } }

cyplot.violin <- function(dataset,combined,type,typename){
  theme_set(theme_bw())
  p1<-ggplot(dataset, aes(type, VLF.real)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p2<-ggplot(dataset, aes(type, VLF.imag)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p3<-ggplot(dataset, aes(type, VMF.real)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p4<-ggplot(dataset, aes(type, VMF.imag)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p5<-ggplot(dataset, aes(type, VHF.real)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p6<-ggplot(dataset, aes(type, VHF.imag)) +
    geom_violin(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2) 
  if(combined){ 
    p7<-ggplot(dataset, aes(type, MEC.ampl)) +
      geom_violin(aes(fill=type)) + 
      labs(x =NULL) +
      theme(axis.text.x = element_text(angle=65, vjust=0.6))
    grid.arrange(p7,ncol=1) }    }

cyplot.boxpoint <- function(dataset,combined,type,typename,varnum){
  theme_set(theme_bw())
  p1<-ggplot(dataset, aes(type, VLF.real)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p2<-ggplot(dataset, aes(type, VLF.imag)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p3<-ggplot(dataset, aes(type, VMF.real)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p4<-ggplot(dataset, aes(type, VMF.imag)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p5<-ggplot(dataset, aes(type, VHF.real)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p6<-ggplot(dataset, aes(type, VHF.imag)) +
    geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
  if(combined){ 
    p7<-ggplot(dataset, aes(type, MEC.ampl)) +
      geom_dotplot(binaxis='y',binwidth = varnum, stackdir='center') +
      labs(x =NULL) +
      theme(axis.text.x = element_text(angle=65, vjust=0.6))
    grid.arrange(p7,ncol=1) }    }

# ------------------------------------------------------------- 

cyplot.wavee <- function(dataset,icell){
  par(mfrow=c(3,2))
  x=as.numeric(dataset[[icell]][,1]);
  y=as.numeric(dataset[[icell]][,2]);
  plot(x,y,xlab="time", ylab="VLF.real")
  y=as.numeric(dataset[[icell]][,3]);
  plot(x,y,xlab="time", ylab="VLF.imag")
  y=as.numeric(dataset[[icell]][,4]);
  plot(x,y,xlab="time", ylab="VMF.real")
  y=as.numeric(dataset[[icell]][,5]);
  plot(x,y,xlab="time", ylab="VMF.imag")
  y=as.numeric(dataset[[icell]][,6]);
  plot(x,y,xlab="time", ylab="VMF.real")
  y=as.numeric(dataset[[icell]][,7]);
  plot(x,y,xlab="time", ylab="VMF.imag") }

# -------------------------------------------------------------

cyplot.wavem <- function(dataset,icell){
  par(mfrow=c(1,1))
  x=as.numeric(dataset[[icell]][,1]);
  y=as.numeric(dataset[[icell]][,2]);
  plot(x,y,xlab="time", ylab="MEC.ampl")  }

# -------------------------------------------------------------
cyplot.spline.coef.e<-function(Xfd.VLF.real,Xfd.VLF.imag,Xfd.VMF.real,Xfd.VMF.imag,
                               Xfd.VHF.real,Xfd.VHF.imag,icell)  {
  par(mfrow=c(3,2))
  plot(Xfd.VLF.real[icell]$coefs)
  plot(Xfd.VLF.imag[icell]$coefs)
  plot(Xfd.VMF.real[icell]$coefs)
  plot(Xfd.VMF.imag[icell]$coefs)
  plot(Xfd.VHF.real[icell]$coefs)
  plot(Xfd.VHF.imag[icell]$coefs)   }

# -------------------------------------

cyplot.spline.coef.m<-function(Xfd.MEC.ampl,icell)  {
  par(mfrow=c(1,1))
  plot(Xfd.MEC.ampl[icell]$coefs)  }

# -------------------------------------

cyplot.spline.time.e<-function(Xfd.VLF.real,Xfd.VLF.imag,Xfd.VMF.real,Xfd.VMF.imag,
                               Xfd.VHF.real,Xfd.VHF.imag,icell)  {
  par(mfrow=c(3,2))
  if(icell!=0) {
    plot(Xfd.VLF.real[icell],xlab="Time", ylab="VLF.real")
    plot(Xfd.VLF.imag[icell],xlab="Time", ylab="VLF.imag")
    plot(Xfd.VMF.real[icell],xlab="Time", ylab="VMF.real")
    plot(Xfd.VMF.imag[icell],xlab="Time", ylab="VMF.imag")
    plot(Xfd.VHF.real[icell],xlab="Time", ylab="VHF.real")
    plot(Xfd.VHF.imag[icell],xlab="Time", ylab="VHF.imag") }
  if(icell==0) {
    plot(Xfd.VLF.real,xlab="Time", ylab="VLF.real")
    plot(Xfd.VLF.imag,xlab="Time", ylab="VLF.imag")
    plot(Xfd.VMF.real,xlab="Time", ylab="VMF.real")
    plot(Xfd.VMF.imag,xlab="Time", ylab="VMF.imag")
    plot(Xfd.VHF.real,xlab="Time", ylab="VHF.real")
    plot(Xfd.VHF.imag,xlab="Time", ylab="VHF.imag") }
}

# ---------------------------------------------------------------------

cyplot.spline.time.m<-function(Xfd.MEC.ampl,icell)  {
  par(mfrow=c(1,1))
  if(icell!=0) {
    plot(Xfd.MEC.ampl[icell],xlab="Time", ylab="MEC.ampl") }
  if(icell==0) {
    plot(Xfd.MEC.ampl,xlab="Time", ylab="MEC.ampl") } }
# ---------------------------------------------------------------------


cyplot.box.scores <- function(dataset,combined,type,typename){
  theme_set(theme_bw())
  p1<-ggplot(dataset, aes(type, VLF.real.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p2<-ggplot(dataset, aes(type, VLF.imag.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p3<-ggplot(dataset, aes(type, VMF.real.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p4<-ggplot(dataset, aes(type, VMF.imag.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  p5<-ggplot(dataset, aes(type, VHF.real.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
  p6<-ggplot(dataset, aes(type, VHF.imag.1)) +
    geom_boxplot(aes(fill=type)) + 
    labs(x =NULL) +
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
  grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2) 
  
  if(combined){ 
    p7<-ggplot(dataset, aes(type, MEC.ampl.1)) +
      geom_boxplot(aes(fill=type)) + 
      labs(x =NULL) +
      theme(axis.text.x = element_text(angle=65, vjust=0.6)) 
    grid.arrange(p7,ncol=1) }    }





