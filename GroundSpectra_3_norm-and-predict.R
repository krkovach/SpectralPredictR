#Kyle R Kovach
#ASD/SED/SIG Spectra Processing
#December 18, 2020

#----Load Libra0ries----
require(doParallel)
require(tools)

#----Initialize----
set.seed(112345)
cores=detectCores()-1
setwd("F:/PostDoc/Internal_Projects/NEON/Database/Data/Ground Spectra")
split_path=function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path),split_path(dirname(path))))
}

#----Read in Data----
fname=file.choose()


#----Process Spectra and Vector Normalize----

  sampledata=read.csv(fname,header=TRUE)
  samplename=file_path_sans_ext(split_path(fname)[1])
  # sampledata$Sample=paste(sampledata$Date,sampledata$Sample,sep=".")
  headcount=ncol(sampledata)-2151
  specstart=headcount+1
  
  cl=makeCluster(cores);registerDoParallel(cl);sampledata_VN=foreach(i=1:nrow(sampledata),
                                                                     .combine=rbind)%dopar%
    {
      a=sampledata[i,]
      b=sqrt(rowSums(a[,specstart:ncol(a)]^2))
      d=a[,-c(1:headcount)]
      e=a[,c(1:headcount)]
      f=d/b
      cbind(e,f)
    }
  stopCluster(cl)
  
  sampledata_wav=sampledata_VN[,-c(1:headcount)]
  sampledata_head=sampledata_VN[,c(1:headcount)]
  sampledata_5nm_wav=sampledata_wav[,seq_len(ncol(sampledata_wav)) %% 5 == 1]
  sampledata_5nm=cbind(sampledata_head,sampledata_5nm_wav)
  sampledata_VN=sampledata_5nm
  
  # write.csv(sampledata_VN,"NEON2019_sampledata_fresh_VN.csv")

#----Predict Traits----
  sampledata_VN_fresh=subset(sampledata_VN,sampledata_VN$Type=="Fresh")
  sampledata_VN_dry=subset(sampledata_VN,sampledata_VN$Type=="Dry")
  
modeldirectory=list.files(path=choose.dir(),
                          pattern=".csv$",
                          recursive=FALSE,
                          full.names = TRUE)

cl=makeCluster(cores);registerDoParallel(cl);finished_output=foreach(model=modeldirectory,
                                                                     .combine=rbind)%dopar%
  {
    modelname=split_path(model)[1]
    modelcoeffs_mat=data.matrix(read.csv(model, header = TRUE, row.names = 1,sep = ","))
    sampledata_VN_tail=sampledata_VN_dry[,-c(1:headcount)]
    sampledata_VN_head=sampledata_VN_dry[,c(1:headcount)]
    sampledata_VN_tail$constant=1
    sampledata_VN_tail=cbind(intercept=sampledata_VN_tail[,ncol(sampledata_VN_tail)],sampledata_VN_tail[,c(1:ncol(sampledata_VN_tail)-1)])
    sampledata_VN_tail_mat=as.matrix(sampledata_VN_tail)
    modeloutput=sampledata_VN_tail_mat %*% t(modelcoeffs_mat)
    modeloutput_bind=cbind(sampledata_VN_head,as.data.frame(modeloutput))
    data.frame(modelname=modelname,
               sampledata_VN_head,
               t_mean=rowMeans(modeloutput),
               t_std=apply(modeloutput,1,sd),
               as.data.frame(modeloutput),
               check.names=FALSE)
  };stopCluster(cl)
write.csv(finished_output,paste0(samplename,"_dry_model_output.csv"))
