#Kyle R Kovach
#SpectralPredictR
#December 1, 2020

#----Load Libraries----
# install.packages("doParallel")
# install.packages("spectrolab")
# install.packages("dplyr")
# install.packages("tidyr")

require(doParallel)
require(spectrolab)
require(dplyr)
require(tidyr)

#----Load Jump Correct Function----
jump_correct_specdal = function(spec_df, splices, reference=2, method="additive") {
  corrected = spec_df
  for (row_idx in seq_len(nrow(spec_df))) {
    series = spec_df[row_idx, ]
    get_sequence_num = function(w) {
      for (i in seq_along(splices)) if (as.numeric(w) <= splices[i]) return(i)
      return(length(splices) + 1)
    }
    groups = split(series, sapply(as.numeric(names(series)), get_sequence_num))
    translate_y = function(ref, mov, right=TRUE) {
      if (method == "additive") {
        diff = if (right) ref[length(ref)] - mov[1] else ref[1] - mov[length(mov)]
        return(mov + diff)
      } else if (method == "multiplicative") {
        ratio = if (right) ref[length(ref)] / mov[1] else ref[1] / mov[length(mov)]
        return(mov * ratio)
      }
    }
    for (i in reference:(length(groups) - 1)) {
      groups[[i + 1]] = translate_y(groups[[i]], groups[[i + 1]], TRUE)
    }
    for (i in seq(reference, 2, by = -1)) {
      groups[[i - 1]] = translate_y(groups[[i]], groups[[i - 1]], FALSE)
    }
    corrected[row_idx, ] = unlist(groups)
  }
  corrected
}

apply_jump_correction = function(spec, instrument, splices=NULL, method="additive") {
  instrument = tolower(instrument)
  # Auto-select splices if not provided
  if (is.null(splices)) {
    if (instrument == "asd") splices = c(1000, 1800)
    if (instrument == "sig") splices = c(990, 1900)
  }
  message(sprintf("%s: %s jump correction.", toupper(instrument), method))
  spec_df = as.data.frame(spec)
  meta_cols = 1
  corrected = jump_correct_specdal(
    spec_df[, -(1:meta_cols)],
    splices=splices,
    reference=2,
    method=method
  )
  spec_df = cbind(spec_df[, 1:meta_cols, drop=FALSE], corrected)
  as_spectra(spec_df, name_idx=1)
}

#----Initialize----
set.seed(112345)
cores=detectCores()-1

specfolder=choose.dir(default="",caption="Please select main folder containing subdirectories of spectral files.")

split_path=function(path) {
  if (dirname(path) %in% c("_", path)) return(basename(path))
  return(c(basename(path),split_path(dirname(path))))
}

for (j in c("asd","sed","sig"))
{
  if (exists('finaloutput')){rm(finaloutput)}
  folderlist=unique(sub("/[^/]+$","",list.files(path=specfolder,
                                                pattern=paste0(".",j,"$"),
                                                recursive=TRUE,
                                                full.names = TRUE)))
  if (length(folderlist)!=0){
    cl=makeCluster(cores)
    registerDoParallel(cl)
    
    finaloutput=foreach(i=folderlist,
                        .combine=rbind,
                        .packages=c('spectrolab'))%dopar%
      {
        spectra=read_spectra(path=i,format=j)
        data.frame(spectra,folder=i,format=j,check.names=FALSE)
        };FST=as_spectra(finaloutput[,1:2152],name_idx=1)
    
    stopCluster(cl)
    
if (j=="asd") {
  FST = apply_jump_correction(FST, instrument="ASD")
}
if (j=="sig") {
  FST = apply_jump_correction(FST, instrument="SVC")
};FSTdf=data.frame(FST,check.names=FALSE)
    if (!exists('FFO')){FFO=FSTdf[0,]}
    FFO=rbind(FFO,FSTdf)
  }};FSTfinal=as_spectra(FFO,name_idx=1);FSTfinal_norm=normalize(FSTfinal)

FST_f=data.frame(FSTfinal)
FST_n=data.frame(FSTfinal_norm)

write.csv(FST_f,"SpectralPredictR_raw_spectra.csv",row.names=FALSE)
write.csv(FST_n,"SpectralPredictR_VN_jumpcorrected_spectra.csv",row.names=FALSE)

#----Process Spectra, Vector Normalize, and Resample----
## This assumes all metadata exists only as the left side columns, and that the spectral dataset is consistant and ends the columns.
headcount=ncol(FST_f)-2151
specstart=headcount+1

cl=makeCluster(cores);registerDoParallel(cl);sampledata_VN=foreach(i=1:nrow(FST_f),
                                                                   .combine=rbind)%dopar%
  {
    a=FST_f[i,]
    b=sqrt(rowSums(a[,specstart:ncol(a)]^2))
    d=a[,-c(1:headcount)]
    sample_name=a[,c(1:headcount)]
    f=d/b
    cbind(sample_name,f)
  }
stopCluster(cl)

sampledata=sampledata_VN

sampledata_wav=sampledata[,-c(1:headcount)]
sampledata_head=sampledata[,c(1:headcount)]
sampledata_5nm_wav=sampledata_wav[,seq_len(ncol(sampledata_wav)) %% 5 == 1]
sampledata_5nm=cbind(sampledata_head,sampledata_5nm_wav)
sampledata=sampledata_5nm
write.csv(sampledata,"SpectralPredictR_5nmresampled_spectra.csv",row.names=FALSE)

#----Predict Traits----
# Choose model directory
modeldirectory=list.files(path=choose.dir(default="",
                                          caption="Please select main folder containing either fresh or dry models (based on spectra being processed)."),
                          pattern=".csv$",
                          recursive=FALSE,
                          full.names = TRUE)

cl=makeCluster(cores);registerDoParallel(cl);finished_output=foreach(model=modeldirectory,
                                                                     .combine=rbind)%dopar%
  {
    modelname=split_path(model)[1]
    modelcoeffs_mat=data.matrix(read.csv(model, header = TRUE, row.names = 1,sep = ","))
    sampledata_tail=sampledata[,-c(1:headcount)]
    sampledata_head=sampledata[,c(1:headcount)]
    sampledata_tail$constant=1
    sampledata_tail=cbind(intercept=sampledata_tail[,ncol(sampledata_tail)],sampledata_tail[,c(1:ncol(sampledata_tail)-1)])
    sampledata_tail_mat=as.matrix(sampledata_tail)
    modeloutput=sampledata_tail_mat %*% t(modelcoeffs_mat)
    modeloutput_bind=cbind(sampledata_head,as.data.frame(modeloutput))
    data.frame(modelname=modelname,
               sampledata_head,
               t_mean=rowMeans(modeloutput),
               t_std=apply(modeloutput,1,sd),
               as.data.frame(modeloutput),
               check.names=FALSE)
  };stopCluster(cl)

predsub=finished_output[,c(1:2,6:7)]

predspread=pivot_wider(data=predsub,id_cols=sample_name,names_from=modelname,values_from=c("t_mean","t_std"))
predunlist=as.data.frame(unnest(predspread))
colnames(predunlist)=sub(".csv", "", colnames(predunlist))
write.csv(predunlist,"SpectralPredictR_trait_output.csv",row.names=FALSE)
