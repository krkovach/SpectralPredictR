#Kyle R Kovach
#ASD/SED/SIG Spectra Processing
#December 1, 2020

#----Load Libraries----
# install.packages("doParallel")
# install.packages("spectrolab")

require(doParallel)
require(spectrolab)

#----Initialize----
set.seed(112345)
cores=detectCores()-1

enspecdrive=choose.dir(default="",caption="Please select main folder containing subdirectories of spectral files.")

#----Import Spectra----
for (j in c("asd","sed","sig"))
{
  folderlist=unique(sub("/[^/]+$","",list.files(path=enspecdrive,
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
    
    if (j=="asd"||j=="sig"){
      splice_val=ifelse(j=="asd",c(1000, 1800),c(990, 1900))
      match_sensors(
        FST,
        splice_at=splice_val
      )
    };FSTdf=data.frame(FST,check.names=FALSE)
    if (!exists('FFO')){FFO=FSTdf[0,]}
    if (!exists('RT')){RT=finaloutput[0,c(1,2153:2154)]}
    FFO=rbind(FFO,FSTdf);RT=rbind(RT,finaloutput[,c(1,2153:2154)])
  }};FSTfinal=as_spectra(FFO,name_idx=1);FSTfinal_norm=normalize(FSTfinal)

#----Data Exploration----
# plot_interactive(FSTfinal_norm)

#----Table Output----
FST_f=data.frame(FSTfinal);FST_f=cbind(RT,FSTfinal)
FST_n=data.frame(FSTfinal_norm);FST_n=cbind(RT,FSTfinal_norm)
FST_f$SpecGrade=NA;FST_f$SpecGrade[FST_f$"450">0.8]=3;FST_f$SpecGrade[FST_f$"X1000"<0.05]=3

write.csv(FST_f,"GroundSpectra-RAW_NEON_2020.csv",row.names=FALSE)
write.csv(FST_n,"GroundSpectra-NORMALIZED_NEON_2020.csv",row.names=FALSE)