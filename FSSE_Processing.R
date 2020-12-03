#Kyle R Kovach
#FieldSpec and Spectral Evolution Processing
#December 1, 2020

#----Load Libraries----
# install.packages("doParallel")
# install.packages("spectrolab")

require(doParallel)
require(spectrolab)

#----Initialize----
set.seed(112345)
cores=detectCores()-1

enspecdrive="Z:/neon/NEON_2020_field_data"

#----Import Spectra----
folderlist=unique(sub("/[^/]+$","",list.files(path=enspecdrive,pattern='.sed$', recursive=TRUE, full.names = TRUE)))

cl=makeCluster(cores[1]-1)
registerDoParallel(cl)

finaloutput=foreach(i=folderlist,
                    .combine=rbind,
                    .packages=c('spectrolab'))%dopar%
  {
    spectra=read_spectra(path=i,format = "sed")
    data.frame(spectra,check.names=FALSE)
  };FST=as_spectra(finaloutput,name_idx=1)

stopCluster(cl)

#----Initial Sweep for WR Outliers----
hist(data.frame(FST)$"X450")
FSTr=FST[-which(data.frame(FST)$"X450">0.75),]
FSTr=FSTr[-which(data.frame(FSTr)$"X1000"<0.05),]
hist(data.frame(FSTr)$"X450")
# write.csv(data.frame(FSTr),"NEON_SED_Spec_2020.csv",row.names=FALSE)


#----Data Exploration----
plot(FSTr, lwd = 0.25, lty = 1, col = "grey50", main="Spectra and 95% Quantile")
plot_quantile(FSTr, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(FSTr, regions = default_spec_regions(), add = TRUE)
plot_interactive(FSTr)

#----Remove Water Bands----