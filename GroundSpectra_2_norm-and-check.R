#Kyle R Kovach
#ASD/SED/SIG Spectra Processing
#December 18, 2020

#----Load Libraries----
# install.packages("ggplot2")
# install.packages("reshape2")

require(ggplot2)
require(reshape2)

#----Initialize----
set.seed(112345)
setwd("F:/PostDoc/FieldSpec_SpecEvo/Output")

#----Read in Data----
fname=list.files(path=getwd(),pattern=glob2rx("*.csv"))
for (i in fname)
{
  j=paste0(substr(i,1,17))
  assign(j,read.csv(i,header=TRUE))
};rm(fname,i,j)

#----Review Spectra by Group----
sampledata=NEON_2020_fshspec
sampledata$Sample=paste(sampledata$Date,sampledata$Sample,sep=".")

tic=0;for (i in unique(sampledata$Sample)){
  a=subset(sampledata[,-c(1:2,5)],Sample==i)
  b=sqrt(apply(a[,3:2151]^2,1,sum))
  d=a[,3:2151]/b
  e=cbind(a[,1:2],d)
  f=melt(e, id.vars=c("Sample","Replicate"))
  f$Replicate=as.factor(f$Replicate)
  print(ggplot(f, aes(x=as.numeric(variable), y=value, group=Replicate,color=Replicate))+
    geom_line()+
    ggtitle(paste("Sample:",f$Sample[1]))+
    ylab("value")+
    theme_bw()+
    scale_x_continuous(name="band",
                       limits=c(0, 2151),
                       breaks=seq(150, 2151,by=500),
                       labels=seq(500,2500, by=500)))
  tic=tic+1;print(paste0(tic,"/",length(unique(sampledata$Sample))))
  readline(prompt="Press [enter] to continue")
};rm(a,b,d,e,f,i,sampledata,tic)