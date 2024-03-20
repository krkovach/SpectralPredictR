setwd("/Users/kylekovach/Downloads")

# install.packages("rjson")
require(rjson)
jsonthing=fromJSON(file="Chemistry_Ground_Models_mean-min.json")

wavelengths=unlist(jsonthing[[1]][[2]]$wavelengths[[1]])

coefftable=setNames(data.frame(matrix(ncol=length(wavelengths)+1,nrow=0)), c("intercept",wavelengths))

coefftable