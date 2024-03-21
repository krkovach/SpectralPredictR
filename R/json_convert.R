# install.packages("rjson")
require(rjson)

HyTraits_model_convert=function(input,output){
    require(rjson)
    jsonthing=fromJSON(file=input)
    wavelengths=unlist(jsonthing[[1]][[2]]$wavelengths[[1]])
    coefftable=setNames(data.frame(matrix(ncol=length(wavelengths)+1,nrow=0)), c("intercept",wavelengths))
    for (i in 1:length(jsonthing[[1]][[2]]$coefficients)){
        j=unlist(jsonthing[[1]][[2]]$coefficients[[i]])
        j=c(1,j)
        coefftable=data.frame(mapply(c,coefftable, j, SIMPLIFY=FALSE))
        }  
    }