spsvmSDR <- function(x, y, H, C, method, lambda=NULL, gamma=3.7, penalty="grSCAD", max.iter=100){
 if(as.character(method) == "pplsvm"){
   obj <- pplssvm(x, y, H, C, lambda, gamma, penalty, max.iter=100)
 }
 if(as.character(method) == "ppals"){
   obj <- ppalssvm(x, y, H, C, lambda, gamma, penalty, max.iter=100)
 }
 if(as.character(method) == "pplr"){
   obj <- pplr(x, y, H, C,  lambda, gamma, penalty, max.iter=100, tol=1.0e-4)
 }
 if(as.character(method) == "ppl2svm"){
   obj <- ppl2svm(x, y, H, C, lambda, gamma, penalty, max.iter=100)
 }
 if(as.character(method) == "wpplr"){
   obj <- wPPLR(x, y, H, C, lambda, gamma, penalty, max.iter=100, tol=1.0e-4)
 }
 return(obj)
}

