cquad_ext <- function(id, yv, X=NULL, be=NULL, w = rep(1,n),Ttol=10){


# Fit a Conditional Logit model using CML (with time varying number of observation)
# the first column of data is the increasing number of the unit with aggregated data
# Y       : vector of binary response variable of dimension rc x 1 (r=units, c=times)
# X       : matrix of covariates of dimension rc x k (k=number of covariates)
# beta    : vector of first guess of the regression parameters of dimension k
# w       : vector of weights for each subject
# betahat : vector of estimated parameters of dimension k
# se      : vector of estimated standard errors of dimension k
# lk      : log-likelihood value at convergence

                                        # preliminaries
    
    ## Sort data by id ###
    input_data = cbind(id,yv,X)
    sorted_data = input_data[order(input_data[,1],decreasing=FALSE),]
    
    id = sorted_data[,1]
    yv = sorted_data[,2]
    X =  sorted_data[,-(1:2)]
    #######################

    
    pid = id
    r = length(pid)
    label = unique(pid)
    n = length(label)
    if(is.null(X)) k=0 else{X = as.matrix(X); k = ncol(X)}
# prepare covariate matrix
    X1 = NULL
    if(k==0){
        for(i in label){
            Ti = sum(id==i)
            tmp = c(rep(0,Ti-1),1)
            X1 = rbind(X1,as.matrix(tmp))
        }
        colnames(X1) = "int"
    }else{
        X = as.matrix(X)
        if(is.null(colnames(X))) colnames(X) = paste("X",1:k,sep="")
        X1 = matrix(0,nrow(X),2*k+1)
        ind = 0
        for(i in label){
            Ti = sum(id==i)
            Xi = matrix(X[id==i,],Ti)
            ind = max(ind)+(1:Ti)
            Tmp = rbind(matrix(0,Ti-1,k+1),c(1,Xi[Ti,]))
            X1[ind,] = cbind(Xi,Tmp)
        }
        names = c(colnames(X),"int")
        for(j in 1:k) names = c(names,paste("diff-",names[j],sep=""))
        colnames(X1) = names
    }
    ind = which((apply(X1,2,max)-apply(X1,2,min))==0)
    if(length(ind)>0) X1 = X1[,-ind]
                                        # use cquad_basic
                                       
    out1 = cquad_basic(id,yv,X1,w=w,dyn=TRUE,Ttol=Ttol)
                                        # adjust output	
    out = list(formula=formula,lk=out1$lk,coefficients=out1$coefficients,vcov=out1$vcov,scv=out1$scv,J=out1$J,se=out1$se,
               ser=out1$ser,Tv=out1$Tv,call=match.call())
    class(out) = c("cquad","panelmodel")	
    return(out)
    
}
