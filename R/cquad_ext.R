cquad_ext <-
function(id, yv, X=NULL, be=NULL, w = rep(1,n)){

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
	pid = id
	r = length(pid)
	label = unique(pid)
	n = length(label)
	if(is.null(X)) k=0 else{X = as.matrix(X); k = ncol(X)}
	c = max(yv)+1 # number of categories
	if(k>0) Xv = X
# chech variabiliy of the covariates
	if(k>0) for(j in 1:k){
		flag = TRUE
		for(i in 1:n){
			il = label[i] 
			if(max(X[pid==il,j])-min(X[pid==il,j])>0) flag = FALSE
		}
		if(flag) stop("at least one covariate without variability within unit")
	}
# variable names	
	varnames = NULL
	if(k>0){
		if(is.null(colnames(X))) for(j in 1:k) varnames = c(varnames,paste("X",j,sep="")) 
		else varnames = colnames(X)
	}
	varnames = c(varnames,"int",varnames,"y_lag")	
#  starting values
	be = rep(0,2*k+2)   
# check for balanced data
	Tv = rep(0,n)
	ind = id
	for(i in 1:n) Tv[i] = sum(ind==label[i])
	balanced = all(Tv==max(Tv))
	Tv = Tv-1
  	if(balanced) cat("Balanced panel data\n") else cat("Unbalanced panel data\n") 
# iterate until convergence    
	Sc = matrix(0,n,1)
	it = 0; lk = -Inf; lk0 = -Inf
	cat("------------|-------------|-------------|\n")
	cat("  iteration |      lk     |    lk-lko   |\n")
	cat("------------|-------------|-------------|\n")
	while(abs(lk-lk0)>10^-6 | it==0){
		it = it+1; lk0 = lk
		scv = matrix(0,n,2*k+2)
		lk = 0; J = 0
		for(cut_point in 1:(c-1)){
			yd = 1*(yv>(cut_point-1))
			for(i in 1:n){
				if(Tv[i]>1){
					il = label[i]
					y_i = yd[pid==il]
					y_i0 = y_i[1]; y_i = y_i[-1]
					sui = sum(y_i)
					if(sui>0 & sui<Tv[i]){
						Z = sq(Tv[i],sui)
						if(k==0) x_i = NULL else x_i = as.matrix(Xv[pid==il,])
						if(k>0) x_i = as.matrix(x_i[-1,])
						if(Tv[i]==2){
							x_i = rbind(c(x_i[1,],0,rep(0,k),0),
									    c(x_i[2,],1,x_i[Tv[i],],0),
									    c(rep(0,k),0,rep(0,k),1))
							Z = cbind(Z,y_i0*Z[,1]+Z[,1]*Z[,2])
						}else{
							x_i = rbind(cbind(x_i[-Tv[i],],rep(0,Tv[i]-1),matrix(0,Tv[i]-1,k),rep(0,Tv[i]-1)),
									    c(x_i[Tv[i],], 1, x_i[Tv[i],], 0),
									    c(rep(0,k), 0, rep(0,k), 1))             	
               							Z = cbind(Z,y_i0*Z[,1]+rowSums(Z[,1:Tv[i]-1]*Z[,2:Tv[i]]))
           				}
           				xb = x_i%*%be
           				den = exp(Z%*%xb)
           				sden = sum(den)
           				y_i = c(y_i,y_i0*y_i[1]+sum(y_i[1:Tv[i]-1]*y_i[2:Tv[i]]))
           				pc_i = exp(y_i%*%xb)/sden
       					Zt = t(Z)
       					lk = lk+w[i]*log(pc_i)
       					pp_i = as.vector(den/sden)
       					e_i = Zt%*%pp_i
       					scv[i,] = scv[i,]+w[i]*(t(y_i-e_i)%*%x_i)
       					V_i = Zt%*%diag(pp_i)%*%Z-e_i%*%t(e_i)
       					J = J-w[i]*(t(x_i)%*%V_i%*%x_i)
       				}
       			}
       		}
		}
		sc = colSums(scv)
		iJ = solve(J)
		be = be-iJ%*%sc
   		cat(sprintf("%11g", c(it,lk,lk-lk0)), "\n", sep = " | ")
	}
	cat("------------|-------------|-------------|\n")
	be = as.vector(be)
	Va = iJ%*%(t(scv)%*%scv)%*%iJ
	se = sqrt(diag(Va))
   	lk = as.vector(lk)
   	names(be) = varnames   	
   	colnames(scv) = varnames
   	rownames(J) = colnames(J) = varnames
   	names(se) = varnames   		
	out = list(lk=lk,be=be,scv=scv,J=J,se=se,Tv=Tv,call=match.call())
	class(out) = "cquad"	
	return(out)

}