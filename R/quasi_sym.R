quasi_sym <- function(eta,s,dyn=FALSE,y0=NULL){

# Compute quasi-symmetric function and its derivatives for given parameters eta
# and total equal to s
# if dyn then the last element of is gamma, in this case y0 must be precised
# preliminaries	

    TT = length(eta)
	if(dyn){
		ga = eta[TT]
		eta = eta[-TT]
		TT = TT-1
		uT1 = c(rep(0,TT),1)
	}
# initialization
	if(dyn){
		g0 = c(1,0)
		g1 = c(0,exp(eta[1]+y0*ga))
		E0 = matrix(0,TT+1,2)
		E1 = matrix(0,TT+1,2)
		E1[1,2] = exp(eta[1]+y0*ga)
		E1[TT+1,2] = y0*exp(eta[1]+y0*ga)
		F0 = array(0,c(TT+1,TT+1,2))
		F1 = array(0,c(TT+1,TT+1,2))
		F1[1,1,2] = exp(eta[1]+y0*ga)
		F1[1,TT+1,2] = F1[TT+1,1,2] = F1[TT+1,TT+1,2] = y0*exp(eta[1]+y0*ga)
	}else{
		f = c(1,exp(eta[1]))
		D1 = t(c(0,exp(eta[1])))
		D2 = array(D1,c(1,1,2))
	}
	if(TT>1) for(t in 2:TT){
		if(dyn){
			g00 = g0; g10 = g1
			E00 = E0; E10 = E1
			F00 = F0; F10 = F1
# function g
			g0 = c(g00+g10,0)
			g1 = c(0,g00*exp(eta[t])+g10*exp(eta[t]+ga))
# first derivative
			E0 = cbind(E00+E10,0)
			E1 = cbind(0,E00*exp(eta[t])+E10*exp(eta[t]+ga))
			E1[t,-1] = E1[t,-1] + g00*exp(eta[t])+g10*exp(eta[t]+ga)
			E1[TT+1,-1] = E1[TT+1,-1] + g10*exp(eta[t]+ga) 
# second derivative
			F0 = array(0,c(TT+1,TT+1,t+1))
			F0[,,1:t] = F00+F10
			F1 = array(0,c(TT+1,TT+1,t+1))
			for(h in 1:t){
				ut = rep(0,TT+1); ut[t] = 1
				ut2 = ut+uT1
				Tmp1 = E00[,h]%o%ut
				Tmp2 = E10[,h]%o%ut2
				F1[,,h+1] = F1[,,h+1]+(F00[,,h]+Tmp1+t(Tmp1)+g00[h]*(ut%o%ut))*exp(eta[t])+
				                      (F10[,,h]+Tmp2+t(Tmp2)+g10[h]*(ut2%o%ut2))*exp(eta[t]+ga)
			}
		}else{
			f0 = f; D10 = D1; D20 = D2
# function f
			f = c(f,0)+c(0,f)*exp(eta[t])
# first derivative
			st = min(s,t)
			s1 = 1:(t-1); s1 = s1[s>=s1]
			ee = exp(eta[t])
			D1 = matrix(0,t,st+1)
			if(length(s1)>0){
				D1[1:(t-1),s1+1] = D10[,s1+1]+D10[,s1]*ee
				D1[t,s1+1] = f0[s1]*ee
			}
			if(s>=t){
				D1[1:(t-1),t+1] = D10[,t]*ee
				D1[t,t+1] = f0[t]*ee
			}
# second derivative
			D2 = array(0,c(t,t,st+1))
			if(length(s1)>0){
				D2[1:(t-1),1:(t-1),s1+1] = D20[,,s1+1]+D20[,,s1]*ee
				D2[1:(t-1),t,s1+1] = D2[t,1:(t-1),s1+1] = D10[,s1]*ee
				D2[t,t,s1+1] = f0[s1]*ee
			}
			if(s>=t){
				D2[1:(t-1),1:(t-1),t+1] = D20[,,t]*ee
				D2[1:(t-1),t,t+1] = D2[t,1:(t-1),t+1] = D10[,t]*ee
				D2[t,t,t+1] = f0[t]*ee
			}
		}
	}
# output
	if(dyn){
		f = g0+g1; f = f[s+1]
		D1 = E0+E1; d1 = D1[,s+1]
		D2 = F0+F1; D2 = D2[,,s+1]
		lf = log(f); dl1 = d1/f; Dl2 = D2/f-dl1%o%dl1
	}else{
		f = f[s+1]; d1 = D1[,s+1]; D2 = D2[,,s+1]
		lf = log(f); dl1 = d1/f; Dl2 = D2/f-dl1%o%dl1
	}
	out = list(f=as.vector(f),d1=as.vector(d1),D2=D2,lf=as.vector(lf),dl1=as.vector(dl1),Dl2=Dl2)
	return(out)
	
}
