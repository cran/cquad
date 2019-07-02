quasi_sym_equ <- function(eta,s,y0=NULL){

# Compute quasi-symmetric function and its derivatives for given parameters eta
# and total equal to s
# if dyn then the last element of is gamma, in this case y0 must be precised
# preliminaries	

    TT = length(eta)
    
    ga = eta[TT]
    eta = eta[-TT]
    TT = TT-1
    uT1 = c(rep(0,TT),1)
    
                                        # initialization
    g0 = c(exp(ga*(y0==0)),0)
    g1 = c(0,exp(eta[1]+ga*(y0==1)))

    E0 = matrix(0,TT+1,2)
    E0[TT+1,1] = (y0==0)*exp(ga*(y0==0)) 
    E1 = matrix(0,TT+1,2)
    E1[1,2] = exp(eta[1]+(y0==1)*ga)
    E1[TT+1,2] = (y0==1)*exp(eta[1]+(y0==1)*ga)

    F0 = array(0,c(TT+1,TT+1,2))
    F0[TT+1,TT+1,1] = (y0==0)*exp(ga*(y0==0))

    F1 = array(0,c(TT+1,TT+1,2))
    F1[1,1,2] = exp(eta[1]+(y0==1)*ga)
    F1[1,TT+1,2] = F1[TT+1,1,2] = F1[TT+1,TT+1,2] = (y0==1)*exp(eta[1]+(y0==1)*ga)
    
    
    if(TT>1) for(t in 2:TT){
                 g00 = g0; g10 = g1
                 E00 = E0; E10 = E1
                 F00 = F0; F10 = F1
                                        # function g
                 g0 = c(g00*exp(ga)+g10,0)
                 g1 = c(0,g00*exp(eta[t])+g10*exp(eta[t]+ga))
                                        # first derivative
                 #browser()
                 E0 = cbind(E00*exp(ga)+E10,0)
                 E0[TT+1,-(ncol(E0))] = E0[TT+1,-(ncol(E0))] + g00*exp(ga)
                 

                 E1 = cbind(0,E00*exp(eta[t])+E10*exp(eta[t]+ga))
                 E1[t,-1] = E1[t,-1] + g00*exp(eta[t])+g10*exp(eta[t]+ga)
                 E1[TT+1,-1] = E1[TT+1,-1] + g10*exp(eta[t]+ga) 
                                        # second derivative
                 F0 = array(0,c(TT+1,TT+1,t+1))
                 F0[,,1:t] = F00*exp(ga)+F10
                 #browser()
                 F0[TT+1,,1:t] = F0[TT+1,,1:t] + E00[,1:t]*exp(ga)
                 
                 F0[,TT+1,1:t] = F0[,TT+1,1:t] + E00[,1:t]*exp(ga)

                  F0[TT+1,TT+1,1:t] = F0[TT+1,TT+1,1:t] + g00*exp(ga) 

                 
                 
                 F1 = array(0,c(TT+1,TT+1,t+1))
                 for(h in 1:t){
                     ut = rep(0,TT+1); ut[t] = 1
                     ut2 = ut+uT1
                     Tmp1 = E00[,h]%o%ut
                     Tmp2 = E10[,h]%o%ut2
                     F1[,,h+1] = F1[,,h+1]+(F00[,,h]+Tmp1+t(Tmp1)+g00[h]*(ut%o%ut))*exp(eta[t])+
                         (F10[,,h]+Tmp2+t(Tmp2)+g10[h]*(ut2%o%ut2))*exp(eta[t]+ga)
                 }
             }
                                        # output
    
    f = g0+g1; f = f[s+1]
    D1 = E0+E1; d1 = D1[,s+1]
    D2 = F0+F1; D2 = D2[,,s+1]
    lf = log(f); dl1 = d1/f; Dl2 = D2/f-dl1%o%dl1
    
    
    out = list(f=as.vector(f),d1=as.vector(d1),D2=D2,lf=as.vector(lf),dl1=as.vector(dl1),Dl2=Dl2)
    return(out)
    
}
