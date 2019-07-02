quasi_sym_pseudo <- function(eta,qi,s,y0=NULL){
    
                                        # Compute quasi-symmetric function and its derivatives for given parameters eta
                                        # and total equal to s
                                        # if dyn then the last element of is gamma, in this case y0 must be precised
                                        # preliminaries
    TT = length(eta)
    
    ga = eta[TT]
    nu = qi*ga
    eta = eta[-TT]
    TT = TT-1
    uT1 = c(rep(0,TT),1)
    
                                        # initialization
    g0 = c(exp(nu[1]*y0),0)
    g1 = c(0,exp(eta[1]+y0*(ga + nu[1])))
    
    E0 = matrix(0,TT+1,2)
    E0[TT+1,1] = y0*(exp(y0*nu[1]))*qi[1]
    
    E1 = matrix(0,TT+1,2)
    E1[1,2] = exp(eta[1]+y0*(ga + nu[1]))
    E1[TT+1,2] = exp(eta[1]+y0*(ga + nu[1]))*y0*(1+qi[1])
    
    F0 = array(0,c(TT+1,TT+1,2))
    F0[TT+1,TT+1,1] = y0*(exp(y0*nu[1]))*qi[1]^2
    
    F1 = array(0,c(TT+1,TT+1,2))
    F1[1,1,2] = exp(eta[1]+y0*(ga + nu[1]))
    F1[1,TT+1,2] = F1[TT+1,1,2] = exp(eta[1]+y0*(ga + nu[1]))*y0*(1+qi[1])
    F1[TT+1,TT+1,2] = exp(eta[1]+y0*(ga + nu[1]))*y0*(1+qi[1])^2




 if(TT>1) for(t in 2:TT){
                 
                 g00 = g0; g10 = g1
                 E00 = E0; E10 = E1
                 F00 = F0; F10 = F1
                 
                                        # function g
                 g0 = c(g00+g10*(exp(nu[t])),0)
                 g1 = c(0,g00*exp(eta[t])+g10*exp(eta[t]+ga+nu[t]))
                 
                                        # first derivative
                  E0 = cbind(E00+E10*exp(nu[t]),0)
                  E0[TT+1,-(t+1)] = E0[TT+1,-(t+1)] + g10*exp(nu[t])*(qi[t]) 
              
                 
                 F0 = array(0,c(TT+1,TT+1,t+1))
                 for(h in 1:t){
                 F0[,,h] = F00[,,h]+F10[,,h]*exp(nu[t])
                 F0[,TT+1,h] = F0[,TT+1,h]+E10[,h]*exp(nu[t])*qi[t]
                 F0[TT+1,,h] = F0[TT+1,,h]+E10[,h]*exp(nu[t])*qi[t]
                 F0[TT+1,TT+1,h] = F0[TT+1,TT+1,h]+g10[h]*exp(nu[t])*qi[t]^2
                 }
              


              
              E1 = cbind(0,E00*exp(eta[t])+E10*exp(eta[t]+ga+nu[t]))
              E1[t,-1] = E1[t,-1] + g00*exp(eta[t]) + g10*exp(eta[t]+ga+nu[t])
              E1[TT+1,-1] = E1[TT+1,-1] + g10*exp(eta[t]+ga+nu[t])*(1+qi[t])
              
              
              
              F1 =  array(0,c(TT+1,TT+1,t+1))
              for(h in 1:t){
                  
                  F1[,,h+1] = F00[,,h]*exp(eta[t])+F10[,,h]*exp(eta[t]+ga+nu[t])
                  F1[,t,h+1] = F1[,t,h+1]+E00[,h]*exp(eta[t])+E10[,h]*exp(eta[t]+ga+nu[t])
                  F1[,TT+1,h+1] = F1[,TT+1,h+1]+E10[,h]*exp(eta[t]+ga+nu[t])*(1+qi[t])
                  
                  
                  F1[t,,h+1] = F1[t,,h+1]+E00[,h]*exp(eta[t])+E10[,h]*exp(eta[t]+ga+nu[t])
                  
                  F1[t,t,h+1] = F1[t,t,h+1]+g00[h]*exp(eta[t])+g10[h]*exp(eta[t]+ga+nu[t])
                
                  F1[t,TT+1,h+1] = F1[t,TT+1,h+1]+g10[h]*exp(eta[t]+ga+nu[t])*(1+qi[t])
                  
                  
                  
                  F1[TT+1,,h+1] = F1[TT+1,,h+1]+E10[,h]*exp(eta[t]+ga+nu[t])*(1+qi[t])
                  
                  F1[TT+1,t,h+1] = F1[TT+1,t,h+1]+g10[h]*exp(eta[t]+ga+nu[t])*(1+qi[t])
                  
                  F1[TT+1,TT+1,h+1] = F1[TT+1,TT+1,h+1]+g10[h]*exp(eta[t]+ga+nu[t])*(1+qi[t])^2
                  
                  
              }
              
              
              
              
                                        # second derivative
          }
             
    
    
                                        # output
    
    f = g0+g1; f = f[s+1]
    D1 = E0+E1; d1 = D1[,s+1]
    D2 = F0+F1; D2 = D2[,,s+1]
    lf = log(f); dl1 = d1/f; Dl2 = D2/f-(dl1%o%dl1)
    
    
    out = list(f=as.vector(f),d1=as.vector(d1),D2=D2,lf=as.vector(lf),dl1=as.vector(dl1),Dl2=Dl2)
    return(out)
    
}
