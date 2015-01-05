summary.cquad <-
function(object, ...){
	
# preliminaries
	out = object
# print output	
	cat("\nCall:\n")
    print(out$call)
    cat("\nLog-likelihood:\n")
    print(round(out$lk,2))
    tstat = out$be/out$se
    pv = 2*(1-pnorm(abs(tstat)))
    Tab = cbind("est."=out$be,"s.e."=out$se,"t-stat"=tstat,"p-value"=pv)
    print(Tab)
    cat("\n")
    
}
