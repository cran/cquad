cquad <- function(formula, data, index=NULL, model=c("basic","equal","extended","pseudo"),
                  w=rep(1,n), dyn=FALSE, Ttol=10){

# INTERFACE FOR CQUAD ACCEPTING A FORMULA IN INPUT
#
# formula : formula with the same syntax as in plm
# data    : data.frame or pdata.frame
# index   : to denote panel structure as in plm
# model   : type of model = "basic", "equal", "extended", "pseudo"
# w       : vector of weights (not used for pseudo version)
# dyn     : for dynamic version (only for basic version)

  # preliminaries
	model = match.arg(model)
	if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index)
	if (!inherits(formula, "Formula")) formula <- Formula(formula)
	data <- model.frame(data, formula)
	X = model.matrix(data)[,-1]
	yv = pmodel.response(data, formula)
	id = attr(data,"index")[,1] 
	n = length(unique(id))
  # call model
  if(model=="basic") out = cquad_basic(id,yv,X,w=w,dyn=dyn,Ttol=Ttol)
  if(model=="equal") out = cquad_equ(id,yv,X,w=w,Ttol=Ttol)
  if(model=="extended") out = cquad_ext(id,yv,X,w=w,Ttol=Ttol)
  if(model=="pseudo") out = cquad_pseudo(id,yv,X,Ttol=Ttol)		

  out$formula = formula
  # return output
  out
}
