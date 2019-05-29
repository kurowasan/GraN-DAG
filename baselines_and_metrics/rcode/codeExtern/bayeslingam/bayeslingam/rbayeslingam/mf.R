# mf:
# returns a negated value of a function, along with negated gradient
# and negated hessian attributes

mf<-function(p,f,...) {
  val<-f(p,...)
  res<-(-1)*value(val)
  attr(res,"gradient")<-(-1)*attr(val,"gradient")
  attr(res,"hessian")<-(-1)*attr(val,"hessian")

  res
}

