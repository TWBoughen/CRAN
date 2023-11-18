source('functions.R')

df = load_data('rpkg_20190129.csv')
df0 = df-1
dat = df0$imports[df0$imports>0]
# -------------------------------------------------------------------------

#distribution functions
dpl_igp = function(x,u,phi,alpha,shape,scale,type=1 ,log=T){
  if(length(x>1)){
    return(dpl_igp.vec(x,u,phi,alpha,shape,scale,type,log))
  }
  if(type==1){
    u=floor(u)
  }else{
    u=ceiling(u)
  }
  if(x<=0){
    out=-Inf
  }
  if(x<=u & x>0){
    out = log(1-phi) - (alpha+1)*log(x) - log(gsl::zeta(alpha+1)-gsl::hzeta(alpha+1, u+1))
  }
  if(x>u){
    out = log(phi) + log(pgpd(x-u,scale=scale,shape=shape) - pgpd(x-u-1,scale=scale,shape=shape))
  }
  if(!log){
    out=exp(out)
  }
  return(out)
}
dpl_igp.vec = Vectorize(dpl_igp,vectorize.args = 'x')
ppl_igp = function(q,u,phi,alpha,shape,scale,type=1,log=T,lower.tail=F){
  if(length(x)>1){
    return(ppl_igp.vec(q,u,phu,alpha,shape,scale,type,log,lower.tail))
  }
  if(type==1){
    u=floor(u)
  }else{
    u=ceiling(u)
  }
  if(q<=0){
    out = 1
  }
  if(q<=u & q>0){
    out = log(1 - (1-phi)* (gsl::zeta(alpha+1)-gsl::hzeta(alpha+1, q+1)) / (gsl::zeta(alpha+1)-gsl::hzeta(alpha+1, u+1)) )
  }
  if(q>u){
    out = log(phi) + log(pgpd(q-u,scale=scale,shape=shape,lower.tail=F ))
  }
  if(lower.tail){
    out=1-out
  }
  if(!log){
    out=exp(out)
  }
  return(out)
}
ppl_igp.vec = Vectorize(ppl_igp,vectorize.args = 'q')

llpl_igp = function(X,u,phi,alpha,shape,scale,type=1){
  if(type==1){
    u=floor(u)
  }else{
    u=ceiling(u)
  }
  n = sum(X<=u)
  if(n==length(X)){
    return(n*log(1-phi) - n*log(gsl::zeta(alpha+1)-gsl::hzeta(alpha+1,u+1)) - (alpha+1)*sum(log(X)))
  }
  if(n==0){
    return(length(X)*log(phi) + sum(log(pgpd(X-u,scale=scale,shape=shape) - pgpd(X-u-1,scale=scale,shape=shape))))
  }
  N=length(X)
  X.pl = X[X<=u]
  X.igp = X[X>u]
  return(n*log(1-phi) - n*log(gsl::zeta(alpha+1)-gsl::hzeta(alpha+1,u+1)) -(alpha+1)*sum(log(X.pl)) +
           (N-n)*log(phi) + sum(log(pgpd(X.igp-u,scale=scale,shape=shape) - pgpd(X.igp-u-1,scale=scale,shape=shape))))
  
}

##priors
u_prior = function(u,type='c',log=T,...){
  args = list(...)
  if(type=='c'){
    out = dgamma(u,args$shape,rate = args$rate,log=log)
  }else if(type=='d'){
    out = -log(args$N)
    if(!log){
      return(exp(out))
    }
    return(out)
  }
}
alpha.prior = function(alpha,log=T,...){
  args = list(...)
  # print(class(alpha))
  return(dgamma(alpha,args$shape,rate=args$rate))
}
shape.prior = function(shape,log=T,...){
  args = c(...)
  return(dnorm(shape,mean=0,sd=args$S,log=log))
}
scale.prior = function(scale,log=T,...){
  args = list(...)
  return(dgamma(scale,args$shape,rate=args$rate))
}

#joint prior
pl_igp.prior = function(u,alpha,shape,scale,pars,u.type='c',log=T){
  
  out = u_prior(u,type=u.type,log=T,N=pars$u$N,shape=pars$u$shape,rate=pars$u$rate)+
    alpha.prior(alpha,log=T,shape=pars$alpha$shape,rate=pars$alpha$rate)+
    scale.prior(scale,log=T,shape=pars$scale$shape,rate=pars$scale$rate)+
    shape.prior(shape,log=T,pars$shape)
  if(!log){
    out=exp(out)
  }
  return(out)
}

###posterior

pl_igp.posterior = function(X,u,alpha,shape,scale,prior.pars,u.type='c',dist.type=1,log=T){
  out = pl_igp.prior(u,alpha,shape,scale,prior.pars,u.type,log=T)+
        llpl_igp(X,u,sum(X>u)/length(X),alpha,shape,scale,dist.type)
  if(!log){
    out=exp(out)
  }
  return(out)
}
pl_igp.posterior.u = Vectorize(pl_igp.posterior, vectorize.args = 'u')
pl_igp.posterior.u.scaled = function(X,u,alpha,shape,scale,prior.pars,u.type='c',dist.type=1){
  log_dens = pl_igp.posterior.u(X,u,alpha,shape,scale,prior.pars,u.type,dist.type)
  return(exp(log_dens-max(log_dens)))
}

#next, the mcmcs 
pl_igp.mcmc = function(n.iter, init, par_cov.init, X, prior.pars, u.type='c', dist.type=1,H=200){
  acc.states = data.frame(u=init$u, alpha=init$alpha, shape=init$shape, scale=init$scale)
  U = 0:100
  for(i in 1:n.iter+1){
    message('iter: ',i,'| accepted: ',nrow(acc.states), '| scale: ',round(mean(acc.states$u)))
    #proposal steps
    if(u.type=='c'){
      if(nrow(acc.states)>H){
        Sig.i = 2.38^2 / 4 * cov(tail(acc.states,H))
      }else{
        Sig.i = par_cov.init
      }
      state.prop = rmvn(1,mu=tail(acc.states,1), Sig.i)
    }else if(u.type=='d'){
      weights = pl_igp.posterior.u.scaled(X,U,tail(acc.states,1)[2][[1]],tail(acc.states,1)[3][[1]],tail(acc.states,1)[4][[1]],prior.pars,u.type='c',dist.type=1)
      u.prop = sample(U,1,prob=weights)
      if(nrow(acc.states)>H){
        Sig.i = 2.38^2 / 3 * cov(tail(acc.states,H)[,-1])
      }else{
        Sig.i = par_cov.init[-1,-1]
      }
      rest.prop = rmvn(1,mu=tail(acc.states,1)[-1], Sig.i)
      state.prop = c(u.prop, rest.prop)
    }
    #preliminary rejection
    if(any(state.prop[-3]<0)){
      next
    }
    #calculting log-acc prob
    A = min(0, pl_igp.posterior(X,state.prop[1], state.prop[2], state.prop[3],state.prop[4],prior.pars,u.type,dist.type)-
              pl_igp.posterior(X,tail(acc.states,1)[1][[1]],tail(acc.states,1)[2][[1]],tail(acc.states,1)[3][[1]],tail(acc.states,1)[4][[1]],prior.pars,u.type,dist.type))
    if(is.nan(A)){
      break
    }
    #accepting / rejecting
    if(log(runif(1))<A){
      acc.states[nrow(acc.states)+1, ] = state.prop
      next
    }else{
      next
    }
  }
  return(acc.states)
}


#post-processing



# -------------------------------------------------------------------------
U = 0:(length(dat)-1)
alpha=1.2
shape=0.9
scale=20
prior.pars = list(
  u=list(N=length(dat), shape=1, rate=0.01), 
  alpha=list(shape=1, rate=0.01), 
  shape=list(S=4,v=16),
  scale=list(shape=1,rate=0.01)
)
init = list(
  u=15,
  alpha=2,
  shape=1,
  scale=10
)
u.type='d'
dist.type=2
par_cov.init = matrix(diag(rep(0.01,4)),nrow=4)

mcmc.out = pl_igp.mcmc(2e4,init, par_cov.init, dat, prior.pars, u.type, dist.type, H=200)
plot(mcmc.out$alpha,type='l')
plot(mcmc.out$shape,type='l')
plot(mcmc.out$scale,type='l')
hist(mcmc.out$u)



# -------------------------------------------------------------------------
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(mcmc.out$u, freq=F, col=c1)
hist(mcmc.out1$u, freq=F,add=T, col=c2)





