source('functions.R')

df = load_data('rpkg_20190129.csv')
df0 = df-1
dat = df0$imports[df0$imports>0]
# -------------------------------------------------------------------------

#distribution functions
dpl_igp = function(x,u,phi,alpha,shape,scale,type=1 ,log=T){
  if(length(x)>1){
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
  if(length(q)>1){
    return(ppl_igp.vec(q,u,phi,alpha,shape,scale,type,log,lower.tail))
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
  U = 0:20
  for(i in 1:n.iter+1){
    message('iter: ',i,'| accepted: ',nrow(acc.states), '| u: ',round(tail(acc.states$u,1),2))
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
U = 0:100
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
u.type='c'
dist.type=2
par_cov.init = matrix(diag(rep(0.01,4)),nrow=4)

mcmc.out = pl_igp.mcmc(2e4,init, par_cov.init, dat, prior.pars, u.type, dist.type, H=200)
plot(mcmc.out$alpha,type='l')
plot(mcmc.out$shape,type='l')
plot(mcmc.out$scale,type='l')
plot(mcmc.out$u, type='l')



# -------------------------------------------------------------------------
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(mcmc.out$u, freq=F, col=c1)
hist(mcmc.out1$u, freq=F,add=T, col=c2)

# -------------------------------------------------------------------------
library(igraph)
sim.degs = degree(barabasi.game(1e5, power=1, m=NULL, zero.appeal = 1), mode='in')
sim.degs=sim.degs[sim.degs>0]

plot(epdf2(sim.degs), log='xy')
plot(epdf2(dat), log='xy')


plot(ecdf2(sim.degs), log='xy', type='l')
plot(ecdf2(dat), log='xy', type='l')



dat.alpha = fit_power_law(dat,xmin=1)$alpha



# -------------------------------------------------------------------------

xmin = 1
k=1:1e3
sim.alpha = fit_power_law(sim.degs[sim.degs>xmin] - xmin,xmin=1)$alpha
plot(ecdf2(sim.degs[sim.degs>xmin] - xmin), log='xy', type='l')
lines(k,pzeta(k,sim.alpha, lower.tail = F))



# -------------------------------------------------------------------------

sim.dat = rzeta(1e6,3)
sim.dat=sim.dat[sim.dat>0]
plot(epdf2(sim.dat), log='xy')
plot(ecdf2(sim.dat), log='xy')
# -------------------------------------------------------------------------

out = zeta.mcmc(1e4, sim.dat, 1,c(1,0.01))
# -------------------------------------------------------------------------

p2p = read.csv('../data/p2pusers.txt',sep = '\t')
p2p.in = table(p2p[,2])
p2p.in = dat

U = 0:100
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
  u=19,
  alpha=2,
  shape=5,
  scale=5
)
u.type='d'
dist.type=1
par_cov.init = matrix(diag(rep(0.005,4)),nrow=4)

mcmc.out = pl_igp.mcmc(1e4,init, par_cov.init, p2p.in, prior.pars, u.type, dist.type, H=1e4)
plot(mcmc.out$alpha,type='l')
plot(mcmc.out$shape,type='l')
plot(mcmc.out$scale,type='l')
plot(mcmc.out$u, type='l')
mcmc.out.burned = mcmc.out[-(1:500),]


q = unique(dat)

means = apply(mcmc.out.burned, 2, mean)


dens.vec.df = data.frame(matrix(nrow=nrow(mcmc.out.burned), ncol=length(q)))
for(i in 1:nrow(mcmc.out.burned)){
  dens.vec.df[i,] = pzigpd(q,sum(p2p.in>mcmc.out.burned[i,1])/length(p2p.in), mcmc.out.burned[i,2] ,mcmc.out.burned[i,1],mcmc.out.burned[i,3],mcmc.out.burned[i,4])
}


dens.mean = apply(dens.vec.df, 2, mean)
dens.upper = apply(dens.vec.df, 2, quantile, probs=0.95)
dens.lower = apply(dens.vec.df, 2, quantile, probs=0.05)
plot(ecdf2(p2p.in), log='xy')
lines(q,dens.upper, lty=1, col='blue')
lines(q,dens.mean, lty=1, col='red')
lines(q,dens.lower, lty=1, col='blue')




# -------------------------------------------------------------------------
u = 33
alpha = 3
shape= 2
scale=13

sim.dat = sample(1:1e4,2e3,  dpl_igp.vec(1:1e4, u, 0.02, alpha, shape, scale, log=F), replace=T)

plot(ecdf2(sim.dat), log='xy')


# -------------------------------------------------------------------------


p2p.in = sim.dat

U = 0:100
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
  u=20,
  alpha=2,
  shape=10,
  scale=15
)
u.type='c'
dist.type=1
par_cov.init = matrix(diag(c(rep(0.1,3),10)),nrow=4)

mcmc.out = pl_igp.mcmc(1e5,init, par_cov.init, p2p.in, prior.pars, u.type, dist.type, H=2e3)
plot(mcmc.out$alpha,type='l')
plot(mcmc.out$shape,type='l')
plot(mcmc.out$scale,type='l')
plot(mcmc.out$u, type='l')
mcmc.out.burned = mcmc.out[-(1:500),]


q = sort(unique(p2p.in), decreasing = F)
# q=1:max(p2p.in)

means = apply(mcmc.out.burned, 2, mean)
meds = apply(mcmc.out.burned,2,median)

dens.vec.df = data.frame(matrix(nrow=nrow(mcmc.out.burned), ncol=length(q)))
for(i in 1:nrow(mcmc.out.burned)){
  message(i)
  dens.vec.df[i,] = pzigpd(q,sum(p2p.in>mcmc.out.burned[i,1])/length(p2p.in), mcmc.out.burned[i,2] ,mcmc.out.burned[i,1],mcmc.out.burned[i,3],mcmc.out.burned[i,4])
}


dens.mean = apply(dens.vec.df, 2, mean)
dens.upper = apply(dens.vec.df, 2, quantile, probs=0.95)
dens.lower = apply(dens.vec.df, 2, quantile, probs=0.05)




par(mar=c(4,4,4,4))
plot(ecdf2(p2p.in), log='xy', type='l')
# polygon(c(q,rev(q)), c(dens.upper,rev(dens.lower)), col=rgb(1, 0.6, 1,0.5), border=NA)
lines(q,dens.upper, lty=2)
lines(q,dens.lower,lty=2)

# u = 50
# alpha = 1
# shape=0.1
# scale=10
# 
# lines(q,ppl_igp.vec(q, u, 0.05, alpha, shape, scale, log=F), lty=2)









