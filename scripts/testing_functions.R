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
  return(dgamma(alpha,args$shape,rate=args$rate, log=log))
}
shape.prior = function(shape,S,log=T){
  return(dnorm(shape,mean=0,sd=S,log=log))
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
    shape.prior(shape,log=T,S=pars$shape$S)
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
pl_igp.mcmc = function(n.iter, init, par_cov.init, X, prior.pars, u.type='c', dist.type=1,H=200, fix.u = NULL, burn.in=0.4, show=T){
  if(!is.null(fix.u)){
    init$u = fix.u
  }
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
    if(!is.null(fix.u)){
      state.prop[1] = fix.u
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
  acc.states = acc.states[-(1:burn.in*nrow(acc.states)),]
  if(show){
    plot.mcmc.pl_igp(acc.states, prior.pars, u.type, X, dist.type)
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
p2p.in = degree(barabasi.game(1e5,directed = T),mode='in')
p2p.in = p2p.in[p2p.in>0]

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

mcmc.out = pl_igp.mcmc(1e4,init, par_cov.init, p2p.in, prior.pars, u.type, dist.type, H=2e2, fix.u=15)
plot(mcmc.out$alpha,type='l')
plot(mcmc.out$shape,type='l')
plot(mcmc.out$scale,type='l')
plot(mcmc.out$u, type='l')
mcmc.out.burned = mcmc.out[-(1:250),]


q = sort(unique(p2p.in))

means = apply(mcmc.out.burned, 2, mean)


dens.vec.df = data.frame(matrix(nrow=nrow(mcmc.out.burned), ncol=length(q)))
for(i in 1:nrow(mcmc.out.burned)){
  dens.vec.df[i,] = pzigpd(q,sum(p2p.in>mcmc.out.burned[i,1])/length(p2p.in), mcmc.out.burned[i,2] ,mcmc.out.burned[i,1],mcmc.out.burned[i,3],mcmc.out.burned[i,4])
}


dens.mean = apply(dens.vec.df, 2, mean)
dens.upper = apply(dens.vec.df, 2, quantile, probs=0.95)
dens.lower = apply(dens.vec.df, 2, quantile, probs=0.05)
plot(ecdf2(p2p.in), log='xy', type='l')
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


p2p.in = degree(barabasi.game(1e5,directed = T),mode='in')
p2p.in = p2p.in[p2p.in>0]
plot(ecdf2(p2p.in), type='l', log='xy')

# -------------------------------------------------------------------------
u.a = 30
shape.a=1
scale.a=70
alpha.a=1
phi.a = 0.01
x = rpl_igp(5e3, alpha.a, phi.a, u.a, shape.a, scale.a)
plot(ecdf2(x), type='l', log='xy')
lines(1:3e3,pzigpd(1:3e3, phi.a, alpha.a,u.a, shape.a, scale.a),lty=2,lwd=2)


# rewriting ---------------------------------------------------------------


pgpd = function(x, shape, scale, threshold=0, lower.tail=T){
  if(!lower.tail){
    return(1-pgpd(x,shape,scale,threshold))
  }
  if(length(x)>1){
    return(pgpd.vec(x,shape,scale,threshold, lower.tail))
  }else{
    if(scale<0){
      stop('error')
    }
    if(threshold<0){
      stop('error')
    }
    if(x<=threshold){
      return(0)
    }
    if(shape<0 & x>=threshold - (scale+shape*threshold)/shape){
      return(1)
    }
    if(shape==0){
      return(1-exp(-(x-threshold)/scale))
    }else{
      A = 1  + ((shape*(x-threshold))/(scale+shape*threshold))
      if(A > 0){
        return(1 - max(0,A^(-1/shape)))
      }else{
        return(0)
      }
    }
  }
}
pgpd.vec = Vectorize(pgpd, vectorize.args = 'x')

dpl = function(x, alpha, upper.lim = NULL){
  if(length(x)>1){
    return(dpl.vec(x,alpha, upper.lim))
  }
  if(is.null(upper.lim)){
    if(alpha<=0){
      stop('ERROR alpha<=0')
    }else{
      return(x^-(alpha+1) / gsl::zeta(alpha+1))
    }
  }else{
    if(x>upper.lim){
      return(0)
    }else{
      return(x^-(alpha+1) / sum((1:upper.lim)^-(alpha+1)))
    }
  }
}
dpl.vec = Vectorize(dpl, vectorize.args = 'x')

ppl = function(x,alpha, upper.lim=NULL, lower.tail=T){
  if(length(x)>1){
    return(ppl.vec(x,alpha,upper.lim, lower.tail))
  }
  if(!lower.tail){
    return(1-ppl(x,alpha, upper.lim))
  }
  if(is.null(upper.lim)){
    return(1 - gsl::hzeta(alpha+1, x+1) / gsl::zeta(alpha+1))
  }else{
    if(x>upper.lim){
      return(1)
    }
    
    return(sum((1:x)^-(alpha+1))/ sum((1:upper.lim)^-(alpha+1)))
  
      
  }
}
ppl.vec = Vectorize(ppl, vectorize.args = 'x')

digp = function(x,shape,scale, threshold=0, type=1){
  #type:1 - modelling floor(H)
  #type:2 - modelling ceiling(H)
  u=threshold
  if(length(x)>1){
    return(digp.vec(x,shape,scale,threshold,type))
  }
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  if(x<=v){
    return(0)
  }
  return(pgpd(x-(v-ceiling(u)),shape,scale, threshold = v) - pgpd(x-1-(v-ceiling(u)), shape, scale, threshold = v))
}
digp.vec = Vectorize(digp, vectorize.args = 'x')

pigp = function(x,shape,scale,threshold, type=1, lower.tail=T){
  u = threshold
  if(length(x)>1){
    return(pigp.vec(x,shape,scale,threshold,type,lower.tail))
  }
  if(!lower.tail){
    return(1-pigp(x,shape,scale,threshold,type))
  }
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  if(x<=v){
    return(0)
  }
  return(pgpd(x-(v-ceiling(threshold)), shape, scale, threshold = v) - pgpd(ceiling(threshold), shape, scale, threshold=v))
}
pigp.vec = Vectorize(pigp, vectorize.args = 'x')
dpli = function(x,phi,alpha,shape,scale,threshold,type=1){
  return((1-phi)*dpl(x,alpha,upper.lim=threshold) + phi*digp(x,shape,scale,threshold,type))
}
ppli = function(x,phi,alpha,shape,scale,threshold,type=1, lower.tail=T){
  if(!lower.tail){
    return(1-ppli(x,phi,alpha,shape,scale,threshold,type))
  }
  return((1-phi)*ppl(x,alpha,upper.lim=threshold) + phi*pigp(x,shape,scale,threshold,type))
}

llpli = function(x,alpha,shape,scale,threshold,type=1){
  u=threshold
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  N = length(x)
  n = sum(x<=v)
  phi = (N-n)/N
  if(n<(N-1)){
    return(n*log(1-phi) + (N-n)*log(phi)- n*log(sum((1:v)^-(alpha+1))) - (alpha+1)*sum(log(x[x<=v]))+
             sum(log(digp(x[x>v],shape,scale,threshold,type))))
  }else{
    return(-Inf)
  }
}

joint.prior = function(phi,alpha,shape,scale,params){
  return(dbeta(phi,params$phi[1], params$phi[2], log = T) +
           dnorm(alpha,0,sd=params$alpha,log=T) +
           dnorm(shape,0,sd=params$shape, log=T)+
           dgamma(scale, 1, rate=params$scale, log=T))
}

joint.posterior = function(x,alpha,shape,scale,threshold,params,type=1){
  u=threshold
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  phi = sum(X[X[,1]<=v,2] / sum(X[,2]))
  return(llpli2(x,alpha,shape,scale,threshold,type) + joint.prior(phi,alpha,shape,scale,params))
}


cov.prop = function(data, prob = 1){
  if(runif(1)<prob){
    return(2.38^2 * cov(data) / ncol(data))
  }else{
    if(ncol(data)==1){
      return(0.01)
    }
    return(diag(rep(0.01, ncol(data))))
  }
}






# -------------------------------------------------------------------------

pli.mcmc = function(n.iter, x, init, prior.params, init.cov, type=1, burn.in=0, cov.period=1e4, fix.u = NULL, backup_period = 200, show=F){
  alpha.acc = data.frame(alpha=init$alpha)
  shapescale.acc = data.frame(shape=init$shape, scale=init$scale)
  thresh.acc = data.frame(threshold=init$threshold)
  alpha.cov = init.cov$alpha
  shapescale.cov = init.cov$shapescale
  thresh.cov = init.cov$thresh
  post.vals = numeric(n.iter+1)
  states = data.frame(alpha = numeric(n.iter+1), shape = numeric(n.iter+1), scale=numeric(n.iter+1), threshold = numeric(n.iter+1))
  states[1,] = init
  post.vals[1] = joint.posterior(x,init$alpha,init$shape,init$scale,init$threshold,prior.params,type)
  thresh.props = numeric(n.iter)
  if(show){
    dev.new()
    plot.new()
    plot.window(c(0,n.iter), c(0,60))
    axis(side=1)
    axis(side=2)
  }
   
  
  x.max = max(x[,1])
  for(i in 1:n.iter+1){
    message('iteration: ',i)
    #alpha steps
    
    if(nrow(alpha.acc)%%cov.period == 0){
      alpha.cov = cov.prop(tail(alpha.acc,cov.period))
    }
    alpha.prop = rnorm(1,tail(alpha.acc,1)[[1]],alpha.cov)
    if(alpha.prop!=0){
      logA = min(0, joint.posterior(x,alpha.prop,tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type)-
                   post.vals[i-1])
      if(log(runif(1))<logA){
        alpha.acc = rbind(alpha.acc,alpha.prop)
      }
    }
    #threshold steps
    
    if(!is.null(fix.u)){
      thresh.prop = fix.u
      thresh.acc = rbind(thresh.acc, thresh.prop)
    }else{
      if(nrow(thresh.acc)%%cov.period==0){
        thresh.cov = cov.prop(tail(thresh.acc,cov.period))
      }
      thresh.prop = rnorm(1,tail(thresh.acc,1)[[1]],thresh.cov)
      if(tail(shapescale.acc$shape,1)[[1]] > - (tail(shapescale.acc$scale,1)[[1]]/ (floor(thresh.prop)-1))){
        logA = min(0,joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],thresh.prop,prior.params,type)-
                     post.vals[i-1])
        if(log(runif(1))<logA){
          thresh.acc = rbind(thresh.acc,thresh.prop)
        }
      }
    }
    thresh.props[i]  =thresh.prop
    #shape and scale steps
    if(nrow(shapescale.acc)%%cov.period==0){
      shapescale.cov = cov.prop(tail(shapescale.acc,cov.period))
    }
    shapescale.prop = rmvnorm(1,unlist(tail(shapescale.acc,1)), shapescale.cov)
    
    if(shapescale.prop[2]>0 & shapescale.prop[1]>(-shapescale.prop[2]/floor(tail(thresh.acc,1)))){
      
      logA = min(0,joint.posterior(x,tail(alpha.acc,1)[[1]],shapescale.prop[1],shapescale.prop[2],tail(thresh.acc,1)[[1]],prior.params,type)-
                   post.vals[i-1])
      if(log(runif(1))<logA){
        shapescale.acc = rbind(shapescale.acc, shapescale.prop)
      }
    }
    states[i,] = c(tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]])
    post.vals[i] = joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type)
    if(i%%backup_period==0){
      mcmc.out.temp <<- list(alpha=alpha.acc, shapescale=shapescale.acc, thresh=thresh.acc, post.vals = post.vals, states = states)
    }
    if(i>1 & show){
      dev.flush()
      plot(1:min(i,500),tail(states$thresh[1:i],500), type='l', main = cor(tail(states$shape[1:i],500), tail(states$scale[1:i],500)))
      abline(h = tail(X[,1],4), lty=2)
      }
  }
  return(list(alpha=alpha.acc, shapescale=shapescale.acc, thresh=thresh.acc, post.vals = post.vals, states = states))
}


# -------------------------------------------------------------------------
plot.pli.mcmc = function(X,mcmc.out.states, burn.in=0, thin.by = 0){
  Y=X
  Y[,2] = cumsum(X[,2])
  Y[,2] = Y[,2]/max(Y[,2])
  
  N = length(mcmc.out$states$alpha)
  mcmc.out$states = tail(mcmc.out$states, N-burn.in)[seq(1,N-burn.in, by=thin.by+1),]
  mcmc.out$states$phi = sapply(mcmc.out$states$threshold, FUN=function(x){return(sum(X[X[,1]>x,2])/sum(X[,2]))})
  k = unique(X[,1])
  
  pmf.df = matrix(nrow=length(mcmc.out$states$shape), ncol=length(k))
  surv.df = matrix(nrow=length(mcmc.out$states$shape), ncol=length(k))
  for(i in 1:nrow(surv.df)){
    print(i)
    surv.df[i,] = ppli(k,mcmc.out$states$phi[i], mcmc.out$states$alpha[i], 
                       mcmc.out$states$shape[i], mcmc.out$states$scale[i], mcmc.out$states$threshold[i], lower.tail=F)
    pmf.df[i,] = dpli(k,mcmc.out$states$phi[i], mcmc.out$states$alpha[i], 
                      mcmc.out$states$shape[i], mcmc.out$states$scale[i], mcmc.out$states$threshold[i])
    }
  
  surv.upper95 = apply(surv.df, 2, quantile, probs = 0.975, na.rm=T)
  surv.lower95 = apply(surv.df, 2, quantile, probs = 1-0.975, na.rm=T)

  pmf.upper95 = apply(pmf.df, 2, quantile, probs = 0.975, na.rm=T)
  pmf.lower95 = apply(pmf.df, 2, quantile, probs = 1-0.975, na.rm=T)
  
  #plots
  p1 = ggplot(data=mcmc.out.states) + geom_density(aes(x=alpha))
  p2 = ggplot(data=mcmc.out.states, aes(x=shape, y=scale)) + geom_density_2d()
  p3 = ggplot(data=mcmc.out.states) + geom_density(aes(x=threshold))
  p4=ggplot(data=Y) + geom_line(aes(x=x, y=1-Freq))+geom_ribbon(aes(x=x,ymin=surv.lower95, ymax=surv.upper95), alpha=0.3) +
    scale_x_log10() +scale_y_log10()+ylab('Surv')
  p5 = ggplot(data=X/sum(X[,2])) +geom_point(aes(x=x, y=Freq))+geom_ribbon(aes(x=x,ymin=pmf.lower95, ymax=pmf.upper95), alpha=0.3)+   scale_x_log10() +scale_y_log10()
  plot.list = list(p1,p2,p3,p4)
  return(marrangeGrob(plot.list, nrow = 2, ncol=2))
}
llpli2 = function(X,alpha,shape,scale,threshold,type=1){
  u=threshold
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  cc = as.numeric(X[,2])
  y = as.numeric(X[,1])
  m=length(y)
  #restrict u to be at most the third largest value
  if(v>y[m-2]){
    return(-Inf)
  }else{
    n=sum(cc[y<=v])
    N = sum(cc)
    phi = (N-n)/N
    return(n*log(1-phi) + (N-n)*log(phi) -(alpha+1)*sum(cc[y<=v]*log(y[y<=v])) -
             n*log(sum((1:v)^-(alpha+1))) + sum(cc[y>v]*log(digp(y[y>v],shape,scale,threshold,type))))
  }
} 




