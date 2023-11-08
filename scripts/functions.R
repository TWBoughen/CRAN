#libraries
{
require(igraph)
require(evd)
require(ismev)
require(extRemes)
require(VGAM)
require(rjags)
require(gsl)
}
#functions
{
##data loading functions
{
rbind_ind = function(df1, df2){
  names1 = names(df1)
  names2 = names(df2)
  names(df2) = names1
  return(rbind(df1, df2))
}
create_graph = function(df) {
  abc = df
  abc = rbind_ind(abc[abc$reverse,c(2,1,3,4)],abc[!abc$reverse,c(1,2,3,4)])
  abc$reverse=F
  abc = abc[!duplicated(abc),]
  abc.g = graph.data.frame(abc)
  E(abc.g)$type = abc$type
  types <<- c('depends', 'suggests','imports','linking to')
  cols = c('darkred', 'orange', 'blue', 'green')
  E(abc.g)$color = cols[match(abc$type,types)]
  deg_size_scale = 1
  V(abc.g)$deg = degree(abc.g)
  V(abc.g)$label = V(abc.g)$name
  return(abc.g)
}
focus_type = function(G,t='depends') {
  abc.g=G
  abc.sg2 = delete.edges(abc.g, E(abc.g)[E(abc.g)$type!=t])
  return(abc.sg2)
}
load_data = function(name=NA,all=T, path='../data/'){
  all=is.na(name)
  if(all){
    df.list = list()
    for(i in 1:length(list.files(path))){
      df.list[[i]] = load_data(name = list.files(path)[i], path=path)
    }
    return(df.list)
  }else{
    df.raw = read.csv(paste0(path,name))
    df.reverse = df.raw[df.raw$reverse,c(2,1,3,4)]
    df = rbind(df.raw[!df.raw$reverse,], df.reverse)
    G = create_graph(df)
    IG = focus_type(G,'imports')
    DG = focus_type(G,'depends')
    degdf = data.frame(imports = degree(IG,mode='in')+1,depends = degree(DG,mode='in')+1)
    return(degdf)
  }
}
}
##misc functions
{
ecdf2 = function(x){
  e = epdf2(x)
  return(data.frame(x=e$x,p=1-cumsum(e$p)))
}
epdf2 = function(x){
  out = table(x)/length(x)
  return(data.frame(x = names(out),p=as.numeric(out)))
}
}
##zeta model functions
{
dzeta  = Vectorize(function(x,alpha){
  return(x^-(alpha+1) / gsl::zeta(alpha+1))
},vectorize.args = 'x')
pzeta = Vectorize(function(x,alpha,lower.tail=T){
  if(lower.tail){
    return(sum((1:x)^-(alpha+1))/gsl::zeta(alpha+1))
  }else{
    return(1-sum((1:x)^-(alpha+1))/gsl::zeta(alpha+1))
  }
},vectorize.args = 'x')
lzeta = function(X,alpha,log=T){
  if(log){
    return(-length(X)*log(gsl::zeta(alpha+1)) - (alpha+1)*sum(log(X)))
  }else{
    return(exp(-length(X)*log(gsl::zeta(alpha+1)) - (alpha+1)*sum(log(X))))
  }
}
alpha.prior = function(alpha,pars){
  return(dgamma(alpha,shape=pars[1],rate=pars[2],log=T))
}
zeta.posterior = function(X,alpha,alpha.pars){
  return(lzeta(X,alpha)+alpha.prior(alpha,alpha.pars))
}
zeta.mcmc = function(n.iter,data,alpha.init,alpha.pars,prop.var.init=0.01,H=200,show=T){
  results = numeric(n.iter+1)
  results[1] = alpha.init
  acc.vec = numeric(n.iter+1)
  acc.points = c()
  prop.var = prop.var.init
  var.vec = c()
  for(i in 1:n.iter+1){
    message(i)
    if(length(acc.points)>H){
      k=length(acc.points)
      prop.var = 2.4^2 * var(acc.points[(k-H+1):k])
    }
    
    a.prop = rnorm(1,mean=results[i-1],sd = sqrt(prop.var))
    var.vec = c(var.vec,prop.var)
    if(a.prop<0){
      results[i]=results[i-1]
      next
    }
    logA = min(0,zeta.posterior(data,a.prop,alpha.pars)-zeta.posterior(data,results[i-1],alpha.pars))
    
    if(log(runif(1))<logA){
      results[i] = a.prop
      acc.vec[i]=1
      acc.points = c(acc.points,a.prop)
      next
    }else{
      results[i] = results[i-1]
      next
    }
  }
  if(show){
    plot.zeta.mcmc(acc.points,data)
  }
  return(list(res = acc.points, prop.var = var.vec))
}
plot.zeta.mcmc = function(mcmc.out,dat){
  layout.mat = t(matrix(c(1,1,1,1,
                          2,2,3,3,
                          2,2,3,3),nrow=4))
  L = layout(layout.mat)
  par(mar=c(2,4,4,0))
  plot(mcmc.out,type='l')
  k=1:2e3
  conf.size=0.95
  plot(epdf2(dat),pch=16,cex=0.5,log='xy')
  lines(k,dzeta(k,quantile(mcmc.out,probs=1-conf.size)), col='blue')
  lines(k,dzeta(k,mean(mcmc.out)), col='red')
  lines(k,dzeta(k,quantile(mcmc.out,probs=conf.size)), col='blue')
  
  plot(ecdf2(dat), type='l', log='xy')
  lines(k,pzeta(k,quantile(mcmc.out,probs=1-conf.size),lower.tail = F), col='blue')
  lines(k,pzeta(k,mean(mcmc.out),lower.tail=F),col='red')
  lines(k,pzeta(k,quantile(mcmc.out,probs=conf.size),lower.tail = F), col='blue')

  
  par(mfrow=c(1,1))
}
}
##zeta-igpd model functions
{
dzigpd = Vectorize(function(x,alpha,u,shape,scale){
  if(x<=0){
    return(0)
  }
  if(x<=u){
    return(dzeta(x,alpha))
  }
  if(x>u){
    return(pzeta(u,alpha,lower.tail=F)*(pgpd(x,u,scale=scale,shape=shape)-pgpd(x-1,u,scale=scale,shape=shape)))
  }
},vectorize.args = 'x')
pzigpd = Vectorize(function(x,alpha,u,shape,scale,lower.tail=T){
  out=NA
  if(x<=u){
    out = pzeta(x,alpha,lower.tail=F)
  }
  if(x>u){
    out = pzeta(u,alpha,lower.tail=F)*pgpd(x,u,scale=scale,shape=shape, lower.tail=F)
  }
  if(lower.tail){
    return(1-out)
  }else{
    return(out)
  }
},vectorize.args = 'x')
lzigpd = function(X,alpha,u,shape,scale){
  n = length(X[X<=u])
  N = length(X)
  
  A = -n*log(gsl::zeta(alpha+1)) 
  B = - (alpha+1)*sum(log(X[X<=u]))
  CC = (N-n)*log(pzeta(u,alpha,lower.tail=F))
  DD = sum(log(pgpd(X[X>u],u,scale=scale,shape=shape)-pgpd(X[X>u]-1,u,scale=scale,shape=shape)))
  

  
  return(A+B+CC+DD)
}
shape.prior = function(shape,pars){
  return(dnorm(shape,mean=0,sd=pars[1],log=T))
}
scale.prior = function(scale,pars){
  return(dgamma(scale,shape=pars[1], rate=pars[2],log=T))
}
zigpd.prior = function(alpha,shape,scale,alpha.pars,shape.pars,scale.pars){
  return(alpha.prior(alpha,alpha.pars) + shape.prior(shape,shape.pars) + scale.prior(scale,scale.pars))
}
zigpd.posterior = function(X,alpha,u,shape,scale,alpha.pars,shape.pars,scale.pars){
  return(lzigpd(X,alpha,u,shape,scale) + zigpd.prior(alpha,shape,scale,alpha.pars,shape.pars,scale.pars))
}
zigpd.mcmc = function(n.iter,data,u,init,pri.pars,prop.cov = diag(c(0.01,0.01,0.01)),H=200, show=T){
  states = data.frame(alpha=numeric(n.iter+1),shape=numeric(n.iter+1),scale=numeric(n.iter+1))
  states[1,] = init
  cov.ls = list()
  cov.ls[[1]] = prop.cov
  acc.states = data.frame(alpha=numeric(),shape=numeric(),scale = numeric())
  state.prop=c(0,0,0)
  for(i in 1:n.iter+1){
    message('iter: ',i,' - accepted:',nrow(acc.states))
    # print(state.prop)
    # print(format(prop.cov))
    if(nrow(acc.states)>H){
      prop.cov =( 2.4^2 / 3 ) * cov(acc.states[(nrow(acc.states)-H+1) : nrow(acc.states),])
    }
    cov.ls[[i]] = prop.cov
    state.prop = rmvn(1,mu=states[i-1,],V=prop.cov)
    if(state.prop[1]<0){
      states[i,] = states[i-1,]
      next
    }
    if(state.prop[3]<0){
      states[i,] = states[i-1,]
      next
    }
    #acceptance
    logA = min(0,(zigpd.posterior(data,state.prop[1],u,state.prop[2],state.prop[3],pri.pars$alpha,pri.pars$shape,pri.pars$scale)-
                 zigpd.posterior(data,states[i-1,1],u,states[i-1,2],states[i-1,3],pri.pars$alpha,pri.pars$shape,pri.pars$scale)))
    
    if(log(runif(1))<logA){
      states[i,] = state.prop
      acc.states = rbind(acc.states,state.prop)
      next
    }else{
      states[i,] = states[i-1,]
      next
    }
  }
  colnames(acc.states) = c('alpha','shape','scale')
  if(show){
    plot.zigpd.mcmc(acc.states,data,u)
  }
  return(list(results=acc.states,vars = cov.ls))
}
plot.zigpd.mcmc = function(mcmc.out,dat,threshold){
  layout.mat = t(matrix(c(rep(1,4),rep(4,3),
                          rep(1,4),rep(4,3),
                          rep(2,4),rep(4,3),
                          rep(2,4),rep(5,3),
                          rep(3,4),rep(5,3),
                          rep(3,4),rep(5,3)),ncol=6))
  L = layout(layout.mat)
  par(mar=c(0,4,0,0))
  
  plot(mcmc.out[-(1:(0.2*nrow(mcmc.out))),1], type='l',xaxt = 'n',xlab='',ylab='alpha')
  plot(mcmc.out[-(1:(0.2*nrow(mcmc.out))),2], type='l',xaxt = 'n',xlab='',ylab='shape')
  plot(mcmc.out[-(1:(0.2*nrow(mcmc.out))),3], type='l',ylab='scale')
  
  
  conf.size = 0.95
  par.means = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,mean)
  par.lower = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,quantile, probs = 1-conf.size)
  par.upper = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,quantile, probs = conf.size)
  
  k=1:2e3
  plot(epdf2(dat), log='xy',pch=16,cex=0.5)
  lines(dzigpd(k,par.means[1], threshold,par.means[2], par.means[3]), col='red')
  lines(dzigpd(k,par.lower[1], threshold,par.upper[2], par.upper[3]), col='blue')
  lines(dzigpd(k,par.upper[1], threshold,par.lower[2], par.lower[3]), col='blue')
  
  plot(ecdf2(dat), log='xy', type='l')
  lines(pzigpd(k,par.means[1], threshold,par.means[2], par.means[3], lower.tail=F), col='red')
  lines(pzigpd(k,par.lower[1], threshold,par.upper[2], par.upper[3], lower.tail=F), col='blue')
  lines(pzigpd(k,par.upper[1], threshold,par.lower[2], par.lower[3], lower.tail=F), col='blue')
  
  par(mfrow=c(1,1))
}
fast.zigpd.mcmc = function(n.iter,dat,threshold){
  mcmc.out.full = zigpd.mcmc(n.iter,dat,threshold,c(3,2,5),list(alpha=c(8,4),shape=10,scale=c(60,2)))
  return(mcmc.out.full)
}
}
##zero-zeta-igpd model functions
{
dzzigpd = function(p,alpha,u,shape,scale){
  if(x<0){
    return(0)
  }
  if(x==1){
    return(p)
  }
  if(x<=u){
    return((1-p)*dzeta(x,alpha))
  }
  if(x>u){
    return( (1-(1-p)*pzeta(u,alpha,lower.tail = F))*(pgpd(x,u,scale,shape)-pgpd(x-1,u,scale,shape)))
  }
}
lzzigpd = function(X,p,alpha,u,shape,scale){
  X.pl = X[1<X & X<=u]
  X.igpd = X[X>u]
  return(
    length(X[X==1])*log(p) + length(X.pl)*log(zeta(alpha+1))-(alpha+1)*sum(log(X.pl))+
      length(X[X>u])*log((1-(1-p)*(pzeta(u,alpha,lower.tail = F)-zeta(alpha+1)^-1))) + sum(log(pgpd(X.igpd,u,scale,shape)-pgpd(X.igpd-1,u,scale,shape)))
  )
}
p.prior = function(p,pars){
  return(dbeta(p,pars[1],pars[2],log=T))
}
zzigpd.prior = function(p,alpha,shape,scale,p.pars,alpha.pars,shape.pars,scale.pars){
  return(p.prior(p,p.pars) + alpha.prior(alpha,alpha.pars) + shape.prior(shape,shape.pars) + scale.prior(scale,scale.pars))
}
zzigpd.posterior = function(X,p,alpha,u,shape,scale,p.pars,alpha.pars,shape.pars,scale.pars){
  return(lzzigpd(X,p,alpha,u,shape,scale) + zzigpd.prior(p,alpha,shape,scale,p.pars,alpha.pars,shape.pars,scale.pars))
}
zzigpd.mcmc = function(n.iter,data,u,init,pri.pars,prop.cov = diag(c(0.000001,0.001,0.0001,0.001)),H=200){
  states = data.frame(p=numeric(n.iter+1),alpha=numeric(n.iter+1),shape=numeric(n.iter+1),scale=numeric(n.iter+1))
  states[1,] = init
  cov.ls = list()
  cov.ls[[1]] = prop.cov
  acc.states = data.frame(p=numeric(),alpha=numeric(),shape=numeric(),scale = numeric())
  state.prop=c(0,0,0,0)
  for(i in 1:n.iter+1){
    message('iter: ',i,' - accepted:',nrow(acc.states))
    # print(state.prop)
    # print(format(prop.cov))
    if(nrow(acc.states)>H){
      prop.cov =( 2.4^2 / 4 ) * cov(acc.states[(nrow(acc.states)-H+1) : nrow(acc.states),])
    }
    cov.ls[[i]] = prop.cov
    state.prop = rmvn(1,mu=states[i-1,],V=prop.cov)
    if(state.prop[2]<0){
      states[i,] = states[i-1,]
      next
    }
    if(state.prop[4]<0){
      states[i,] = states[i-1,]
      next
    }
    if(state.prop[1]<0 | state.prop[1]>1){
      print('p')
      states[i,] = states[i-1,]
      next
    }
    #acceptance
    logA = min(0,(zzigpd.posterior(data,state.prop[1],state.prop[2],u,state.prop[3],state.prop[4],pri.pars$p,pri.pars$alpha,pri.pars$shape,pri.pars$scale)-
                    zzigpd.posterior(data,states[i-1,1],states[i-1,2],u,states[i-1,3],states[i-1,4],pri.pars$p,pri.pars$alpha,pri.pars$shape,pri.pars$scale)))
    # print(logA)
    if(log(runif(1))<logA){
      states[i,] = state.prop
      acc.states = rbind(acc.states,state.prop)
      next
    }else{
      print(logA)
      states[i,] = states[i-1,]
      next
    }
  }
  colnames(acc.states) = c('p','alpha','shape','scale')
  return(list(results=acc.states,vars = cov.ls))
}
}
}











