#libraries
{
  require(igraph)
  require(evd)
  require(ismev)
  require(extRemes)
  require(VGAM)
  require(gsl)
  require(MASS)
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
    alpha.prior0 = function(alpha,pars,log=T){
      return(dgamma(alpha,shape=pars[1],rate=pars[2],log=log))
    }
    zeta.posterior = function(X,alpha,alpha.pars){
      return(lzeta(X,alpha)+alpha.prior0(alpha,alpha.pars))
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
    plot.zeta.mcmc = function(mcmc.out,dat,burn.in = 0.2*length(mcmc.out)){
      mcmc.out = mcmc.out[-(1:burn.in)]
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
    dzigpd = Vectorize(function(x,phi,alpha,u,shape,scale){
      if(x<=0){
        return(0)
      }
      if(x<=u){
        return((1-phi)*x^-(alpha+1) / (gsl::hzeta(alpha+1,1)-gsl::hzeta(alpha+1,u+1)))
      }
      if(x>u){
        return(phi*(pgpd(x-u,scale=scale,shape=shape)-pgpd(x-u-1,scale=scale,shape=shape)))
      }
    },vectorize.args = 'x')
    pzigpd = Vectorize(function(x,phi,alpha,u,shape,scale,lower.tail=F){
      if(x<=0){
        return(0)
      }
      if(x<=u){
        return( 1 -   (1-phi)*(gsl::zeta(alpha+1) - gsl::hzeta(alpha+1,x+1))/(gsl::zeta(alpha+1) - gsl::hzeta(alpha+1,u+1))   )
      }
      if(x>u){
        return(phi*pgpd(x-u,shape=shape,scale=scale,lower.tail = F))
      }
    },vectorize.args = 'x')
    lzigpd = function(X,phi,alpha,u,shape,scale){
      N = length(X)
      n = length(X[X<=u])
      X.pl = X[X<=u]
      X.igp = X[X>u]
      return(n*log(1-phi) + (N-n)*log(phi) - n*log(gsl::zeta(alpha+1) - gsl::hzeta(alpha+1,u+1)) - (alpha+1)*sum(log(X.pl)) + 
               sum(log(pgpd(X.igp-u,scale=scale,shape=shape)- pgpd(X.igp-u-1,shape = shape,scale=scale))))
    }
    
    shape.prior0 = function(shape,pars,log=T){
      return(dnorm(shape,mean=0,sd=pars[1],log=log))
    }
    scale.prior0 = function(scale,pars,log=T){
      return(dgamma(scale,shape=pars[1], rate=pars[2],log=log))
    }
    zigpd.prior = function(alpha,shape,scale,alpha.pars,shape.pars,scale.pars){
      return(alpha.prior0(alpha,alpha.pars) + shape.prior0(shape,shape.pars) + scale.prior0(scale,scale.pars))
    }
    zigpd.posterior = function(X,phi,alpha,u,shape,scale,alpha.pars,shape.pars,scale.pars){
      return(lzigpd(X,phi,alpha,u,shape,scale) + zigpd.prior(alpha,shape,scale,alpha.pars,shape.pars,scale.pars))
    }
    zigpd.mcmc = function(n.iter,data,phi,u,init,pri.pars,prop.cov = diag(c(0.01,0.01,0.01)),H=200, show=T){
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
        logA = min(0,(zigpd.posterior(data,phi,state.prop[1],u,state.prop[2],state.prop[3],pri.pars$alpha,pri.pars$shape,pri.pars$scale)-
                        zigpd.posterior(data,phi,states[i-1,1],u,states[i-1,2],states[i-1,3],pri.pars$alpha,pri.pars$shape,pri.pars$scale)))
        
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
        plot.zigpd.mcmc(acc.states,data,phi,u,pri.pars)
      }
      return(list(results=acc.states,vars = cov.ls))
    }
    plot.zigpd.mcmc = function(mcmc.out,dat,phi,threshold,pri.pars){
      layout.mat = t(matrix(c(rep(1,4),rep(4,3),
                              rep(1,4),rep(4,3),
                              rep(2,4),rep(4,3),
                              rep(2,4),rep(5,3),
                              rep(3,4),rep(5,3),
                              rep(3,4),rep(5,3)),ncol=6))
      L = layout(layout.mat)
      par(mar=c(2,4,0,1))
      
      conf.size = 0.95
      par.means = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,mean)
      par.lower = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,quantile, probs = 1-conf.size)
      par.upper = apply(mcmc.out[-(1:(0.2*nrow(mcmc.out))),],2,quantile, probs = conf.size)
      
      a=seq(0.01,2*par.upper[1],length.out=1e5)
      
      plot(density(mcmc.out[-(1:(0.2*nrow(mcmc.out))),1]), type='l',xlab='alpha',main='')
      lines(a,alpha.prior0(a,pri.pars$alpha,log=F), lty=2)
      legend('topright',legend=c('Prior','Posterior'),lty=c(2,1))
      s = seq(par.lower[2]-1,par.upper[2]+1,length.out=1e5)
      plot(density(mcmc.out[-(1:(0.2*nrow(mcmc.out))),2]), type='l',xlab='shape',main='')
      lines(s,shape.prior0(s,pri.pars$shape,log=F), lty=2)
      
      s = seq(par.lower[3]/2,par.upper[3]*2,length.out=1e5)
      plot(density(mcmc.out[-(1:(0.2*nrow(mcmc.out))),3]), type='l',xlab='scale',main='')
      lines(s,scale.prior0(s,pri.pars$scale,log=F), lty=2)
      
      
      
      k=1:2e3
      plot(epdf2(dat), log='xy',pch=16,cex=0.5)
      lines(dzigpd(k,phi,par.means[1], threshold,par.means[2], par.means[3]), col='red')
      lines(dzigpd(k,phi,par.lower[1], threshold,par.upper[2], par.upper[3]), col='blue')
      lines(dzigpd(k,phi,par.upper[1], threshold,par.lower[2], par.lower[3]), col='blue')
      
      plot(ecdf2(dat), log='xy', type='l')
      lines(pzigpd(k,phi,par.means[1], threshold,par.means[2], par.means[3], lower.tail=F), col='red')
      lines(pzigpd(k,phi,par.lower[1], threshold,par.upper[2], par.upper[3], lower.tail=F), col='blue')
      lines(pzigpd(k,phi,par.upper[1], threshold,par.lower[2], par.lower[3], lower.tail=F), col='blue')
      legend('bottomleft',legend = c(paste('Threshold = ', threshold),
                                     paste('# <= ',threshold,': ',sum(dat<=threshold)),
                                     paste('# > ',threshold,': ',sum(dat>threshold)),
                                     paste('Alpha: ', round(par.means[1],2)),
                                     paste('Shape: ', round(par.means[2],2)),
                                     paste('Scale: ', round(par.means[3],2))))
      
      par(mfrow=c(1,1))
    }
    fast.zigpd.mcmc = function(n.iter,dat,threshold,H=100){
      mcmc.out.full = zigpd.mcmc(n.iter,phi=sum(dat>threshold)/length(dat),dat,threshold,c(3,2,5),list(alpha=c(1,0.01),shape=3,scale=c(1,0.01)),H=H)
      return(mcmc.out.full)
    }
  }
  #zeta-changepoint model
  
  dzc = Vectorize(function(x,phi,u,alpha,beta){
    if(x<0){
      return(0)
    }
    if(x<=u){
      return((1-phi)*x^-(alpha+1) / (gsl::hzeta(alpha+1,1)-gsl::hzeta(alpha+1,u+1)))
    }
    if(x>u){
      return(phi*x^-(beta+1) / gsl::hzeta(beta+1,u+1))
    }
  },vectorize.args = 'x')
  pzc = Vectorize(function(x,phi,u,alpha,beta){
    if(x<0){
      return(1)
    }
    if(x<=u){
      return(1 - (1-phi)*((gsl::hzeta(alpha+1,1)-gsl::hzeta(alpha+1,x+1))/(gsl::hzeta(alpha+1,1)-gsl::hzeta(alpha+1,u+1))))
    }
    if(x>u){
      return(phi* gsl::hzeta(beta+1,x+1)/gsl::hzeta(beta+1,u+1))
    }
  },vectorize.args = 'x')
  lzc = function(X,phi,u,alpha,beta){
    X1 = X[X<=u]
    X2 = X[X>u]
    n = length(X1)
    N = length(X)
    
    return(n*log(1-phi) + (N-n)*log(phi) - n*log((gsl::hzeta(alpha+1,1)-gsl::hzeta(alpha+1,u+1))) + (n-N)*log(gsl::hzeta(beta+1,u+1)) - 
             (alpha+1)*sum(log(X1)) - (beta+1)*sum(log(X2)))
  }
  
  zc.pos = function(X,phi,u,alpha,beta,pri.pars){
    return(lzc(X,phi,u,alpha,beta) + alpha.prior0(alpha,pri.pars$alpha) + alpha.prior0(beta,pri.pars$alpha))
  }
  
  zc.mcmc = function(n.iter,dat,phi,u,init,pri.pars,prop.cov = diag(c(0.01,0.01)),H=200){
    states = data.frame(alpha = numeric(n.iter+1),beta=numeric(n.iter+1))
    states[1,] = init
    acc.states = data.frame(alpha=numeric(),beta=numeric())
    state.prop = c(0,0)
    
    for(i in 1:n.iter+1){
      message('iter: ',i,' - accepted: ',nrow(acc.states))
      if(nrow(acc.states)>H){
        prop.cov = (2.38^2 / 2)*cov(acc.states[(nrow(acc.states)-H+1):nrow(acc.states),])
      }
      
      state.prop = rmvn(1,mu=states[i-1,],V=prop.cov)
      
      if(any(state.prop<=0)){
        states[i,] = states[i-1,]
        next
      }
      
      logA = min(0,zc.pos(dat,phi,u,state.prop[1],state.prop[2],pri.pars) - zc.pos(dat,phi,u,states[i-1,1],states[i-1,2],pri.pars) )
      # message(logA)
      if(log(runif(1))<logA){
        states[i,] = state.prop
        acc.states[nrow(acc.states)+1,] = state.prop
        next
        
      }else{
        states[i,] = states[i-1,]
        next
      }
      
    }
    return(acc.states)
  }
  
  fast_zc_mcmc <- function(n.iter,dat, p ) {
    u = quantile(dat,p)
    phi = sum(dat>u)/length(dat)
    pri.pars = list(alpha=c(1,0.01))
    init = c(2,2)
    mcmc.out = zc.mcmc(n.iter,dat,phi,u,init,pri.pars,H=100)
    plot_zc_mcmc(mcmc.out,pri.pars,dat,phi,u)
    return(mcmc.out)
  }
  
  alpha_beta_prior = function(a,b,pri.pars){
    return(alpha.prior0(a,pri.pars,log=F)*alpha.prior0(b,pri.pars,log=F))
  }
  
  plot_zc_mcmc <- function(mcmc.out, pri.pars, dat, phi, u,burn.in = 0.2*nrow(mcmc.out)) {
    
    mcmc.out = mcmc.out[-(1:burn.in),]
    a = seq(0.1,max(mcmc.out),length.out=2e3)
    pri.dens = outer(a,a,FUN=alpha_beta_prior,pri.pars=pri.pars$alpha)
    
    xlim = range(mcmc.out$alpha)
    ylim = range(mcmc.out$beta)
    
    layout_mat = t(matrix(c(4,2,2,2,
                            3,1,1,1,
                            3,1,1,1,
                            3,1,1,1),nrow=4))
    layout(layout_mat)
    par(mar=c(4,0,0,0))
    contour(kde2d(mcmc.out[-(1:burn.in),]$alpha,mcmc.out[-(1:burn.in),]$beta),xlim=xlim,ylim=ylim,xlab='alpha')
    contour(a,a,pri.dens,add=T,lty=2,nlevels=15)
    par(mar=c(0,0,0,0))
    plot(density(mcmc.out[-(1:burn.in),]$alpha),main='',ylab='',xlab = 'alpha',xlim=xlim,xaxt='n',yaxt='n')
    lines(a,alpha.prior0(a,pri.pars$alpha,log=F),lty=2)
    legend('topright',legend=c('Prior','Posterior'),lty=c(2,1))
    b.dens = density(mcmc.out[-(1:burn.in),]$beta)
    par(mar=c(4,0,0,0))
    plot(b.dens$y,b.dens$x,type='l',xlim=c(max(b.dens$y),0),ylim=ylim,xaxt='n',xlab='',ylab='beta')
    lines(alpha.prior0(a,pri.pars$alpha,log=F),a,lty=2)
    
    par(mfrow=c(1,1),mar=c(3,3,3,3))
    
    abline(a=0,b=1,lty=2)
    means = apply(mcmc.out,2,mean)
    uppers = apply(mcmc.out,2,quantile,probs=0.95)
    lowers = apply(mcmc.out,2,quantile,probs=0.05)
    k=1:3e3
    plot(ecdf2(dat),type='l',log='xy')
    
    lines(k,pzc(k,phi,u,lowers[1],lowers[2]),lty=2,col='blue')
    lines(k,pzc(k,phi,u,means[1],means[2]),col='red')
    lines(k,pzc(k,phi,u,uppers[1],uppers[2]),lty=2,col='blue')
    
    legend('bottomleft',legend = c(paste('Threshold = ', u),
                                   paste('# <= ',u,': ',sum(dat<=u)),
                                   paste('# > ',u,': ',sum(dat>u)),
                                   paste('Phi: ',signif(phi,2)),
                                   paste('Alpha: ', round(means[1],2)),
                                   paste('Beta: ', round(means[2],2))))
    
    
  }
  
}





# new functions -----------------------------------------------------------


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
  acc.states = tail(acc.states, (1-burn.in)*nrow(acc.states))
  if(show){
    plot.mcmc.pl_igp(acc.states, prior.pars, u.type, X, dist.type)
  }
  return(acc.states)
}

plot.mcmc.pl_igp <- function(mcmc.out.burned, prior.pars, u.type, p2p.in, dist.type) {
  layout.mat = t(matrix(c(1,2,3,4,
                          5,5,6,6,
                          5,5,6,6), nrow=4))
  
  L = layout(layout.mat)
  par(mar = c(0,2,2,2))
  k=seq(range(mcmc.out.burned$alpha)[1]*0.75,range(mcmc.out.burned$alpha)[2]*1.25, length.out=1e3)
  pri.dens = alpha.prior(k,log=F,shape=prior.pars$alpha$shape, rate = prior.pars$alpha$rate)
  post.dens = density(mcmc.out.burned$alpha)
  plot(post.dens, ylim=c(0,max(c(pri.dens, post.dens$y))), xlab='',main='alpha')
  lines(k,pri.dens,lty=2)
  
  par(mar = c(0,0,2,2))
  k=seq(range(mcmc.out.burned$shape)[1]*0.75,range(mcmc.out.burned$shape)[2]*1.25, length.out=1e3)
  pri.dens = shape.prior(k,log=F,S=prior.pars$shape$S)
  post.dens = density(mcmc.out.burned$shape)
  plot(post.dens, ylim=c(0,max(c(pri.dens, post.dens$y))),main='shape', xlab='', ylab='')
  lines(k,pri.dens,lty=2)
  
  k=seq(range(mcmc.out.burned$scale)[1]*0.75,range(mcmc.out.burned$scale)[2]*1.25, length.out=1e3)
  pri.dens = scale.prior(k,log=F,shape=prior.pars$scale$shape, rate=prior.pars$scale$rate)
  post.dens = density(mcmc.out.burned$scale)
  plot(post.dens, ylim=c(0,max(c(pri.dens, post.dens$y))),main='scale', xlab='', ylab='')
  lines(k,pri.dens,lty=2)
  par(mar = c(0,0,2,0))
  k=seq(range(mcmc.out.burned$u)[1]*0.75,range(mcmc.out.burned$u)[2]*1.25, length.out=1e3)
  pri.dens = u_prior(k,type=u.type,log=F,N=prior.pars$u$N,shape=prior.pars$u$shape,rate=prior.pars$u$rate)
  post.dens = density(mcmc.out.burned$u)
  plot(post.dens, ylim=c(0,max(c(pri.dens, post.dens$y))),main='threshold', xlab='', ylab='')
  lines(k,pri.dens,lty=2)
  polygon(c(max(p2p.in),max(post.dens$x)*1.25,max(post.dens$x)*1.25,max(p2p.in)),c(-1,-1,max(post.dens$y)*1.25,max(post.dens$y)*1.25), density=5)
  
  q = sort(unique(p2p.in), decreasing = F)
  # q=1:max(p2p.in)
  
  means = apply(mcmc.out.burned, 2, mean)
  meds = apply(mcmc.out.burned,2,median)
  
  post.vec = numeric(nrow(mcmc.out.burned))
  for(i in 1:length(post.vec)){
    post.vec[i] = pl_igp.posterior(p2p.in,mcmc.out.burned[i,1],mcmc.out.burned[i,2], mcmc.out.burned[i,3],mcmc.out.burned[i,4],prior.pars,u.type,dist.type)
  }
  #survival#######################################################################
  dens.vec.df = data.frame(matrix(nrow=nrow(mcmc.out.burned), ncol=length(q)))
  for(i in 1:nrow(mcmc.out.burned)){
    message(i)
    dens.vec.df[i,] = pzigpd(q,sum(p2p.in>mcmc.out.burned[i,1])/length(p2p.in), mcmc.out.burned[i,2] ,mcmc.out.burned[i,1],mcmc.out.burned[i,3],mcmc.out.burned[i,4])
  }
  
  dens.vec.df[is.na(dens.vec.df)]=0
  dens.mean = apply(dens.vec.df, 2, mean)
  dens.upper95 = apply(dens.vec.df, 2, quantile, probs=0.975)
  dens.lower95 = apply(dens.vec.df, 2, quantile, probs=0.025)
  dens.upper99 = apply(dens.vec.df, 2, quantile, probs=0.995)
  dens.lower99 = apply(dens.vec.df, 2, quantile, probs=0.005)
  
  par(mar=c(4,4,2,0))
  plot(ecdf2(p2p.in), log='xy', type='l', col='red', lwd =1, ylab='S(x)')
  polygon(c(q,rev(q)), c(dens.upper95,rev(dens.lower95)), col=rgb(0.6, 0.6, 0.7,0.2), border=NA, lty=2)
  polygon(c(q,rev(q)), c(dens.upper99,rev(dens.lower99)), col=rgb(0.6, 0.6, 0.7,0.2), border=NA, lty=2)
  lines(ecdf2(p2p.in)$x, ecdf2(p2p.in)$p, col='red', lwd =1)
  legend('bottomleft', legend = 'Empirical survival', col='red', lty=1, lwd=1)
  #pmf############################################################################
  dens.vec.df = data.frame(matrix(nrow=nrow(mcmc.out.burned), ncol=length(q)))
  for(i in 1:nrow(mcmc.out.burned)){
    message(i)
    dens.vec.df[i,] = dzigpd(q,sum(p2p.in>mcmc.out.burned[i,1])/length(p2p.in), mcmc.out.burned[i,2] ,mcmc.out.burned[i,1],mcmc.out.burned[i,3],mcmc.out.burned[i,4])
  }
  
  dens.vec.df[is.na(dens.vec.df)]=0
  dens.mean = apply(dens.vec.df, 2, mean)
  dens.upper95 = apply(dens.vec.df, 2, quantile, probs=0.975)
  dens.lower95 = apply(dens.vec.df, 2, quantile, probs=0.025)
  dens.upper99 = apply(dens.vec.df, 2, quantile, probs=0.995)
  dens.lower99 = apply(dens.vec.df, 2, quantile, probs=0.005)
  
  par(mar=c(4,4,2,0))
  plot(epdf2(p2p.in), log='xy', col='red', ylab='f(x)', pch=19, cex = 0.5)
  polygon(c(q,rev(q)), c(dens.upper95,rev(dens.lower95)), col=rgb(0.6, 0.6, 0.7,0.2), border=NA, lty=2)
  polygon(c(q,rev(q)), c(dens.upper99,rev(dens.lower99)), col=rgb(0.6, 0.6, 0.7,0.2), border=NA, lty=2)
  
  
  par(mfrow=c(1,1))
}


rpl_igp = function(n,alpha,phi,threshold,shape,scale, deg_max = 2e3){
  if(n==1){
    u=runif(1)
    if(u>=phi){
      k=1:threshold
      probs = k^-(alpha+1) / (gsl::zeta(alpha+1) - gsl::hzeta(alpha+1, threshold+1))
      x = sample(k,1, prob=probs)
      return(x)
    }
    k = (threshold+1) : deg_max
    probs = pgpd(k-threshold, shape=shape, scale=scale) - pgpd(k-threshold-1, shape=shape, scale=scale)
    x = sample(k,1,prob=probs)
    return(x)
  }else{
    return(replicate(n,rpl_igp(1,alpha,phi,threshold,shape,scale,deg_max)))
  }
  
}






























