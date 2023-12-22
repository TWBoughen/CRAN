library(ggplot2)
library(gridExtra)
library(networkdata)
pgpd = function(x, shape, scale, threshold=0, lower.tail=T){
  if(!lower.tail){
    return(1-pgpd(x,shape,scale,threshold))
  }
  if(length(x)>1){
    return(pgpd.vec(x,shape,scale,threshold, lower.tail))
  }else{
    if(scale<0){
      stop('error:scale<0')
    }
    if(threshold<0){
      stop('error:threshold<0')
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

pli.joint.prior = function(phi,alpha,shape,scale,params){
  return(dbeta(phi,params$phi[1], params$phi[2], log = T) +
           dnorm(alpha,0,sd=params$alpha,log=T) +
           dnorm(shape,0,sd=params$shape, log=T)+
           dgamma(scale, 1, rate=params$scale, log=T))
}

pli.joint.posterior = function(X,alpha,shape,scale,threshold,params,type=1){
  u=threshold
  if(type==1){
    v = ceiling(u)
  }else if(type==2){
    v = floor(u)
  }
  phi = sum(X[X[,1]<=v,2] / sum(X[,2]))
  return(llpli2(X,alpha,shape,scale,threshold,type) + pli.joint.prior(phi,alpha,shape,scale,params))
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

pli.mcmc = function(n.iter, x, init, prior.params, init.cov, type=1, burn.in=0, cov.period=1e4, fix.u = NULL, backup_period = 1e7, show=F, logging=T, log_id = NULL){
  alpha.acc = data.frame(alpha=init$alpha)
  shapescale.acc = data.frame(shape=init$shape, scale=init$scale)
  thresh.acc = data.frame(threshold=init$threshold)
  alpha.cov = init.cov$alpha
  shapescale.cov = init.cov$shapescale
  thresh.cov = init.cov$thresh
  post.vals = numeric(n.iter+1)
  states = data.frame(alpha = numeric(n.iter+1), shape = numeric(n.iter+1), scale=numeric(n.iter+1), threshold = numeric(n.iter+1))
  states[1,] = init
  post.vals[1] = pli.joint.posterior(x,init$alpha,init$shape,init$scale,init$threshold,prior.params,type)
  thresh.props = numeric(n.iter)
  if(show){
    dev.new()
    plot.new()
    plot.window(c(0,n.iter), c(0,60))
    axis(side=1)
    axis(side=2)
  }
  if(logging){
    if(is.null(log_id)){
      log_id = as.character(signif(abs(post.vals[1]),min(log10(abs(post.vals[1])),5)))
    }
    writeLines(c(""), paste0('logs/log-',log_id,'.txt'))
  }
  
  x.max = max(x[,1])
  for(i in 1:n.iter+1){
    if(logging){
      cat('Iteration: ', i, "\n", file=paste0('logs/log-',log_id,'.txt'))
    }else{
      message('Iteration: ', i)
    }
    # 
    #alpha steps
    
    if(nrow(alpha.acc)%%cov.period == 0){
      alpha.cov = cov.prop(tail(alpha.acc,cov.period))
    }
    alpha.prop = rnorm(1,tail(alpha.acc,1)[[1]],alpha.cov)
    if(alpha.prop!=0){
      logA = min(0, pli.joint.posterior(x,alpha.prop,tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type)-
                   pli.joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type))
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
        logA = min(0,pli.joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],thresh.prop,prior.params,type)-
                     pli.joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type))
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
      
      logA = min(0,pli.joint.posterior(x,tail(alpha.acc,1)[[1]],shapescale.prop[1],shapescale.prop[2],tail(thresh.acc,1)[[1]],prior.params,type)-
                   pli.joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type))
      if(log(runif(1))<logA){
        shapescale.acc = rbind(shapescale.acc, shapescale.prop)
      }
    }
    states[i,] = c(tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]])
    post.vals[i] = pli.joint.posterior(x,tail(alpha.acc,1)[[1]],tail(shapescale.acc$shape,1)[[1]],tail(shapescale.acc$scale,1)[[1]],tail(thresh.acc,1)[[1]],prior.params,type)
    if(i%%backup_period==0){
      mcmc.out.temp <<- list(alpha=alpha.acc, shapescale=shapescale.acc, thresh=thresh.acc, post.vals = post.vals, states = states)
    }
    if(i>1 && show){
      dev.flush()
      plot(max(1,i-499):i,tail(states$thresh[1:i],500), type='l', 
           main = paste(round(tail(states[1:i,],1),3),sep=','))
      legend('topleft', legend = paste0('Posterior Log-Likelihood: ',round(tail(post.vals[1:i],1),3)))
    }
  }
  return(list(alpha=alpha.acc, shapescale=shapescale.acc, thresh=thresh.acc, post.vals = post.vals, states = states))
}


# -------------------------------------------------------------------------
plot.pli.mcmc = function(X,mcmc.out.states, burn.in=0, thin.by = 0, top='', show=F, by=1){
  Y=X
  Y[,2] = cumsum(X[,2])
  Y[,2] = Y[,2]/max(Y[,2])
  
  N = length(mcmc.out.states$alpha)
  mcmc.out.states = tail(mcmc.out.states, N-burn.in)[seq(1,N-burn.in, by=thin.by+1),]
  mcmc.out.states$phi = sapply(mcmc.out.states$threshold, FUN=function(x){return(sum(X[X[,1]>x,2])/sum(X[,2]))})
  k = seq(1,max(X[,1])+by, by=by)
  
  pmf.df = matrix(nrow=length(mcmc.out.states$shape), ncol=length(k))
  surv.df = matrix(nrow=length(mcmc.out.states$shape), ncol=length(k))
  for(i in 1:nrow(surv.df)){
    if(show){
      message(i,'/',nrow(surv.df))
    }
    
    surv.df[i,] = ppli(k,mcmc.out.states$phi[i], mcmc.out.states$alpha[i], 
                       mcmc.out.states$shape[i], mcmc.out.states$scale[i], mcmc.out.states$threshold[i], lower.tail=F)
    # pmf.df[i,] = dpli(k,mcmc.out.states$phi[i], mcmc.out.states$alpha[i], 
    #                   mcmc.out.states$shape[i], mcmc.out.states$scale[i], mcmc.out.states$threshold[i])
  }
  
  surv.upper95 = apply(surv.df, 2, quantile, probs = 0.975, na.rm=T)
  surv.lower95 = apply(surv.df, 2, quantile, probs = 1-0.975, na.rm=T)
  surv.mean = apply(surv.df, 2, mean, na.rm=T)
  # pmf.upper95 = apply(pmf.df, 2, quantile, probs = 0.975, na.rm=T)
  # pmf.lower95 = apply(pmf.df, 2, quantile, probs = 1-0.975, na.rm=T)
  
  #plots
  p1 = ggplot(data=mcmc.out.states) + geom_density(aes(x=alpha), linewidth=0.1)+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  p2 = ggplot(data=mcmc.out.states, aes(x=shape, y=scale)) + geom_density_2d(linewidth=0.2)+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  p3 = ggplot(data = mcmc.out.states) + geom_histogram(aes(x = ceiling(threshold)))+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  p4=ggplot() + geom_line(data=Y,aes(x=x, y=1-Freq),linewidth=0.1)+geom_ribbon(aes(x=k,ymin=surv.lower95, ymax=surv.upper95), alpha=0.3) + geom_line(aes(x=k, y=surv.mean),colour='red')+
    scale_x_log10() +scale_y_log10()+ylab('Surv')+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  # p5 = ggplot(data=X/sum(X[,2])) +geom_point(aes(x=x, y=Freq))+geom_ribbon(aes(x=x,ymin=pmf.lower95, ymax=pmf.upper95), alpha=0.3)+   scale_x_log10() +scale_y_log10()+theme(text = element_text(size=5))
  plot.list = list(p1,p2,p3,p4)
  return(arrangeGrob(p1,p2,p3,p4, nrow = 2, ncol=2,top=top))
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

# Power Law --------------------------------------------------------------

llpl = function(X, alpha){
  N = sum(X[,2])
  return(-N*log(gsl::zeta(alpha+1)) - (alpha+1)*sum(X[,2]*log(X[,1])))
}
pl.prior = function(alpha, param){
  return(dgamma(alpha, 1, param, log=T))
}
pl.posterior = function(X,alpha, prior.param){
  return(llpl(X,alpha) + pl.prior(alpha, prior.param))
}

pl.mcmc = function(n.iter,X, init, prior.param,prop.var.init = 0.1, S=1e3, show=F){
  alpha.acc = c(init)
  states = numeric(n.iter+1)
  states[1] = init
  prop.var = prop.var.init
  var.vec = numeric(S)
  for(i in 1:n.iter+1){
    if(show){
      message(i-1)
    }
    
    var.vec[i] = prop.var
    alpha.prop = rnorm(1, tail(alpha.acc, 1), prop.var)
    if(alpha.prop>0){
      logA = min(0,pl.posterior(X,alpha.prop, prior.param) - pl.posterior(X,tail(alpha.acc, 1), prior.param))
      if(log(runif(1))< logA){
        alpha.acc = c(alpha.acc, alpha.prop)
        if(i-1 <=S){
          prop.var = prop.var*(1+ 0.03/sqrt(i-1))
        }
      }else if(i-1 <=S){
        prop.var = prop.var*(1 - 0.01/sqrt(i-1))
      }
    }else if(i-1<=S){
      prop.var = prop.var*(1 - 0.01/sqrt(i-1))
    }
  }
  return(list(states = states, acc = alpha.acc, var = var.vec))
}

pl.mcmc.plot = function(X,mcmc.out.states, burn.in=1e3, thin.by=0, top=''){
  mcmc.out.states = tail(mcmc.out.states, length(mcmc.out.states)-burn.in)
  mcmc.out.states = mcmc.out.states[seq(1,length(mcmc.out.states), by=thin.by+1)]
  k = unique(X[,1])
  surv.df = data.frame(matrix(nrow=length(mcmc.out.states), ncol=length(k)))
  for(i in 1:nrow(surv.df)){
    # print(i)
    surv.df[i,] = ppl(k,mcmc.out.states[i], lower.tail = F)
  }
  Y=X
  Y[,2] = cumsum(X[,2])
  Y[,2] = Y[,2]/max(Y[,2])
  surv.upper95 = apply(surv.df, 2, quantile, probs = 0.975, na.rm=T)
  surv.lower95 = apply(surv.df, 2, quantile, probs = 1-0.975, na.rm=T)
  surv.mean = apply(surv.df,2,mean, na.rm=T)
  
  p1 = ggplot(data=Y) + geom_line(aes(x=x, y=1-Freq))+geom_ribbon(aes(x=x,ymin=surv.lower95, ymax=surv.upper95), alpha=0.3) + geom_line(aes(x=x, y=surv.mean),colour='red')+
    scale_x_log10() +scale_y_log10()+ylab('Surv') +ggtitle('Survival Function')+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  p2 = ggplot() + geom_density(aes(x=mcmc.out.states)) +xlab(latex2exp::TeX('\\alpha')) + ggtitle(latex2exp::TeX('Posterior Density of \\alpha'))+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  plot.list = list(p2,p1)
  return(arrangeGrob(p2,p1, nrow = 1, ncol=2,top=top))
}

#Power Law Mixture####################################

pplpl = function(x,alpha, beta, u, phi, lower.tail=F){
  if(length(x)>1){
    return(pplpl.vec(x,alpha,beta,u,phi,lower.tail))
  }
  if(lower.tail){
    return(1-pplpl(x,alpha,beta,u,phi))
  }
  if(x<=ceiling(u)){
    return(1-(1-phi)* sum((1:x)^-(alpha+1))/sum((1:ceiling(u))^-(alpha+1)))
  }
  if(x>ceiling(u)){
    return(phi*(1 - sum(  ( (ceiling(u)+1) :x) ^-(beta+1)   ) / gsl::hzeta(beta+1, ceiling(u)+1)    ))
  }
}
pplpl.vec = Vectorize(pplpl, vectorize.args = 'x')




llplpl = function(X,alpha, beta, u){
  v=ceiling(u)
  N = sum(X[,2])
  n = sum(X[X[,1]<=v,2])
  phi = (N-n) / N
  if(n< 3|(N-n)<3){
    return(-Inf)
  }else{
    return(n*log(1-phi) + (N-n)*log(phi) - n*log(sum((1:v)^-(alpha+1))) + 
             (n-N)*log(gsl::hzeta(beta+1, v+1) )- (alpha+1)*sum(X[X[,1]<=v,2]*log(X[X[,1]<=v,1])) -
                         (beta+1)*sum(X[X[,1]>v,2]*log(X[X[,1]>v,1])))
  }
}
plpl.joint.prior = function(alpha,beta,phi, params){
  return(
    dnorm(alpha, 0, sd  =sqrt(1/params$alpha), log=T)+
      dgamma(beta, 1, params$beta, log=T)+
      dbeta(phi, 1, 1, log=T)
  )
}
plpl.joint.posterior = function(X,alpha,beta,u,prior.params){
  phi = sum(X[X[,1]>ceiling(u),2])/sum(X[,2])
  return(llplpl(X,alpha,beta,u) + plpl.joint.prior(alpha, beta, phi, prior.params))
}

plpl.mcmc = function(n.iter, X, init, prior.params, cov.init, S=1e3, H=1e3, log_id=NULL){
  threshold.cov = cov.init$threshold
  ab.cov = cov.init$ab
  states = data.frame(alpha=numeric(n.iter),beta = numeric(n.iter), threshold=numeric(n.iter), phi = numeric(n.iter))
  ab.acc = data.frame(init$alpha,init$beta)
  threshold.acc = c(init$threshold)
  if(!is.null(log_id)){
    writeLines(c(""), paste0('logs/log-',log_id,'.txt'))
  }
  for(i in 1:n.iter){
    #threshold steps
    # print(i)
    if(!is.null(log_id)){
      cat('Iteration: ', i, "\n", file=paste0('logs/log-',log_id,'.txt'))
    }
    threshold.prop = rnorm(1,tail(threshold.acc,1),sqrt(threshold.cov))
    if(threshold.prop>0){
      logA = min(0, plpl.joint.posterior(X, tail(ab.acc,1)[1][[1]], tail(ab.acc, 1)[2][[1]], threshold.prop, prior.params))-
                    plpl.joint.posterior(X, tail(ab.acc,1)[1][[1]], tail(ab.acc, 1)[2][[1]],tail(threshold.acc,1), prior.params)
      if(log(runif(1))<logA){
        threshold.acc = c(threshold.acc,threshold.prop)
        if(i<=S){
          threshold.cov = threshold.cov*(1+0.03/sqrt(i))
        }
      }else if(i<=S){
        threshold.cov = threshold.cov*(1-0.01/sqrt(i))
      }
    }
    #alpha and beta steps
    if(length(ab.acc)>H){
      ab.cov = 0.5*2.38^2* cov(states[1:2,]) 
    }
    ab.prop = rmvn(1,tail(ab.acc,1), ab.cov)
    
    if(ab.prop[2]>0){
      logA = min(0,plpl.joint.posterior(X, ab.prop[1], ab.prop[2],tail(threshold.acc,1), prior.params) -
                   plpl.joint.posterior(X, tail(ab.acc,1)[1][[1]], tail(ab.acc, 1)[2][[1]],tail(threshold.acc,1), prior.params))
      if(log(runif(1))<logA){
        ab.acc  = rbind(ab.acc, ab.prop)
      }
    }
    u = tail(threshold.acc,1)
    phi.curr = sum(X[X[,1]>ceiling(u),2])/sum(X[,2])
   states[i,] = c(tail(ab.acc,1), tail(threshold.acc,1), phi.curr)
  }
  return(states)
}

plot.plpl.mcmc = function(X,mcmc.out.states, burn.in=0, thin.by=0,top='', show=F){
  mcmc.out.states = mcmc.out.states[seq(burn.in, nrow(mcmc.out.states), by=thin.by+1),]
  k = unique(X[,1])
  surv.df = data.frame(matrix(ncol=length(k), nrow=nrow(mcmc.out.states)))
  
  for(i in 1:nrow(surv.df)){
    if(show){
      message(i,'/',nrow(surv.df))
    }
    surv.df[i,] = pplpl(k, mcmc.out.states$alpha[i], mcmc.out.states$beta[i],ceiling(mcmc.out.states$threshold[i]),mcmc.out.states$phi[i], lower.tail=F)
  }
  
  Y=X
  Y[,2] = cumsum(X[,2])
  Y[,2] = Y[,2]/max(Y[,2])
  surv.upper95 = apply(surv.df, 2, quantile, probs = 0.975, na.rm=T)
  surv.lower95 = apply(surv.df, 2, quantile, probs = 1-0.975, na.rm=T)
  surv.mean = apply(surv.df,2,mean, na.rm=T)
  
  
  p1 = ggplot(data=Y) + geom_line(aes(x=x, y=1-Freq),linewidth=0.1)+geom_ribbon(aes(x=x,ymin=surv.lower95, ymax=surv.upper95), alpha=0.3) + geom_line(aes(x=x, y=surv.mean),colour='red')+
    scale_x_log10() +scale_y_log10()+ylab('Surv') +ggtitle('Survival Function')+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  
  p2 = ggplot(data=mcmc.out.states) + geom_density_2d(aes(x=alpha, y=beta),linewidth=0.1)+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  p3 = ggplot(data = mcmc.out.states) + geom_histogram(aes(x = ceiling(threshold)))+theme(text = element_text(size=20))+theme(plot.title = element_text(hjust = 0.5))
  plot.list = list(p2,p3,p1)
  
  return(arrangeGrob(p2,p3,p1, nrow=2, ncol=2, top=top))
}

# misc --------------------------------------------------------------------
plot.surv = function(X){
  Y=X
  Y[,2] = cumsum(X[,2])
  Y[,2] = Y[,2]/max(Y[,2])
  p1 = ggplot(data=Y) + geom_line(aes(x=x, y=1-Freq))+scale_x_log10() + scale_y_log10() + ylab('Survival') +xlab('Degree')
  return(p1)
}