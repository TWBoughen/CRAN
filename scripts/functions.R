# libs --------------------------------------------------------------------
require(igraph)
require(ggraph)
require(evd)
require(ismev)
require(extRemes)
require(VGAM)
require(rjags)
# funs --------------------------------------------------------------------

rbind_ind = function(df1, df2){
  names1 = names(df1)
  names2 = names(df2)
  names(df2) = names1
  return(rbind(df1, df2))
}

# reading data ------------------------------------------------------------

load_depend_data <- function() {
  df.raw = read.csv('../data/rpkg_20190129.csv')
  df.local = df.raw
  df.reverse = df.local[df.local$reverse,c(2,1,3,4)]
  df1 = rbind(df.local[!df.local$reverse,], df.reverse)
  message("data loaded to variable called 'dat'")
  dat<<-df1
}

# -------------------------------------------------------------------------
create_graph <- function() {
  df=dat
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
  message("graph created called 'G' and a backup called 'OG' ")
  G<<-abc.g
  OG<<-abc.g
}

reduce_graph <- function(degCutoff = 1000,degCutoff_neighbors = 50,neighbor_limit=1) {
  abc.g = G
  vs_close_to =  unique(unlist(ego(abc.g, order = neighbor_limit, nodes = V(abc.g)[V(abc.g)$deg>=degCutoff])))
  abc.sg = subgraph(abc.g,vs_close_to[V(abc.g)[vs_close_to]$deg>=degCutoff_neighbors])
  G<<-abc.sg
}

focus_type <- function(t='depends') {
  abc.g=G
  abc.sg2 = delete.edges(abc.g, E(abc.g)[E(abc.g)$type!=t])
  abc.sg2 = delete.vertices(abc.sg2,V(abc.sg2)[degree(abc.sg2)==0])
  G<<-abc.sg2
}


reset = function(){
  G<<-OG
}


showG = function(){
  plot(G, vertex.label = V(G)$name, vertex.size = 4,vertex.label.font = 11,vertex.label.color = 'blue',vertex.color = 'white',
       edge.arrow.size = 0.2, layout = layout_nicely(G))
}


toggle_labels = function(){
  if(all(is.na(V(G)$name))){
    V(G)$name<<-V(G)$label
  }else{
    V(G)$name <<- NA
  }
}

label_hubs = function(cutoffPc = 0.1){
  hub = degree(G)/max(degree(G)) >= cutoffPc
  V(G)$name[hub] <<- V(G)$label[hub]
  V(G)$name[!hub] <<-NA
  
}

plot_graph = function(){
  ggraph(G, layout = 'auto')+
    geom_edge_link(arrow = arrow(length = unit(2, 'mm')),end_cap = circle(2, 'mm'), color = 'darkgrey')+
    geom_node_point(aes(size=degree(G)), color='black', fill='darkgrey', pch=21)+
    geom_node_text(aes(label=V(G)$name),color='black', repel=T,fontface='bold')
  
}

deg_exp = function(mode='all', min = 1, max = NA){
  degs = degree.distribution(G, mode=mode)
  if(is.na(max)){
    max = length(degs)+1 
  }
  degs = degs[min:max-1]
  return(fit_power_law(degs))
}

overview.type = function(type, threshold=NA, probs = c(0.5,0.999)){
  reset()
  focus_type(type)
  degs = degree(G, mode='in')
  if(is.na(threshold)){
    threshold = max(degs)
  }
  kmax=threshold
  degs.dist = degree.distribution(G, mode='in')
  degs = which(degs.dist>0)
  ps = degs.dist[degs.dist>0]
  
  degs.raw = degree(G, mode='in')
  degs.tofit = degs.raw[degs.raw<=kmax]
  fit =power.law.fit(degs.tofit+1)
  a = fit$alpha
  k = seq(0,1000)
  pk = k^-a / zeta(a)
  
  par(mfrow = c(2,2))
  
  plot(degs, ps, log='xy', ylab='Density', xlab = 'Degree', main = paste("Type: ", type ))
  lines(k, pk)
  abline(v=kmax, lty=2)
  legend('topright', legend = c(paste('Power law coefficient: ', round(a,3))))
  
  degrees = degree(G, mode='in')
  thresh.limits = quantile(degrees,probs )
  
  evd::mrlplot(degrees,tlim=thresh.limits)
  tcplot(degs, tlim = thresh.limits, model='gpd', vci=F, lty = c(2,1,2))
  reset()
  par(mfrow = c(1,1))
}


pl.mcmc <- function(x,a.init = 1.2, scale=0.02,N=1e4) {
  B.init = log(a.init-1)
  B.vec = numeric(N+1)
  B.vec[1] = B.init
  
  for(i in 2:N+1){
    B.prop = rnorm(1, mean=B.vec[i-1], sd=scale)
    a.prop = exp(B.prop)+1
    a.cur = exp(B.vec[i-1])+1
    A = exp(min(0,pl.ll(x+1,a.prop) - pl.ll(x+1,a.cur)))
    
    if(runif(1)<A){
      B.vec[i] = B.prop
    }else{
      B.vec[i] = B.vec[i-1]
    }
  }
  rej.rate = sum(duplicated(B.vec))/N
  B.vec = B.vec[!duplicated(B.vec)]
  a.vec = exp(B.vec)+1
  return(list(results = a.vec, acc.rate = 1-rej.rate))
  
}


dpl = function(x,a){
  return(x^-a / zeta(a))
}
pl.ll = function(X,a){
  return(sum(log(dpl(X,a))))
}
spl = function(x,a){
  return(zeta(a,shift=x)/zeta(a))
}
dpl.v = Vectorize(dpl, 'a')








