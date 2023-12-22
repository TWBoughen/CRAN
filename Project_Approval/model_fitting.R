source('../scripts/functions.R')
source('../scripts/new_functions.R')
library(networkdata)
library(mvtnorm)
#loading data############################
x=degree(atp[[54]],mode='in')
x=x[x>0]
tennis = as.data.frame(table(x))
tennis[,1] = as.numeric(as.character(tennis[,1]))

harvard = read.csv('../data/harvard.txt')
colnames(harvard) = c('x', 'Freq')

data("protein", package='networkdata')
x = degree(protein)
protein.dat = as.data.frame(table(x[x>0]))
protein.dat[,1] = as.numeric(as.character(protein.dat[,1]))/2
colnames(protein.dat) = c('x', 'Freq')


df = load_data('../data/rpkg_20190129.csv')
df0 = df-1
x = df0$depends[df0$depends>0]
depends = as.data.frame(table(x))
depends[,1] = as.numeric(as.character(depends[,1]))

set.seed(123)
G = barabasi.game(3e4)
x = degree(G, mode='in')
x=x[x>0]
sim = as.data.frame(table(x))
sim[,1] = as.numeric(as.character(sim[,1]))

#PLIGP######################################
library(foreach)
library(doParallel)
library(doSNOW)
library(futile.logger)

data.list = list(tennis, harvard, depends, sim, protein.dat)
data.names = list('tennis', 'harvard', 'depends', 'sim', 'protein')
mcmc.out.list = list()

init = list(alpha=1, shape=1, scale=5, threshold=10)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=0.1)

#multitail -N 1 -m 1 ./* -w --no-mark-change -d

cl = makePSOCKcluster(length(data.list))
registerDoSNOW(cl)
N=2e5
results <- foreach(i=1:length(data.list)) %dopar% {
  source('../scripts/functions.R')
  source('../scripts/new_functions.R')
  library(networkdata)
  library(mvtnorm)
  out = pli.mcmc(N,data.list[[i]],init,prior.params,init.cov, cov.period = 1e7, type=1, show=F, logging = T, log_id = data.names[[i]])
  saveRDS(out, file=paste0('mcmc.outputs/pli-',data.names[[i]], '.rds'))
  return(out)
}
stopCluster(cl)

# PLPL --------------------------------------------------------------------

init = list(alpha=1, beta=1, threshold=5)
prior.params = list(alpha=0.01, beta=0.01, threshold = 0.01)
cov.init = list(ab = diag(c(0.1,0.1)), threshold=0.5)

cl = makePSOCKcluster(length(data.list))
registerDoSNOW(cl)
N=2e5
results <- foreach(i=1:length(data.list)) %dopar% {
  source('../scripts/functions.R')
  source('../scripts/new_functions.R')
  library(networkdata)
  library(mvtnorm)
  out = plpl.mcmc(N, data.list[[i]],init ,
                  prior.params,
                  cov.init, log_id=data.names[[i]])
  saveRDS(out, file=paste0('mcmc.outputs/plpl-',data.names[[i]], '.rds'))
  return(out)
}
stopCluster(cl)
# PL ----------------------------------------------------------------------



# -------------------------------------------------------------------------



tennis.pli = readRDS('mcmc.outputs/pli-tennis.rds')
tennis.pli.plot = plot.pli.mcmc(tennis,tennis.pli$states,burn.in=1e4, thin.by=20, top = grid::textGrob('Tennis', gp=grid::gpar(fontsize=24)), show=T)

harvard.pli = readRDS('mcmc.outputs/pli-harvard.rds')
harvard.pli.plot = plot.pli.mcmc(harvard,harvard.pli$states,burn.in=1e4, thin.by=20, top = grid::textGrob('Harvard', gp=grid::gpar(fontsize=24)), show=T)

protein.pli = readRDS('mcmc.outputs/pli-protein.rds')
protein.pli.plot = plot.pli.mcmc(protein.dat,protein.pli$states,burn.in=1e4, thin.by=20, top = grid::textGrob('Protein', gp=grid::gpar(fontsize=24)), show=T)

depends.pli = readRDS('mcmc.outputs/pli-depends.rds')
depends.pli.plot = plot.pli.mcmc(depends,depends.pli$states,burn.in=1e4, thin.by=20, top = grid::textGrob('CRAN', gp=grid::gpar(fontsize=24)), show=T)

sim.pli = readRDS('mcmc.outputs/pli-sim.rds')
sim.pli.plot = plot.pli.mcmc(sim,sim.pli$states,burn.in=1e4, thin.by=20, top = grid::textGrob('BA', gp=grid::gpar(fontsize=24)), show=T)



saveRDS(marrangeGrob(list(tennis.pli.plot,harvard.pli.plot,protein.pli.plot,depends.pli.plot,sim.pli.plot), nrow=2, ncol=3, top=''), file='plotRDS/pli-plot.rds')

ggsave(marrangeGrob(list(tennis.pli.plot,harvard.pli.plot,protein.pli.plot,depends.pli.plot,sim.pli.plot),layout_matrix = layout_mat, top=''), file='../poster3/pli-plot.svg', device='svg', width=26, height=11)

# -------------------------------------------------------------------------


tennis.plpl = readRDS('mcmc.outputs/plpl-tennis.rds')
tennis.plpl.plot = plot.plpl.mcmc(tennis,tennis.plpl,burn.in=1e4, thin.by=20,  top = grid::textGrob('Tennis', gp=grid::gpar(fontsize=24)), show=T)

harvard.plpl = readRDS('mcmc.outputs/plpl-harvard.rds')
harvard.plpl.plot = plot.plpl.mcmc(harvard,harvard.plpl,burn.in=1e4, thin.by=20,  top = grid::textGrob('Harvard', gp=grid::gpar(fontsize=24)), show=T)

protein.plpl = readRDS('mcmc.outputs/plpl-protein.rds')
protein.plpl.plot = plot.plpl.mcmc(protein.dat,protein.plpl,burn.in=1e4, thin.by=20,  top = grid::textGrob('Protein', gp=grid::gpar(fontsize=24)), show=T)

depends.plpl = readRDS('mcmc.outputs/plpl-depends.rds')
depends.plpl.plot = plot.plpl.mcmc(depends,depends.plpl,burn.in=1e4, thin.by=20,  top = grid::textGrob('CRAN', gp=grid::gpar(fontsize=24)), show=T)

sim.plpl = readRDS('mcmc.outputs/plpl-sim.rds')
sim.plpl.plot = plot.plpl.mcmc(sim,sim.plpl,burn.in=1e4, thin.by=20, top = grid::textGrob('BA', gp=grid::gpar(fontsize=24)), show=T)


layout_mat = t(matrix(c(1,1,2,2,3,3,
                        1,1,2,2,3,3,
                        NA,4,4,5,5,NA,
                        NA,4,4,5,5,NA),nrow=6,ncol=4))

ggsave(marrangeGrob(list(tennis.plpl.plot,harvard.plpl.plot,protein.plpl.plot,depends.plpl.plot,sim.plpl.plot),layout_matrix = layout_mat, top=''), file='../poster3/plpl-plot.svg', device='svg', width=26, height=13)


saveRDS(marrangeGrob(list(tennis.plpl.plot,harvard.plpl.plot,protein.plpl.plot,depends.plpl.plot,sim.plpl.plot),layout_matrix = layout_mat, top=''), file='plotRDS/plpl-plot.rds')

# -------------------------------------------------------------------------

tennis.surv= plot.surv(tennis) +ggtitle('Tennis')+theme(text = element_text(size=35))+theme(plot.title = element_text(hjust = 0.5))
harvard.surv = plot.surv(harvard) + ggtitle('Harvard')+theme(text = element_text(size=35))+theme(plot.title = element_text(hjust = 0.5))
protein.surv = plot.surv(protein.dat) + ggtitle('Protein')+theme(text = element_text(size=35))+theme(plot.title = element_text(hjust = 0.5))
depends.surv = plot.surv(depends)+ggtitle('CRAN')+theme(text = element_text(size=35))+theme(plot.title = element_text(hjust = 0.5))

plot.list = list(tennis.surv,
                 harvard.surv,
                 protein.surv,
                 depends.surv)
surv.plots = marrangeGrob(plot.list, nrow=2, ncol=2, top='')
ggsave(file='../poster3/surplots.svg', surv.plots, device='svg', width=23, height=15)


surv.plots

# -------------------------------------------------------------------------

sim.surv = plot.surv(sim)
ggsave(file='../poster3/simsurv.svg', sim.surv, device='svg', width=4.5, height=2.5)

# -------------------------------------------------------------------------
tennis.pl.mcmc = pl.mcmc(1e5, tennis, 1, 0.01, S=1e4)
harvard.pl.mcmc = pl.mcmc(1e5, harvard, 1, 0.01, S=1e4)
protein.pl.mcmc = pl.mcmc(1e5, protein.dat, 1, 0.01, S=1e4)
depends.pl.mcmc = pl.mcmc(1e5, depends, 1, 0.01, S=1e4)
sim.pl.mcmc = pl.mcmc(1e5, sim, 1, 0.01, S=1e4)

tennis.pl.plot = pl.mcmc.plot(tennis, tennis.pl.mcmc$acc, burn.in=500, thin.by=5,top = grid::textGrob('Tennis', gp=grid::gpar(fontsize=45)))
harvard.pl.plot = pl.mcmc.plot(harvard, harvard.pl.mcmc$acc, burn.in=500, thin.by=5, top = grid::textGrob('Harvard', gp=grid::gpar(fontsize=45)))
protein.pl.plot = pl.mcmc.plot(protein.dat, protein.pl.mcmc$acc, burn.in=500, thin.by=5, top = grid::textGrob('Protein', gp=grid::gpar(fontsize=45)))
depends.pl.plot = pl.mcmc.plot(depends, depends.pl.mcmc$acc, burn.in=500, thin.by=5, top = grid::textGrob('CRAN', gp=grid::gpar(fontsize=45)))
sim.pl.plot = pl.mcmc.plot(sim, sim.pl.mcmc$acc, burn.in=500, thin.by=5, top = grid::textGrob('BA', gp=grid::gpar(fontsize=45)))

layout_mat = t(matrix(c(1,1,2,2,
                        1,1,2,2,
                        3,3,4,4,
                        3,3,4,4,
                        NA,5,5,NA,
                        NA,5,5,NA), nrow=4, ncol=6))

pl.plots = marrangeGrob(
  list(
    tennis.pl.plot,
    harvard.pl.plot,
    protein.pl.plot,
    depends.pl.plot,
    sim.pl.plot),layout_matrix = layout_mat, top='')
ggsave(file='../poster3/plsurv.svg', pl.plots, device='svg', width=25, height=27)
