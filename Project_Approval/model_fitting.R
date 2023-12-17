source('../scripts/functions.R')
source('../scripts/new_functions.R')
library(networkdata)
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


G = barabasi.game(3e4)
x = degree(G, mode='in')
x=x[x>0]
sim = as.data.frame(table(x))
sim[,1] = as.numeric(as.character(sim[,1]))

#PLIGP######################################
{
#tennis data


init = list(alpha=1, shape=1, scale=5, threshold=15)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=0.5)

tennis.out = pli.mcmc(1e4,tennis,init,prior.params,init.cov, cov.period = 1e7, type=1, show=F)
#harvard data


init = list(alpha=1, shape=1, scale=5, threshold=15)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=1)

harvard.out = pli.mcmc(1e5,harvard,init,prior.params,init.cov, cov.period = 1e7, type=1, show=F)
#protein data


init = list(alpha=1, shape=1, scale=5, threshold=15)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=1)

protein.out = pli.mcmc(1e4,protein.dat,init,prior.params,init.cov, cov.period = 1e7, type=1, show=F)
#CRAN data



init = list(alpha=1, shape=1, scale=5, threshold=15)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=1)

depends.out = pli.mcmc(1e4,depends,init,prior.params,init.cov, cov.period = 1e7, type=1, show=T)

#simulated data

# init = list(alpha=1, shape=1, scale=5, threshold=15)
# prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
# init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
#                                                -0.0003,1),ncol=2), threshold=1)
# 
# sim.out = pli.mcmc(1e4,sim,init,prior.params,init.cov, cov.period = 1e7, type=1, show=F)


# -------------------------------------------------------------------------

tennis.plot = plot.pli.mcmc(tennis, tennis.out$states)
harvard.plot = plot.pli.mcmc(harvard, harvard.out$states)
protein.plot = plot.pli.mcmc(protein.dat, protein.out$states)
depends.plot = plot.pli.mcmc(depends, depends.out$states)

}
#Power Law




tennis.out = pl.mcmc(1e5, tennis, 1, 0.01, S=5e3)
plot(tennis.out$acc, type='l')
plot(tennis.out$var, type='l')


sim.out = pl.mcmc(1e5, sim, 0.1, 0.01, S=1e3)
plot(sim.out$acc, type='l')
plot(sim.out$var, type='l')



protein.out = pl.mcmc(1e5, protein.dat, 1, 0.01, S=1e3)
plot(protein.out$acc, type='l')
plot(protein.out$var, type='l')
# -------------------------------------------------------------------------

init = list(alpha=1, beta=1, threshold=5)
prior.params = list(alpha=0.01, beta=0.01, threshold = 0.01)
cov.init = list(ab = diag(c(0.1,0.1)), threshold=0.5)

tennis.out = plpl.mcmc(1e4, tennis,init ,
                       prior.params,
                       cov.init)



plot(as.mcmc(tennis.out))


X = tennis
mcmc.out.states = tennis.out





plot.plpl.mcmc(tennis, tennis.out)


 












