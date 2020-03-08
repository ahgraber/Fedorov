here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/")) 
source(paste(here,"Fedorov_cholesky.R",sep="/")) 
source(paste(here,"FedorovGA.R",sep="/")) 
# source(paste(here,"FedorovGAparallel.R",sep="/"))

# -- Initialize DesignMatrix -----------------------------
# how many cards?
n <- 50

### create DesignMatrix object
dm <- design(n=n)

# create attributes
dm$add_attribute(name="age", levels=3, dist=c(25,50,25))
dm$add_attribute(name="gender", levels=2, dist=c(50,50))
dm$add_attribute(name="bmi", levels=3, dist=c(20,40,40))
dm$add_attribute(name="race", levels=4, dist=c(40,20,20,20))
dm$add_attribute(name="dia", levels=2, dist=c(50,50))
dm$add_attribute(name="strk", levels=2, dist=c(75,25))
dm$add_attribute(name="angn", levels=3, dist=c(50,25,25))
#dm$add_attribute(name="ldl", levels=3, dist=c(25,50,25))
#dm$add_attribute(name="bp", levels=3, dist=c(25,50,25))
dm$add_attribute(name="a1c", levels=4, dist=c(25,25,25,25))
#dm$add_attribute(name="ren", levels=3, dist=c(50,25,25))
#dm$add_attribute(name="srcl", levels=3, dist=c(50,25,25))
#dm$add_attribute(name="uacr", levels=3, dist=c(33,33,34))
#dm$add_attribute(name="ptx", levels=4, dist=c(25,25,25,25))
#dm$add_attribute(name="hist", levels=3, dist=c(50,25,25))
#dm$add_attribute(name="smk", levels=3, dist=c(50,25,25))

# generate
dm$generate_design()

# add interaction
dm$add_interaction(n1="dia", n2="a1c", l1=0, l2=3, eq=F) # nondiabetics can't have high a1c
dm$add_interaction(n1="dia", n2="a1c", l1=1, l2=0, eq=F) # diabetics can't have low a1c

# how to penalize slacks?
lmda = 1

# view attributes of dm
dm$names
dm$values
dm$islacks
dm$dslacks
dm$X

# duplicate design matrix for testing
dm_ga <- dm$copy(shallow=T)

# test initial optimality
doptimality(dm, lambda=lmda, how='det')
doptimality(dm, lambda=lmda, how='chol')
sumfisherz(dm, lambda=lmda)

## -- GENETIC + CHOLESKY -----------------------------------
### GENETIC
system.time({
  ga_DM <- gen_alg(dm_ga, pop=16, gens=1000, test='doptimality', lambda=lmda)
})
# ga_DM
ga_DM$X
doptimality(ga_DM, lambda=lmda, how='det')
doptimality(ga_DM, lambda=lmda, how='chol')
sumfisherz(ga_DM, lambda=lmda)

# test fisherz
system.time(
  ga_DMf <- gen_alg(ga_DM, pop=100, gens=1000, test="sumfisherz", lambda=lmda)
)
# ga_DMf
ga_DMf$X
doptimality(ga_DMf, lambda=lmda, how='det')
doptimality(ga_DMf, lambda=lmda, how='chol')
sumfisherz(ga_DMf, lambda=lmda)


# notes:
# time doubles as you double population
# time doubles as you double generations