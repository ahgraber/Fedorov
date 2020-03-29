here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/")) 
source(paste(here,"Fedorov_cholesky.R",sep="/")) 
source(paste(here,"FedorovGA.R",sep="/")) 
# source(paste(here,"FedorovGAparallel.R",sep="/"))

# -- Initialize DesignMatrix -----------------------------
# how many cards?
n <- 40

### create DesignMatrix object
dm <- DesignMatrix(n=n, cholesky=T)

# create attributes
dm$add_attribute(name="age", levels=3, dist=c(25,50,25))
dm$add_attribute(name="gender", levels=2, dist=c(50,50))
dm$add_attribute(name="bmi", levels=3, dist=c(20,40,40))
dm$add_attribute(name="race", levels=4, dist=c(40,20,20,20))
dm$add_attribute(name="dia", levels=2, dist=c(50,50))
dm$add_attribute(name="angn", levels=3, dist=c(50,25,25))
dm$add_attribute(name="ldl", levels=3, dist=c(25,50,25))
dm$add_attribute(name="bp", levels=3, dist=c(25,50,25))
# dm$add_attribute(name="a1c", levels=4, dist=c(25,25,25,25))
# dm$add_attribute(name="ren", levels=3, dist=c(50,25,25))
# dm$add_attribute(name="srcl", levels=3, dist=c(50,25,25))
# dm$add_attribute(name="uacr", levels=3, dist=c(33,33,34))
# dm$add_attribute(name="strk", levels=2, dist=c(75,25))
# dm$add_attribute(name="ptx", levels=4, dist=c(25,25,25,25))
# dm$add_attribute(name="hist", levels=3, dist=c(50,25,25))
# dm$add_attribute(name="smk", levels=3, dist=c(50,25,25))

# enable/disable cholesky functionality
dm$set_cholesky(cholesky=T)

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
dm_fed <- dm$copy(shallow=T)
dm_fed$set_cholesky(cholesky=F)

dm_chol <- dm$copy(shallow=T)
dm_ga <- dm$copy(shallow=T)

# check cholesky enable/disable
dm_fed$cholesky
dm_chol$cholesky
dm_ga$cholesky

# test initial optimality
doptimality(dm, lambda=lmda, how='det')
doptimality(dm, lambda=lmda, how='chol')
sumfisherz(dm, lambda=lmda)

## -- CANDIDATE SET -------------------------------------
library(AlgDesign)
system.time(
  candidate_set <- gen.factorial(dm$levels, factors="all")
  # 12-attribute candidate_set generation takes ~3 secconds on fast computer
  # full 16-attribute candidate_set generation takes ~3-5 minutes on fast computer
)

# convert to numeric matrix
indx <- sapply(candidate_set, is.factor)
candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
candidate_set <- as.matrix(candidate_set)
# zero base
candidate_set <- candidate_set-1

## -- FEDOROV --------------------------------------------
f_time <- system.time(
  f_DM <- fedorov(dm_fed, candidate_set, n, lambda=lmda) # iter=20
)
# f_DM
f_DM$X
doptimality(f_DM, lambda=lmda, how='det')
# force update L to confirm doptimality w/ cholesky method
f_DM$L <- try(t(chol(t(f_DM$X)%*%f_DM$X)))
doptimality(f_DM, lambda=lmda, how='chol')
sumfisherz(f_DM, lambda=lmda)up

## -- FEDOROV + CHOLESKY -----------------------------------
fc_time <- system.time(
  fc_DM <- fedorov(dm_chol, candidate_set, n, lambda=lmda) # iter=20
)
# fc_DM
fc_DM$X
doptimality(fc_DM, lambda=lmda, how='det')
doptimality(fc_DM, lambda=lmda, how='chol')
sumfisherz(fc_DM, lambda=lmda)

## -- GA -----------------------------------
ga_time <- system.time({
  ga_DM <- gen_alg(dm_ga, pop=16, gens=1000, test='doptimality', lambda=lmda)
})
# ga_DM
ga_DM$X
doptimality(ga_DM, lambda=lmda, how='det')
doptimality(ga_DM, lambda=lmda, how='chol')
sumfisherz(ga_DM, lambda=lmda)

## -- SUMMARY -----------------------------------
print('Fedorov')
f_time
doptimality(f_DM, lambda=lmda, how='det')
doptimality(f_DM, lambda=lmda, how='chol')
sumfisherz(f_DM, lambda=lmda)

print('Fedorov - Cholesky')
fc_time
doptimality(fc_DM, lambda=lmda, how='det')
doptimality(fc_DM, lambda=lmda, how='chol')
sumfisherz(fc_DM, lambda=lmda)

print('Fedorov - GA')
ga_time
doptimality(ga_DM, lambda=lmda, how='det')
doptimality(ga_DM, lambda=lmda, how='chol')
sumfisherz(ga_DM, lambda=lmda)