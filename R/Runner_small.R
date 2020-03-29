here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/")) 
source(paste(here,"Fedorov_cholesky.R",sep="/")) 
source(paste(here,"FedorovGA.R",sep="/")) 
#source(paste(here,"FedorovGAparallel.R",sep="/")) 

# -- Initialize DesignMatrix -----------------------------
# how many cards?
n <- 8

# create DesignMatrix object
dm <- DesignMatrix(n=n, cholesky=T)

# create attributes
dm$add_attribute(name="body", levels=3, dist=c(25,50,25))
dm$add_attribute(name="engine", levels=2, dist=c(50,50))
dm$add_attribute(name="seats", levels=4, dist=c(25,25,25,25))

# generate
dm$generate_design()

# add interacts
dm$add_interaction(n1="body", n2="seats", l1=0, l2=3, eq=F)
dm$add_interaction("body","seats",1,0,F)

# set penalty
lmda <- 1

# view attributes of dm
dm$names
dm$islacks
dm$dslacks
dm$X

### duplicate design matrix for testing
# fedorov - deterimant
dm_fed <- dm$copy(shallow=T)
dm_fed$set_cholesky(cholesky=F)
# fedorov - cholesky
dm_chol <- dm$copy(shallow=T)
dm_chol$set_cholesky(cholesky=T)
# genetic - determinant
dm_ga <- dm$copy(shallow=T)
dm_ga$set_cholesky(cholesky=F)
# genetic - cholesky
dm_gac <- dm$copy(shallow=T)
dm_gac$set_cholesky(cholesky=T)

# check cholesky enable/disable
dm_fed$cholesky
dm_chol$cholesky
dm_ga$cholesky
dm_gac$cholesky

### test initial optimality
doptimality(dm, lambda=lmda, how='det')
doptimality(dm, lambda=lmda, how='chol')
sumfisherz(dm, lambda=lmda)


### generate candidate set 
library(AlgDesign)
system.time(
  candidate_set <- gen.factorial(dm$levels, factors="all")
  # full 16-attribute candidate_set generation takes ~3 minutes on fast computer
)

# convert to numeric matrix
indx <- sapply(candidate_set, is.factor)
candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
candidate_set <- as.matrix(candidate_set)
# zero base
candidate_set <- candidate_set-1


## -- FEDOROV --------------------------------------------
f_time <- system.time({
  f_DM <- fedorov(dm_fed, candidate_set, n, lambda=lmda, iter=20)
})
# f_DM
f_DM$X
doptimality(f_DM, lambda=lmda, how='det')

# force update L to confirm doptimality w/ cholesky method
f_DM$L <- try(t(chol(t(f_DM$X)%*%f_DM$X)))
doptimality(f_DM, lambda=lmda, how='chol')
sumfisherz(f_DM, lambda=lmda)

## -- FEDOROV + CHOLESKY -----------------------------------
### generate candidate set 
fc_time <- system.time({
  fc_DM <- fedorov_chol(dm_chol, candidate_set, n, lambda=lmda, iter=20)
})
# fc_DM
doptimality(fc_DM, lambda=lmda, how='det')

doptimality(fc_DM, lambda=lmda, how='chol')
sumfisherz(fc_DM, lambda=lmda)

## -- GENETIC -----------------------------------
ga_time <- system.time({
  ga_DM <- gen_alg(dm_ga, pop=16, gens=1000, test='doptimality', lambda=lmda, how='det')
})
# ga_DM
ga_DM$X
doptimality(ga_DM, lambda=lmda, how='det')

# force update L to confirm doptimality w/ cholesky method
ga_DM$L <- try(t(chol(t(ga_DM$X)%*%ga_DM$X)))
doptimality(ga_DM, lambda=lmda, how='chol')

sumfisherz(ga_DM, lambda=lmda)


## -- GENETIC + CHOLESKY -----------------------------------
gac_time <- system.time({
  gac_DM <- gen_alg(dm_gac, pop=16, gens=1000, test='doptimality', lambda=lmda, how='chol')
})
# ga_DM
gac_DM$X
doptimality(gac_DM, lambda=lmda, how='det')
doptimality(gac_DM, lambda=lmda, how='chol')
sumfisherz(gac_DM, lambda=lmda)


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

print('Fedorov - GA - Cholesky')
gac_time
doptimality(gac_DM, lambda=lmda, how='det')
doptimality(gac_DM, lambda=lmda, how='chol')
sumfisherz(gac_DM, lambda=lmda)
