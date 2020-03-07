here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov_cholesky.R",sep="/")) 
source(paste(here,"FedorovGA.R",sep="/")) 
#source(paste(here,"FedorovGAparallel.R",sep="/")) 

# -- Initialize DesignMatrix -----------------------------
# how many cards?
n <- 8

# create DesignMatrix object
dm <- design(n=n)

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
dm$values
dm$islacks
dm$dslacks
dm$X

# duplicate design matrix for testing
dm_ga <- dm$copy()
dm_fed <- dm$copy()
dm_chol <- dm$copy()

# test optimality
doptimality(dm, lambda=lmda)

## -- FEDOROV --------------------------------------------
### generate candidate set 
system.time(
  candidate_set <- AlgDesign::gen.factorial(dm$levels, factors="all")
  # full 16-attribute candidate_set generation takes ~3 minutes on fast computer

  # convert to numeric matrix
  indx <- sapply(candidate_set, is.factor)
  candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
  candidate_set <- as.matrix(candidate_set)
  # zero base
  candidate_set <- candidate_set-1

  ### test optimality
  f_DM <- fedorov(dm_fed, candidate_set, n, lambda=lmda)
)
f_DM
f_DM$X
doptimality(dm_fed, lambda=lmda)


## -- FEDOROV + CHOLESKY -----------------------------------
### generate candidate set 
system.time(
  candidate_set <- AlgDesign::gen.factorial(dm$levels, factors="all")
  # full 16-attribute candidate_set generation takes ~3 minutes on fast computer

  # convert to numeric matrix
  indx <- sapply(candidate_set, is.factor)
  candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
  candidate_set <- as.matrix(candidate_set)
  # zero base
  candidate_set <- candidate_set-1

  ### test optimality
  fc_DM <- fedorov_chol(dm_chol, candidate_set, n, lambda=lmda)
)
fc_DM
fc_DM$X
doptimality(dm_chol, lambda=lmda)


## -- GENETIC + CHOLESKY -----------------------------------
system.time(
  ga_DM <- gen_alg(dm_ga, pop=50, gens=1000, test='doptimality', lambda=lmda)
)
ga_DM
ga_DM$X
doptimality(dm_ga, lambda=lmda)

