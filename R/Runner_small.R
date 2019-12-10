here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/")) 
source(paste(here,"FedorovGA.R",sep="/")) 
#source(paste(here,"FedorovGAparallel.R",sep="/")) 


# how many patients?
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

# preserve design matrix for testing
oldX <- dm$X

# test optimality
doptimality(dm, dm$X, lambda=lmda)


## -- GENETIC --------------------------------------------
system.time(
  gaX <- gen_alg(dm, pop=50, gens=1000, test='doptimality', lambda=lmda)
)
gaX

## -- FEDOROV --------------------------------------------
dm$X <- oldX
doptimality(dm, dm$X, lambda=lmda)

### generate candidate set 
system.time(
  candidate_set <- AlgDesign::gen.factorial(dm$levels, factors="all")
  # full 16-attribute candidate_set generation takes ~3 minutes on fast computer
)
# bkup <- candidate_set

# convert to numeric matrix
indx <- sapply(candidate_set, is.factor)
candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
candidate_set <- as.matrix(candidate_set)
# zero base
candidate_set <- candidate_set-1

### test optimality
system.time(
  fX <- fedorov(dm, candidate_set, n, lambda=lmda)
)
fX
doptimality(dm, fX, lambda=lmda)
