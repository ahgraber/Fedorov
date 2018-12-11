here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/"))
# source(paste(here,"Fedorovparallel.R",sep="/"))

# how many patients?
n <- 8 # toy problem
# n <- 50

# how to penalize slacks?
lambda = 1

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
dm$add_attribute(name="ldl", levels=3, dist=c(25,50,25))
dm$add_attribute(name="bp", levels=3, dist=c(25,50,25))
dm$add_attribute(name="a1c", levels=4, dist=c(25,25,25,25))
dm$add_attribute(name="ren", levels=3, dist=c(50,25,25))
dm$add_attribute(name="srcl", levels=3, dist=c(50,25,25))
dm$add_attribute(name="uacr", levels=3, dist=c(33,33,34))
dm$add_attribute(name="ptx", levels=4, dist=c(25,25,25,25))
dm$add_attribute(name="hist", levels=3, dist=c(50,25,25))
dm$add_attribute(name="smk", levels=3, dist=c(50,25,25))

# generate
dm$generate()

# add interaction
dm$add_interaction(n1="dia", n2="a1c", l1=0, l2=3, eq=F) # nondiabetics can't have high a1c
dm$add_interaction(n1="dia", n2="a1c", l1=1, l2=0, eq=F) # diabetics can't have low a1c
dm$X
dm$islacks

### generate candidate set 
library(AlgDesign)
system.time(
  candidate_set <- gen.factorial(dm$levels, factors="all")
  # full 16-attribute candidate_set generation takes ~3 minutes on fast computer
)
# bkup <- candidate_set

# convert to numeric matrix
indx <- sapply(candidate_set, is.factor)
candidate_set[,indx] <- lapply(candidate_set[,indx], function(x) as.numeric(as.character(x)))
candidate_set <- as.matrix(candidate_set)
# zero base
candidate_set <- candidate_set-1

# # randomly choose 1% of rows for proof of concept
# indx <- as.numeric(sample(candidate_set, nrow(candidate_set)*.01)) 
# candidate_set <- candidate_set[indx,]
# rm(indx)
#  

# initial optimality
doptimality(dm, dm$X, lambda=lambda)

### test optimality
system.time(
X <- fedorov(dm, candidate_set, n, lambda=lambda)
)
X
doptimality(dm, X, lambda=lambda)
sumfisherz(dm, X)

# 
# ### test optimality (parallelized)
# system.time(
#   X <- fedorovpar(dm, candidate_set, n, lambda=lambda)
# )
# # notes: throws error trying to pass the DesignMatrix object.  Might be error in foreach package?
# X
# doptimality(dm, X, lambda=lambda)
# sumfisherz(dm, X)