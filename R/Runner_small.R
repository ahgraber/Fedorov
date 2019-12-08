here <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(here,"FedorovDesignClass.R",sep="/")) # ignore error
source(paste(here,"Fedorov.R",sep="/")) 
# source(paste(here,"FedorovGA.R",sep="/")) 
source(paste(here,"FedorovGAparallel.R",sep="/")) 


# how many patients?
n <- 8

# create DesignMatrix object
dm <- design(n=n)

# create attributes
dm$add_attribute(name="age", levels=3, dist=c(25,50,25))
dm$add_attribute(name="sex", levels=2, dist=c(50,50))
dm$add_attribute(name="bmi", levels=3, dist=c(33,33,34))

# generate
dm$generate()

# add interacts
dm$add_interaction(n1="age", n2="sex", l1=1, l2=1, kind=T)
dm$add_interaction("age","bmi",0,2,F)


# view attributes of dm
dm$names
dm$values
dm$islacks
dm$dslacks
dm$X

# test optimality
doptimality(dm, dm$X)


X <- gen_alg(dm, pop=50, gens=1000, test='doptimality')
