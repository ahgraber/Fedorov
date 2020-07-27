#-- supporting functions ---------------------------------------
breed <- function(p1, p2){
  # breeding function for genetic algorithm; approx 50/50 mix of parents
  # params:
    # p1: DesignClass object (parent 1)
    # p2: DesignClass object (parent 2)
  # returns:
    # list(c1,c2) where c1,c2 are child DesignClass objects of p1,p2
  
  if (nrow(p1$X) != nrow(p2$X)) {stop('Parents do not have equivalent rows')}
  if (ncol(p1$X) != ncol(p2$X)) {stop('Parents do not have equivalent columns')}
  
  # DEPRECATED: use cholesky functionality instead
  # probsA <- runif(nrow(p1$X),0,1)
  # ax <- p1$X[probsA > 0.5,]  # chromosomes from parent X
  # ay <- p2$X[probsA <= 0.5,] # chromosomes from parent Y
  # A <- rbind(ax,ay)

  # probsB <- runif(nrow(X),0,1)
  # bx <- p2$X[probsB > 0.5,]  # chromosomes from parent X
  # by <- p2$Y[probsB <= 0.5,] # chromosomes from parent Y
  # # cholesky update
  # B <- rbind(bx,by)

  # cholesky update:
  # find which is longer, copy that parent
  # then add new rows from other parent 
  # and delete rows not inherited from first
  probs_c1 <- runif(nrow(p1$X),0,1)
  p1_rowidx <- which(probs_c1 > 0.5)
  p2_rowidx <- which(probs_c1 <= 0.5)
  if (length(p1_rowidx) > length(p2_rowidx)) {
    # if p1 inheritence > p2 inheritence, copy p1
    c1 <- p1$copy() 
    for (i in p2_rowidx) {
      # add rows inherited from p2
      c1$add_row(p2$X[i,])
      # remove rows NOT inherited from p1
      c1$del_row(i)
    }
  } else {
    # if p1 inheritence < p2 inheritence, copy p2
    c1 <- p2$copy()
    for (i in p1_rowidx) {
      # add rows inherited from p2
      c1$add_row(p1$X[i,])
      # remove rows NOT inherited from p1
      c1$del_row(i)
    }
  }

  probs_c2 <- runif(nrow(p2$X),0,1)
  p1_rowidx <- which(probs_c2 > 0.5)
  p2_rowidx <- which(probs_c2 <= 0.5)
  if (length(p1_rowidx) > length(p2_rowidx)) {
    # if p1 inheritence > p2 inheritence, copy p1
    c2 <- p1$copy() 
    for (i in p2_rowidx) {
      # add rows inherited from p2
      c2$add_row(p2$X[i,])
      # remove rows NOT inherited from p1
      c2$del_row(i)
    }
  } else {
    # if p1 inheritence < p2 inheritence, copy p2
    c2 <- p2$copy()
    for (i in p1_rowidx) {
      # add rows inherited from p2
      c2$add_row(p1$X[i,])
      # remove rows NOT inherited from p1
      c2$del_row(i)
    }
  }

  return(list(c1, c2))
} # end breed

mutate <- function(dm, alpha) {
  # mutation function for genetic algorithm
  # params:
  # dm: Design Matrix object
  # alpha: num 0-1 indicating likelihood for mutation (lower increases mutation)
  # returns: mutated matrix X
  
  X <- dm$X
  row_test <- runif(nrow(X), 0, 1)
  row_mask <- matrix(rep(row_test>=alpha, ncol(X)), 
                     nrow(X), 
                     ncol(X), 
                     byrow=F)
  
  cell_test <- matrix(runif(length(X), 0, 1), 
                      nrow(X), 
                      ncol(X), 
                      byrow=T)
  cell_mask <- cell_test >= alpha
  
  mutation_mask <- row_mask & cell_mask

  # generate possible mutation for every cell
  bulk_mutations <- sapply(X=dm$levels, 
                           FUN=function (y) {(sample(c(1:y-1), nrow(X), replace=T))} 
                          )

  # apply mutations to edit matrix
  X[mutation_mask] <- bulk_mutations[mutation_mask]

  # update DesignMatrix
  row_idx <- which(row_mask[,1])
  for (i in row_idx) {
    # add mutation as new row
    dm$add_row(X[i,])
    # remove unmuated row
    dm$del_row(i)
  }

  return(dm)
  # return(list(X, score)) # return updated matrix AND updated fitness??
  
} # end mutate

cull <- function(elite, stock, children, pop, dir) {
  # function to reduce population back down to pop
  # params:
    # elite: list of elite design matrices to be preserved
    # stock: list of non-elite parents design matrices
    # children: list of new design matrices
    # pop: int, population size to achieve
    # dir: direction to sort
  # returns: list (length pop) of best design matrices
  
  # combine stock and child lists
  fill <- append(stock, children) 
  fill <- sorter(fill, dir)
  
  # find number of herd to fill after elites are kept
  nfill <- pop-length(elite) 
  
  # generate new herd with elites and best of rest
  herd <- append(elite, head(fill, nfill))

  return(sorter(herd, dir))
} # end cull

sorter <- function(herd, dir) {
  # function to sort the herd based on objective function
  # params:
  # herd: list of (dval, matrix) tuples to be sorted
  # dir: direction to sort
  # returns: sorted herd
  
  if (dir=="min") {
    # value <- sapply(herd, function(x) x[[1]])
    # herd[order(value, decreasing=T)]
    return(herd[order(sapply(herd, function(x) x[[1]]), decreasing=F)])
  } else if (dir=="max") {
    # value <- sapply(herd, function(x) x[[1]])
    # herd[order(value, decreasing=T)]
    return(herd[order(sapply(herd, function(x) x[[1]]), decreasing=T)])
  } else {
    stop("Direction for objective function not defined")
  }
} # end sorter

#-- genetic algorithm  ---------------------------------------
gen_alg <- function(dm, pop, gens, test, alpha=0.25, lambda=0, how='chol', return_iter=FALSE, debug=FALSE) {
  # genetic algorithm to find d-optimal design
  # params:
    # dm: Design Matrix object
    # pop: population size
    # gens: int, maximum number of generations
    # test: objective function to use
    # lambda: weight for slack penalties
    # how: use cholesky updates
    # return_iter: whether to include total number of iterations in run in output
    # debug: whether to print debugging logs
  # returns: optimal design
  
  if (pop < 16){ stop('Suggested population minimum of 16') }
  if (gens < 100){ stop('Suggested generation minimum of 100') }
  
  if (test=="doptimality"){
    objfun <- doptimality
    dir <- "max"
  } else if (test=="sumfisherz") {
    objfun <- sumfisherz
    dir <- "min"
  } else {
    stop("Test value not in c('doptimality','sumfisherz')")
  }

  alpha <- 0.17 # probability threshold; lower increases variation/mutation

  ### create herd (list of (dval, DesignMatrix object) tuples)
  if (debug) {
    print('Creating Generation 0')
  }
  herd <- list()
  for (p in 1:pop) {
    new <- dm$copy(shallow=TRUE)
    new$generate_design()
    herd[[p]] <- list(objfun(new, lambda, how), new)
  }
  herd <- sorter(herd, dir)
  top <- herd[[1]]
  
  ### pick elite designs to leave unchanged/unculled
  if (pop < 32) { 
    nelite <- 2 
  } else { 
    nelite <- 4
  }
  
  g <- 1  # initialize iterator
  converge <- 0 # initialize convergence criteria counter
  while ((g < gens) && (converge < log2(gens)+10)) {
    # stop if reach maximum generations OR 
    # if difference between top designs remains small for some number of generations
    
    gen_time <- system.time({
      if ((pop %% 2) != 0) { nelite <- nelite-1} # adjust for odd population
      elite <- head(herd, nelite)
      stock <- tail(herd, -nelite)
      
      ### breed randomized pairs
      x <- sample(c(1:length(stock)))
      y <- sample(c(1:length(stock)))
      children <- list()
      for (i in 1:length(stock)) { 
        if (x[i] != y[i]) { # no self-replication
          if (runif(1,0,1) >= alpha){
            # if test passed, breed & save children
            kids <- breed(
              stock[[ x[i] ]][[2]], 
              stock[[ y[i] ]][[2]]
            ) 
            children[[length(children)+1]] <- kids[[1]]
            children[[length(children)+1]] <- kids[[2]]
          }
        }
      } # end for i (breed)
      
      ### mutation
      for (j in 1:length(children)) { 
        if (runif(1,0,1) >= alpha) {
          # if test passed, mutate child
          children[[j]] <- mutate(children[[j]], alpha)
        }
      } # end for j (mutate)
  
      ### assess fitness
      for (k in 1:length(children)) {
        children[[k]] <- list(objfun(children[[k]], lambda, how), children[[k]])
      }
      
      ### cull
      herd <- cull(elite, herd, children, pop, dir)
  
      ### updates
      g <- g+1
      if (dir == "max") {
        if ((herd[[1]][[1]]-top[[1]]) < 10e-6 ) {
          # maximizing, so change should be positive
          # if the change in objval (new - old) 0 and small pos number, system is converging
            # smallest possible change is 0 (i.e., same best design) b/c preserving elites
          # if converging, count as a converge step to potentially break out of while loop
          top <- herd[[1]][[1]]
          converge <- converge+1
        } else {
          # if not a converge step, reset converge coutner to 0
          converge <- 0
        }
      } else if (dir == "min") {
        if ((herd[[1]][[1]]-top[[1]]) > -10e-6 ) {
          # minimizing, so change should be negative
          # if the change in objval (new - old) small neg number and 0, system is converging
            # largest possible change is 0 (i.e., same best design) b/c preserving elites
          # if converging, count as a converge step to potentially break out of while loop
          top <- herd[[1]][[1]]
          converge <- converge+1
        } else {
          # if not a converge step, reset converge coutner to 0
          converge <- 0
        }
      }
  
      # top <- herd[[1]][[1]]
    })
    if (debug) {
      print(paste(paste("Generation", g, "in", round(gen_time[3],4), "seconds", sep=" "), round(top,5), sep=" | "))
    }
  } # end while
  
  print(paste("Convergence achieved in ", g," iterations"))
  
  if (return_iter) {
    return(list(herd[[1]][[2]], g))
  } else {
    return(herd[[1]][[2]])
  }
  
} # end gen_alg
