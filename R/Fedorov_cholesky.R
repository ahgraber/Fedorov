#-- Fedorov with Choleksy --------------------------
fedorov_chol <- function(dm, candidate_set, n, lambda=0, iter=100) {
  # fedorov algorithm find d-optimal design
  # params:
    # dm: Design Matrix object
    # candidate_set: matrix of all possible combinations of attributes
    # n: number of rows
    # lambda: weight for slack penalties
  # returns: optimal design
  
  # set inital values of algorithm
  n_iter <- 1
  obj_delta_best <- .0001
  obj <- doptimality(dm, lambda, how='chol')
  
  # iterate until the improvement in D-optimality is minimal or 100 iterations is reached
  while ((obj_delta_best > 10e-6) && (n_iter < iter)) {

    dm_best <- NULL
    obj_delta_best <- 0

    iter_time <- system.time({
      for (i in 1:n) {                      # iterate through rows in design matrix
        for (j in 1:nrow(candidate_set)) {  # iterate through rows in candidate set 
  
          # create test case
          dm_test <- dm$copy(shallow=TRUE)
          dm_test$add_row(candidate_set[j,])
          dm_test$del_row(i)
  
          # calculate change from the row swap
          # obj_test <- det_chol(dm_test$L) - penalty(dm_test, lambda)
          obj_test <- doptimality(dm_test, lambda, how='chol')
          obj_delta <- obj_test - obj
  
          # if that is greater than the best candidate so far, make it the new best pair
          if (obj_delta > obj_delta_best) {
            dm_best <- dm_test$copy(shallow=TRUE)
            obj_delta_best <- obj_delta
            obj_best <- obj_test
  
            print(paste(paste("Iteration", n_iter, "candidate_set swap", i, sep=" "), obj_delta_best, obj_best, sep=" | "))
          } else {
            next
          }
          try(rm(dm_test), silent=T)
          try(rm(obj_test), silent=T)
          try(rm(obj_delta), silent=T)
  
        } # end for j
      } # end for i
  
      ### updates
      if (is.null(dm_best)) {
        # no better swaps found
        break
      } else {
        # retain swap
        dm <- dm_best$copy(shallow=TRUE)
        
        # update doptimality following the swap
        obj <- obj_best
        
        n_iter <- n_iter+1
      }
    })
  print(paste("Iteration", n_iter-1, "in", iter_time, "seconds", sep=" "))    
  } # end while
  
  if (n_iter == iter) {
    print("Algorithm stopped due to reaching iteration limit")
  } else {
    print(paste("Convergence achieved in ",n_iter," iterations"))
  }
  print(paste("DEBUG: internal objfun - ",obj))
  return(dm)
} # end fedorov

# asdf