#-- Helper functions -----------------------------------------------------------
var_est <- function(D, i, j) {
  # variance estimator for swapping rows
  # params:
  # D: current design matrix
  # i: leaving row
  # j: entering row
  
  # attempting protection from singular matrix inversions
  X <- tryCatch(
    # expr = { solve( t(D)%*%D ) },
    # error = function(e) { return( MASS::ginv( t(D)%*%D ) ) }
    expr = { solve( crossprod(D) ) },
    error = function(e) { return( MASS::ginv( crossprod(D) ) ) }
  )
  est <- j %*% X %*% i

  return(est)
  # return(j %*% solve( t(as.matrix(D))%*%as.matrix(D) ) %*% t(as.matrix(i)))
} # end var_est

delta_var <- function(D, i, j, lambda) {
  # calculates the difference in variance for swapping rows
  # params:
  # D: current design matrix
  # i: leaving row
  # j: entering row
  # lambda: penalty for slacks
  # returns: variance estimator
  Dii <- var_est(D,i,i)
  Djj <- var_est(D,j,j)
  Dij <- var_est(D,i,j)

  est <- Djj - ( (Djj * Dii) - Dij^2 ) - Dii

  return(est)
} # end delta_var

update_obj <- function(dm_old, dm_new, lambda, det, dvar) {
  # calculates the % change in the objective function for the row swap
  # params:
  # dm_old: DesignMatrix object with before swap
  # dm_new: DesignMatrix object with rows swapped
  # lambda: penalty for slacks
  # det: determinant of dm_old
  # dvar: delta variance estimator to update det
  # returns: updated doptimality metric

  # calculate penalties
  old_p <- penalty(dm, lambda)
  new_p <- penalty(dm_new, lambda)

  dopt <- (det*(1+dvar)-new_p) / (det-old_p) - 1 
  if (is.na(dopt)) {
    dopt <- 0
  }
  return(dopt)
} # end update_obj

#-- Fedorov --------------------------
fedorov <- function(dm, candidate_set, n, lambda=0, iter=100, return_iter=FALSE, debug=FALSE) {
  # fedorov algorithm find d-optimal design
  # params:
    # dm: Design Matrix object
    # candidate_set: matrix of all possible combinations of attributes
    # n: number of rows
    # lambda: weight for slack penalties
    # iter: maximal number of iterations
    # return_iter: whether to include total number of iterations in run in output
    # debug: whether to print debugging logs
  # returns: optimal design
  
  # set inital values of algorithm
  n_iter <- 1
  obj_delta_best <- .0001
  obj <- doptimality(dm, lambda, how='det')

  # preserve deterimant for use in estimator
  det <- as.numeric(determinant(crossprod(dm$X) )$modulus)
  # det <- as.numeric(determinant(t(dm$X)%*%dm$X)$modulus)
  
  # iterate until the improvement in D-optimality is minimal or 100 iterations is reached
  while ((obj_delta_best > 10e-6) && (n_iter < iter)) {

    iter_time <- NULL
    dm_best <- NULL
    obj_delta_best <- 0

    iter_time <- system.time({
      for (i in 1:n) {                      # iterate through rows in design matrix
        for (j in 1:nrow(candidate_set)) {  # iterate through rows in candidate set
  
          # create test case
          dm_test <- dm$copy(shallow=TRUE)
          dm_test$add_row(candidate_set[j,])
          dm_test$del_row(i)
  
          # calculate the potential % improvement in D-optimality from the row swap
          dvar <- delta_var(dm$X, dm$X[i,], candidate_set[j,])
          obj_delta <- update_obj(dm, dm_test, lambda, det, dvar)
  
          # calculate the objective function using variance estimators
          pnlt <- penalty(dm_test, lambda)
          # det_est <- det*(1+dvar)
          det_est <- det*(1+dvar) + pnlt
          obj_test <- (100 * det_est^(1/ncol(dm_test$X)))/ nrow(dm_test$X) - pnlt
          # print(paste('Debug: actual objfun: ', doptimality(dm_test, lambda, how='chol')))
          # print(paste('Debug: estimated objfun: ', obj_test))
  
          # if swapped design is better than the best candidate so far, make it the new best design
          if (obj_delta > obj_delta_best) {
            dm_best <- dm_test$copy(shallow=TRUE)
            dvar_best <- dvar
            obj_delta_best <- obj_delta
            obj_best <- obj_test
            
            if (debug) {
              # print(paste(paste("Iteration", n_iter, "candidate_set swap", i, sep=" "), dvar_best, obj_delta_best, obj_best, sep=" | "))
              print(
                paste(
                  paste("Iteration", n_iter, "design", i, "candidate", j, sep=" "), 
                  round(obj_delta_best,5), 
                  round(obj_best,5), 
                  sep=" | "
                )
              )              
            }

          } else {
            next
          }
          try(rm(dm_test), silent=T)
          try(rm(dvar), silent=T)
          try(rm(det_est), silent=T)
          try(rm(obj_test), silent=T)
          try(rm(obj_delta), silent=T)
        } # end for j
      } # end for i
    }) # end system.time
    
    if (debug) {
      print(paste("Iteration", n_iter, "in", round(iter_time[3],4), "seconds", sep=" "))
    }

    ### updates
    if (is.null(dm_best)) {
      # no better swaps found
      break
    } else {
      # retain swap
      dm <- dm_best$copy(shallow=TRUE)

      # update the determinant following the swap
      # det <- det*(1+dvar_best)
      det <- det*(1+dvar_best) + penalty(dm_best, lambda)
      obj <- obj_best
      
      n_iter <- n_iter+1
    } # end updates

  } # end while
  
  if (n_iter == iter) {
    print("Algorithm stopped due to reaching iteration limit")
  } else {
    print(paste("Convergence achieved in ",n_iter," iterations"))
  }
  if (debug) {
    print(paste("DEBUG: internal objfun - ",obj))
  }
  
  if (return_iter) {
    return(list(dm, n_iter))
  }
  else {
    return(dm)
  }
} # end fedorov

# asdf