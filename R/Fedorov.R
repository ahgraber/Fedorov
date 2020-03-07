#-- Helper functions -----------------------------------------------------------
# d <- function(D, row) {
#   # function to calculate variance estimator
#   # params:
#   # D: current design matrix
#   # row: row to evaluate wrt. design
#   estimator <- row %*% solve( t(D)%*%D ) %*% t(row)
#   return(estimator)
# }

var_est <- function(D,i_row,j_row) {
  # variance estimator for swapping rows
  # params:
  # D: current design matrix
  # i_row: leaving row
  # j_row: entering row
  
  # attempting protection from singular matrix inversions
  X <- tryCatch(
    expr = { solve( t(D)%*%D ) },
    error = function(e) { return( MASS::ginv( t(D)%*%D ) ) }
  )
  est <- j_row %*% X %*% i_row
  return(est)
  # return(j %*% solve( t(as.matrix(D))%*%as.matrix(D) ) %*% t(as.matrix(i)))
} # end var_est

delta_var <- function(D, i, j, lambda) {
  # variance estimator for swapping rows
  # params:
  # D: current design matrix
  # i: leaving row
  # j: entering row
  # lambda: penalty for slacks
  # returns: variance estimator
  
  est <- var_est(D,j,j) - ( var_est(D,j,j)*var_est(D,i,i)-var_est(D,i,j)^2 ) - var_est(D,i,i)
  return(est)
} # end delta_var

update_obj <- function(D, i, j_row, lambda, det, dvar) {
  # calculates the % increase in the objective function for the row swap
  # params:
  # D: current design matrix
  # i: leaving row *index*
  # j: entering row
  # lambda: penalty for slacks
  # det: determinant
  # dvar: delta variance estimator
  # returns: updated doptimality metric
  
  # calculate penalties
  old_p <- penalty(dm, D, lambda)
  i_row <- D[i,]
  dm$del_row(i)
  dm$add_row(j_row)
  dm$update_slacks()
  new_p <- penalty(dm, dm$X, lambda)
  
  # revert dm$X
  dm$X <- D 
  return( (det*(1+dvar)-new_p) / (det-old_p) - 1 )
} # end update_obj

#-- Fedorov --------------------------
fedorov <- function(dm, candidate_set, n, lambda=0) {
  # fedorov algorithm find d-optimal design
  # params:
    # dm: Design Matrix object
    # candidate_set: matrix of all possible combinations of attributes
    # n: number of rows
    # lambda: weight for slack penalties
  # returns: optimal design
  
  # set inital values of algorithm
  iter <- 1
  obj_delta_best <- .0001
  det <- as.numeric(determinant(t(dm$X)%*%dm$X)$modulus)
  
  # iterate until the improvement in D-optimality is minimal or 100 iterations is reached
  while ((obj_delta_best > 10e-6) && (iter < 100)) {

    i_best <- NULL
    j_best <- NULL
    obj_delta_best <- 0
    
    for (i in 1:n) {                      # iterate through rows in design matrix
      for (j in 1:nrow(candidate_set)) {  # iterate through rows in candidate set
        
        dvar <- NULL
        obj_delta <- NULL    

        # calculate the potential improvement in D-optimality by replacing the
        # current row in the design with the current row in the candidate set
        dvar <- delta_var(dm$X, dm$X[i,], candidate_set[j,])
        # delta <- delta_var(dm$X, t(dm$X[i,]), as.matrix(candidate_set[j,]))
        obj_delta <- update_obj(dm$X, i, candidate_set[j,], lambda, det, dvar)

        # if that is greater than the best candidate so far, make it the new best pair
        if (obj_delta > obj_delta_best) {
          obj_delta_best <- obj_delta
          dvar_best <- dvar
          i_best <- i
          j_best <- j
          print(paste(paste("iteration", iter, sep=" "), obj_delta_best, dvar_best, sep=" | "))
        } else {
          next
        }
      } # end for j
    } # end for i

    ### updates
    if (is.null(i_best)) {
      # no better swaps found
      break
    } else {
      # perform swap
      dm$del_row(i_best)
      dm$add_row(candidate_set[j_best,])
      dm$update_slacks()
      
      # update the determinate following the swap
      det <- det*(1+dvar_best)
      
      iter <- iter+1
    }
    
  } # end while
  
  if (iter == 100) {
    print("Algorithm stopped due to reaching iteration limit")
  } else {
    print(paste("Convergence achieved in ",iter," iterations"))
  }
  return(dm$X)
} # end fedorov

# asdf