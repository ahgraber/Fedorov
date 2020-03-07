#-- Objective functions -----------------------------------------------------------
doptimality <- function(dm, lambda=0) {
  # calculates doptimality of design (and optionally penalizes distribution constraints)
  # params:
  # dm: DesignMatrix object containing attribute & constraint information
  # design: design matrix where columns are attributes and rows are patients
  # lambda: weight to penalize constraints.  lambda=0 means no distribution constraints
  # returns: d-efficiency metric
  
  # calculate slacks for the design
  dm$update_slacks()

  objective <- (100 * det_chol(dm$L) )^(1/ncol(design)) / nrow(design)
  # objective <- det( t(design)%*%design ) / nrow(design)
  penalty <- lambda*( sum(abs(unlist(dm$dslacks))) + lambda*(sum(abs(unlist(dm$islacks)))) )
  # this double-penalizes islacks b/c we really don't want impossible interactions
  return(objective - penalty)
}

sumfisherz <- function(dm, design, lambda=0) {
  # calculates the sum of the fisher z score of the absolute values of the correlation matrix
  # minimization objective function
  # params
  # dm: DesignMatrix object containing attribute & constraint information
  # design: design matrix where columns are attributes and rows are patients
  # lambda: weight to penalize constraints.  lambda=0 means no distribution constraints
  # design: design matrix where columns are attributes and rows are patients
  # returns correlation score
  
  # calculate slacks for the design
  dm$X <- design
  dm$update_slacks()

  r <- abs(cor(design))
  z <- .5*(log(1+r)/(1-r))
  objective <- sum(z[is.finite(z)])
  penalty <- lambda*( sum(abs(unlist(dm$dslacks))) + lambda*(sum(abs(unlist(dm$islacks)))) )
  return(objective + penalty)
}

#-- Helper functions -----------------------------------------------------------
# d <- function(D, row) {
#   # function to calculate variance estimator
#   # params:
#   # D: current design matrix
#   # row: row to evaluate wrt. design
#   estimator <- row %*% solve( t(D)%*%D ) %*% t(row)
#   return(estimator)
# }

# var_est <- function(D,i_row,j_row) {
#   # variance estimator for swapping rows
#   # params:
#   # D: current design matrix
#   # i_row: leaving row
#   # j_row: entering row
#   
#   # attempting protection from singular matrix inversions
#   X <- tryCatch(
#     expr = { solve( t(D)%*%D ) },
#     error = function(e) { return( MASS::ginv( t(D)%*%D ) ) }
#   )
#   est <- j_row %*% X %*% i_row
#   return(est)
#   # return(j %*% solve( t(as.matrix(D))%*%as.matrix(D) ) %*% t(as.matrix(i)))
# }

penalty <- function(dm, lambda) {
  # penatly calculator
  # params:
  # dm: DesignMatrix object with attributes and constraints
  # X: design matrix
  # lambda: penalty for slacks
  # returns: penalty
  
  # calculate slacks for the current design
  dm$update_slacks()
  
  penalty <- lambda*( sum(abs(unlist(dm$dslacks))) + lambda*(sum(abs(unlist(dm$islacks)))) )
  return(penalty)
}

# delta_var <- function(D, i, j, lambda) {
#   # variance estimator for swapping rows
#   # params:
#   # D: current design matrix
#   # i: leaving row
#   # j: entering row
#   # lambda: penalty for slacks
#   # returns: variance estimator
#   
#   est <- var_est(D,j,j) - ( var_est(D,j,j)*var_est(D,i,i)-var_est(D,i,j)^2 ) - var_est(D,i,i)
#   return(est)
# }

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
}

update_chol <- function(L, x) {
  # updates the Cholesky with given row addition
  n <- length(x)
  for (k in 1:n) {
    r <- sqrt(L[k, k]^2 + x[k]^2)
    c <- r / L[k, k]
    s <- x[k] / L[k, k]
    L[k, k] <- r
    if (k < n) {
      L[(k+1):n, k] <- (L[(k+1):n, k] + s * x[(k+1):n]) / c;
      x[(k+1):n] <- c * x[(k+1):n] - s * L[(k+1):n, k];
    }
  }
  return(L)
}

downdate_chol <- function(L, x) {
  # "downdates"" the Cholesky with given row removal
  n <- length(x)
  for (k in 1:n) {
    r <- sqrt(L[k, k]^2 - x[k]^2)
    c <- r / L[k, k]
    s <- x[k] / L[k, k]
    L[k, k] <- r
    if (k < n) {
      L[(k+1):n, k] <- (L[(k+1):n, k] - s * x[(k+1):n]) / c;
      x[(k+1):n] <- c * x[(k+1):n] - s * L[(k+1):n, k];
    }
  }
  return(L)
}

det_chol <- function(L) {
  # returns the determinant of the cholesky
  return(prod(diag(L))^2)
}

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
  obj <- det_chol(dm$L) - penalty(dm, lambda)
  
  # iterate until the improvement in D-optimality is minimal or 100 iterations is reached
  while ((obj_delta_best > 10e-6) && (iter < 100)) {

    dm_best <- NULL
    obj_delta_best <- 0
    
    for (i in 1:n) {                      # iterate through rows in design matrix
      for (j in 1:nrow(candidate_set)) {  # iterate through rows in candidate set
        
        dm_test <- dm$copy(shallow=TRUE)
        dm_test$add_row(candidate_set[j,])
        dm_test$del_row(i)
        obj_test <- det_chol(dm_test$L) - penalty(dm_test, lambda)
        obj_delta <- obj_test - obj

        # calculate the potential improvement in D-optimality by replacing the
        # current row in the design with the current row in the candidate set
        #dvar <- delta_var(dm$X, dm$X[i,], candidate_set[j,])
        # delta <- delta_var(dm$X, t(dm$X[i,]), as.matrix(candidate_set[j,]))
        #obj_delta <- update_obj(dm$X, i, candidate_set[j,], lambda, det, dvar)

        # if that is greater than the best candidate so far, make it the new best pair
        if (obj_delta > obj_delta_best) {
          obj_delta_best <- obj_delta
          obj_best <- obj_test
          dm_best <- dm_test$copy(shallow=TRUE)
          print(paste(paste("iteration", iter, sep=" "), obj_delta_best, dvar_best, sep=" | "))
        } else {
          next
        }
        rm(dm_test)
      } # end for j
    } # end for i

    ### updates
    if (is.null(dm_best)) {
      # no better swaps found
      break
    } else {
      # perform swap
      dm_best$update_slacks()
      dm <- dm_best$copy(shallow=TRUE)
      
      # update the determinant following the swap
      obj <- obj_best
      
      iter <- iter+1
    }
    
  } # end while
  
  if (iter == 100) {
    print("Algorithm stopped due to reaching iteration limit")
  } else {
    print(paste("Convergence achieved in ",iter," iterations"))
  }
} # end fedorov

# asdf