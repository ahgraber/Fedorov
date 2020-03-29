# create attribute object
DesignMatrix <- setRefClass("DesignMatrix", 
  fields = list(n = "numeric", names = "character", levels = "numeric", dist = "list", 
                values = "list", interacts = "list", X = "matrix", L = "matrix",
                dslacks = "list", islacks = "numeric", cholesky = "logical")
) # end attributeClass class definition

DesignMatrix$methods(
  # initialize = function(n) {
  #   .self$n <- n
  #   .self$islacks <- list()
  # }, # end initialize
  
  set_cholesky = function(cholesky) {
    .self$cholesky = cholesky  
  }, 
  
  add_attribute = function(name, levels, dist) {
    # adds an attribute & params to design
    # calculates initial attribute values and slacks

    if ( (sum(dist) != 1 ) && (sum(dist) != 100) )  {
      stop("Error: distribution does not sum to 1 or 100")
      # return(NULL)
    } else { 
      if (sum(dist) == 100) { dist <- dist/100 }
    }
    
    if (levels != length(dist)) {
      stop("Error: number levels & distribution list do not match")
      # .self$values <- NULL  
    } else {
      # save params to object
      .self$names <- c(.self$names, name)
      .self$levels <- c(.self$levels, levels)
      .self$dist <- c(.self$dist, list(dist))

    }
  }, # end add_attribute
  
  add_interaction = function(n1, n2, l1, l2, eq) {
    # adds information for interaction constraints (i.e., if male, cannot be pregnant)
    # in example above: "gender","pregnant", 0, 1, F --> gender(0) != pregnant(1)
    # note that this is one direction: from n1 operating on n2.#
    # if mutual, need to add reverse interaction constraint
    if (is.logical(eq)) {
      .self$interacts <- append(.self$interacts, list(list(n1, n2, l1, l2, eq)))
      .self$update_islacks()
    } else {
      stop("Error: 'eq' must be T or F")
    }  
  }, # end add_interaction

  generate_design = function() {
    # generates matrix based on distributions
    
    .self$X <- matrix()

    if (.self$cholesky) { .self$L <- matrix() }
    
    # create column for each attribute
    for (j in 1:length(.self$names)) {
      # find approximately accurate distribution of values
      column <- c()                 # create empty column
      reps <- round(.self$dist[[j]]*.self$n)         # calculate number per level
      for (i in 1:.self$levels[[j]]) {         # fill column
        fill <- rep(i, reps[i])
        column <- append(column, fill)
      } # end for i
      
      # check column is appropriate length
      if (length(column) < .self$n) { 
        # if short, add most frequent level fill times for correct length
        remaining <- .self$n - length(column)
        column <- append(column, rep(which.max(reps), remaining))
      } else { 
        # if long, drop most frequent level until correct length
        while(length(column) > .self$n){
          drop_index <- which.max(table(column))
          column <- column[-drop_index]
        }
      } 
      
      # save values to object
      vals <- sample(column-1) # convert to 0-based, randomize order
      .self$values[[j]] <- vals
      
      # calculate slacks
      .self$dslacks <- list()
      vals <- c()
      for (k in 1:.self$levels[j]-1) {
        vals <- append(vals, length(which(.self$values[[j]]==k)))
      }
      .self$dslacks[[j]] <- .self$dist[[j]]*.self$n - vals
      
    } # end for j
    
    # save design
    .self$X <- matrix(unlist(.self$values), nrow=.self$n, ncol=length(.self$names))
    
    if (.self$cholesky) {
      # save cholesky
      A <- t(.self$X)%*%.self$X
      .self$L <- try(t(chol(A)))      
    }


  }, # end generate
  
  add_row = function(row){
    if (length(row) != ncol(.self$X)) {
      # ensure row is appropriate length
      stop("Row is wrong dimension")
    # } else if ( all(row > .self$levels-1) ) {
    } else if ( any(row > .self$levels-1) ) {
      # ensure row[j] is valid value for its column[,j]
      stop("Row values exceed allowed levels")        
    } else {
      #add row
      .self$X <- rbind(.self$X, row)
      row.names(.self$X) <- NULL
      
      # recalculate values
      if (.self$cholesky) { .self$update_chol(row) }
      .self$update_values()
      .self$update_slacks()

    }
  }, # end add_row
  
  del_row = function(i) {
    if (i <= 0) {
      stop("Row index is < 0")
    } else if (i > nrow(.self$X)) {
      stop("Row index > number rows in design matrix")
    } else {
      row <- .self$X[i,]
      
      # delete row
      .self$X <- .self$X[-i,] 
     
      # recalculate values
      if (.self$cholesky) { .self$downdate_chol(row) }
      .self$update_values()
      .self$update_slacks()

    }
  }, # end del_row
  
  update_values = function() {
    # reset .self$values from matrix X
    .self$values <- lapply(seq_len(ncol(.self$X)), function(i) {.self$X[,i]})
  }, # end update_values
  
  update_chol = function(row) {
    # updates the Cholesky with given row addition
    n <- length(row)
    for (k in 1:n) {
      r <- sqrt(.self$L[k, k]^2 + row[k]^2)
      c <- r / .self$L[k, k]
      s <- row[k] / .self$L[k, k]
      .self$L[k, k] <- r
      if (k < n) {
        .self$L[(k+1):n, k] <- (.self$L[(k+1):n, k] + s * row[(k+1):n]) / c
        row[(k+1):n] <- c * row[(k+1):n] - s * .self$L[(k+1):n, k]
      }
    }
    # return(.self$L)
  }, # end update_chol

  downdate_chol = function(row) {
    # "downdates"" the Cholesky with given row removal
    n <- length(row)
    for (k in 1:n) {
      r <- sqrt(.self$L[k, k]^2 - row[k]^2)
      c <- r / .self$L[k, k]
      s <- row[k] /.self$ L[k, k]
      .self$L[k, k] <- r
      if (k < n) {
        .self$L[(k+1):n, k] <- (.self$L[(k+1):n, k] - s * row[(k+1):n]) / c
        row[(k+1):n] <- c * row[(k+1):n] - s * .self$L[(k+1):n, k]
      }
    }
    # return(.self$L)
  }, # end downdate_chol

  # det_chol = function() {
  #   # returns the determinant of the cholesky
  #   return(prod(diag(.self$L))^2)
  # }, # end det_chol

  update_dslacks = function() {
    # updates all slacks from each attribute's distribution constraints
    .self$dslacks <- lapply(seq_len(ncol(dm$X)), function(i) {
      prop.table(
        apply(dm$X, MARGIN=2, FUN=table)[[i]]
      )
    })
  }, # end update_dslacks
  
  update_islacks = function() {
    # updates all slacks from interaction constraints
    .self$islacks <- 0
    
    if (length(interacts) > 0) {
      for (i in 1:length(interacts)) {
        indx1 <- which(.self$names == .self$interacts[[i]][[1]])
        indx2 <- which(.self$names == .self$interacts[[i]][[2]])

        testA <- .self$X[,indx1] == .self$interacts[[i]][[3]]
        testB <- .self$X[,indx2] == .self$interacts[[i]][[4]]

        if (.self$interacts[[i]][[5]]) {
          # if A must == B, count all where !=
          .self$islacks <- .self$islacks + (sum(testA) - sum(testA & testB))
        } else {
          # if A must != B, count all where ==
          .self$islacks <- .self$islacks + (sum(testA & testB))
        }
      } # end for 
    } # end if
  },  # end update_islacks
  
  update_slacks = function () {
    # wrapper for dslacks and islacks
    .self$update_dslacks()
    .self$update_islacks()
  }
) # end methods

#-- Objective functions -----------------------------------------------------------
penalty <- function(dm, lambda) {
  # penatly calculator
  # params:
  # dm: DesignMatrix object with attributes and constraints
  # X: design matrix
  # lambda: penalty for slacks 
  # returns: penalty
  
  # calculate slacks for the current design
  dm$update_slacks()
  
  penalty <- lambda*( sum(abs(unlist(dm$dslacks))) + 2*lambda*(sum(abs(unlist(dm$islacks)))) )
  return(penalty)
}

doptimality <- function(dm, lambda=0, how='chol') {
  # calculates doptimality of DesignMatrix object 
  # (and optionally penalizes distribution constraints)
  # params:
  # dm: DesignMatrix object containing attribute & constraint information
  # lambda: weight to penalize constraints.  lambda=0 means no distribution constraints
  # how: using standard determinant 'det' or cholesky update matrix 'chol' for calculation
  # returns: d-efficiency metric

  if (how == 'det') {
    obj <- objective(dm)
  } else if (how == 'chol') {
    obj <- objective_chol(dm)
  } else {
    stop('Error: "how" not in c("det","chol")')
  }

  pen <- penalty(dm, lambda)
  # this double-penalizes islacks b/c we really don't want impossible interactions

  return(obj - pen)
}

objective <- function(dm) {
  obj <- (100 * det( t(dm$X)%*%dm$X )^(1/ncol(dm$X)))/ nrow(dm$X)
  # obj <- det( t(dm$X)%*%dm$X ) / nrow(dm$X)
  return(obj)
}

det_chol <- function(L) {
  # returns the determinant of the cholesky
  return( prod(diag(L))^2 )
}

objective_chol <- function(dm) {
  obj <- (100 * det_chol(dm$L)^(1/ncol(dm$X))) / nrow(dm$X)
  # obj <- (100 * det_chol(dm$L) )
  return(obj)
}

sumfisherz <- function(dm, lambda=0) {
  # calculates the sum of the fisher z score of the absolute values of the correlation matrix
  # minimization objective function
  # params
  # dm: DesignMatrix object containing attribute & constraint information
  # lambda: weight to penalize constraints.  lambda=0 means no distribution constraints
  # design: design matrix where columns are attributes and rows are patients
  # returns correlation score
  
  # calculate slacks for the design
  dm$update_slacks()

  r <- abs(cor(dm$X))
  z <- .5*(log(1+r)/(1-r))
  obj <- sum(z[is.finite(z)])
  pen <- penalty(dm, lambda)
  return(obj + pen)
}
