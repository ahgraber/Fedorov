# create attribute object
design <- setRefClass("DesignMatrix", 
  fields = list(n = "numeric", names = "character", levels = "numeric", dist = "list", 
                values = "list", interacts = "list", X = "matrix", 
                dslacks = "list", islacks = "numeric")
) # end attributeClass class definition

design$methods(
  # initialize = function(n) {
  #   .self$n <- n
  #   .self$islacks <- list()
  # }, # end initialize
  
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

  generate = function() {
    # generates matrix based on distributions
    
    .self$X <- matrix()
    
    # create column for each attribute
    for (j in 1:length(.self$names)) {
      # find approximately accurate distribution of values
      column <- c()                 # create empty column
      reps <- round(.self$dist[[j]]*.self$n)         # calculate number per level
      for (i in 1:.self$levels[[j]]) {         # fill column
        x <- rep(i, reps[i])
        column <- append(column, x)
      } # end for i
      
      # check column is appropriate length
      if (length(column) < .self$n) { 
        # if short, add most frequent level fill times for correct length
        fill <- .self$n - length(column)
        column <- append(column, rep(which.max(reps), fill))
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
    
  }, # end generate
  
  add_row = function(row){
    if (length(row) != ncol(.self$X)) {
      # ensure row is appropriate length
      stop("Row is wrong dimension")
    } else if ( all(row > X$levels-1) ) {
      # ensure row[j] is valid value for its column[,j]
      stop("Row values exceed allowed levels")        
    } else {
      #add row
      .self$X <- rbind(.self$X, row)
      row.names(.self$X) <- NULL
      # recalculate values
      .self$update_values()
      # # recalculate slacks
      # .self$update_dslacks()
      # .self$update_islacks()
    }
  }, # end add_row
  
  del_row = function(j) {
    if (j <= 0) {
      stop("Row index is < 0")
    } else if (j > nrow(.self$X)) {
      stop("Row index > number rows in design matrix")
    } else {
      # delete row
      .self$X <- .self$X[-j,] 
      # recalculate values
      .self$update_values()
      # # recalculate slacks
      # .self$update_dslacks()
      # .self$update_islacks()
    }
  }, # end del_row
  
  update_values = function() {
    .self$values <- lapply(seq_len(ncol(.self$X)), function(i) {.self$X[,i]})
  }, # end update_values
  
  update_dslacks = function() {
    # updates all slacks from each attribute's distribution constraints
    .self$dslacks <- list()
    for (i in 1:length(.self$names)) {
      vals <- c()
      for (j in 1:.self$levels[i]-1) {
        vals <- append(vals, length(which(.self$values[[i]]==j)))
      }
      .self$dslacks[[i]] <- .self$dist[[i]]*.self$n - vals
    }
    # .self$dslacks <- mapply(function(d,v) {(d*.self$n) - (as.numeric(table(v)))},
    #                         d=.self$dist, v=.self$values)
  }, # end update_dslacks
  
  update_islacks = function() {
    # updates all slacks from interaction constraints
    .self$islacks <- 0
    
    if (length(interacts) > 0) {
      for (i in 1:length(interacts)) {
        indx1 <- which(.self$names == .self$interacts[[i]][[1]])
        indx2 <- which(.self$names == .self$interacts[[i]][[2]])
        if (.self$interacts[[i]][[5]]) {
          # if A must == B, count all where !=
          testA <- .self$X[,indx1] == .self$interacts[[i]][[3]]
          testB <- .self$X[,indx2] == .self$interacts[[i]][[4]]
          .self$islacks <- .self$islacks + (sum(testA) - sum(testA & testB))
        } else {
          # if A must != B, count all where ==
          testA <- .self$X[,indx1] == .self$interacts[[i]][[3]]
          testB <- .self$X[,indx2] == .self$interacts[[i]][[4]]
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
