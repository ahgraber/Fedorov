import numpy as np

#%% objective functions(s)
# TODO: add to DesignMatrix class ???
def doptimality (dm, X, lmda=0):
    """ calculates doptimality of design (and optionally penalizes
    distribution constraints)

    params:
    dm: DesignMatrix object containing attribute & constraint information
    X: current design matrix for which to calculate doptimality
    lmda: weight to penalize constraints.  lmda=0 means no distribution constraints

    returns: d-efficiency metric
    """

    nrow, ncol = X.shape
    dm.X = X
    dm.update_slacks()
    

    obj = (100 * np.linalg.det(X.T @ X)**(1/ncol)) / nrow

    if lmda = 0:
        pnlty = 0
    else:
        pnlty = penalty(dm, X, lmda)

    return obj + pnlty



#%% Helper functions
def penalty(dm, X ,lmda=0):
    dm.update_slacks()

    return lmda*sum([np.sum(np.abs(d)) for d in dm.dslacks]) + lmda*sum([np.sum(np.abs(i)) for i in dm.islacks])
    
def var_est(X, i_row, j_row):
    """variance estimator for swapping rows
    params:
    X: current design matrix
    i_row: leaving row
    j_row: entering row
    """

    # attempting protection from singular matrix inversions
    try:
        Xinv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        Xinv = np.linalg.pinv(X.T @ X)

    est = j_row @ Xinv @ i_row
    return est

def delta_var(X, i, j, lmda):
    """variance estimator for swapping rows

    params:
    X: current design matrix
    i: leaving row
    j: entering row
    lmda: penalty for slacks

    returns: variance estimator
    """

    d_var = var_est(D,j,i) - ( var_est(D,j,j)*var_est(D,i,i)-var_est(D,i,j)**2 ) - var_est(D,i,i)
    return d_var

def update_obj (X, i, j_row, lmda, det, d_var):
    """# calculates the % increase in the objective function for the row swap
    
    # params:
    # D: current design matrix
    # i: leaving row *index*
    # j: entering row
    # lambda: penalty for slacks
    # det: determinant
    # dvar: delta variance estimator

    # returns: updated doptimality metric
    """

    # calculate penalties
    old_p <- penalty(dm, X, lambda)
    i_row <- X[i,]
    dm.del_row(i)
    dm.add_row(j_row)
    dm.update_slacks()
    new_p <- penalty(dm, dm.X, lambda)

    # revert dm.X
    dm.X <- X
    return( (det*(1+dvar)-new_p) / (det-old_p) - 1 )


#%% Fedorov
def fedorov (dm, candidate_set, n, lmda=0):
    """# fedorov algorithm find d-optimal design
    
    # params:
    # dm: Design Matrix object
    # candidate_set: matrix of all possible combinations of attributes
    # n: number of rows
    # lambda: weight for slack penalties
    
    # returns: optimal design
    """

    # set initial values
    n_iter = 1
    obj_delta_best = 0.0001
    det = np.linalg.det(X.T @ X) #abs()?

    # iterate until the improvement in D-optimality is minimal or 100 iterations is reached
    while ((obj_delta_best > 10e-6) && (n_iter < 100)):
        i_best = np.nan
        j_best = np.nan
        obj_delta_best = 0

        for i in range(n):                         # iterate through rows in design
            for j in range(len(candidate_set)):    # iterate through rows in candidate set
                dvar = np.nan
                obj_delta = np.nan

                # calculate the potential improvement in D-optimality by replacing the
                # current row in the design with the current row in the candidate set
                dvar = delta_var(dm.X, dm.X[i,], candidate_set[j,])
                obj_delta = update_obj(dm.X, i, candidate_set[j,], lmda, det, dvar)

                # if that is greater than the best candidate so far, make it the new best pair
                if obj_delta > obj_delta_best:
                    obj_delta_best = obj_delta
                    dvar_best = dvar
                    i_best = i
                    j_best = j
                    print(f'Iteration {n_iter}: {obj_delta_best} | {dvar_best}')
                else:
                    continue 

        # updates
        if np.isnan(i_best):
            # no better swaps found
            break
        else:
            dm.del_row(i_best)
            dm.add_row(candidate_set[j_best,])
            dm.update_slacks()

            # update the determinate following the swap
            det = det*(1+dvar_best)
            
            n_iter += 1

    # end while

    if n_iter == 100:
        print('Algorithm stopped due to reaching iteration limit of 100')
    else:
        print(f'Convergence achieved in {n_iter} iterations')
    
    return dm.X