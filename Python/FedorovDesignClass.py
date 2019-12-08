#%%
import numpy as np
import random

#%%
class DesignMatrix():
    """Object that contains a design matrix for optimization"""
    
    # attributes
    def __init__(self, n=0):
        # self.nrow = nrow
        # self.ncol = ncol
        self.n = n
        self.names = np.zeros(0)
        self.levels = np.zeros(0)
        self.dist = [] # list of np.array
        # self.values = np.zeros(0)
        self.interacts = [] # list of np.array
        self.X = np.zeros(0)
        self.dslacks = [] # list of np.array
        self.islacks = []

    # properties
    def shape(self):
        return self.X.shape

    # TODO: include optimality as class attributes & functions ???
        # TODO: Add lmda to constructor
        # TODO: Add X update function
        # TODO: Add optimality calculation to update funcitons; update self
        
    # def determinant(self):
    #     return np.linalg.det(self.X)
    
    # def det_XtX(self):
    #     return np.linalg.det(self.X.T @ self.X)

    # def doptimality(self, lmda=0):
    #     """ calculates doptimality of design (and optionally penalizes
    #     distribution constraints)

    #     params:
    #     dm: DesignMatrix object containing attribute & constraint information
    #     lmda: weight to penalize constraints.  lmda=0 means no distribution constraints

    #     returns: d-efficiency metric
    #     """

    #     self.update_slacks()
    #     nrow, ncol = self.shape()
        

    #     objective = (100 * det_XtX**(1/ncol)) / nrow

    #     penalty = lmda*sum([np.sum(np.abs(d)) for d in self.dslacks]) + lmda*sum([np.sum(np.abs(i)) for i in self.islacks])

    #     return objective + penalty

    # methods
    def set_size(self, n):
        self.n = n

    def add_attribute(self, name, levels, dist):
        """Adds and attribute & params to design
        Calculates initial attribute values & slacks"""
        dist = np.array(dist)
        if (np.sum(dist) !=1) & (np.sum(dist) !=100):
            raise SystemExit('Error: Distribution does not sum to 1 or 100') 
        else:
            if sum(dist) == 100: # convert to %
                dist = dist/100

        if levels != len(dist):
            raise SystemExit('Error: Number of levels & distributions do not match')
        else:
            self.names = np.append(self.names, name)
            self.levels = np.append(self.levels, levels).astype(int)
            self.dist.append(dist)

    def add_interaction(self, n1, n2, l1, l2, eq):
        """Adds information for interaction constraints (i.e., if male, cannot be pregnant)
        In example above: "gender","pregnant", 0, 1, F --> gender(0) != pregnant(1)
        Note that this is one direction: from n1 operating on n2.#
        If mutual, need to add reverse interaction constraint"""
        if isinstance(eq, bool):
            self.interacts.append(np.array([n1, n2, l1, l2, eq]))
        else:
            raise SystemExit('Error: "eq" must be "True" or "False"')

    def generate_design(self):
        """Generates a design matrix based on provided distributions"""
        ncol = len(self.names)
        nrow = self.n
        self.X = np.zeros((nrow,ncol))

        # find approx accurate distribution of values
        for j in range(len(self.names)):
            column = []
            reps = np.round(self.dist[j] * self.n) # calc # per level
            for i in range(self.levels[j]):
                fill = np.repeat(i, reps[i])
                column.extend(fill)

            # check column is appropriate length
            if len(column) < self.n:
                # if short, fill with most frequent level
                remaining = self.n - len(column)
                column.extend(np.repeat(np.argmax(reps), remaining))
            else:
                # if long, remove most frequent level
                while len(column) > self.n:
                    elements, counts = np.unique(column, return_counts=True)
                    drop_index = np.where(column == elements[np.argmax(counts)])[0][0]
                    column.pop(drop_index)

            # save values
            random.shuffle(column) #inplace
            self.X[:,j] = column

        # calculate slacks
        self.update_slacks()

    def add_row(self, row):
        """manually add row to design matrix"""
        if len(row) != self.X.shape[1]:
            # ensure row is proper length
            raise SystemExit('Error: Row length does not match matrix shape')
        elif all(row > self.levels-1):
            # ensure row[j] is valid for its column[,j]
            raise SystemExit('Error: Row values exceed defined levels')
        else:
            # add row
            self.X = np.append(self.X, [row], axis=0)
            # self.X.nrow, self.X.ncol = X.shape
            # recalculate values
            # self.update_values()
            # # recalculate slacks
            # self.update_dslacks()
            # self.update_islacks()

    def del_row(self, i):
        """delete row at index i"""
        if i < 0:
            raise SystemExit('Error: Row index is < 0.')
        elif i > self.X.shape[0]:
            raise SystemExit('Error: Row index > number of rows in design matrix.')
        else:
            # delete row
            self.X = np.delete(self.X, i, axis=0)
            # self.X.nrow, self.X.ncol = X.shape
            # # recalcualte values
            # self.update_values()
            # # recalculate slacks
            # self.update_dslacks()
            # self.update_islacks()

    # def update_values(self):
    #     self.values = np.apply(ncol(self.X), lambda i: self.X[,i])

    def update_dslacks(self):
        """update all slacks from each attribute's distribution constraints"""
        self.dslacks = []
        # iterate through attributes
        for j in range(self.X.shape[1]):
            vals = []
            # iterate through attribute levels
            for lvl in range(self.levels[j]):
                # count the number of rows with the current attrb level
                vals.append(sum(self.X[:,j]==lvl))
            
            # slacks are the difference between the intended count of profiles
            # and the actual count
            self.dslacks.append(self.dist[j]*self.n - vals)

    def update_islacks(self):
        """update all slacks from each attribute's interaction constraints"""
        self.islacks = 0
        if len(self.interacts) > 0:
            for i in range(len(self.interacts)):
                att1 = np.where(self.names == self.interacts[i][0])
                att2 = np.where(self.names == self.interacts[i][1])

                testA = self.X[:,att1] == int(self.interacts[i][2])
                testB = self.X[:,att2] == int(self.interacts[i][3])

                if self.interacts[i][5]:
                    # if A == B, count all where !=
                    self.islacks = self.islacks + (np.sum(testA) - np.sum(testA & testB))
                else:
                    # if A != B, count all where ==
                    self.islacks = self.islacks + (np.sum(testA & testB))
    
    def update_slacks(self):
        """update all slacks"""
        self.update_dslacks()
        self.update_islacks()
 


# %%
