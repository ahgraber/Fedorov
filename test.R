# create test design
dm <- design(n=8)

# test add_attribute
dm$add_attribute(name="age", levels=3, dist=c(25,50,25))
dm$add_attribute(name="sex", levels=2, dist=c(50,50))
dm$add_attribute(name="bmi", levels=3, dist=c(33,33,34))
dm$generate()

dm$values
dm$dslacks
dm$islacks
dm$X

# test add_row
dm$add_row(c(3,3,3)) # note: throws an error because 3 is not in any of the levels

dm$values
dm$dslacks
dm$X

# try again with valid vals
dm$add_row(c(1,1,1)) 

dm$values
dm$dslacks
dm$X

# test del_row
dm$del_row(1)
dm$update_dslacks()

dm$values
dm$dslacks
dm$X
doptimality(dm, dm$X)


