## R CMD build sInverse
## R CMD INSTALL sInverse_0.1.0.tar.gz

library(sInverse)
library(SparseM)
## create sample matrix to test
B = diag(8,8)
B[2,1] = 1
B[1,2] = 1
B[3,2] = 1
B[2,3] = 1
B[5,1] = 1
B[6,2] = 1
B[1,5] = 1
B[2,6] = 1
testM = B
testM = as.matrix.csr(testM)
colptr = testM@ia
rowind = testM@ja
nzvals = testM@ra
nnz = length(nzvals)
nnodes = nrow(testM)

## run select inserve
selinv2r(nnodes, nnz, colptr, rowind, nzvals)
