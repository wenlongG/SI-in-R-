library(SparseM)

    B = diag(1,8)
    B[2,1] = 1
    B[1,2] = 1
    B[3,2] = 1
    B[2,3] = 1
    B[5,1]=1
    B[6,2]=1
    B[1,5]=1
    B[2,6]=1
    testM=B*4
    testM = as.matrix.csr(testM)
    colptr = testM@ja
    rowind = testM@ia
    nzvals = testM@ra
    nnz = length(nzvals)
    nnodes = nrow(testM)

selinv2julia(nnodes, nnz, colptr, rowind, nzvals)
