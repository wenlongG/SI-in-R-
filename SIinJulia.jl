#A function to get selinv of matrix corresponding to certain sparsity spB
function SIplusD(spA,spB)
  spM=spA+1.0e-30*spB
  return SelInvD(spM)
end


# Selected inversion of sparse symmetrix matrix corresponding to nonzero entries.
# Output the selected inverse matrix and logdet of the input matrix
function SelInvD(testM_sparse)
nnodes=size(testM_sparse,1)
colptr = testM_sparse.colptr;
rowind = testM_sparse.rowval;
nzvals = testM_sparse.nzval;
nnz = length(nzvals);

nnodes = convert(Int32,nnodes);
nnz = convert(Int32,nnz);
colptr = convert(Array{Int32,1},colptr);
rowind = convert(Array{Int32,1},rowind);

Lnnz = [0]
Lnnz = convert(Array{Int32,1},Lnnz);
nnzlplus = [0]
nnzlplus = convert(Array{Int32,1},nnzlplus);
permout=zeros(Int32,nnodes);
#permout = convert(Array{Int32,1},permout);
DIAG=zeros(nnodes);
LDL_D = zeros(nnodes);

V = ccall((:selinv2julia,
        "/home/grad/wgong/Spatial/newSelInv3/SelInvCopy/EXAMPLES/julia.so"), Ptr{Cdouble}, (Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble},
	Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble}),
       nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout,nnzlplus,DIAG,LDL_D);
v2 = pointer_to_array(V, (3*Lnnz[1]+nnzlplus[1]),true)# for julia-0.4.5
#Libc.free(V)
#v2=unsafe_wrap(Array, V, (3*Lnnz[1]+nnzlplus[1]),true);# was pointer_to_array, now deprecated in julia-0.5
#  This unsafe_wrap function is labelled "unsafe" because it will crash if `pointer` is not a valid memory address to data of the requested length.
nnzInv = Lnnz[1];
iInv = v2[1:nnzInv];
jInv = v2[(nnzInv+1):(2*nnzInv)];
iInv =  convert(Array{Int32,1},iInv);
jInv = convert(Array{Int32,1},jInv);
LDL_L =  v2[(2*nnzInv+1):(3*nnzInv)];
invElement =  v2[(3*nnzInv+1):(3*nnzInv+nnzlplus[1])];

inv_pos = sparse(iInv,jInv,1);

k=Int64[1]
newLNZ = zeros(nnzInv);
for i=1:nnzInv
    if i == 1
       k = convert(Int64,DIAG[jInv[i]]);
    elseif jInv[i]>jInv[i-1]
       k = convert(Int64,DIAG[jInv[i]]);
    end
    newLNZ[i] = invElement[k];
    k = k + 1;
end



for i=1:nnzInv
    iInv[i] = permout[iInv[i]];
    jInv[i] = permout[jInv[i]];
end


## get the LDL component
#LDL_L_matrix = sparse(iInv,jInv,LDL_L);
#for i=1:nnodes
#    LDL_L_matrix[i,i] = 1;
#end
#LDL_D_matrix = sparse([1:1:nnodes],[1:1:nnodes],LDL_D);

##log determinent
logDet=sum(log(LDL_D))
#logDet =sum(log(diag((LDL_D_matrix))))

## reconstruct testM
#testM_reconstruct = LDL_L_matrix*LDL_D_matrix*transpose(LDL_L_matrix);#LDLT
#testM_reconstruct-testM_sparse

testM_inv = sparse(iInv,jInv,newLNZ);
testM_inv = testM_inv + transpose(testM_inv);
for i=1:nnodes
    testM_inv[i,i] = testM_inv[i,i]/2
end



## cholesky decomposition from LDLT
#cholLDL=Any
#cholLDL=(LDL_D_matrix.^0.5)*LDL_L_matrix'

  return(testM_inv,
          logDet
         )

end
