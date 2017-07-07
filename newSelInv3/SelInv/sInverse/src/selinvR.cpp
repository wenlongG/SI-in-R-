#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

#define nzvals(i) nzvals[(i)-1]
#define rowind(i) rowind[(i)-1]
#define colptr(i) colptr[(i)-1]
#define xsol(i)   xsol[(i)-1]
#define rhs(i)    rhs[(i)-1]
#define diag(i)   diag[(i)-1]
#define diag2(i)  diag2[(i)-1]

extern "C" void ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
extern "C" void ldlt_fact__(int *, int *, int *, double *,int*);
extern "C" void ldlt_solve__(int *, double *, double *);
extern "C" void ldlt_free__(int *);
extern "C" void ldlt_selinv__(int *, double *, int*,double*,int*,
                         int*,int*,int, double*, double*);

#ifdef TIMING
extern double getime(void);
#endif

// Input
// nnodes: nrows
// nnz: number of non-zero elemet
// colptr_, col pointer of the sparse matrix stucture
// rowind_, row index

// [[Rcpp::export]]
NumericVector selinv2r(int nnodes, int nnz,
                 IntegerVector colptr_, IntegerVector rowind_,
                 NumericVector nzval_)
{
    // convert R data to C pointer
    int *colptr, *rowind;
    colptr = (int*) malloc(sizeof(int)*colptr_.size());
    for(int i=0; i<colptr_.size(); i++)
        colptr[i] = colptr_[i];
    rowind = (int*) malloc(sizeof(int)*rowind_.size());
    for(int i=0; i<rowind_.size(); i++)
        rowind[i] = rowind_[i];
    double* nzvals;
    nzvals = (double*) malloc(sizeof(double)*nzval_.size());
    for(int i=0; i<nzval_.size(); i++)
        nzvals[i] = nzval_[i];


    int* nnzlplus;
    nnzlplus = (int*) malloc(sizeof(int));
    double* LDL_D;
    LDL_D = (double*) malloc(nnodes*sizeof(double));
    int* permout;
    permout = (int*) malloc(nnodes*sizeof(int));


    int i;
    int *perm;
    double  *diag2;
    int token, order=0;
    int Lnnz;

    token = 0;
    perm  = (int*)malloc((nnodes)*sizeof(int));
    for(i=0;i<nnodes;i++){
        perm[i]=i+1;
    }
    int nndess = nnodes;
    ldlt_preprocess__(&token, &nndess, colptr, rowind, &Lnnz, &order, perm);
    return nzval_;

    ldlt_fact__(&token, colptr, rowind, nzvals, nnzlplus);

    double *LNZXX, *LDL_L;
    int *RowXX,*ColXX;
    double *V;
    int dumpL=1;

    RowXX =(int *) malloc(Lnnz*sizeof(int));
    ColXX =(int *) malloc(Lnnz*sizeof(int));
    LDL_L =(double *) malloc(Lnnz*sizeof(double));
    LNZXX = (double *)malloc(*nnzlplus*sizeof(double));


    V=(double*)calloc(3*Lnnz+*nnzlplus,sizeof(double));
    /* selected inversion */
    diag2 = (double*)calloc(nnodes,sizeof(double));
    ldlt_selinv__(&token, diag2, &dumpL, LNZXX, RowXX, ColXX,
                  permout,*nnzlplus,LDL_L,LDL_D);

    int LnnzOutput;
    NumericVector diagOutput(nnodes);


    for(i=0;i<Lnnz;i++)
    {
        V[i] = (double) RowXX[i];
        V[i+Lnnz] = (double) ColXX[i];
        V[i+2*Lnnz] = LDL_L[i];
    }
    for (i=0;i<*nnzlplus;i++){
        V[i+3*Lnnz] = LNZXX[i];
    }


    ldlt_free__(&token);

    if (order == 0) free(perm);
    free(diag2);

    for(i=0;i<nnodes;i++)
        diagOutput[i] = diag2[i];
    LnnzOutput = Lnnz;


    return diagOutput;
}
