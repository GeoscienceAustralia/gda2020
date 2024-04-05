#include <stdio.h>

  /*
   *-----------------------------------------------------------------------------
   *	internal matrix structure
   *-----------------------------------------------------------------------------
   */
  typedef struct {
    int row;
    int col;
  } MATHEAD;

  typedef struct {
    MATHEAD head;
    /*
     * only the starting address of the following will be
     * returned to the C programmer, like malloc() concept
     */
    double *matrix;
  } MATBODY;

  typedef double **MATRIX;

#define	Mathead(a)	((MATHEAD *)((MATHEAD *)(a) - 1))
#define MatRow(a)	(Mathead(a)->row)
#define	MatCol(a)	(Mathead(a)->col)

  /*
   *----------------------------------------------------------------------------
   *	mat_errors definitions
   *----------------------------------------------------------------------------
   */
#define	MAT_MALLOC	1
#define MAT_FNOTOPEN	2
#define	MAT_FNOTGETMAT	3

  /*
   *----------------------------------------------------------------------------
   *	matrice types
   *----------------------------------------------------------------------------
   */
#define UNDEFINED	-1
#define ZERO_MATRIX	0
#define	UNIT_MATRIX	1


  MATRIX mat_error(int errno);

  MATRIX _mat_creat(int row, int col);

  MATRIX mat_creat(int row, int col, int type);

  MATRIX mat_fill(MATRIX A, int type);

  int mat_free(MATRIX A);

  MATRIX mat_copy(MATRIX A);

  MATRIX mat_colcopy1(MATRIX A, MATRIX B, int cola, int colb);

  int fgetmat(MATRIX A, FILE *fp);

  MATRIX mat_dumpf(MATRIX A, const char *s);

  MATRIX mat_dump(MATRIX A);

  MATRIX mat_fdump(MATRIX A, FILE *fp);

  MATRIX mat_fdumpf(MATRIX A, const char *s, FILE *fp);



  MATRIX mat_add(MATRIX A, MATRIX B);

  MATRIX mat_sub(MATRIX A, MATRIX B);

  MATRIX mat_mul(MATRIX A, MATRIX B);

  double mat_diagmul(MATRIX A);

  MATRIX mat_tran(MATRIX A);

  MATRIX mat_inv(MATRIX a);

  MATRIX mat_SymToeplz(MATRIX R);



  int mat_lu(MATRIX A, MATRIX P);

  MATRIX mat_backsubs1(MATRIX A, MATRIX B, MATRIX X, MATRIX P, int xcol);

  MATRIX mat_lsolve(MATRIX a, MATRIX b);



  MATRIX mat_submat(MATRIX A, int i, int j);

  double mat_cofact(MATRIX A, int i, int j);

  double mat_det(MATRIX a);

  double mat_minor(MATRIX A, int i, int j);



  MATRIX mat_durbin(MATRIX R);

  MATRIX mat_lsolve_durbin(MATRIX A, MATRIX B);


