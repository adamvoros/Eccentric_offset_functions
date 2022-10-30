/*****************************************************************************/
/* lmfit.h																     */
/*****************************************************************************/

#ifndef	__LMFIT_H_INCLUDED
#define	__LMFIT_H_INCLUDED	1

/* matrix_alloc():
   Allocates a 2 dimensional array of 'double' types with a size of 'd'*'d'. */
double **matrix_alloc(int d);
/* matrix_alloc_gen():
   Allocates a 2 dimensional array of 'double' types with a size of 'x'*'y'. */
double **matrix_alloc_gen(int x,int y);
/* matrix_free():
   Releases the array allocated by matrix_alloc() or matrix_alloc_gen().     */
void	matrix_free(double **m);

/* vector_alloc():
   Allocates a 1 dimensional array of 'double' types with a size of 'd'.     */
double	*vector_alloc(int d);
/* vector_free():
   Releases the array allocated by vector_alloc().			     */
void	vector_free(double *v);

/* solve_gauss():
   Solves the linear equation 'a'*x='b', where 'a' should be an 'n'*'n' matrix
   (maybe allocated by matrix_alloc()), and 'b' is a vector with the size of 
   'n'. The result is stored in 'b' and the matrix 'a' is destroyed.	     */
int	solve_gauss(double **a,double *b,int n);

/* invert_gauss():
   Inverts the matrix 'a' (where 'a' is an 'n' times 'n' matrix, maybe
   alloceted by matrix_alloc()). The result is stored in 'a' itself.	     */
int	invert_gauss(double **a,int n);

/* lin_fit(), nlm_fit_base():
   These two routines fit the function 'funct' to the values 'y'. 
   In the linear fit routine, lin_fit(), the function 'funct' should be
   linear in its parameters 'a' and should have 'n' independent parameters.
   The fitted parameters are stored in the array 'a', and no previous
   assumption for these parameters is needed (so, the array a[] might be
   uninitialized after a dynamical allocation). The array 'y' should 
   contain 'ndata' real numbers defined on the arbitrary kind of points x[]. 
   The optional argument 'weight' might contain the weights (also 'ndata' real 
   numbers) for each data point. If 'weight' is NULL, all points do have the 
   same weight during the fit.
   The function 'funct' must have the following prototype: 

      void funct(void *parameter,double *a,double *y,double *dyda,void *param).

 */


int	lin_fit(	void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void *),
			int n,int ndata,void *param,
			double *err);

double	nlm_fit_base(	void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void*),
			int n,int ndata,void *param,
			double lam,double lam_mpy);

typedef struct
 {	double	**cmatrix;	/* constraint matrix, [0..nc-1][0..n-1]	*/
	double	*cvector;	/* constraint vector, [0..nc-1]		*/
	int	nc;		/* number of linear constraints		*/
 } constraint;

/* lin_fit_con(), nlm_fit_base_con():
   The same as lin_fit() and nlm_fit_base(), but optional linear constraints 
   also can be taken into account, specified by the constraint data 'cc'. 
   If 'cc' is NULL, the call falls back to lin_fit().			     */

int	lin_fit_con(	void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void *),
			int n,int ndata,void *param,constraint *cc,
			double *err);

double	nlm_fit_base_con(void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void*),
			int n,int ndata,void *param,constraint *cc,
			double lam,double lam_mpy);

/* nlm_fit_nmdf_con(): 
   The same as nlm_fit_base_con() but the partial derivatives of the function 
   'funct' are calculated numerically, assuming a symmetric finite difference 
   defined by the values of 'da'. 	    				     */

double	nlm_fit_nmdf_con(void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void*),
			int n,int ndata,void *param,constraint *cc,
			double lam,double lam_mpy,double *da);

/* nm_fit_nmdf():
   Same as above but omitting the correlations 'cc' (equivalent to the call
   of nlm_fit_nmdf_con() with setting cc=NULL).				     */

double	nlm_fit_nmdf(	void **x,double *y,double *a,double *weight,
			void (*funct)(void *,double *,double *,double *,void*),
			int n,int ndata,void *param,
			double lam,double lam_mpy,double *da);


#endif
                                                                         
