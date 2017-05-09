/* function pointer to objective function for an empirical cost function.  
   in this case, no analytical solution is available for the cost function.
	inputs:
		float *params: parameters required for evaluating cost function
		void **extras: extra things needed to compute cost function.
	returns:
		value of cost function 
*/
typedef float (*obj_fun_ptr)(float*, void**);


/* function: performs brute force optimization given arbitrary cost function 
	inputs:
		(obj_fun_ptr) objective_function_ptr: function pointer
		(int) ndim: number of dimensions in optimization problem
		(float*) pl: array containing low values of parameters
		(float*) ph: array containing max values of parameters
		(float*) pd: array containing increment of parameters
		(void**) extras: time series for empirical comparisons
		(int*) ne : length of each array in extras
		(float*) full_output : pass NULL if not requested. allocation for full_output will happen inside fmin_brute
	returns:
		float* fmin_brute: array containing best fitting parameters.
*/
float* fmin_brute( obj_fun_ptr objective_function_ptr, 
						float **p, int ndim, int np, void** extras, float* full_output );
