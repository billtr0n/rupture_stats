/* function pointer to objective function for an empirical cost function.  
   in this case, no analytical solution is available for the cost function.
	inputs:
		float *params: parameters required for evaluating cost function
		float **extras: extra things needed to compute cost function.
		int n : length of the time-series
	returns:
		value of cost function 
*/
typedef float (*obj_fun_ptr)(float*, float**);


/* function: performs brute force optimization given arbitrary cost function 
	inputs:
		emp_obj_fun_ptr objective_function_ptr: function pointer
		int ndim: number of dimensions in optimization problem
		float* pl: array containing low values of parameters
		float* ph: array containing max values of parameters
		float* pd: array containing increment of parameters
		float** extras: time series for empirical comparisons
		float* full_output : pass NULL if not requested. allocation for full_output will happen inside fmin_brute
	returns:
		float* fmin_brute: array containing best fitting parameters.
*/
float* fmin_brute( obj_fun_ptr objective_function_ptr, int ndmin, float *pl, float *ph, float *pd, float** extras, float* full_output );