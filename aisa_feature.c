/*prof.dr.selamibeyhan
 *
 *aisa_feature.c - This extension is used in order to implement the feature extraction algorithm for AISA. 
 *This extension eventually returns an array with best individuals from the current population. 
 *
 *Most of variables are used with memory allocation in the end they make free the allocated memories.
 *  
*/
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
static  PyObject*  aisa_feature(PyObject *self, PyObject *args) { // beginning of main function
	PyObject* _populations;
	PyObject* _output;
	PyObject* item;
	PyObject* res1;
	PyObject* result;
	int N,D,i = 0;
	double* frame = NULL;
	double* yd0 = NULL;   

	/*pass along the arguments and convert them to C/pyObject objects*/
	if (!PyArg_ParseTuple(args, "iiOO", &N, &D, &_populations, &_output)) {
		PyErr_SetString(PyExc_RuntimeError, "arguments were not received properly from python");
	}

	if (!PyList_Check(_populations) || !PyList_Check(_output)) {
		PyErr_SetString(PyExc_RuntimeError, "observations and/or initial yd are not of type list - terminate");
	}
	/*allocate memory and convert observations and initial clusters array from python-types*/
	yd0  = calloc(N, sizeof(double));
	if (yd0 == NULL) PyErr_SetString(PyExc_MemoryError, "Memory Allocation Error");
	
	frame = calloc(N*D, sizeof(double));
	if (frame == NULL) {
		free(yd0);
		PyErr_SetString(PyExc_MemoryError, "Memory Allocation Error");
	}

	for (i = 0; i < N; i++) {
		item = PyList_GetItem(_output, i); /*extract values from the centroid array*/
		if (!PyFloat_Check(item)) continue;
		yd0[i] = PyFloat_AsDouble(item); /*w is the initial observation assigned to current i'th cluster*/
	}

	for (i = 0; i < N*D; i++) {
		item = PyList_GetItem(_populations, i); /*convert value at i'th index in the _obs python list to C-type*/
		if (!PyFloat_Check(item)) continue; /*assert that observation coordinate if of type float*/
		frame[i] = PyFloat_AsDouble(item);
	}
        
        /* MAIN FEATURE EXTRACTION ALGORITHM */
        /* This function can be written in another .c file in case of more functions needed for the module. */
        
	/* Variable Definitions */
	int a, b, dimbase, r;
	float lamda = 0.99, delta = 1e-3;
	a = N; b = D; dimbase = 3*D;
	float min_val, max_val, yd[a], x[a][b], best_params[D], res_prod; 
	//float xbase[a], yhat[a], u[dimbase], Pu[dimbase], Ka[dimbase], Par[dimbase];
        //float xn[a][b], Par2[D][3], P[dimbase][dimbase], Pd[dimbase][dimbase], Phi[a][dimbase], Psi[a][3], yhat_new[b][a];	
	
	/* Vector memory allocation */
	float *xbase = NULL, *yhat = NULL, *u = NULL, *Pu = NULL, *Ka = NULL,*Par = NULL;      
        xbase = calloc(sizeof(float),a);
        yhat  = calloc(sizeof(float),a);
     	u     = calloc(sizeof(float),dimbase); 
     	Pu    = calloc(sizeof(float),dimbase);    
     	Ka    = calloc(sizeof(float),dimbase);
     	Par   = calloc(sizeof(float),dimbase);
     	
     	float **xn, **Par2, **P, **Pd, **Phi, **Psi, **yhat_new;
        /* Matrix memory allocation */
        xn  = calloc(a, sizeof(float*)); 
        Phi = calloc(a, sizeof(float*)); 
        Psi = calloc(a, sizeof(float*)); 

    	for(int i = 0;i<a;i++) {
    	    xn[i]  = calloc(b, sizeof(float));
    	    Phi[i] = calloc(dimbase, sizeof(float)); 
    	    Psi[i] = calloc(3, sizeof(float)); }
    	    
        P   = calloc(dimbase, sizeof(float*)); 
        Pd  = calloc(dimbase, sizeof(float*)); 
    	for(int i = 0;i<dimbase;i++) {
    	    P[i]  = calloc(dimbase, sizeof(float));
    	    Pd[i] = calloc(dimbase, sizeof(float));      } 
    	       	    
        Par2 = calloc(b, sizeof(float*)); 
        yhat_new = calloc(b, sizeof(float*));
        for(int i = 0;i<b;i++) {
            yhat_new[i] = calloc(a,sizeof(float));
    	    Par2[i]  = calloc(3, sizeof(float));    }    
   
        /* Data read from python objects in a proper matrix form */
	for (int i = 0 ;i < a; i++) {
	         yd[i] = *(yd0+i);
	        for (int j = 0 ;j < b; j++) {
	        x[i][j] = *(frame+i*b+j);   }	  	 
	   }
       	free(frame); 
	free(yd0);
	
	 /* Normalization of parameters columnwise normalization*/
	 for(int n=0;n<b;n++){
		 for(int r=0;r<a;r++){
		 xbase[r] = x[r][n]; }
		 
		  /* Find min max values of columns*/
		  min_val = xbase[0]; 
		  max_val = xbase[0];
		  
		 for (int ii = 1; ii < a; ii++){
		    if (xbase[ii] < min_val){
		     min_val = xbase[ii];   }
		    if (xbase[ii] > max_val){
		    		     max_val = xbase[ii];    
		    		     }}
	  	for(int i = 0; i < a; i++) {
    	  	xn[i][n] = 2*(xbase[i]-min_val)/(max_val-min_val)-1;   }
	 	 	     }// End of normalization  
	free(xbase); 	 	     
	/* Chebyshev basis of parameters */
     	for(int i=0;i<a;i++){
    		 int jj = 0;
    	 	for(int j=0;j<b;j++){
    		 Phi[i][jj]   = xn[i][j];
    		 Phi[i][jj+1] = 2*xn[i][j]*xn[i][j]-1;
    		 Phi[i][jj+2] = 4*xn[i][j]*xn[i][j]*xn[i][j]-3*xn[i][j];
    		 jj = jj+3;  	 }
     	 	 	 	 	   }// End of basis matrix
     	 	 	 	 	   
     /**********************RLSE ALGORITHM********************************/
     /* Initialization of P and Pd matrix, Par vector*/
     for(int r = 0; r < dimbase; r++) {
    	 Par[r] = 0.0;
            for(int c = 0; c < dimbase; c++)    {
                if(r == c){P[r][c] = delta*1.0;}
                else {P[r][c] = 0.0;}
                Pd[r][c] = 0.0;		      } }
                
     /* BEGINNING OF RLSE LOOP */
     /* Main loop of RLSE with index n */
     for(int n = 0; n < a; n++) {
    	 /* zero initial values*/
    	 yhat[n] = 0;
    	 res_prod = 0;
         for(int i=0;i<dimbase;i++)
         { u[i] = 0; Ka[i] = 0; Pu[i] = 0;}
    	 /* Rows of Phi matrix*/
    	 for (int j=0;j<dimbase;j++) {
             u[j] = Phi[n][j]; }
    	 /* Pu product */
    	 for (int i=0;i<dimbase;i++) {
             for (int j=0;j<dimbase;j++) {
    		  Pu[i] += P[i][j]*u[j]; }   }
    	 /* uPu product */
    	 for (int j=0;j<dimbase;j++) {
    	    		  res_prod += u[j]*Pu[j];}
    	 /* Kalman Gain and yhat */
     	 for(int j=0;j<dimbase;j++) {
                      Ka[j]    = Pu[j]/(lamda+res_prod);
     	              yhat[n] += Par[j]*u[j];}
     	 /*Parameter update*/
     	 for(int j=0;j<dimbase;j++) {
         	              Par[j] = Par[j]+Ka[j]*((yd[n])-yhat[n]); }
         /* K*P*u Product or Pd calculation*/
         for(int i=0;i<dimbase;i++) {
        	 for (int jj=0;jj<dimbase;jj++) {
                  Pd[jj][i] = Pu[jj]*Ka[i];
                                        }}
         /* P calculation */
         for(int i=0;i<dimbase;i++) {
                 for (int j=0;j<dimbase;j++)     {
                	 P[i][j] = (P[i][j]-Pd[i][j])/lamda; }
                  	 	 	 	 	 } }
     /********************************END OF RLSE****************************/
     
     /* Padded parameter matrix */
     for(int i=0;i<D;i++) {
    	 int j = 0;
    	 r = i*3; /* shifting the index of parameter vector */
    	 for(int jj=r;jj<r+3;jj++)   {
    		 Par2[i][j] = Par[jj];
    	     j++;	 	 	 	 	 } }
    	     
     /* Selection of best properties */
     for(int i = 0; i < b; i++) {
    	 for(int j = 0; j < a; j++) {
		 yhat_new[i][j] = 0.0; }}
		 
     for(int n = 0; n < b; n++) {  
		 /* Psi matrix */
		for(int i = 0; i < a; i++) {
                		Psi[i][0] = xn[i][n];
				Psi[i][1] = 2.0*xn[i][n]*xn[i][n]-1;
				Psi[i][2] = 4.0*xn[i][n]*xn[i][n]*xn[i][n]-3.0*xn[i][n];
				            }
               	for(int i = 0; i < a; i++) {
			for(int j = 0; j < 3; j++) { 
			    yhat_new[n][i] += Psi[i][j]*Par2[n][j];
						   }   
			        	   }					     
		 /* Sorting yhat_new */
		int min_index = 0;
		for (int ii = 1; ii < a; ii++){
		    if (yhat_new[n][ii] < yhat_new[n][min_index]){
		      min_index = ii; 
		   				 }}
		  best_params[n] = x[min_index][n];
		  		}  // end of property selection

	/*Construct PyObect for return*/
	result = PyList_New(D);
	for (int j = 0 ;j < D; j++) {
		res1 = Py_BuildValue("d", best_params[j]);
		PyList_SetItem(result, j, res1);
				    }  
	/* MAKE MEMORY FREE */
	free(yhat);
        free(Ka);
        free(Pu);
        free(Par);
        free(u);
        for (int i = 0; i < a; ++i) {
               free(xn[i]);
               free(Phi[i]);
               free(Psi[i]);        }
        for (int i = 0; i < dimbase; ++i) {
               free(P[i]);
               free(Pd[i]); } 
         for (int i = 0; i < b; ++i) {
               free(Par2[i]); 
               free(yhat_new[i]);    }                              
        free(xn);
        free(Phi);
        free(Psi); 
        free(P);
        free(Pd);
        free(Par2);
        free(yhat_new);
          
/* RETURN YOUR RESULT*/
return result;     
}  // end of main function
        

/*THIS PART IS TO TRANSFER INFORMATION BETWEEN PYTHON AND C EXTENSION*/
static PyMethodDef capiMethods[] = {
    {"aisa_feature",                    /* the Python method name that will be used */
      (PyCFunction) aisa_feature,       /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,                     /* flags indicating parameters
                                           accepted for this function */
      PyDoc_STR(" ")},                  /* The docstring for the function */
    {NULL, NULL, 0, NULL}               /* The last entry must be all NULL as shown to act as a
                                           sentinel. Python looks for this entry to know that all
                                           of the functions for the module have been defined. */
};
static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "ChebyshevFeatureSelection",   /* name of module */
    NULL,         /* module documentation, may be NULL */
    -1,           /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods   /* the PyMethodDef array from before containing the methods of the extension */
};
/*return the module to python*/
PyMODINIT_FUNC
PyInit_ChebyshevFeatureSelection(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef); 
    if (!m) {
        return NULL;
    }
    return m;
}


