//=============================================================================
//
//  paradmm.h
//
// == Authors ==
// Ning Hao
// Nate Derbinsky
//
//=============================================================================

#ifndef __PARADMM_H
#define __PARADMM_H

#ifdef NEEDS_EXTERN_C
#define MY_EXTERN_C extern "C"
#else
#define MY_EXTERN_C
#endif

#include <stddef.h>

typedef struct params_matrix
{
	int flag_1stiter;
	int rows;
	double* A;
	double* B;
	double* R;
	double* ATB;
	double* ATA;
} mparams;

typedef struct graphs
{
	int num_opts;  // total number of operators in the graph
	int num_vars;  // total number of gloval variables in the graph
	int num_dims;	// dimension of each variable

	int num_edges; // total number of edges in the graph

	void (**list_of_operators_functions_pointers)(double*, double*, double*, int, int, void*); // pointer of pointers pointing to each operator
	void** list_of_operators_param_pointers; // pointer of pointers pointing to parameters of each operator
	int* list_of_operators_param_sizes; // this is the size of each parameter

	int* operatorix_to_numedges;	// pointer pointing to the number of edges of each operator
	int** operatorix_to_edgelistix;// pointer of pointers poing to the list of edges of each operator

	int* varix_to_numedges; // pointer pointing to the number of edges of each variable
	int** varix_to_edgelistix; // pointer of pointers pointing to the list of edges of each variable

	int* edgeix_to_varlistix; //pointer of pointers pointing to the list of variables of each edg

	double* X; // X is of size num_dims * num_edges
	double* M; // M is of size num_dims * num_edges
	double* N; // N is of size num_dims * num_edges
	double* U; // U is of size num_dims * num_edges

	double* Z; // Z is of size num_dims * num_vars

	double* RHOS; // one RHO per edge. RHOS is of size num_edges
	double* ALPHAS; // one ALPHA per edge. ALPHA is of size num_edges

} graph;

//////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////

MY_EXTERN_C double abs_d(double x);
MY_EXTERN_C void quicksort(int* array, int n);
MY_EXTERN_C double get_random(double lowerbound, double upperbound);

MY_EXTERN_C void *realloc2 (void* ptr, size_t sizesbefore, size_t sizesnow);

//////////////////////////////////////////////////
//////////////////////////////////////////////////

MY_EXTERN_C void startG(graph* G, int numdims);
MY_EXTERN_C void closeG(graph* G);
MY_EXTERN_C void printG(graph* G);

MY_EXTERN_C void initialize_RHOS_APHAS(graph* G, double rho0, double alpha0);
MY_EXTERN_C void initialize_X_N_Z_M_U_zeros(graph* G);
MY_EXTERN_C void initialize_X_N_Z_M_U_rand(graph* G, double X_l, double X_u, double N_l, double N_u,double Z_l, double Z_u, double M_l, double M_u, double U_l, double U_u);
MY_EXTERN_C void initialize_X_N_Z_M_U_optimal(graph* G, double Zinit[], double Uinit[], int maxiter);

MY_EXTERN_C void addNode(graph* G, void (*op)(double*, double*, double*, int, int, void*), void* operator_para, int opt_numVars, int* listvarix);
MY_EXTERN_C void addNodeWithParamSize(graph* G, void (*Operator)(double*, double*, double*, int, int, void*), void* operator_para, int operator_para_size, int opt_numvars, int* listvarix);

MY_EXTERN_C void update_outgoing_rhos(graph* G);

MY_EXTERN_C void updateXM(graph* G);
MY_EXTERN_C void updateX(graph* G);
MY_EXTERN_C void updateM(graph* G);

MY_EXTERN_C void updateZ(graph* G);

MY_EXTERN_C void updateUN(graph* G);
MY_EXTERN_C void updateU(graph* G);
MY_EXTERN_C void updateN(graph* G);

MY_EXTERN_C void updateG(graph* G);

//////////////////////////////////////////////////
// OpenMP support
//////////////////////////////////////////////////

#ifdef OPENMP
MY_EXTERN_C void updateXOpenMP(graph* G);
MY_EXTERN_C void updateMOpenMP(graph* G);
MY_EXTERN_C void updateZOpenMP(graph* G);
MY_EXTERN_C void updateUOpenMP(graph* G);
MY_EXTERN_C void updateNOpenMP(graph* G);
MY_EXTERN_C void updateGOpenMP(graph* G, int maxiter, int numthreads);
#endif

#endif
