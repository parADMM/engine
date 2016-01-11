//=============================================================================
//
//  paradmm.c
//
// == Authors ==
// Ning Hao
// Nate Derbinsky
//
//=============================================================================

#include "paradmm.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

//////////////////////////////////////////////////
// Utility methods
//////////////////////////////////////////////////

double abs_d(double x)
{
	if (x >=0) return x;
	return -x;
}

// implement quicksort, is used to sort variable index list for each operator passed by user
void quicksort(int* array, int n)
{
	if (n > 1)
	{
		int pivot_ix;
		pivot_ix = rand() % n;
		int pivot_val;
		pivot_val = array[pivot_ix];
		int left;
		int right;
		left = 0;
		right = n - 1;
		while (left <= right)
		{
				while (array[left] < pivot_val)
				{
					left = left + 1;
				}
				while (array[right] > pivot_val)
				{
					right = right - 1;
				}
				if (left <= right)
				{
					int tmp; tmp = array[left]; array[left] = array[right]; array[right] = tmp;
					left = left + 1;
					right = right - 1;
				}
		}
		quicksort(&array[0], right + 1);
		quicksort(&array[left], n - left);
	}
}

// get random number with upper and lower bounds
double get_random(double lowerbound, double upperbound)
{
	return ((double)rand()/RAND_MAX) * (upperbound - lowerbound) + lowerbound;
}

// this modify realloc guarantees that, if we call realloc starting from 0 elements we basically do a malloc
void *realloc2 (void* ptr, size_t sizesbefore, size_t sizesnow)
{
	if (sizesbefore == 0) return malloc(sizesnow);
	else return realloc(ptr, sizesnow);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

// initialized graph G
void startG (graph* G, int numdims)
{
	(*G).num_dims = numdims;
	(*G).num_edges = 0;
	(*G).num_opts = 0;
	(*G).num_vars = 0;
}

// it frees all the double arrays with edges and nodes IDs and the lists of operators
// and also the list of pointers for parameters
// but this function does not free the data that the list of pointers to parameters is pointing to
// that is of the responsibility of the user
void closeG(graph* G)
{
	free((*G).X);
	free((*G).M);
	free((*G).N);
	free((*G).U);
	free((*G).Z);
	free((*G).RHOS);
	free((*G).ALPHAS);

	// free all operatorix to ... variables
	for (int i=0; i<(*G).num_opts; i++) free(((*G).operatorix_to_edgelistix)[i]);
	free((*G).list_of_operators_param_pointers);
	free(((*G).operatorix_to_edgelistix));

	// free lists of sizes
	free((*G).list_of_operators_param_sizes);
	free((*G).operatorix_to_numedges);
	free((*G).varix_to_numedges);

	// free the list of function pointers. we don't free the functions themselves
	free((*G).list_of_operators_functions_pointers);

	// free all varix to ... variables
	for (int i=0; i<(*G).num_vars; i++) free(((*G).varix_to_edgelistix)[i]);
	free(((*G).varix_to_edgelistix));

	// free all edgeix to ... variables
	free(((*G).edgeix_to_varlistix));
}

// print graph G
void printG(graph* G)
{
	printf("number of operator in G is %d \n", (*G).num_opts);
	printf("\n");

	printf("number of global variables in G is %d \n", (*G).num_vars);
	printf("\n");

	printf("Dimension of each variable in G is %d \n", (*G).num_dims);
	printf("\n");

	printf("number of edges in G is %d \n", (*G).num_edges);
	printf("\n");


	printf("RHOS \n");
	for (int i=0; i<(*G).num_edges; i++) printf("%f \n", (*G).RHOS[i]);
	printf("\n");


	printf("ALPHAS \n");
	for(int i = 0 ; i < (*G).num_edges; i++) printf("%f\n", (*G).ALPHAS[i]);
	printf("\n");

	printf("List of operators_pointers \n");
	for (int i=0; i<(*G).num_opts; i++) printf("%p\n",(((*G).list_of_operators_functions_pointers)[i]));
	printf("\n");

	printf("List of operatorix_to_numedges \n");
	for (int i=0; i<(*G).num_opts; i++) printf("%d\n", (((*G).operatorix_to_numedges)[i]));
	printf("\n");

	printf("List of operatorix_to_edgelistin \n");
	for (int i=0; i<(*G).num_opts; i++)
	{
		printf("edge list of operator %d : ", i);
		int* temp_opt_edge = ((*G).operatorix_to_edgelistix)[i];
		for (int j=0; j<((*G).operatorix_to_numedges)[i]; j++) printf("\t%d", temp_opt_edge[j]);
		printf("\n");
	}
	printf("\n");

	printf("List of varix_to_numedges \n");
	for (int i=0; i<(*G).num_vars; i++) printf("%d\n",(((*G).varix_to_numedges)[i]));
	printf("\n");

	printf("List of varix_to_edgelistix \n");
	for (int i=0; i<(*G).num_vars; i++)
	{
		printf("edge list of variable %d : ", i);
		int* temp_var_edge = ((*G).varix_to_edgelistix)[i];
		for (int j=0; j<((*G).varix_to_numedges)[i]; j++) printf("\t%d", (temp_var_edge[j]));
		printf("\n");
	}
	printf("\n");

	printf("List of edgeix_to_varlistix: \n");
	for (int j=0; j<(*G).num_edges; j++) printf("%d ",(*G).edgeix_to_varlistix[j]);
	printf("\n");
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

// initialize RHOs and ALPHAS, here we use the same value of rho and alpha for all edges
void initialize_RHOS_APHAS(graph* G, double rho0, double alpha0)
{
	size_t size = (*G).num_edges* sizeof(double);

	(*G).RHOS = (double*) malloc(size);

	(*G).ALPHAS = (double*) malloc(size);

	if ((*G).RHOS == NULL || (*G).ALPHAS == NULL)
	{
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}

	for (int i=0; i<(*G).num_edges; i++)
	{
		(*G).RHOS[i] = rho0;
		(*G).ALPHAS[i] = alpha0;
	}
}

// initialize X, N, Z, M,U with zeros
void initialize_X_N_Z_M_U_zeros(graph* G)
{
	int length1 = (*G).num_dims * (*G).num_edges;
	int length2 = (*G).num_dims * (*G).num_vars;

	(*G).X = (double *) calloc(length1, sizeof(double));
	(*G).N = (double *) calloc(length1, sizeof(double));
	(*G).M = (double *) calloc(length1, sizeof(double));
	(*G).U = (double *) calloc(length1, sizeof(double));
	(*G).Z = (double *) calloc(length2, sizeof(double));

	if ((*G).X == NULL || (*G).N == NULL || (*G).M == NULL || (*G).U == NULL || (*G).Z == NULL)
	{
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}
}

// initialize X, N, Z, M,U with random numbers use get_random
void initialize_X_N_Z_M_U_rand(graph* G, double X_l, double X_u, double N_l, double N_u,double Z_l, double Z_u, double M_l, double M_u, double U_l, double U_u)
{

	size_t size1 = (*G).num_dims * (*G).num_edges * sizeof(double); // size of X,N,U,M
	size_t size2 = (*G).num_dims * (*G).num_vars* sizeof(double);	// size of Z

	(*G).X = (double*) malloc(size1);
	(*G).N = (double*) malloc(size1);
	(*G).M = (double*) malloc(size1);
	(*G).U = (double*) malloc(size1);
	(*G).Z = (double*) malloc(size2);

	if ((*G).X == NULL || (*G).N == NULL || (*G).M == NULL || (*G).U == NULL || (*G).Z == NULL)
	{
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}

	srand((unsigned)time(NULL));

	for(int i = 0; i<(*G).num_dims * (*G).num_edges; i++)
	{
		(*G).X[i] = get_random(X_l, X_u); //generate a random real number between X_l and X_u. (*G).X[i] = (double) (((double)rand())/RAND_MAX);
		(*G).N[i] = get_random(N_l, N_u); //generate a random real number between N_l and N_u.(*G).N[i] = (double) (((double)rand())/RAND_MAX);
		(*G).M[i] = get_random(M_l, M_u); //generate a random real number between M_l and M_u.(*G).M[i] = (double) (((double)rand())/RAND_MAX);
		(*G).U[i] = get_random(U_l, U_u); //generate a random real number between U_l and U_u.(*G).U[i] = (double) (((double)rand())/RAND_MAX);
	}

	for (int i = 0; i<((*G).num_dims) * ((*G).num_vars); i++)
	{
		(*G).Z[i] = get_random(Z_l, Z_u); //generate a random real number between Z_l and Z_u.(*G).Z[i] = (double) (((double)rand())/RAND_MAX);
	}

}

// initialize X, N, Z, M,U with optimal choice of N,Z,M,U with the initialized Z
// rhos should be initialized before calling this function.
void initialize_X_N_Z_M_U_optimal(graph* G, double Zinit[], double Uinit[], int maxiter)
{
	int num_dims = (*G).num_dims;
	int length1 = num_dims * (*G).num_edges;
	int length2 = num_dims * (*G).num_vars;

	size_t size1 = length1 * sizeof(double); // size of X,N,U,M
	size_t size2 = length2 * sizeof(double); // size of Z

	(*G).X = (double*) malloc(size1);
	(*G).N = (double*) malloc(size1);
	(*G).M = (double*) malloc(size1);
	(*G).U = (double*) malloc(size1);
	(*G).Z = (double*) malloc(size2);

	if ((*G).X == NULL || (*G).N == NULL || (*G).M == NULL || (*G).U == NULL || (*G).Z == NULL)
	{
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}

	// initialize Z
	for (int i=0; i<length2; i++)
	{
		(*G).Z[i] = Zinit[i];
	}

	// initialize X from Z
	for (int i=0; i<(*G).num_vars; i++)
	{
		int edge_num;
		edge_num = (*G).varix_to_numedges[i];

		for (int k=0; k<edge_num; k++)
		{
				int edge_ix = (*G).varix_to_edgelistix[i][k];

				for (int j=0; j<num_dims; j++)
				{
					(*G).X[edge_ix * num_dims + j] = Zinit[i * num_dims + j];
				}
		}
	}

	double *xminusu;
	xminusu = (double*) malloc(size1);
	for (int i=0; i<length1; i++)
	{
		xminusu[i] = (*G).X[i] - Uinit[i];
	}

	double *Xiter;
	Xiter = (double*) malloc(size1);
	for (int i=0; i<length1; i++)
	{
		Xiter[i] = (*G).X[i];
	}

	double* tempRHOS;
	tempRHOS = (double*) malloc((*G).num_edges * sizeof(double));
	for (int i=0; i<(*G).num_edges; i++){
		tempRHOS[i] = (*G).RHOS[i];
	}

	// initialize U
	for (int opt_ix=0; opt_ix<(*G).num_opts; opt_ix++)
	{
		//get each operator
		void (*curr_opt)(double*, double*, double*, int, int, void*);
		curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];

		// get number of edges (number of variables) of each operator
		int num_edges;
		num_edges = (*G).operatorix_to_numedges[opt_ix];

		// get parameters of each operator
		void * params;
		params = (*G).list_of_operators_param_pointers[opt_ix];

		// starting position of x for each operator
		int start_edge_ix;
		start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];

		for (int i=0; i<maxiter; i++)
		{
				curr_opt(&(Xiter[num_dims * start_edge_ix]), &(xminusu[num_dims * start_edge_ix]), &(tempRHOS[start_edge_ix]), num_edges, num_dims, params);
				for(int j=0; j<(num_dims * num_edges); j++)
				{
					(*G).U[num_dims * start_edge_ix + j] = - xminusu[num_dims * start_edge_ix + j] + Xiter[num_dims * start_edge_ix + j];
					xminusu[num_dims * start_edge_ix + j] = (*G).X[num_dims * start_edge_ix + j] - (*G).U[num_dims * start_edge_ix + j];
				}
		}
	}

	// initiazlie M
	for (int i=0; i<length1; i++)
	{
		(*G).M[i] = (*G).X[i] + (*G).U[i];
	}

	// initialize N
	for (int i=0; i<length1; i++)
	{
		(*G).N[i] = xminusu[i];
	}

	free(Xiter);
	free(tempRHOS);
	free(xminusu);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

// add nodes to graph G, update corresponding fields
void addNode(graph* G, void (*operator)(double*, double*, double*, int, int, void*), void* operator_para, int opt_numvars, int* listvarix)
{
	quicksort(listvarix, opt_numvars);

	(*G).num_opts += 1;

	// find maximum variable index in listvarix
	int max_index = 0;
	for (int i=0; i<opt_numvars; i++)
	{
		if (max_index < listvarix[i])
		{
				max_index = listvarix[i];
		}
	}
	// update the maximum number of vairiables used so far
	int prev_num_vars = (*G).num_vars; // we keep track of the previous number of variables
	if ((*G).num_vars < max_index + 1)
	{
		(*G).num_vars = max_index + 1;
	}

	// create the space so that we can increament the list of function pointers and respective parameters
	void (*tmp)(double*, double*, double*, int, int, void*);
	(*G).list_of_operators_functions_pointers = realloc2((*G).list_of_operators_functions_pointers, (-1 + (*G).num_opts) * sizeof(tmp), ((*G).num_opts) * sizeof(tmp));
	(*G).list_of_operators_param_pointers = (void**) realloc2((*G).list_of_operators_param_pointers, (-1 + (*G).num_opts ) * sizeof(void*), ((*G).num_opts) * sizeof(void*));

	// add the new operator and parameters to the list
	((*G).list_of_operators_functions_pointers)[((*G).num_opts) - 1] = operator;
	((*G).list_of_operators_param_pointers)[((*G).num_opts) - 1] = operator_para;

	// extend the list with the number of edges incident in each operator
	(*G).operatorix_to_numedges = (int*) realloc2((*G).operatorix_to_numedges, (-1 + (*G).num_opts ) * sizeof(int), ((*G).num_opts) * sizeof(int));

	// store the number of edges incident on the function node
	((*G).operatorix_to_numedges)[((*G).num_opts) - 1] = opt_numvars;

	// extend the list of lists associated to the edges associated to the operators
	(*G).operatorix_to_edgelistix = (int**) realloc2((*G).operatorix_to_edgelistix, (-1 + (*G).num_opts ) * sizeof(int*), ((*G).num_opts) * sizeof(int*));

	// extend the list of edges associated to the operator at hand
	(*G).operatorix_to_edgelistix[-1 + (*G).num_opts] = (int*) realloc2((*G).operatorix_to_edgelistix[-1 + (*G).num_opts], 0 * sizeof(int) ,(opt_numvars) * sizeof(int));

	for (int i=0; i<opt_numvars; i++)
	{
		(*G).operatorix_to_edgelistix[-1 + (*G).num_opts][i] = (*G).num_edges + i;
	}

	// extend the lists associated with the variable nodes
	(*G).varix_to_numedges = (int*) realloc2((*G).varix_to_numedges, prev_num_vars * sizeof(int) ,((*G).num_vars) * sizeof(int));

	// THIS IS VERY IMPORTANT ! BECAUSE WE ARE INCREMENTING THIS VARIABLE ONE BY ONE WE NEED TO MAKE SURE WE START AT ZERO. SIMILAR BUGS MIGHT EXIST IN OTHER, SIMILAR, PARTS OF THE CODE
	if (prev_num_vars < (*G).num_vars){
		for (int i=prev_num_vars; i<((*G).num_vars); i++)
		{
				(*G).varix_to_numedges[i] = 0;
		}
	}

	(*G).varix_to_edgelistix = (int**) realloc2((*G).varix_to_edgelistix, prev_num_vars * sizeof(int*) ,((*G).num_vars) * sizeof(int*));

	for (int i=0; i<opt_numvars; i++)
	{

		int tmp_var_ix = listvarix[i]; //current variable node being proccess
		int tmp_edge_ix = ((*G).operatorix_to_edgelistix)[(*G).num_opts - 1][i]; //edge associated to the current variable node being processed

		(*G).varix_to_numedges[tmp_var_ix] += 1;

		// allocate space in the list of edges seens from the current variable node
		((*G).varix_to_edgelistix)[tmp_var_ix] = (int*) realloc2( ((*G).varix_to_edgelistix)[tmp_var_ix], (-1 + (*G).varix_to_numedges[tmp_var_ix]) * sizeof(int),((*G).varix_to_numedges[tmp_var_ix]) * sizeof(int));
		((*G).varix_to_edgelistix)[tmp_var_ix][-1 + (*G).varix_to_numedges[tmp_var_ix]] = tmp_edge_ix;
	}

	// we increment the number of edges used which also corresponds to the last edge index we can use
	(*G).num_edges += opt_numvars;

	int num_edges = (*G).num_edges;
	// extend the list of variables associated to each edge
	(*G).edgeix_to_varlistix = (int*) realloc2((*G).edgeix_to_varlistix , (num_edges - opt_numvars) * sizeof(int) ,(num_edges) * sizeof(int));

	for (int i=0; i<opt_numvars; i++)
	{
		(*G).edgeix_to_varlistix[num_edges - opt_numvars + i] = listvarix[i];
	}
}

// This function is very similar to the addNode of the C-framework but also takes into account the size of each parameter associated to each operator
// this size is necessary for when we copy the CPU-graph to the GPU-graph.
void addNodeWithParamSize(graph* G, void (*operator)(double *, double *, double *, int, int, void *), void *operator_para, int operator_para_size, int opt_numvars, int *listvarix)
{
	quicksort(listvarix, opt_numvars);

	(*G).num_opts += 1;

	// find maximum variable index in listvarix
	int max_index = 0;
	for (int i=0; i<opt_numvars; i++)
	{
		if (max_index < listvarix[i])
		{
				max_index = listvarix[i];
		}
	}
	// update the maximum number of vairiables used so far
	int prev_num_vars = (*G).num_vars; // we keep track of the previous number of variables
	if ((*G).num_vars < max_index + 1)
	{
		(*G).num_vars = max_index + 1;
	}

	// create the space so that we can increament the list of function pointers and respective parameters
	void (*tmp)(double*, double*, double*, int, int, void*);
	(*G).list_of_operators_functions_pointers =  (void (**)(double*, double*, double*, int, int, void*)) realloc2((*G).list_of_operators_functions_pointers, (-1 + (*G).num_opts) * sizeof(tmp), ((*G).num_opts) * sizeof(tmp));
	(*G).list_of_operators_param_pointers = (void**) realloc2((*G).list_of_operators_param_pointers, (-1 + (*G).num_opts ) * sizeof(void*), ((*G).num_opts) * sizeof(void*));
	(*G).list_of_operators_param_sizes = (int*) realloc2((*G).list_of_operators_param_sizes, (-1 + (*G).num_opts ) * sizeof(int), ((*G).num_opts) * sizeof(int));

	// add the new operator and parameters to the list
	((*G).list_of_operators_functions_pointers)[((*G).num_opts) - 1] = operator;
	((*G).list_of_operators_param_pointers)[((*G).num_opts) - 1] = operator_para;
	((*G).list_of_operators_param_sizes)[((*G).num_opts) - 1] = operator_para_size;

	// extend the list with the number of edges incident in each operator
	(*G).operatorix_to_numedges = (int *) realloc2((*G).operatorix_to_numedges, (-1 + (*G).num_opts ) * sizeof(int), ((*G).num_opts) * sizeof(int));

	// store the number of edges incident on the function node
	((*G).operatorix_to_numedges)[((*G).num_opts) - 1] = opt_numvars;

	// extend the list of lists associated to the edges associated to the operators
	(*G).operatorix_to_edgelistix = (int**) realloc2((*G).operatorix_to_edgelistix, (-1 + (*G).num_opts ) * sizeof(int*), ((*G).num_opts) * sizeof(int*));

	// extend the list of edges associated to the operator at hand
	(*G).operatorix_to_edgelistix[-1 + (*G).num_opts] = (int*) realloc2((*G).operatorix_to_edgelistix[-1 + (*G).num_opts], 0 * sizeof(int), (opt_numvars) * sizeof(int) );

	for (int i=0; i<opt_numvars; i++)
	{
		(*G).operatorix_to_edgelistix[-1 + (*G).num_opts][i] = (*G).num_edges + i;
	}

	// extend the lists associated with the variable nodes
	(*G).varix_to_numedges = (int*) realloc2((*G).varix_to_numedges, prev_num_vars * sizeof(int) ,((*G).num_vars) * sizeof(int));

	// THIS IS VERY IMPORTANT ! BECAUSE WE ARE INCREMENTING THIS VARIABLE ONE BY ONE WE NEED TO MAKE SURE WE START AT ZERO. SIMILAR BUGS MIGHT EXIST IN OTHER, SIMILAR, PARTS OF THE CODE
	if (prev_num_vars < (*G).num_vars){
		for (int i=prev_num_vars; i<((*G).num_vars); i++)
		{
				(*G).varix_to_numedges[i] = 0;
		}
	}

	(*G).varix_to_edgelistix = (int **) realloc2( (*G).varix_to_edgelistix, prev_num_vars * sizeof(int *) ,((*G).num_vars) * sizeof(int *)  );

	for (int i=0; i<opt_numvars; i++)
	{
		int tmp_var_ix = listvarix[i]; //current variable node being proccess
		int tmp_edge_ix = ((*G).operatorix_to_edgelistix)[(*G).num_opts - 1][i]; //edge associated to the current variable node being processed

		(*G).varix_to_numedges[tmp_var_ix] += 1;

		// allocate space in the list of edges seens from the current variable node
		((*G).varix_to_edgelistix)[tmp_var_ix] = (int *) realloc2( ((*G).varix_to_edgelistix)[tmp_var_ix], (-1 + (*G).varix_to_numedges[tmp_var_ix]) * sizeof(int),((*G).varix_to_numedges[tmp_var_ix]) * sizeof(int)  );

		((*G).varix_to_edgelistix)[tmp_var_ix][-1 + (*G).varix_to_numedges[tmp_var_ix]] = tmp_edge_ix;
	}

	// we increment the number of edges used which also corresponds to the last edge index we can use
	(*G).num_edges += opt_numvars;

	int num_edges = (*G).num_edges;
	// extend the list of variables associated to each edge
	(*G).edgeix_to_varlistix = (int*) realloc2((*G).edgeix_to_varlistix , (num_edges - opt_numvars) * sizeof(int) ,(num_edges) * sizeof(int));

	for (int i=0; i<opt_numvars; i++)
	{
		(*G).edgeix_to_varlistix[num_edges - opt_numvars + i] = listvarix[i];
	}

}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

// placeholder
void update_outgoing_rhos(graph* G)
{
}

// update X and M
void updateXM(graph* G)
{
	int num_dims;
	num_dims = (*G).num_dims;

	for(int opt_ix=0; opt_ix<(*G).num_opts; opt_ix++)
	{

		void (*curr_opt)(double*, double*, double*, int, int, void*);
		curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];

		int num_vars;
		num_vars = (*G).operatorix_to_numedges[opt_ix];

		void * params;
		params = (*G).list_of_operators_param_pointers[opt_ix];

		// starting position of x for each operator
		int start_edge_ix;
		start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];

		curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);
	}

	for (int i=0; i<((*G).num_edges) * (num_dims); i++)
	{
		(*G).M[i] = (*G).X[i] + (*G).U[i];
	}
}


void updateX(graph* G)
{
	int num_dims;
	num_dims = (*G).num_dims;

	for(int opt_ix=0; opt_ix<(*G).num_opts; opt_ix++)
	{

		void (*curr_opt)(double *, double *, double *, int, int, void *);
		curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];

		int num_vars;
		num_vars = (*G).operatorix_to_numedges[opt_ix];

		void* params;
		params = (*G).list_of_operators_param_pointers[opt_ix];

		// starting position of x for each operator
		int start_edge_ix;
		start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];

		curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);
	}
}

void updateM(graph* G)
{
	int num_dims;
	num_dims = (*G).num_dims;

	for (int i=0; i<((*G).num_edges) * (num_dims); i++)
	{
		(*G).M[i] = (*G).X[i] + (*G).U[i];
	}
}

// update Z
void updateZ(graph* G)
{
	int num_dims = (*G).num_dims;

	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		// computing the denominator
		double denomrhos;
		denomrhos = 0;

		for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
		{
				int edge_ix;
				edge_ix = (*G).varix_to_edgelistix[var_ix][i];
				denomrhos += (*G).RHOS[edge_ix];
		}

		for (int j=0; j<num_dims; j++)
		{
				double tempZ;
				tempZ = 0;

				for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
				{
					int edge_ix;
					edge_ix = (*G).varix_to_edgelistix[var_ix][i];

					tempZ = tempZ + (*G).M[num_dims * edge_ix + j] * ((*G).RHOS[edge_ix]);
				}

				(*G).Z[num_dims * var_ix + j] = tempZ /denomrhos;
		}
	}

	update_outgoing_rhos(G);
}

// update U and N
void updateUN(graph* G)
{
	int num_dims = (*G).num_dims;
	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
		{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[var_ix][j];

				double alpha = (*G).ALPHAS[edge_ix];

				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).U[edge_ix * num_dims + dim_ix] += alpha * ((*G).X[num_dims * edge_ix + dim_ix] - (*G).Z[num_dims * var_ix + dim_ix]);
					(*G).N[edge_ix * num_dims + dim_ix] = (*G).Z[num_dims * var_ix + dim_ix] - (*G).U[edge_ix * num_dims + dim_ix];
				}
		}
	}
}

void updateU(graph* G)
{
	int num_dims = (*G).num_dims;
	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
		{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[var_ix][j];

				double alpha = (*G).ALPHAS[edge_ix];

				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).U[edge_ix * num_dims + dim_ix] += alpha * ((*G).X[num_dims * edge_ix + dim_ix] - (*G).Z[num_dims * var_ix + dim_ix]);

				}
		}
	}
}


void updateN(graph* G)
{
	int num_dims = (*G).num_dims;
	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
		{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[var_ix][j];

				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).N[edge_ix * num_dims + dim_ix] = (*G).Z[num_dims * var_ix + dim_ix] - (*G).U[edge_ix * num_dims + dim_ix];
				}
		}
	}
}

void updateG(graph* G)
{
	updateXM(G);
	updateZ(G);
	updateUN(G);
}

//////////////////////////////////////////////////
// Matlab support
//////////////////////////////////////////////////

#ifdef MATLAB

void matlab_proximal_operator(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
{
	matparams* matparams = params;
	Engine* eng = matparams->eng;
	mxArray *n_mat = NULL, *rhos_mat = NULL, *params_mat = NULL, *output_mat = NULL, *outrhos_mat = NULL;

	n_mat = mxCreateDoubleMatrix(numdims, numvars, mxREAL);

	for (int i=0; i<numdims * numvars; i++)
	{
		mxGetPr(n_mat)[i] = (double) n[i];
	}

	rhos_mat = mxCreateDoubleMatrix(numvars, 1, mxREAL);
	for (int i=0; i<numvars; i++)
	{
		mxGetPr(rhos_mat)[i] = (double) rhos[i];
	}

	int paramsize = matparams->opt_params_size;
	params_mat = mxCreateDoubleMatrix(paramsize, 1, mxREAL);
	for (int i=0; i<paramsize; i++)
	{
		mxGetPr(params_mat)[i] = (double) (*matparams).params[i];
	}

	engPutVariable(eng, "N", n_mat);
	engPutVariable(eng, "RHOS", rhos_mat);
	engPutVariable(eng, "PARAMS", params_mat);

	char* optname = matparams->opt_name;
	char* leftstr = "[OUTPUT, OUTRHOS] = ";
	char* rightstr = "(N, RHOS, PARAMS)";
	char* callerstring = malloc((20 + 18 + strlen(optname)) * sizeof(char));
	strcpy(callerstring, leftstr);
	strcat(callerstring, optname);
	strcat(callerstring, rightstr);

	engEvalString(eng, callerstring);

	output_mat = engGetVariable(eng, "OUTPUT");
	for (int i=0; i<numvars*numdims; i++)
	{
		x[i] = (double) mxGetPr(output_mat)[i];
	}

	outrhos_mat = engGetVariable(eng, "OUTRHOS");
	for (int i=0; i<numvars; i++)
	{
		rhos[i] = (double) mxGetPr(outrhos_mat)[i];
	}
}

#endif

//////////////////////////////////////////////////
// OpenMP support
//////////////////////////////////////////////////

#ifdef OPENMP

#include <omp.h>

void updateXOpenMP(graph* G)
{
	int num_dims;
	num_dims = (*G).num_dims;

#pragma omp parallel for
	for (int opt_ix=0; opt_ix<(*G).num_opts; opt_ix++)
	{
		void (*curr_opt)(double*, double*, double*, int, int, void*);
		curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];

		int num_vars;
		num_vars = (*G).operatorix_to_numedges[opt_ix];

		void * params;
		params = (*G).list_of_operators_param_pointers[opt_ix];

		// starting position of x for each operator
		int start_edge_ix;
		start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];

		curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);

	}

}

void updateMOpenMP(graph* G)
{

	int num_dims;
	num_dims = (*G).num_dims;

#pragma omp parallel for
	for (int i=0; i<((*G).num_edges) * (num_dims); i++)
	{
		(*G).M[i] = (*G).X[i] + (*G).U[i];
	}

}

void updateZOpenMP(graph* G)
{
	int num_dims = (*G).num_dims;

#pragma omp parallel for
	for(int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		// computing the denominator
		double denomrhos;
		denomrhos = 0;

		for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
		{
				int edge_ix;
				edge_ix = (*G).varix_to_edgelistix[var_ix][i];
				denomrhos += (*G).RHOS[edge_ix];
		}


		for (int j=0; j<num_dims; j++)
		{

				double tempZ;
				tempZ = 0;

				for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
				{
					int edge_ix;
					edge_ix = (*G).varix_to_edgelistix[var_ix][i];

					tempZ = tempZ + (*G).M[num_dims * edge_ix + j] * ((*G).RHOS[edge_ix]);

				}

				(*G).Z[num_dims * var_ix + j] = tempZ /denomrhos;

		}

	}

	update_outgoing_rhos(G);
}

void updateUOpenMP(graph* G)
{
	int num_dims = (*G).num_dims;

#pragma omp parallel for
	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
		{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[var_ix][j];

				double alpha = (*G).ALPHAS[edge_ix];

				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).U[edge_ix * num_dims + dim_ix] += alpha * ((*G).X[num_dims * edge_ix + dim_ix] - (*G).Z[num_dims * var_ix + dim_ix]);

				}

		}

	}

}


void updateNOpenMP(graph* G)
{
	int num_dims = (*G).num_dims;

#pragma omp parallel for
	for (int var_ix=0; var_ix<(*G).num_vars; var_ix++)
	{
		for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
		{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[var_ix][j];


				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{

					(*G).N[edge_ix * num_dims + dim_ix] = (*G).Z[num_dims * var_ix + dim_ix] - (*G).U[edge_ix * num_dims + dim_ix];
				}

		}

	}
}



void updateGOpenMP(graph* G, int maxiter, int numthreads)
{
	omp_set_num_threads(numthreads);

	int num_dims;
	num_dims = (*G).num_dims;

#pragma omp parallel
	{
		int id = omp_get_thread_num();

		for(int iterIX=0; iterIX<maxiter; iterIX++)
		{
				// update the variable X
				int N = (*G).num_opts;
				int istart = id * N / numthreads;
				int iend = (id+1) * N / numthreads;
				if (id == numthreads - 1) iend = N;
				for (int opt_ix=istart; opt_ix<iend; opt_ix++){

					void (*curr_opt)(double *, double *, double *, int, int, void *);
					curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];

					int num_vars;
					num_vars = (*G).operatorix_to_numedges[opt_ix];

					void * params;
					params = (*G).list_of_operators_param_pointers[opt_ix];

					// starting position of x for each operator
					int start_edge_ix;
					start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];

					curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);
				}
#pragma omp barrier
				// now everything is syncronized. This also syncronizes memory so we can keep going


				// update the variable M
				N = ((*G).num_edges) * (num_dims);
				istart = id * N / numthreads;
				iend = (id+1) * N / numthreads;
				if (id == numthreads - 1) iend = N;
				for (int i=istart; i<iend; i++)
				{
					(*G).M[i] = (*G).X[i] + (*G).U[i];
				}
#pragma omp barrier
				// now everything is syncronized. This also syncronizes memory so we can keep going

				// update the variable Z
				N = (*G).num_vars;
				istart = id * N / numthreads;
				iend = (id+1) * N / numthreads;
				if (id == numthreads - 1) iend = N;
				for (int var_ix=istart; var_ix<iend; var_ix++){

					// computing the denominator
					double denomrhos;
					denomrhos = 0;

					for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
					{
						int edge_ix;
						edge_ix = (*G).varix_to_edgelistix[var_ix][i];
						denomrhos += (*G).RHOS[edge_ix];
					}

					for (int j=0; j<num_dims; j++)
					{

						double tempZ;
						tempZ = 0;

						for (int i=0; i<(*G).varix_to_numedges[var_ix]; i++)
						{
								int edge_ix;
								edge_ix = (*G).varix_to_edgelistix[var_ix][i];

								tempZ = tempZ + (*G).M[num_dims * edge_ix + j] * ((*G).RHOS[edge_ix]);

						}

						(*G).Z[num_dims * var_ix + j] = tempZ /denomrhos;

					}

				}
				// we are skipping the update the rhos line. We have not parallized that line yet !
#pragma omp barrier
				// now everything is syncronized. This also syncronizes memory so we can keep going

				// update the variable U
				N = (*G).num_vars;
				istart = id * N / numthreads;
				iend = (id+1) * N / numthreads;
				if (id == numthreads - 1) iend = N;
				for (int var_ix = istart; var_ix<iend; var_ix++)
				{
					for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
					{
						int edge_ix;
						edge_ix= (*G).varix_to_edgelistix[var_ix][j];

						double alpha = (*G).ALPHAS[edge_ix];

						for (int dim_ix = 0; dim_ix <(*G).num_dims; dim_ix++)
						{
								(*G).U[edge_ix * num_dims + dim_ix] += alpha * ((*G).X[num_dims * edge_ix + dim_ix] - (*G).Z[num_dims * var_ix + dim_ix]);

						}

					}
				}
#pragma omp barrier // I don't think we need this sink here.
				// now everything is syncronized. This also syncronizes memory so we can keep going

				// update the variable N
				N = (*G).num_vars;
				istart = id * N / numthreads;
				iend = (id+1) * N / numthreads;
				if (id == numthreads - 1) iend = N;
				for (int var_ix=istart; var_ix<iend; var_ix++)
				{
					for (int j=0; j<(*G).varix_to_numedges[var_ix]; j++)
					{
						int edge_ix;
						edge_ix= (*G).varix_to_edgelistix[var_ix][j];

						for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
						{
								(*G).N[edge_ix * num_dims + dim_ix] = (*G).Z[num_dims * var_ix + dim_ix] - (*G).U[edge_ix * num_dims + dim_ix];
						}

					}
				}
#pragma omp barrier
				// now everything is syncronized. This also syncronizes memory so we can keep going
		}

	}

}

#endif

//////////////////////////////////////////////////
//////////////////////////////////////////////////
