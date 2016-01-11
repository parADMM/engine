//=============================================================================
//
//  paradmm-gpu.cu
//
// == Authors ==
// Ning Hao
// Nate Derbinsky
//
//=============================================================================

#include <stdio.h>
#include "paradmm-gpu.h"

__device__ double d_abs_d(double x)
{
	if (x >=0){
		return x;
	}

	return -x;
}

__global__ void d_updateX(graph *G)
{
	int opt_ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	
	if (opt_ix < (*G).num_opts) //this makes the code work even if we accidentally launch more threads than necessary
	{
		void (*curr_opt)(double *, double *, double *, int, int, void *);
		curr_opt = (*G).list_of_operators_functions_pointers[opt_ix];
		
		int num_vars;
		num_vars = (*G).operatorix_to_numedges[opt_ix];
		
		int num_dims;
		num_dims = (*G).num_dims;
		
		void * params;
		params = (*G).list_of_operators_param_pointers[opt_ix];
		
		// starting position of x for each operator
		int start_edge_ix;
		start_edge_ix = (*G).operatorix_to_edgelistix[opt_ix][0];
		
		curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);
	}
}

__global__ void d_updateM(graph *G) 
{	
	int ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	
	int num_dims;
	num_dims = (*G).num_dims;
	
	//compute the M messages, each component goes to one thread;
	if(ix < ((*G).num_edges) * (num_dims))
	{
		(*G).M[ix] = (*G).X[ix] + (*G).U[ix];
	}
}

__global__ void d_updateZ(graph* G) 
{	
	int num_dims = (*G).num_dims;
	//int var_ix = blockIdx.x*blockDim.x + threadIdx.x;
	int var_ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
		
	if (var_ix < (*G).num_vars)
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
		
		
		for (int j=0; j < num_dims; j++)
		{
			double tempZ;
			tempZ = 0;
			
			for(int i = 0; i < (*G).varix_to_numedges[var_ix]; i++)
			{
				int edge_ix;
				edge_ix = (*G).varix_to_edgelistix[var_ix][i];
				
				tempZ = tempZ + (*G).M[num_dims * edge_ix + j] * ((*G).RHOS[edge_ix]);
				
			}
			
			(*G).Z[num_dims * var_ix + j] = tempZ /denomrhos;
		}
	}
}

__global__ void d_updateU_per_edge(graph *G)
{
	int num_dims = (*G).num_dims;
	int edge_ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	
	if (edge_ix < (*G).num_edges)
	{
		double alpha = (*G).ALPHAS[edge_ix];
		
		for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
		{
			(*G).U[edge_ix * num_dims + dim_ix] += alpha * ((*G).X[num_dims * edge_ix + dim_ix] - (*G).Z[num_dims * ((*G).edgeix_to_varlistix[edge_ix]) + dim_ix]);		
		}
	}	
}

__global__ void d_updateU(graph *G)
{
	int num_dims = (*G).num_dims;
	int var_ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	
	if (var_ix < (*G).num_vars)
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

__global__ void d_updateN_per_edge(graph *G) {

		int ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

		int num_dims;
		num_dims = (*G).num_dims;

		//compute the N messages, each component goes to one thread;
		if (ix < ((*G).num_edges))
		{
			for (int i=0; i<num_dims; i++)
			{
				(*G).N[ix*num_dims + i] = (*G).Z[  num_dims*((*G).edgeix_to_varlistix[ix]) + i  ] - (*G).U[ix*num_dims + i];
			}			
		}
}

__global__ void d_updateN(graph *G)
{
	int num_dims = (*G).num_dims;
	int var_ix = blockIdx.y*gridDim.x*blockDim.x*blockDim.y + blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
		
	if (var_ix < (*G).num_vars)
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

__global__ void d_updateG_XMZUN_single_block(graph *G, int numiterations)
{
	int ix = blockIdx.x * blockDim.x + threadIdx.x;
	
	for (int iterix=1; iterix<=numiterations; iterix++)
	{
		if (ix < (*G).num_opts) //this makes the code work even if we accidentally launch more threads than necessary
		{
			void (*curr_opt)(double *, double *, double *, int, int, void *);
			curr_opt = (*G).list_of_operators_functions_pointers[ix];
			
			int num_vars;
			num_vars = (*G).operatorix_to_numedges[ix];
			
			int num_dims;
			num_dims = (*G).num_dims;
			
			void * params;
			params = (*G).list_of_operators_param_pointers[ix];
			
			// starting position of x for each operator
			int start_edge_ix;
			start_edge_ix = (*G).operatorix_to_edgelistix[ix][0];
			
			curr_opt(&((*G).X[num_dims * start_edge_ix]), &((*G).N[num_dims * start_edge_ix]), &((*G).RHOS[start_edge_ix]), num_vars, num_dims, params);
		}
		
		__syncthreads();
		
		if (ix < ((*G).num_edges) * ((*G).num_dims))
		{
			(*G).M[ix] = (*G).X[ix] + (*G).U[ix];
		}
		
		__syncthreads();
		
		if (ix < (*G).num_vars)
		{
			// computing the denominator
			double denomrhos;
			denomrhos = 0;
			
			for (int i=0; i<(*G).varix_to_numedges[ix]; i++)
			{
				int edge_ix;
				edge_ix = (*G).varix_to_edgelistix[ix][i];
				denomrhos += (*G).RHOS[edge_ix];
			}
			
			
			for (int j=0; j<(*G).num_dims; j++)
			{
				double tempZ;
				tempZ = 0;
				
				for (int i=0; i<(*G).varix_to_numedges[ix]; i++)
				{
					int edge_ix;
					edge_ix = (*G).varix_to_edgelistix[ix][i];
					
					tempZ = tempZ + (*G).M[(*G).num_dims * edge_ix + j] * ((*G).RHOS[edge_ix]);
				}
				
				(*G).Z[(*G).num_dims * ix + j] = tempZ /denomrhos;
			}
		}
		
		__syncthreads();
		
		if (ix < (*G).num_vars)
		{
			for (int j=0; j<(*G).varix_to_numedges[ix]; j++)
			{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[ix][j];
				
				double alpha = (*G).ALPHAS[edge_ix];
				
				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).U[edge_ix * (*G).num_dims + dim_ix] += alpha * ((*G).X[(*G).num_dims * edge_ix + dim_ix] - (*G).Z[(*G).num_dims * ix + dim_ix]);
				}
			}
		}
		
		__syncthreads();
		
		if (ix < (*G).num_vars)
		{
			for (int j=0; j<(*G).varix_to_numedges[ix]; j++)
			{
				int edge_ix;
				edge_ix= (*G).varix_to_edgelistix[ix][j];
				
				for (int dim_ix=0; dim_ix<(*G).num_dims; dim_ix++)
				{
					(*G).N[edge_ix * (*G).num_dims + dim_ix] = (*G).Z[(*G).num_dims * ix + dim_ix] - (*G).U[edge_ix * (*G).num_dims + dim_ix];
				}
			}
		}
		
		__syncthreads();
		
	}
}

__global__ void printKernel(graph* G)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	if(ix == 0)
	{
		printf("Number of operators = %d \n" , (*G).num_opts);
		printf("Number of edges = %d \n" , (*G).num_edges);
		printf("Number of variables = %d \n" , (*G).num_vars);
		printf("Number of dimensions = %d \n" , (*G).num_dims);
		
		printf("X \n");
		for (int i=0 ; i<(*G).num_dims * (*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).X[i]);
		}
		
		printf("M \n");
		for (int i=0; i<(*G).num_dims * (*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).M[i]);
			
		}
		
		printf("N \n");
		for (int i=0; i<(*G).num_dims * (*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).N[i]);
		}
		
		printf("U \n");
		for (int i=0; i<(*G).num_dims * (*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).U[i]);
		}
		
		printf("Z \n");
		for (int i=0; i<(*G).num_dims * (*G).num_vars; i++)
		{
			printf("%lf \n" , (*G).Z[i]);
		}
		
		printf("RHOS \n");
		for (int i=0; i<(*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).RHOS[i]);
		}
		
		printf("Alpha \n");
		for (int i=0; i<(*G).num_edges; i++)
		{
			printf("%lf \n" , (*G).ALPHAS[i]);
		}
		
		printf("Parameter pointers and param Sizes \n");
		for (int i=0; i<(*G).num_opts; i++)
		{
			void* ptr = (*G).list_of_operators_param_pointers[i];
			printf("Param ptr %p of operator %d has size %d \n", ptr, i, (*G).list_of_operators_param_sizes[i] );
		}
		
		printf("Operator Edge Numbers \n");
		for (int i=0; i<(*G).num_opts; i++)
		{
			printf("%d \n" , (*G).operatorix_to_numedges[i]);
		}
		
		printf("Var Edge Numbers \n");
		for (int i=0; i<(*G).num_vars; i++)
		{
			printf("%d \n" , (*G).varix_to_numedges[i]);
		}
		
		for (int j=0; j<(*G).num_opts; j++)
		{
			printf("Parameters of operator %d \n", j);
			
			double* paramstemp = (double*) (*G).list_of_operators_param_pointers[j];
			for (int k=0; k<(*G).list_of_operators_param_sizes[j]/sizeof(double); k++)
			{
				printf("%lf \n ",paramstemp[k]);
			}
		}
		
		for (int j=0; j<(*G).num_opts; j++)
		{
			printf("Edge list of operator %d \n", j);
			
			int* opedgetemp = (int*) (*G).operatorix_to_edgelistix[j];
			for (int k = 0; k < (*G).operatorix_to_numedges[j]; k++)
			{
				printf("%d\n ",opedgetemp[k]);
			}
			
			printf("%%%%");
		}
		
		
		for (int j=0; j<(*G).num_vars; j++)
		{
			printf("Edge list of variables %d \n", j);
			
			int* varedgetemp = (int*) (*G).varix_to_edgelistix[j];
			for (int k=0; k<(*G).varix_to_numedges[j]; k++)
			{
				printf("%d\n ",varedgetemp[k]);
			}
		}
		
		printf("List of edgeix_to_varlistix: \n");
		for (int j=0; j<(*G).num_edges; j++)
		{
			printf("%d ",(*G).edgeix_to_varlistix[j]);
		}
		printf("\n");
		
		double Ntemp[] = { 0.0, 0.0, 0.75, 0.75, 1.0, 0.0, 0.75, 0.75 }; // x1, y1, r1, x2, y2, r2
		double rhostmp[] = { 10000.0, 1.0, 10000.0, 1.0 }; // x1, r1, x2, r2
		
		double Xtemp[] = { -1, -2, -3, -4, -5, -6, -7, -8 };
		
		for (int i=10; i<13; i++)
		{
			void* partmp = (*G).list_of_operators_param_pointers[i];
			(*G).list_of_operators_functions_pointers[i](Xtemp,Ntemp,rhostmp,4,2,partmp);	
			printf("TEST OPT OUTPUT = %f %f %f %f %f %f %f %f\n",Xtemp[0],Xtemp[1],Xtemp[2],Xtemp[3],Xtemp[4],Xtemp[5],Xtemp[6],Xtemp[7]);
		}	
		
		for (int i=0; i<(*G).num_opts; i++) 
		{
			double* a = (double*) (*G).list_of_operators_param_pointers[i];
			printf("Param pointer = %p\n", a);
			
			void* opt_ptr = (void*) (*G).list_of_operators_functions_pointers[i];
			printf("X oprt %d with pointer %p \n", i, opt_ptr );
		}
	}
}

// this function copies a grapm from the CPU to a graph in the GPU where ALL the variables that the GPU graph points to have the correct values
// that are the same as the values that the CPU graph points to. Because of this, we need to first create space in GPU, the copy stuff there from the
// CPU and then make the GPU-graph point to these locations
void copyGraphFromCPUtoGPU(graph c_graph, graph *h_d_graph)
{
	// copy all the integers
	(*h_d_graph).num_opts = c_graph.num_opts;
	(*h_d_graph).num_vars = c_graph.num_vars;
	(*h_d_graph).num_dims = c_graph.num_dims;
	(*h_d_graph).num_edges = c_graph.num_edges;

	// here we copy all arrays of doubles and integers
	cudaMalloc((void**)&(*h_d_graph).X, (c_graph.num_dims * c_graph.num_edges )*sizeof(double));
	cudaMemcpy((*h_d_graph).X, c_graph.X, (c_graph.num_dims * c_graph.num_edges )*sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).M, (c_graph.num_dims * c_graph.num_edges )*sizeof(double));
	cudaMemcpy((*h_d_graph).M, c_graph.M, (c_graph.num_dims * c_graph.num_edges )*sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).N, (c_graph.num_dims * c_graph.num_edges )*sizeof(double));
	cudaMemcpy((*h_d_graph).N, c_graph.N, (c_graph.num_dims * c_graph.num_edges )*sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).U, (c_graph.num_dims * c_graph.num_edges )*sizeof(double));
	cudaMemcpy((*h_d_graph).U, c_graph.U, (c_graph.num_dims * c_graph.num_edges )*sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).Z, (c_graph.num_dims * c_graph.num_vars )*sizeof(double));
	cudaMemcpy((*h_d_graph).Z, c_graph.Z, (c_graph.num_dims * c_graph.num_vars )*sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).RHOS, c_graph.num_edges * sizeof(double));
	cudaMemcpy((*h_d_graph).RHOS, c_graph.RHOS, c_graph.num_edges * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).ALPHAS, c_graph.num_edges * sizeof(double));
	cudaMemcpy((*h_d_graph).ALPHAS, c_graph.ALPHAS,  c_graph.num_edges * sizeof(double), cudaMemcpyHostToDevice);


	cudaMalloc((void**)&(*h_d_graph).operatorix_to_numedges, c_graph.num_opts * sizeof(int));
	cudaMemcpy((*h_d_graph).operatorix_to_numedges, c_graph.operatorix_to_numedges,   c_graph.num_opts * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(*h_d_graph).varix_to_numedges, c_graph.num_vars * sizeof(int));
	cudaMemcpy((*h_d_graph).varix_to_numedges, c_graph.varix_to_numedges,  c_graph.num_vars * sizeof(int), cudaMemcpyHostToDevice);


	cudaMalloc((void**)&(*h_d_graph).list_of_operators_param_sizes, c_graph.num_opts * sizeof(int));
	cudaMemcpy((*h_d_graph).list_of_operators_param_sizes, c_graph.list_of_operators_param_sizes,  c_graph.num_opts * sizeof(int), cudaMemcpyHostToDevice);

	// here we copy the list of operator pointers
	void (** d_list_of_operators_functions_pointers)(double *, double *, double *, int, int, void *);
	cudaMalloc((void**)&d_list_of_operators_functions_pointers, c_graph.num_opts * sizeof(void (*)(double *, double *, double *, int, int, void *)));

	cudaMemcpy(d_list_of_operators_functions_pointers, c_graph.list_of_operators_functions_pointers, c_graph.num_opts * sizeof(void (*)(double *, double *, double *, int, int, void *)), cudaMemcpyHostToDevice);

	(*h_d_graph).list_of_operators_functions_pointers = d_list_of_operators_functions_pointers;


	// here we copy all the double arrays. we create new GPU-param pointers in the GPU and copy whatever the CPU-params pointers are pointing to also the GPU
	void** d_list_of_operators_param_pointers;
	cudaMalloc((void**)&d_list_of_operators_param_pointers, c_graph.num_opts*sizeof(void*) );

	(*h_d_graph).list_of_operators_param_pointers = (void**) malloc( c_graph.num_opts * sizeof(void *));

	for (int i=0; i<c_graph.num_opts; i++)
	{
		void* d_paramtemp;
		cudaMalloc((void **)&d_paramtemp, c_graph.list_of_operators_param_sizes[i] );
		cudaMemcpy(d_paramtemp,c_graph.list_of_operators_param_pointers[i], c_graph.list_of_operators_param_sizes[i] , cudaMemcpyHostToDevice);
		(*h_d_graph).list_of_operators_param_pointers[i] = d_paramtemp;
	}

	cudaMemcpy(d_list_of_operators_param_pointers,(*h_d_graph).list_of_operators_param_pointers,c_graph.num_opts * sizeof(void *), cudaMemcpyHostToDevice);

	free((*h_d_graph).list_of_operators_param_pointers); //it is very important to free all the memory that we allocated using malloc
	(*h_d_graph).list_of_operators_param_pointers = d_list_of_operators_param_pointers;

	int** d_operatorix_to_edgelistix;
	cudaMalloc((int**)&d_operatorix_to_edgelistix, c_graph.num_opts *sizeof(int*) );

	(*h_d_graph).operatorix_to_edgelistix = (int**) malloc( c_graph.num_opts * sizeof(int*));

	for (int i=0; i<c_graph.num_opts; i++)
	{
		int* d_opedgelisttemp;
		cudaMalloc((int **)&d_opedgelisttemp, c_graph.operatorix_to_numedges[i]*sizeof(int) );
		cudaMemcpy(d_opedgelisttemp,c_graph.operatorix_to_edgelistix[i], c_graph.operatorix_to_numedges[i]*sizeof(int) , cudaMemcpyHostToDevice);
		(*h_d_graph).operatorix_to_edgelistix[i] = d_opedgelisttemp;
	}

	cudaMemcpy(d_operatorix_to_edgelistix,(*h_d_graph).operatorix_to_edgelistix,c_graph.num_opts * sizeof(int *), cudaMemcpyHostToDevice);

	free((*h_d_graph).operatorix_to_edgelistix); //it is very important to free all the memory that we allocated using malloc
	(*h_d_graph).operatorix_to_edgelistix = d_operatorix_to_edgelistix;

	int** d_varix_to_edgelistix;
	cudaMalloc((int**)&d_varix_to_edgelistix, c_graph.num_vars *sizeof(int*));

	(*h_d_graph).varix_to_edgelistix = (int**) malloc( c_graph.num_vars * sizeof(int*));

	for (int i=0; i<c_graph.num_vars; i++)
	{
		int* d_varedgelisttemp;
		cudaMalloc((int**)&d_varedgelisttemp, c_graph.varix_to_numedges[i]*sizeof(int));
		cudaMemcpy(d_varedgelisttemp,c_graph.varix_to_edgelistix[i], c_graph.varix_to_numedges[i]*sizeof(int), cudaMemcpyHostToDevice);
		(*h_d_graph).varix_to_edgelistix[i] = d_varedgelisttemp;
	}

	cudaMemcpy(d_varix_to_edgelistix,(*h_d_graph).varix_to_edgelistix,c_graph.num_vars * sizeof(int *), cudaMemcpyHostToDevice);

	free((*h_d_graph).varix_to_edgelistix); //it is very important to free all the memory that we allocated using malloc
	(*h_d_graph).varix_to_edgelistix = d_varix_to_edgelistix;

	// copy the list of edge_ix to varlistix
	cudaMalloc((void**)&(*h_d_graph).edgeix_to_varlistix, c_graph.num_edges * sizeof(int));
	cudaMemcpy((*h_d_graph).edgeix_to_varlistix, c_graph.edgeix_to_varlistix,   c_graph.num_edges * sizeof(int), cudaMemcpyHostToDevice);
}
