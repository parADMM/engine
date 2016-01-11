#include <stdio.h>
#include "paradmm-gpu.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

__device__ double d_max_alpha_radius(double n, double rho, double alpha)
{
	const int newton_iterations = 10;
	double r = d_abs_d( n / 2.0 );

	for ( int i=0; i<newton_iterations; i++ )
	{
		const double num = -pow( r, alpha-1 ) + ( rho / alpha )*( r - n );
		const double denom = -( alpha - 1 ) * pow( r, alpha-2 ) + ( rho / alpha );

		r -= num / denom;

		if ( r < 0 )
		{
			return 0.0;
		}
	}

	return r;
}

__device__ void d_maximize_sum_of_squared_disk_radii_alpha(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
{
	double* paramsCasted = (double*) params;

	for ( int i=0; i<numvars; i++ )
	{
		double rho_i = rhos[i];
		double rad_i = d_max_alpha_radius(n[ i*numdims ], rho_i, paramsCasted[0]);

		for ( int j=0; j<numdims; j++ )
		{
			x[ i*numdims + j ] = rad_i;
		}
	}
}
__device__ void (*d_ptr_maximize_sum_of_squared_disk_radii_alpha)(double *x, double *n, double *rhos, int numvars, int numdims, void *params) = d_maximize_sum_of_squared_disk_radii_alpha;

__device__ void d_disk_plane_collision_variable_radius(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
{
	double* paramsCasted = (double*) params;

	double* v = paramsCasted;
	double* p = ( paramsCasted + numdims );
	double inner = 0;

	for ( int i=0; i<numdims; i++ )
	{
		inner += v[i]*( n[i] - p[i] );
	}

	if ( inner + n[numdims] > 0 )
	{
		double r = ( -inner*rhos[0] + n[numdims]*rhos[1] ) / ( rhos[0] + rhos[1] );
		if ( r < 0 )
		{
			r = 0;
		}

		for ( int i=0; i<numdims; i++ )
		{
			x[i] = n[i] - ( inner*v[i] ) + ( -v[i]*r );
			x[ i+numdims ] = r;
		}
	}
	else
	{
		for ( int i=0; i<numdims*numvars; i++ )
		{
			x[i] = n[i];
		}
	}
}
__device__ void (*d_ptr_disk_plane_collision_variable_radius)(double *x, double *n, double *rhos, int numvars, int numdims, void *params) = d_disk_plane_collision_variable_radius;

__device__ void d_disk_collision_variable_radius(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
{
	double nr1 = n[ numdims ];
	double nr2 = n[ 3*numdims ];

	double sumOfRadiiSquared = ( nr1 + nr2 );
	sumOfRadiiSquared *= sumOfRadiiSquared;

	double distanceSquared = 0.0;
	for ( int i=0; i<numdims; i++ )
	{
		double diff = ( n[i] - n[ i + 2*numdims ] );
		distanceSquared += ( diff*diff );
	}

	for ( int i=0; i<numdims*numvars; i++ )
	{
		x[i] = n[i];
	}

	if ( distanceSquared < sumOfRadiiSquared )
	{
		double nx1 = 0;
		double nx2 = sqrt( distanceSquared );

		double rx1 = rhos[0];
		double rr1 = rhos[1];
		double rx2 = rhos[2];
		double rr2 = rhos[3 ];

		double denom = (rr1*rr2*rx1 + rr2*rx1*rx2 + rr1*(rr2 + rx1)*rx2);

		double x1lin = ((-nr1 - nr2 + nx2)*rr1*rr2*rx2 + nx1*rx1*(rr2*rx2 + rr1*(rr2 + rx2)))/denom;
		double x2lin = ((nr1 + nr2 + nx1)*rr1*rr2*rx1 + nx2*(rr2*rx1 + rr1*(rr2 + rx1))*rx2)/denom;
		double r1lin = ((-nr2 - nx1 + nx2)*rr2*rx1*rx2 + nr1*rr1*(rx1*rx2 + rr2*(rx1 + rx2)))/denom;
		double r2lin = ((-nr1 - nx1 + nx2)*rr1*rx1*rx2 + nr2*rr2*(rx1*rx2 + rr1*(rx1 + rx2)))/denom;

		for ( int i=0; i<numdims; i++ )
		{
			double directi = ( n[ i + 2*numdims ] - n[ i ] ) / ( nx2 );
			x[ i + numdims ] = r1lin;
			x[ i + 3*numdims ] = r2lin;

			x[i] += directi*x1lin;
			x[ i+2*numdims ] += directi * ( x2lin - nx2 );
		}

	}

}
__device__ void (*d_ptr_disk_collision_variable_radius)(double *x, double *n, double *rhos, int numvars, int numdims, void *params) = d_disk_collision_variable_radius;

void circle_packing_triangle_test(int num_balls, int num_dims, int numiterations, int factor_X, int factor_M, int factor_Z, int factor_U, int factor_N)
{
	cudaSetDevice(0);
	cudaDeviceReset(); //it is important to make this to free resources from the GPU

	//set the values of dimensions, operators, variables, etc
	int num_walls = 3;

	int numedges = 4*0.5*num_balls*(num_balls-1) + 1*num_balls + 2*num_balls*num_walls;
	int numvars = 2*num_balls;
	int numproxop = 0.5*num_balls*(num_balls-1) + num_balls + num_balls*num_walls;

	dim3 numthreads_X(1, 1024/factor_X); // we arrange threads in 2d arrays
	dim3 numblocks_X(  factor_X , ceil(((double) numproxop)/((double) 1024.0))  );// we arrange blocks in 2d arrays

	dim3 numthreads_M(1, 1024/factor_M); // we arrange threads in 2d arrays
	dim3 numblocks_M(num_dims*factor_M,    ceil(((double) numedges)/((double) 1024.0))       );// we arrange blocks in 2d arrays

	dim3 numthreads_Z(1, 1024/factor_Z); // we arrange threads in 2d arrays
	dim3 numblocks_Z(factor_Z,  ceil(((double) numvars)/((double) 1024.0))       );// we arrange blocks in 2d arrays

	dim3 numthreads_U(1, 1024/factor_U); // we arrange threads in 2d arrays
	dim3 numblocks_U(factor_U,  ceil(((double) numedges)/((double) 1024.0))   );// we arrange blocks in 2d arrays

	dim3 numthreads_N(1, 1024/factor_N); // we arrange threads in 2d arrays
	dim3 numblocks_N(factor_N, ceil(((double) numedges)/((double) 1024.0))     );// we arrange blocks in 2d arrays

	double theta1 = 60.0*PI/180.;
	double theta2 = 60.0*PI/180.;

	// create the walls in the CPU
	double walls[3][4] = { // num_walls, 2*num_dims = direction, center
		{ 0., -1., 0., 0. },
		{ -sin(theta1), cos(theta1), 0., 0. },
		{ sin(theta2), cos(theta2), 1., 0. }
	};

	float time;
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0); //Time to form the CPU graph

	void (* opt1)(double *, double *, double *, int, int, void *);
	cudaMemcpyFromSymbol(&opt1, d_ptr_disk_plane_collision_variable_radius,sizeof(opt1));

	void (* opt2)(double *, double *, double *, int, int, void *);
	cudaMemcpyFromSymbol(&opt2, d_ptr_maximize_sum_of_squared_disk_radii_alpha,sizeof(opt2));

	void (* opt3)(double *, double *, double *, int, int, void *);
	cudaMemcpyFromSymbol(&opt3, d_ptr_disk_collision_variable_radius,sizeof(opt3));

	int opt1_num_vars = 2;
	int opt2_num_vars = 1;
	int opt3_num_vars = 4;

	graph c_graph;
	startG(&c_graph, num_dims);

	int opt1_listvarix[] = { 0, 0 };

	// wall proximal operators
	for ( int i=0; i<num_walls; i++ )
	{
		for ( int j=0; j<num_balls; j++ )
		{
			opt1_listvarix[0] = 2*j;
			opt1_listvarix[1] = ( 2*j + 1 );

			addNodeWithParamSize(&c_graph, opt1, (void *) walls[i], 4*sizeof(double), opt1_num_vars, opt1_listvarix);
		}
	}

	// maximum radius proximal operator
	double alpha = 2.0;

	int* opt2_listvarix = (int*) malloc(num_balls * sizeof(int));
	for ( int i=0; i<num_balls; i++ )
	{
		opt2_listvarix[i] = ( 2*i + 1 );
		addNodeWithParamSize(&c_graph, opt2, (void *) &alpha, 1*sizeof(double), opt2_num_vars, &opt2_listvarix[i]);
	}

	// ball collision proximal operators
	int opt3_listvarix[] = { 0, 0, 0, 0 };
	for (int i = 0; i < (num_balls -1); i++)
	{
		opt3_listvarix[0] = ( 2*i );
		opt3_listvarix[1] = ( 2*i + 1 );

		for (int j = i + 1; j < num_balls; j++)
		{
			opt3_listvarix[2] = ( 2*j );
			opt3_listvarix[3] = ( 2*j + 1 );

			addNodeWithParamSize(&c_graph, opt3, (void *) NULL, 0*sizeof(double), opt3_num_vars, opt3_listvarix);
		}
	}

	initialize_RHOS_APHAS(&c_graph, 10.0, 0.1);

	initialize_X_N_Z_M_U_rand(&c_graph, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
	for ( int i=0; i<num_balls; i++ )
	{
		c_graph.RHOS[i] = 5.0;
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	if (!SILENCE) printf("Time elapsed to form c_graph: %f ms \n",time);

	cudaEventRecord(start,0); //Time to copy CPU graph to the GPU

	graph h_d_graph ;

	copyGraphFromCPUtoGPU(c_graph, &h_d_graph);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	if (!SILENCE) printf("Time elapsed to copyGraphFromCPUtoGPU: %f ms \n",time);

	graph *d_graph;
	cudaMalloc((void **)&d_graph, sizeof(graph));

	cudaEventRecord(start,0);

	cudaMemcpy(d_graph, &h_d_graph, sizeof(graph), cudaMemcpyHostToDevice);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	if (!SILENCE) printf("Time elapsed to copy from h_d_graph to d_graph: %f ms \n",time);

	cudaEventRecord(start,0);

	for (int i = 0; i< numiterations; i++)
	{
		d_updateX<<<numblocks_X,numthreads_X>>>(d_graph);
		d_updateM<<<numblocks_M,numthreads_M>>>(d_graph);
		d_updateZ<<<numblocks_Z,numthreads_Z>>>(d_graph);
		d_updateU_per_edge<<<numblocks_U,numthreads_U>>>(d_graph);
		d_updateN_per_edge<<<numblocks_N,numthreads_N>>>(d_graph);
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	if (!SILENCE) printf("Time elapsed to run %d iterations: %f ms \n",numiterations,time);

	cudaError_t cudaerr = cudaGetLastError();
	if (!SILENCE) printf("Last error message \"%s\".\n",cudaGetErrorString(cudaerr));

	cudaDeviceSynchronize();

	cudaEventRecord(start,0); //Time to Copy Z from GPU to CPU
	cudaMemcpy(c_graph.Z, ((h_d_graph).Z), c_graph.num_dims * c_graph.num_vars * sizeof(double),cudaMemcpyDeviceToHost);
	cudaEventRecord(stop,0);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	if (!SILENCE) printf("Time elapsed to copy Z from GPU to CPUs: %f ms \n",time);

	fflush(stdout);
	cudaFree(d_graph); //this needs to be done properly because right now we're not freeing all the memory in the GPU. We need a "closegraph" functions
}

int main()
{
	// NOT WAITING FOR CONVERGENCE - TO DO SO, WAIT MANY ITERATIONS OR START NEAR SOLUTION VALUES
	// 5000 circles, 2 dimensions, 100 iterations
	// 1024 / 32 threads-per-thread-block for X-update
	// 1024 / 32 threads-per-thread-block for M-update
	// 1024 / 64 threads-per-thread-block for Z-update
	// 1024 / 32 threads-per-thread-block for U-update
	// 1024 / 32 threads-per-thread-block for N-update
	// ~ 21 sec on Tesla K40
	// For 5000 circles, 100 iterations: CPU code on a single core takes ~339 sec
	// Factor Graph -> GPU copy time ~ 7 minutes, which is negligable compared to number of iterations for convergence
	circle_packing_triangle_test(5000, 2, 100, 32, 32, 64, 32, 32);
	return EXIT_SUCCESS;
}
