#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "paradmm.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

double max_alpha_radius(double n, double rho, double alpha)
{
    const int newton_iterations = 10;
    double r = abs_d( n / 2.0 );

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

void disk_plane_collision_variable_radius(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
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

void maximize_sum_of_squared_disk_radii_alpha(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
{
    double* paramsCasted = (double*) params;

    for ( int i=0; i<numvars; i++ )
    {
        double rho_i = rhos[i];
        double rad_i = max_alpha_radius(n[ i*numdims ], rho_i, paramsCasted[0]);

        for ( int j=0; j<numdims; j++ )
        {
            x[ i*numdims + j ] = rad_i;
        }
    }
}

void disk_collision_variable_radius(double *x, double *n, double *rhos, int numvars, int numdims, void *params)
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
        double rr2 = rhos[3];

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

#ifdef OPENMP
#include <omp.h>

void test_circlepacking_variableradius_triangle_auto_with_OpenMP(int num_balls, int maxiter, int nu_threads)
{
	int num_dims = 2;

	double theta1 = 60.0*PI/180.;
	double theta2 = 60.0*PI/180.;

	int num_walls = 3;
	double walls[3][4] = { // num_walls, 2*num_dims = direction, center
		{ 0., -1., 0., 0. },
		{ -sin(theta1), cos(theta1), 0., 0. },
		{ sin(theta2), cos(theta2), 1., 0. }
	};


	double** wallptr = malloc( num_walls*sizeof(double*) );
	for ( int i=0; i<num_walls; i++ )
	{
		wallptr[i] = walls[i];
	}

	// Initialize graph G
	graph G;
	startG(&G, num_dims);

	void* opt1 = &disk_plane_collision_variable_radius;
	void* opt2 = &maximize_sum_of_squared_disk_radii_alpha;
	void* opt3 = &disk_collision_variable_radius;

	int opt1_num_vars = 2;
	int opt2_num_vars = num_balls;
	int opt3_num_vars = 4;

	double *opt2_params = (double*) malloc( num_balls * sizeof(double) );
	for ( int i=0; i<num_balls; i++ )
	{
		opt2_params[i] = 0.5;
	}
	opt2_params[0] = 4;
	opt2_params[0] = 2.0; // for alpha

	int opt1_listvarix[] = { 0, 0 };
	int* opt2_listvarix = (int*) malloc( num_balls * sizeof(int));
	for ( int i=0; i<num_balls; i++ )
	{
		opt2_listvarix[i] = ( 2*i + 1 );
		addNode(&G, (void *) opt2, (void *) opt2_params, 1, &opt2_listvarix[i]);
	}

	// wall proximal operators
	int opt3_listvarix[] = { 0, 0, 0, 0 };
	for ( int i=0; i<num_walls; i++ )
	{
		for ( int j=0; j<num_balls; j++ )
		{
			opt1_listvarix[0] = 2*j;
			opt1_listvarix[1] = ( 2*j+1 );

			addNode(&G, (void *) opt1, (void *) walls[i], opt1_num_vars, opt1_listvarix);
		}
	}

	// ball collision proximal operators
	for (int i = 0; i < (num_balls -1); i++)
	{

		opt3_listvarix[0] = ( 2*i );
		opt3_listvarix[1] = ( 2*i + 1 );

		for (int j = i + 1; j < num_balls; j++)
		{

			opt3_listvarix[2] = ( 2*j );
			opt3_listvarix[3] = ( 2*j + 1 );

			addNode(&G, (void *) opt3, (void *) NULL, opt3_num_vars, opt3_listvarix);

		}

	}

	initialize_X_N_Z_M_U_rand(&G, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2);
	initialize_RHOS_APHAS(&G, 10.0, 0.1); // rho=1.5

	for ( int i=0; i<num_balls; i++ )
	{
		G.RHOS[i] = 5.0;
	}


	 	double start;
	 	double end;

	 	omp_set_num_threads(nu_threads);

	struct timeval startsystime, endsystime;
	 	double delta;

	gettimeofday(&startsystime, NULL);

	for(int i = 0; i < maxiter; i++)
	{
		updateXOpenMP(&G);
		updateMOpenMP(&G);
		updateZOpenMP(&G);
		updateUOpenMP(&G);
		updateNOpenMP(&G);
	}

	gettimeofday(&endsystime, NULL);
	delta = ((endsystime.tv_sec  - startsystime.tv_sec) * 1000000u + endsystime.tv_usec - startsystime.tv_usec) / 1.e6;
	if (!SILENCE) printf("total = %f\n",delta);

	free( opt2_params );
	free( opt2_listvarix );
	free(wallptr);
}
#endif

void test_circlepacking_variableradius_triangle_auto(int num_balls, int maxiter)
{
	int num_dims = 2;

	double theta1 = 60.0*PI/180.;
	double theta2 = 60.0*PI/180.;

	int num_walls = 3;
	double walls[3][4] = { // num_walls, 2*num_dims = direction, center
		{ 0., -1., 0., 0. },
		{ -sin(theta1), cos(theta1), 0., 0. },
		{ sin(theta2), cos(theta2), 1., 0. }
	};


	double** wallptr = malloc( num_walls*sizeof(double*) );
	for ( int i=0; i<num_walls; i++ )
	{
		wallptr[i] = walls[i];
	}

	// Initialize graph G
	graph G;
	startG(&G, num_dims);

	void* opt1 = &disk_plane_collision_variable_radius;
	void* opt2 = &maximize_sum_of_squared_disk_radii_alpha;
	void* opt3 = &disk_collision_variable_radius;

	int opt1_num_vars = 2;
	int opt2_num_vars = num_balls;
	int opt3_num_vars = 4;

	double *opt2_params = (double*) malloc( num_balls * sizeof(double) );
	for ( int i=0; i<num_balls; i++ )
	{
		opt2_params[i] = 0.5;
	}
	opt2_params[0] = 4;
	opt2_params[0] = 2.0; // for alpha

	int opt1_listvarix[] = { 0, 0 };
	int* opt2_listvarix = (int*) malloc( num_balls * sizeof(int));
	for ( int i=0; i<num_balls; i++ )
	{
		opt2_listvarix[i] = ( 2*i + 1 );
		addNode(&G, (void *) opt2, (void *) opt2_params, 1, &opt2_listvarix[i]);
	}

	// wall proximal operators
	int opt3_listvarix[] = { 0, 0, 0, 0 };
	for ( int i=0; i<num_walls; i++ )
	{
		for ( int j=0; j<num_balls; j++ )
		{
			opt1_listvarix[0] = 2*j;
			opt1_listvarix[1] = ( 2*j+1 );

			addNode(&G, (void *) opt1, (void *) walls[i], opt1_num_vars, opt1_listvarix);
		}
	}

	// ball collision proximal operators
	for (int i = 0; i < (num_balls -1); i++)
	{

		opt3_listvarix[0] = ( 2*i );
		opt3_listvarix[1] = ( 2*i + 1 );

		for (int j = i + 1; j < num_balls; j++)
		{

			opt3_listvarix[2] = ( 2*j );
			opt3_listvarix[3] = ( 2*j + 1 );

			addNode(&G, (void *) opt3, (void *) NULL, opt3_num_vars, opt3_listvarix);

		}

	}

	initialize_X_N_Z_M_U_rand(&G, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2);
	initialize_RHOS_APHAS(&G, 10.0, 0.1); // rho=1.5

	for ( int i=0; i<num_balls; i++ )
	{
		G.RHOS[i] = 5.0;
	}


	 	clock_t start;
	 	clock_t startclock;
	 	clock_t end;
	clock_t endclock;

 	double total;

	startclock = clock();

	for(int i = 0; i < maxiter; i++)
	{
		updateXM(&G);
		updateZ(&G);
		updateUN(&G);
	}

	endclock = clock();
	if (!SILENCE) printf("total = %f\n",((float)(endclock - startclock) / CLOCKS_PER_SEC));

	free( opt2_params );
	free( opt2_listvarix );
	free(wallptr);
}

int main()
{
   // 2500 circles, 100 iterations ~ 81 sec on 2.8GHz AMD CPU
	test_circlepacking_variableradius_triangle_auto(2500, 100);
#ifdef OPENMP
   // 2500 circles, 10k iterations, 32 cores (higher number of iterations to mask OpenMP overhead) ~ 13 sec on 32-cores of 2.8GHz AMD CPU
	test_circlepacking_variableradius_triangle_auto_with_OpenMP(2500, 10000, 32);
#endif
}
