
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Domain size.
const double A = 4.0;
const double B = 4.0;
const double EPSILON = 0.00001;

int NX, NY;							        // the number of internal points on axes (ox) and (oy).
double hx, hy;                              // mesh steps on (0x) and (0y) axes

// Buffers for exchanging with other processors
double * sendLeftBuf,  * recLeftBuf;
double * sendRightBuf, * recRightBuf;
double * sendUpBuf,    * recUpBuf;
double * sendDownBuf,  * recDownBuf;

// ----------------------
#define TRUE  ((int) 1)
#define FALSE ((int) 0)
// ----------------------

//#define Test
//#define Print
#define Step 10

#define Max(A,B) ((A)>(B)?(A):(B))
#define R2(x,y) ((x)*(x)+(y)*(y))
#define Cube(x) ((x)*(x)*(x))



#define x(i) ((i)*hx - 2.0)
#define y(j) ((j)*hy - 2.0)

#define LeftPart(P,i,j)\
((-(P[NX*(j)+i+1]-P[NX*(j)+i])/hx+(P[NX*(j)+i]-P[NX*(j)+i-1])/hx)/hx+\
 (-(P[NX*(j+1)+i]-P[NX*(j)+i])/hy+(P[NX*(j)+i]-P[NX*(j-1)+i])/hy)/hy)

int IsPower(int Number)
{
    unsigned int M;
    int p;

    if(Number <= 0)
        return(-1);

    M = Number; p = 0;
    while(M % 2 == 0)
    {
        ++p;
        M = M >> 1;
    }
    if((M >> 1) != 0)
        return(-1);
    else
        return(p);

}

int SplitFunction(int NX, int NY, int p)
{
    float n0, n1;
    int p0, i;

    n0 = (float) NX; n1 = (float) NY;
    p0 = 0;

    for(i = 0; i < p; i++)
        if(n0 > n1)
        {
            n0 = n0 / 2.0;
            ++p0;
        }
        else
            n1 = n1 / 2.0;

    return(p0);
}

double Solution(double x,double y)
{    
    return exp(1-(x+y)*(x+y));
}


double BoundaryValue(double x, double y)
{
	return Solution(x,y);
}

int RightPart(double * rhs)
{
    int i, j;
    //double kappa2 = (16.0/R2(A,B));

    #ifdef Test
        for(j=0; j<NY; j++)
            for(i=0; i<NX; i++)
                rhs[j*NX+i] = 4.0 * ( 1 - 2.0 * (x(i)+y(j))*(x(i)+y(j)) )* exp( 1- (x(i)+y(j))*(x(i)+y(j)) ); 
                //rhs[j*NX+i] = 4.0*sqrt(R2(A,B))*kappa2*(1.0-kappa2*R2(x(i)-A/2,y(j)-B/2))\
                //           /Cube(1.0 + kappa2*R2(x(i)-A/2,y(j)-B/2));
        return 0;
    #else
        memset(rhs,0,NX*NY*sizeof(double));
    
    
        return 0;
    #endif
}

void down_data_exchange(int down, int first_index, int last_index, int n0, int grid_coord_down, double *data) {
    int size = last_index - first_index + 1;
    double *recv_solv_vect_buf = (double *)malloc(size*sizeof(double));
    double *send_solv_vect_buf = (double *)malloc(size*sizeof(double));
    MPI_Request recv_request[1];
    MPI_Irecv(recv_solv_vect_buf, size, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, recv_request);
    int index;
    for (index = first_index; index <= last_index; index++) 
    {
        send_solv_vect_buf[index - first_index] = data[n0 * grid_coord_down + index];
    }
    MPI_Request send_request[1];
    MPI_Isend(send_solv_vect_buf, size, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, send_request);
    MPI_Waitall(1, recv_request, MPI_STATUS_IGNORE);
    for (index = first_index; index <= last_index; index++) 
    {
        data[n0 * (grid_coord_down - 1) + index] = recv_solv_vect_buf[index - first_index];
    }
    free(send_solv_vect_buf);
    free(recv_solv_vect_buf);
}

void left_data_exchange(int left, int first_index, int last_index, int n0, int grid_coord_left, double *data) {
    int size = last_index - first_index + 1;
    double *recv_solv_vect_buf = (double *)malloc(size*sizeof(double));
    double *send_solv_vect_buf = (double *)malloc(size*sizeof(double));
    MPI_Request recv_request[1];
    MPI_Irecv(recv_solv_vect_buf, size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, recv_request);
    int index;
    for (index = first_index; index <= last_index; index++) 
    {
        send_solv_vect_buf[index - first_index] = data[n0 * index + grid_coord_left];
    }
    MPI_Request send_request[1];
    MPI_Isend(send_solv_vect_buf, size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, send_request);
    MPI_Waitall(1, recv_request, MPI_STATUS_IGNORE);
    for (index = first_index; index <= last_index; index++) 
    {
        data[n0 * index + grid_coord_left - 1] = recv_solv_vect_buf[index - first_index];
    }
    free(send_solv_vect_buf);
    free(recv_solv_vect_buf);
}

void right_data_exchange(int right, int first_index, int last_index, int n0, int grid_coord_right, double *data) {
    int size = last_index - first_index + 1;
    double *recv_solv_vect_buf = (double *)malloc(size*sizeof(double));
    double *send_solv_vect_buf = (double *)malloc(size*sizeof(double));
    MPI_Request recv_request[1];
    MPI_Irecv(recv_solv_vect_buf, size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, recv_request);
    int index;
    for (index = first_index; index <= last_index; index++) 
    {
        send_solv_vect_buf[index - first_index] = data[n0 * index + grid_coord_right];
    }
    MPI_Request send_request[1];
    MPI_Isend(send_solv_vect_buf, size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, send_request);
    MPI_Waitall(1, recv_request, MPI_STATUS_IGNORE);
    for (index = first_index; index <= last_index; index++) 
    {
        data[n0 * index + grid_coord_right + 1] = recv_solv_vect_buf[index - first_index];
    }
    free(send_solv_vect_buf);
    free(recv_solv_vect_buf);
}

void up_data_exchange(int up, int first_index, int last_index, int n0, int grid_coord_up, double *data) {
    int size = last_index - first_index + 1;
    double *recv_solv_vect_buf = (double *)malloc(size*sizeof(double));
    double *send_solv_vect_buf = (double *)malloc(size*sizeof(double));
    MPI_Request recv_request[1];
    MPI_Irecv(recv_solv_vect_buf, size, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, recv_request);
    int index;
    for (index = first_index; index <= last_index; index++) 
    {
        send_solv_vect_buf[index - first_index] = data[n0 * grid_coord_up + index];
    }
    MPI_Request send_request[1];
    MPI_Isend(send_solv_vect_buf, size, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, send_request);
    MPI_Waitall(1, recv_request, MPI_STATUS_IGNORE);
    for (index = first_index; index <= last_index; index++)
    {
        data[n0 * (grid_coord_up + 1) + index] = recv_solv_vect_buf[index - first_index];
    }
    free(send_solv_vect_buf);
    free(recv_solv_vect_buf);
}



void boundaryExchange(double * vect, int * Coords, int * dims, int n0, int n1, int x1, int x2, int y1, int y2, MPI_Comm Grid_Comm, int left, int right, int up, int down)
{
    MPI_Request requests[8];
    MPI_Status statuses[2];

    if (Coords[0] > 0) 
    {
        for (int j = 0; j<n1; j++) 
        {
            left_data_exchange(left, y2, y1, n0, x1, vect);
        }
    }

    if (Coords[0] < dims[0]-1) 
    {
        for (int j = 0; j<n1; j++)
        {
            right_data_exchange(right, y2, y1, n0, x2, vect);
        }
    }
    if (Coords[1] > 0) 
    {
        for (int i = 0; i<n0; i++) 
        {
            up_data_exchange(up, x1, x2, n0, y1, vect);
        }
    }
    if (Coords[1] < dims[1]-1) 
    {
        for (int i = 0; i<n0; i++) 
        {
            down_data_exchange(down, x1, x2, n0, y2, vect);
        }
    }
}


int main(int argc, char * argv[])
{
    int NX, NY;                     // Mesh has NX x NY nodes.
    int ProcNum, rank, rank1;   // the number of processes and rank in communicator.
    int power, p0, p1;              // ProcNum = 2^(power), power splits into sum p0 + p1.
    int dims[2];                    // dims[0] = 2^p0, dims[1] = 2^p1 (--> M = dims[0]*dims[1]).
    int n0,n1, k0,k1;               // NX = n0*dims[0] + k0, NY = n1*dims[1] + k1.
    int Coords[2];                  // the process coordinates in the cartesian topology created for mesh.
    int x1, x2, y1, y2;             // coords of the field for which this processor is responsible for
    int fict;
    int k, l, curY;

    MPI_Comm Grid_Comm;             // this is a handler of a new communicator.
    MPI_Status status;
    const int ndims = 2;            // the number of a process topology dimensions.
    int periods[2] = {0,0};         // it is used for creating processes topology.
    int Coords2[2];
    int left, right, up, down;      // the neighbours of the process.

	double * SolVect;					// the solution array.
	double * ResVect;					// the residual array.
	double * BasisVect;					// the vector of A-orthogonal system in CGM.
	double * RHS_Vect;					// the right hand side of Puasson equation.
	double sp, sumSp, alpha, sumAlpha, tau, sumTau, NewValue, err;			// auxiliary values.
	int SDINum, CGMNum;					// the number of steep descent and CGM iterations.
	int iteration;                        // the current iteration number.
	double s1;
    double ex_time;
    double t_start, t_finish;
    double sumErr;
    double * bufD;
    int bufI[3];

	int i,j;
	char str[127];
	FILE * fp;

	switch (argc)
	{
	case 4:{
				SDINum = 1;
				CGMNum = atoi(argv[3]);
				break;
		   }
	case 5:{
				SDINum = Max(atoi(argv[3]),1);		// SDINum >= 1
				CGMNum = atoi(argv[4]);
				break;
		   }
	default:{
				printf("Wrong number of parameters in command line.\nUsage: <ProgName> "
                       "<Nodes number on (0x) axis> <Nodes number on (0y) axis> "
                       "[the number of steep descent iterations] "
                       "<the number of conjugate gragient iterations>\nFinishing...\n");
				return(-1);
			}
	}


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank == 0) 
    {
        t_start =  MPI_Wtime();
    }

	NX = atoi(argv[1]); NY = atoi(argv[2]);
    hx = A / (NX-1);    hy = B / (NY-1);

	if((power = IsPower(ProcNum)) < 0)
    {
        if(rank == 0)
            printf("The number of procs must be a power of 2.\n");
        MPI_Finalize();
        return(3);
    }
    p0 = SplitFunction(NX, NY, power);
    p1 = power - p0;
    dims[0] = (unsigned int) 1 << p0;   dims[1] = (unsigned int) 1 << p1;
    n0 = NX >> p0;                      n1 = NY >> p1;
    k0 = NX - dims[0]*n0;               k1 = NY - dims[1]*n1;
    if(rank == 0)
    {
        printf("The number of processes ProcNum = 2^%d. It is split into %d x %d processes.\n"
               "The number of nodes NX = %d, NY = %d. Blocks B(i,j) have size:\n", power, dims[0],dims[1], NX,NY);
        if((k0 > 0)&&(k1 > 0))
            printf("-->\t %d x %d iff i = 0 .. %d, j = 0 .. %d;\n", n0+1,n1+1, k0-1,k1-1);
        if(k1 > 0)
            printf("-->\t %d x %d iff i = %d .. %d, j = 0 .. %d;\n", n0,n1+1, k0,dims[0]-1, k1-1);
        if(k0 > 0)
            printf("-->\t %d x %d iff i = 0 .. %d, j = %d .. %d;\n", n0+1,n1, k0-1, k1,dims[1]-1);
        printf("-->\t %d x %d iff i = %d .. %d, j = %d .. %d.\n", n0,n1, k0,dims[0]-1, k1,dims[1]-1);
    }
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, TRUE, &Grid_Comm);
    MPI_Comm_rank(Grid_Comm, &rank);
    MPI_Cart_coords(Grid_Comm, rank, ndims, Coords);
    if(Coords[0] < k0)
        ++n0;
    if(Coords[1] < k1)
        ++n1;
    MPI_Cart_shift(Grid_Comm, 0, 1, &left, &right);
    MPI_Cart_shift(Grid_Comm, 1, 1, &up, &down);
    if(rank == 0) 
    {
        sprintf(str,"PuassonSerial_ECGM_%dx%d.log", NX, NY);
        fp = fopen(str,"w");
        fprintf(fp,"The Domain: [0,%f]x[0,%f], number of points: N[0,A] = %d, N[0,B] = %d;\n"
                   "The steep descent iterations number: %d\n"
                   "The conjugate gradient iterations number: %d\n",
                    A,B, NX,NY, SDINum,CGMNum);
	}
	SolVect   = (double *)malloc((n0+1)*(n1+1)*sizeof(double));
	ResVect   = (double *)malloc((n0+1)*(n1+1)*sizeof(double));
	RHS_Vect  = (double *)malloc((n0+1)*(n1+1)*sizeof(double));
    if (Coords[0] < k0) 
    {
        x1 = Coords[0]*n0+1;
        x2 = x1 + n0 - 1;
    } else 
    {
        x1 = k0*(n0+1) + (Coords[0] - k0)*n0 + 1;
        x2 = x1 + n0 - 1;
    }
    if (Coords[1] < k1) 
    {
        y1 = Coords[1]*n1+1;
        y2 = y1 + n1 - 1;
    } else 
    {
        y1 = k1*(n1+1) + (Coords[1] - k1)*n1 + 1;
        y2 = y1 + n1 - 1;
    }

    sendLeftBuf = (double *)malloc(n1*sizeof(double));
    recLeftBuf = (double *)malloc(n1*sizeof(double));
    sendRightBuf = (double *)malloc(n1*sizeof(double));
    recRightBuf = (double *)malloc(n1*sizeof(double));
    sendUpBuf = (double *)malloc(n0*sizeof(double));
    recUpBuf = (double *)malloc(n0*sizeof(double));
    sendDownBuf = (double *)malloc(n0*sizeof(double));
    recDownBuf = (double *)malloc(n0*sizeof(double));

	memset(ResVect,0,(n0+1)*(n1+1)*sizeof(double));
	memset(SolVect,0,(n0+1)*(n1+1)*sizeof(double));
	RightPart(RHS_Vect);

	for(i=0; i<NX; i++)
	{
		SolVect[i] = BoundaryValue(x(i),0.0);
		SolVect[n0*(n1-1)+i] = BoundaryValue(x(i),B);
	}
	for(j=0; j<NY; j++)
	{
		SolVect[n0*j] = BoundaryValue(0.0,y(j));
		SolVect[n0*j+(n0-1)] = BoundaryValue(A,y(j));
	}

	#ifdef Print
        if(rank == 0) 
        {
            printf("\nSteep descent iterations begin ...\n");
        }
	#endif

	for(iteration=1; iteration<=SDINum; iteration++)
	{
	    if (iteration > 1) 
	    {
            boundaryExchange(SolVect, Coords, dims, n0, n1, x1, x2, y1, y2, Grid_Comm, left, right, up, down);
	    }
// The residual vector r(k) = Ax(k)-f is calculating ...
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				ResVect[n0*j+i] = LeftPart(SolVect,i,j)-RHS_Vect[n0*j+i];

// The value of product (r(k),r(k)) is calculating ...
		sp = 0.0;
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				sp += ResVect[n0*j+i]*ResVect[n0*j+i]*hx*hy;

        MPI_Allreduce(&sp, &sumSp, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        sp = sumSp;
		tau = sp;

        boundaryExchange(ResVect, Coords, dims, n0, n1, x1, x2, y1, y2, Grid_Comm, left, right, up, down);


// The value of product sp = (Ar(k),r(k)) is calculating ...
		sp = 0.0;
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				sp += LeftPart(ResVect,i,j)*ResVect[n0*j+i]*hx*hy;

        MPI_Allreduce(&sp, &sumSp, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        sp = sumSp;
		tau = tau/sp;

// The x(k+1) is calculating ...
        err = 0.0;
        for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++) 
            {
                NewValue = SolVect[n0*j+i]-tau*ResVect[n0*j+i];
                s1 = NewValue-SolVect[n0*j+i];
                err += s1*s1*hx*hy;
                SolVect[n0*j+i] = NewValue;
            }
        // err summarizing
        MPI_Allreduce(&err, &sumErr, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        // err has been summirized

    }
// the end of steep descent iteration.

	BasisVect = ResVect;    // g(0) = r(k-1).
	ResVect = (double *)malloc(n0*n1*sizeof(double));
	memset(ResVect,0,n0*n1*sizeof(double));

// CGM iterations begin ...
// sp == (Ar(k-1),r(k-1)) == (Ag(0),g(0)), k=1.
	#ifdef Print
        if(rank == 0) 
        {
            printf("\nCGM iterations begin ...\n");
        }
	#endif

	for(iteration=0; iteration<CGMNum; iteration++)
	{
	    // SolVect exchange
        boundaryExchange(SolVect, Coords, dims, n0, n1, x1, x2, y1, y2, Grid_Comm, left, right, up, down);

        // SolVect has been changed
	// The residual vector r(k) is calculating ...
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				ResVect[n0*j+i] = LeftPart(SolVect,i,j)-RHS_Vect[n0*j+i];

        // ResVect exchange
        boundaryExchange(ResVect, Coords, dims, n0, n1, x1, x2, y1, y2, Grid_Comm, left, right, up, down);

        // ResVect has been changed

	// The value of product (Ar(k),g(k-1)) is calculating ...
		alpha = 0.0;
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				alpha += LeftPart(ResVect,i,j)*BasisVect[n0*j+i]*hx*hy;

        MPI_Allreduce(&alpha, &sumAlpha, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        alpha = sumAlpha;
		alpha = alpha/sp;

	// The new basis vector g(k) is being calculated ...
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				BasisVect[n0*j+i] = ResVect[n0*j+i]-alpha*BasisVect[n0*j+i];

	// The value of product (r(k),g(k)) is being calculated ...
		tau = 0.0;
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				tau += ResVect[n0*j+i]*BasisVect[n0*j+i]*hx*hy;

        MPI_Allreduce(&tau, &sumTau, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        tau = sumTau;
	// The value of product sp = (Ag(k),g(k)) is being calculated ...

        // BasisVect exchange
        boundaryExchange(BasisVect, Coords, dims, n0, n1, x1, x2, y1, y2, Grid_Comm, left, right, up, down);
        
        // BasisVect has been changed

		sp = 0.0;
		for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++)
				sp += LeftPart(BasisVect,i,j)*BasisVect[n0*j+i]*hx*hy;


        MPI_Allreduce(&sp, &sumSp, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        sp = sumSp;
		tau = tau/sp;

	// The x(k+1) is being calculated ...
        err = 0.0;
        for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++) 
            {
                NewValue = SolVect[n0*j+i]-tau*BasisVect[n0*j+i];
                s1 = NewValue-SolVect[n0*j+i];
                err = Max(err, fabs(NewValue - SolVect[n0 * j + i]));
                SolVect[n0*j+i] = NewValue;

            }
        // err summarizing

        MPI_Allreduce(&err, &sumErr, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
        // err has been summirized
        if(rank == 0) 
        {
            if(iteration%Step == 0)
            {
                //printf("The %d iteration of CGM method has been carried out.\n", counter);

    #ifdef Print
                fprintf(fp,"\nThe iteration %d of conjugate gradient method has been finished.\n"
                           "The value of \\alpha(k) = %f, \\tau(k) = %f. The difference value is %f.\n",\
                            iteration, alpha, tau, sumErr);
    #endif

            }
        }
		if (sumErr < EPSILON) 
		{
            break;
		}
	}
// the end of CGM iterations.



    if(rank == 0) 
    {
        t_finish =  MPI_Wtime();
        ex_time = t_finish - t_start;
    }
    err = 0.0;
        for(i=x1; i <= x2; i++)
            for(j=y1; j <= y2; j++) 
            {
                s1 = Solution(x(i),y(i))-SolVect[NX*j+i];
                err += s1*s1*hx*hy;
            }
    // err summarizing

    MPI_Allreduce(&err, &sumErr, 1, MPI_DOUBLE, MPI_SUM, Grid_Comm);
    // err has been summirized
    if(rank == 0) 
    {
        printf("\nThe %d iterations are carried out. The error of iterations is estimated by %.12f.\n",
                    SDINum+iteration, sumErr);
        printf("\nThe executing time is %f.\nNumber of processors = %d\n",
                    ex_time, ProcNum);
    // printing some results ...
        fprintf(fp,"\nThe %d iterations are carried out. The error of iterations is estimated by %.12f.\n",
                    SDINum+iteration, sumErr);
        fclose(fp);
    }


	free(SolVect); free(ResVect); free(BasisVect); free(RHS_Vect);
	free(sendLeftBuf); free(recLeftBuf); free(sendRightBuf); free(recRightBuf);
	free(sendUpBuf); free(recUpBuf); free(sendDownBuf); free(recDownBuf);

	MPI_Finalize();
    // The end of MPI session ...

	return(0);
}
