#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"

#define DAT double
#define maxname 20

int main(void)
{
	MPI_Init(NULL, NULL);
	
    // Get the number of processes
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	printf("Heat Diffusion in 1D (processor: %d)\n\n", world_rank);
	
	
	DAT xs, xe, dx, Ts, Te, K, dt, time;
	int i, t, nt, dimx, dimx_loc;
	double sendbuf, recievebuf;
	
	dimx = 102;
	dimx_loc = 100;
	
	if(world_rank == 0 && world_size > 1)
	{
		xs = 0;
		xe = 4;
	}
	else if(world_rank == 1 && world_size > 1)
	{
		xs = 4;
		xe = 8;
	}
	else
	{
		xs = 0;
		xe = 8;
		dimx -= 2;
	}
	
	// Define geometry parameters
	DAT x[dimx], T[dimx], dT1[dimx-1], dT2[dimx];

	Ts = 0;
	Te = 100;
	
	K = 1;
	nt = 50;
	dx = (xe-xs)/(dimx_loc-1);
	dt = (dx*dx)/K/2.1;  // Magic stabil criterium
	
	// Create coordinates and initial conditions
	x[0] = xs;      T[0] = 0;
	x[dimx-1] = xe; T[dimx-1] = Te;
	for(i=1; i<dimx-1; i++)
	{
		x[i] = xs+i*dx;
		if(x[i]<4)
		{
			T[i] = 0;
		}
		else
		{
			T[i] = Te;
		}	
	}
	
	// Time-loop
	for(t=0; t<nt; t++)
	{
		// MPI Communication
		if(world_rank == 0 && world_size > 1)
		{
			sendbuf = T[dimx-2];
			MPI_Send(&sendbuf,    1, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
			MPI_Recv(&recievebuf, 1, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			T[dimx-1] = recievebuf;
		}
		else if(world_rank == 1 && world_size > 1)
		{
			sendbuf = T[1];
			MPI_Send(&sendbuf,    1, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(&recievebuf, 1, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			T[0] = recievebuf;
		}
		
		// 2 step FD - stable
		for(i=0; i<dimx-1; i++)
		{
			dT1[i] = -K * (T[i+1] - T[i])/dx;
		}
		
		for(i=1; i<dimx-1; i++)
		{
			dT2[i] = -((dT1[i] - dT1[i-1])/dx);
		}
		
		/*
		// Using central FD - super unstable
		for(i=1; i<dimx-1; i++)
		{
			dT2[i] = -K * (T[i+1] - 2*T[i] + T[i-1])/(dx*dx);
		}
		*/
		
		for(i=0; i<dimx; i++)
		{
			T[i] += dT2[i] * dt;
		}
		
		time += dt;
		
		if(world_rank == 0)
		{
			printf("Time = %.5f\n",time);
		}
		
		FILE *fb;
		char fname[maxname];
		snprintf(fname,maxname,"T_%d_p%d.dat",t+1,world_rank);
		fb = fopen(fname,"w");
		fwrite(&T,sizeof(DAT),dimx,fb);
		fclose(fb);
	}

	    MPI_Finalize();
}