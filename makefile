CFLAGS= -lm

all: pointDistributor mpi_Distribute_Points

pointDistributor: pointDistributor.c
	gcc  -o pointDistributor pointDistributor.c $(CFLAGS) 


mpi_Distribute_Points: mpi_Distribute_Points.c
	mpicc -o mpi_Distribute_Points  mpi_Distribute_Points.c  $(CFLAGS) 
