CFLAGS= -lm

all: pointDistributor mpi_Distribute_Points

pointDistributor: pointDistributor.o
	gcc $(CFLAGS) -o pointDistributor pointDistributor.o

pointDistributor.o: pointDistributor.c   
	gcc $(CFLAGS) -c pointDistributor.o pointDistributor.c

mpi_Distribute_Points: mpi_Distribute_Points.o
	mpicc $(CFLAGS) -o mpi_Distribute_Points mpi_Distribute_Points.o

mpi_Distribute_Points.o: mpi_Distribute_Points.c   
	mpicc $(CFLAGS) -c mpi_Distribute_Points.o mpi_Distribute_Points.c

clean: 
	$(RM) count *.o *~