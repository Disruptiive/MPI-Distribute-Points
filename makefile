CFLAGS= -lm

all: pointDistributor mpi_Distribute_Points

pointDistributor: pointDistributor.o
	gcc -o pointDistributor pointDistributor.o $(CFLAGS) 

pointDistributor.o: pointDistributor.c   
	gcc -c pointDistributor.o pointDistributor.c $(CFLAGS)

mpi_Distribute_Points: mpi_Distribute_Points.o
	mpicc -o mpi_Distribute_Points mpi_Distribute_Points.o  $(CFLAGS) 

mpi_Distribute_Points.o: mpi_Distribute_Points.c   
	mpicc -c mpi_Distribute_Points.o mpi_Distribute_Points.c $(CFLAGS) 

clean: 
	$(RM) count *.o *~
