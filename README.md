# Παράλληλα και Διανεμημένα Συστήματα

## 2η Εργασία

Για να τρέξετε την εργασία αφού κάνετε make, τρέξτε το script της Julia μετά το pointDistributor και τέλος το mpi_Distribute_Points.

Αν δεν δουλεύει το makefile οι εντολές για το compile:
```
gcc -lm pointDistributor.c -o pointDistributor   
mpicc -lm mpi_Distribute_Points.c -o mpi_Distribute_Points
```
