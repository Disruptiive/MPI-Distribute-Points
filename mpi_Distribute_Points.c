#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h> 
#include <math.h>
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include "mpi.h"

//exchange locally in every procedure sent points with the ones that they received
double **finishTransfers(double *points_received,double median,double **points,int size,int nums_to_transfer,int dimensions,int taskid, double *distance_array,int pointsNum,int low){
    int i,j;
    int count = 0;

    if(taskid>=size/2+low){
        for(i=0;i<pointsNum;i++){
            if(distance_array[i]<median){
                memcpy(points[i],points_received+count,dimensions*sizeof(double));
                count += dimensions;
            }
        }
    }
    else{
        for(i=0;i<pointsNum;i++){
            if(distance_array[i]>median){
                memcpy(points[i],points_received+count,dimensions*sizeof(double));
                count += dimensions;
            }
        }
    }
    return points;
}

//execute transfers based on the transfer table and return the points after the transfers are done
double **executeTransfers(double median,double*leader, int **transfer_array,double *pts_to_send,double **points,int size,int low,int high,int taskid,int nums_to_transfer,int dimensions,double *distance_array,int pointsNum){
    int i,j,count,swaps,destination,sender;
    MPI_Request request;
    //allocate 1d array to hold all points received (i used 1d array to be sure that the memory is contiguous when transfering them) 
    double *pts_received = malloc(nums_to_transfer*dimensions*sizeof(double));

    //traversing transfer_array (row-wise of column-wise depending on taskid) execute Sendrecv order between the two procedures and go on continue until every procedure has sent and received all the points needed
    count=0;
    if(taskid-low<size/2){
        for(i=0;i<size/2;i++){
            if (transfer_array[taskid-low][i] != 0){
                swaps = transfer_array[taskid-low][i];
                sender = taskid;
                destination = low+size/2+i;
                //printf("1)TASKID: %d SENDING TO: %d COUNT: %d SWAPS: %d NUMS_to_transfer: %d\n",sender,destination,count,swaps,nums_to_transfer);
                MPI_Sendrecv(&pts_to_send[count],dimensions*swaps,MPI_DOUBLE,destination,0,&pts_received[count],dimensions*swaps,MPI_DOUBLE,destination,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                count += swaps*dimensions;
                
            }
        }
    }
    else{
        for(i=0;i<size/2;i++){
            if (transfer_array[i][taskid-low-size/2] != 0){
                swaps = transfer_array[i][taskid-low-size/2];
                sender = taskid;
                destination = low+i;
                //printf("2)TASKID: %d SENDING TO: %d COUNT: %d SWAPS: %d NUMS_to_transfer: %d\n",sender,destination,count,swaps,nums_to_transfer);
                MPI_Sendrecv(&pts_to_send[count],dimensions*swaps,MPI_DOUBLE,destination,0,&pts_received[count],dimensions*swaps,MPI_DOUBLE,destination,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                count += swaps*dimensions;
                
            }
        }
    }
    //replace locally in every procedure the sent points with the received ones
    points = finishTransfers(pts_received, median, points,size, nums_to_transfer, dimensions, taskid, distance_array, pointsNum, low);
    return points;
}
//creates a size/2 x size/2 table with the amount of points that need to be exchanged between procedures
int **createTransferArray(int *swaps_array, int size){
    int i,j;
    int **transferArray;
    
    transferArray = calloc(size/2,sizeof(int*));
    for (i=0;i<size/2;i++){
        transferArray[i] = calloc(size/2,sizeof(int));
    }
    i=0;
    j=size-1;
    //create the transfer table 
    while(i<size/2 && j>=size/2){
        if(swaps_array[i] > swaps_array[j]){
            transferArray[i][j-size/2] = swaps_array[j];
            swaps_array[i] -= swaps_array[j];
            swaps_array[j] = 0;
            j--;
        }
        else if(swaps_array[i] < swaps_array[j]){
            transferArray[i][j-size/2] = swaps_array[i];
            swaps_array[j] -= swaps_array[i];
            swaps_array[i] = 0;
            i++;
        }
        else{
            transferArray[i][j-size/2] = swaps_array[i];
            swaps_array[i] = 0;
            swaps_array[j] = 0;
            i++;
            j--;
        }
    }
    return transferArray;
}

//prepare procedures for point swap. Each procedure counts how many points they need to swap and groups them in an array
double *prepareProcedures(int *nums_to_transfer,double median,int low, int high, int pointsNum,double **points,double *distance_array,int taskid,int dimensions){
    int i,count,size;
    size = high-low+1;
    count = 0;
    //calculate how many elements are larger or smaller than median depending on taskid and the range [low,high] and update nums_to_transfer
    if(taskid>=size/2+low){
        for(i=0;i<pointsNum;i++){
            if(distance_array[i]<median){
                (*nums_to_transfer)++;
            }
        }
    }
    else{
        for(i=0;i<pointsNum;i++){
            if(distance_array[i]>median){
                (*nums_to_transfer)++;
            }
        }
    }

    //group all points that need to be swapped in 1d array and return it
    double *pts_to_swap = malloc(*nums_to_transfer*dimensions*sizeof(double));
   
    if(taskid>=size/2+low){
        for(i=0;i<pointsNum;i++){
            if(distance_array[i]<median){
                memcpy(pts_to_swap+count*dimensions, points[i],dimensions*sizeof(double));
                count += 1;
            }
        }
    }
    else{
        for(i= 0;i<pointsNum;i++){
            if(distance_array[i]>median){
                memcpy(pts_to_swap+count*dimensions, points[i],dimensions*sizeof(double));
                count+= 1;
            }
        }
    }
    
    return pts_to_swap;
}

void swap(double *arr, int a, int b)
{
    double tmp = arr[a];
    arr[a] = arr[b];
    arr[b] = tmp;
}

int partition (double *arr, int left, int right)
{
    double pivot = arr[right]; // pivot
    int i = left,j;
    for (j = left; j <= (right - 1); j++)
    {
        // If current element is smaller than the pivot
        if (arr[j] < pivot)
        {
            swap(arr,i,j);
            i++;
        }
    }
    swap(arr,i,right);
    return i;
}
//standard quickselect algorithm
double quickselect(double *arr,int left,int right, int k){
    int index = partition(arr,left,right);
    if (index - left == k - 1)
        return arr[index];
    if (index - left > k - 1)
        return quickselect(arr, left, index - 1, k);
    return quickselect(arr, index + 1, right,k - index + left - 1);
}

//find median of an even array 
double findMedian(double *arr,int size){
    //find n/2 smallest and n/2+1 smallest number of the array sum them and divide by 2
    double median1 = quickselect(arr,0,size-1,size/2);
    double median2 = quickselect(arr,0,size-1,size/2+1);
    return (median1+median2)/2;
}

//calculate distance between 2 double multidimensional points
double calculateDistance(double *point1, double *point2,int dimensions){
    int i;
    double result = 0;
    for (i=0;i<dimensions;i++){
        result += pow(point2[i]-point1[i],2);
    }
    result = sqrt(result);
    return result;
}

double** distributeByMedian(int low,int high,double *leader,int taskid,int pointsNum,int dimensions,double **points){
    MPI_Status status;
    MPI_Comm new_comm;

    double median,*all_distances;
    int i,color,size,nums_to_transfer,*swaps_array, **transfer_array;
        
    size = high-low+1; //size of range [low,high]
    color = ((taskid<=high) && (taskid>=low)) ? 1:MPI_UNDEFINED; //color for new communicator
   
    //create a new communicator only for procedures in range [low,high]
    MPI_Comm_split(MPI_COMM_WORLD, color, taskid, &new_comm);

    //distribute points only in [low,high] procedure range
    if(taskid<=high && taskid>=low){
        //procedure with taskid == low acts as master procedure
        if (taskid==low){
            //create and allocate array to store all distances
            all_distances = malloc(size*pointsNum*sizeof(double));
         }

        //all procedures inside [low,high] range calculate distance to leader from every point they own
        double *distance_array = malloc(pointsNum*sizeof(double));
            
        for(i=0;i<pointsNum;i++){
            distance_array[i] = calculateDistance(leader,points[i],dimensions);
        }
        
        //send all distance arrays to master 
        MPI_Gather(distance_array,pointsNum,MPI_DOUBLE,all_distances,pointsNum,MPI_DOUBLE,0,new_comm);   

        //master finds median and broadcasts it to the other procedures involved
        if (taskid==low){
            int arr_size = size*pointsNum;
            median = findMedian(all_distances,arr_size);
        }

        MPI_Bcast(&median,1,MPI_DOUBLE,0,new_comm);
        nums_to_transfer = 0;
        
        //find how many points to swap and pack them in an array
        double *pts_to_send = prepareProcedures(&nums_to_transfer,median,low,high,pointsNum,points,distance_array,taskid,dimensions);
        //all procedures share how many swaps they need 
        swaps_array = malloc(size*sizeof(int));
        MPI_Allgather(&nums_to_transfer,1,MPI_INT,swaps_array,1,MPI_INT,new_comm);   
        //and calculate a transfer table
        transfer_array = createTransferArray(swaps_array,size);
        //free(swaps_array);
        
        //execute the transfers as calculated on the transferArray
        points = executeTransfers(median,leader,transfer_array,pts_to_send,points,size,low,high,taskid,nums_to_transfer,dimensions,distance_array,pointsNum);
        MPI_Comm_free(&new_comm);

        if (taskid==low){
            free(all_distances);   
        }
        free(distance_array);
        free(swaps_array);
        free(transfer_array);
        free(pts_to_send);
    
        return points;
    }
    else{
        return points;
    }
    
}
//Every process reads its points
double** readPoints(int task,int* points_per_procedure,int* dimensions){
    int i,j;
    char name[15];
    FILE *f;
    sprintf(name,"points-%d",task);
    
    f = fopen(name,"rb"); //open binary file
    fread(points_per_procedure,sizeof(int),1,f); //read amount of points per procedure
    fread(dimensions,sizeof(int),1,f); //read dimensions of points
    double **pts = malloc(*points_per_procedure*sizeof(double*));
    for(i=0;i<*points_per_procedure;i++){
        pts[i]=malloc(*dimensions*sizeof(double));
    }
    //read the points
    for(i=0;i<*points_per_procedure;i++){
        for(j=0;j<*dimensions;j++){
            fread(&pts[i][j],sizeof(double),1,f);
        }
    }
    fclose(f);
    return pts;
}

//distributePoints() calls distributeByMedian() recursively
double** distributePoints(int low,int high,double *leader,int taskid,int pointsNum,int dimensions,double **points){
    int mid;
    //call distributeByMedian for current range [low,high] of procedures
    points = distributeByMedian(low,high,leader,taskid,pointsNum,dimensions,points);
    //MPI_Barrier(MPI_COMM_WORLD); //wait for everything to finish
    //call it recursively for upper and bottom half of procedures
    if (high - low > 1){
        mid = (high-low+1)/2 + low;   
        points = distributePoints(low,mid-1,leader,taskid,pointsNum,dimensions,points);
        points= distributePoints(mid,high,leader,taskid,pointsNum,dimensions,points);
    }   
    
    return points;
}

void validate(int taskid,double **points,int low,int high,int dimensions,double *leader,int pointsNum,int proceduresNum,double exec_time){
    int i;
    int size = high-low+1;
    double *distance_array = malloc(pointsNum*sizeof(double));
    double min = DBL_MAX;
    double max = DBL_MIN; 
    double *mins,*maxs;
    for(i=0;i<pointsNum;i++){
        distance_array[i] = calculateDistance(leader,points[i],dimensions);
    }
    for(i=0;i<pointsNum;i++){
        if (distance_array[i] > max)
            max = distance_array[i];
        if (distance_array[i] < min)
            min = distance_array[i];
    }
    if (taskid==0){
        mins = malloc(proceduresNum*sizeof(double));
        maxs = malloc(proceduresNum*sizeof(double));
    }
    MPI_Gather(&min,1,MPI_DOUBLE,&mins[taskid],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&max,1,MPI_DOUBLE,&maxs[taskid],1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
    if (taskid==0){
        printf("\tRESULTS\n");
        for(i=0;i<proceduresNum;i++){
            printf("TASK) %d MIN: %lf-MAX: %lf\n",i,mins[i],maxs[i]);
            
        }
        for(i=1;i<proceduresNum;i++){
            if(mins[i-1]>maxs[i]){
                printf("----ERROR----\n");
                printf("TASK) %d min: %lf TASK %d) max: %lf",i-1,mins[i-1],i,maxs[i]);
            }
            
        }
        printf("Execution Time: %.3lf s\n",exec_time/1000);
    }    
    free(distance_array);
    if (taskid==0){
        free(mins);
        free(maxs);
    }
}

int main (int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    struct timeval t_start,t_end;
    int  numprocedures, taskid,points_power_of_2,points_per_procedure,dimensions,i,low,high;
    double *leader,exec_time;
    MPI_Status status;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocedures);
    
    points_power_of_2=16; //pick 2^points_power_of_2 points from all availlable
    //range of procedures to sort points
    low = 0;
    high = numprocedures-1;
    
    double **pts = readPoints(taskid,&points_per_procedure,&dimensions);
    
    leader = malloc(dimensions*sizeof(double));
    if (taskid==0){
        //pick a random point from root as leader
        int rnd_num = rand()%points_per_procedure;
        leader = pts[rnd_num];
    }
    MPI_Bcast(&leader[0],dimensions,MPI_DOUBLE,0,MPI_COMM_WORLD); //send it to the other procedures
    gettimeofday(&t_start,NULL);
    pts = distributePoints(low,high,leader,taskid,points_per_procedure,dimensions,pts);
    gettimeofday(&t_end,NULL);
    exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;
    exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0;
    validate(taskid,pts, low, high, dimensions, leader,points_per_procedure,numprocedures,exec_time);

    //free(pts);
    
    MPI_Finalize();
}