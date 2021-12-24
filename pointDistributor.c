#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

void shuffle(int *array, size_t n) {    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand48(usec);
    if (n > 1) {
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48()*(i+1));
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}
void createPoints(char* filename,int points_power_of_2,int procedures_num){
    int64_t n, dims;
    int points_num = pow(2,points_power_of_2);
    int idx,i,j,dimensions;
    FILE *f;
    double **points_arr;
    f = fopen(filename,"rb"); //open binary file
    fread(&dims,sizeof(int64_t),1,f); //read dimensions of points
    fread(&n,sizeof(int64_t),1,f); //read amount of points
    points_arr = malloc(n*sizeof(double*)); 
    for (i=0;i<n;i++){
        points_arr[i] = malloc(dims*sizeof(double));
    }
    dimensions = (int)dims;
    //create a matrix with all the points
    for(i=0;i<n;i++){
        for (j=0;j<dims;j++){
            fread(&points_arr[i][j],sizeof(double),1,f);
        }
    }
    fclose(f);
    
    //create an array with the indexes of the points
    
    int *points_index_arr = malloc(n*sizeof(int));
    for(i=0;i<points_num;i++){
        points_index_arr[i] = i;
    }
    //shuffle points indexes
    shuffle(points_index_arr,(size_t)points_num);
    int points_per_procedure = points_num/procedures_num;
    for (i=0;i<procedures_num;i++){
        double **pts = malloc(points_per_procedure*sizeof(double*));
        for(j=0;j<points_per_procedure;j++){
            idx = points_index_arr[i*points_per_procedure+j];
            pts[j]=points_arr[idx];
        }
        char new_file[15];
        sprintf(new_file,"points-%d",i);
        f = fopen(new_file,"wb");
        fwrite(&points_per_procedure,sizeof(int),1,f);
        fwrite(&dimensions,sizeof(int),1,f);
        for(j=0;j<points_per_procedure;j++){
            fwrite(pts[j],dimensions*sizeof(double),1,f);
        }
        fclose(f);
    } 
}

int main(int argc,char *argv[]){
    int points_power_of_2 = 16;
    int procedures;
    if(argc<2){
        procedures=16;
    }
    else{
        procedures = atoi(argv[1]);
    }
    createPoints("points.bin",points_power_of_2,procedures);
    return 0;
}