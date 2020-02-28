#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

float tolerance = 0.00001;

float my_random(){
   float val = (double)rand()/(double)RAND_MAX;
   return val;
}

void read_data(float* array,char str[],int rows,int cols){
   FILE *file_pointer;
   file_pointer = fopen(str,"r");
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         fscanf(file_pointer,"%f ",&array[i*cols + j]);
      }
   }

   fclose(file_pointer);
}

void write_data(char file_name[],float* matrix,int rows,int cols){
   FILE *file_pointer;
   file_pointer = fopen(file_name,"w");
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         fprintf(file_pointer, "%f",matrix[i*cols + j]);
         if(j!=cols-1){
            fprintf(file_pointer, " ");
         }else if(i!=rows-1){
            fprintf(file_pointer, "\n");
         }
      }
   }
   fclose(file_pointer);
   return;
}

void generate_data(float* arr,int rows,int cols){
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         arr[i*cols + j] = my_random();
      }
   }
}

void print(float* array,int rows,int cols){
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         printf("%f ",array[i*cols + j]);
      }
      printf("\n");
   }
   printf("**********THE END**********\n");
}

void transpose(float* array,float* arrayT,int rows,int cols){
   for(int i = 0;i<cols;i++){
      for(int j = 0;j<rows;j++){
         arrayT[i*rows+j] = array[j*cols + i];
      }
   }
}

void multiply(float* A,float* B,float* prod,int n,int m){
   int i,j,k;
   float val = 0;

   for(i = 0;i<n;i++){
      for(j = 0;j<n;j++){
         val = 0;
         for(k = 0;k<m;k++){
            val+=(A[i*m + k] * B[k*n + j]);
         }
         prod[i*n + j] = val;
      }
   }
}

void IsEqual(float* A,float* B,int n){
   for(int i = 0;i<n;i++){
      if(fabs(A[i] - B[i])>tolerance){
         printf("%d,%lf,%lf\n",i,A[i],B[i]);
         printf("Incorrect\n");
         return;
      }
   }
   printf("Correct\n");
}

int main(int argc, char* argv[]){
   int m = 32;
   int n = atoi(argv[1]);
   char* file_name1 = argv[2];
   char* file_name2 = argv[3];
   int num_workers = atoi(argv[4])-1;

   float* A = (float*)malloc(n*m*sizeof(float));
   float* B = (float*)malloc(n*m*sizeof(float));
   float* C = (float*)malloc(n*m*sizeof(float));
   float* prodM = (float*)malloc(n*n*sizeof(float));
   float* prodS = (float*)malloc(n*n*sizeof(float));

   // Initialize the MPI Environment
   MPI_Init(&argc,&argv);
   clock_t t;
   int rank,idx;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int load,extra,num_rows,num_cols;

   if(rank==0){
      // read_data(A,file_name1,n,m);
      // read_data(C,file_name2,m,n);

      idx = 0;
      load = n/num_workers;

      generate_data(A,n,m);
      generate_data(C,m,n);

      t = clock();
      transpose(C,B,m,n);

      for(int i = 1;i<=num_workers;i++){
         if(i==num_workers && num_workers>1){
            load = n%num_workers;
         }
         MPI_Send(&n,1,MPI_INT,i,0,MPI_COMM_WORLD);
         MPI_Send(&m,1,MPI_INT,i,0,MPI_COMM_WORLD);
         MPI_Send(&load,1,MPI_INT,i,0,MPI_COMM_WORLD);
         MPI_Send(&A[idx*m],load*m,MPI_FLOAT,i,0,MPI_COMM_WORLD);
         MPI_Send(&B[0],n*m,MPI_FLOAT,i,0,MPI_COMM_WORLD);
         idx+=load;
      }

      idx = 0;
      load = n/num_workers;
      for(int i = 1;i<=num_workers;i++){
         if(i==num_workers && num_workers>1){
            load = n%num_workers;
         }
         MPI_Recv(&prodM[idx*m],load*n,MPI_FLOAT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         idx+=load;
      }

      t = clock()-t;
      double time_taken = ((double)t)/CLOCKS_PER_SEC;
      printf("Time taken for MPI P2P Blocking call:- %lf\n",time_taken);

      t = clock();
      multiply(A,C,prodS,n,m);
      t = clock()-t;
      time_taken = ((double)t)/CLOCKS_PER_SEC;
      printf("Time taken for sequential call:- %lf\n",time_taken);

      // print(prodM,n,n);
      // print(prodS,n,n);
      IsEqual(prodS,prodM,n*n);

      char file_name[30] = "product_P2PB_";
      sprintf(file_name+13,"%d",n);
      strcat(file_name,".txt");
      write_data(file_name,prodM,n,n);
      
   }else{
      MPI_Recv(&num_rows,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&num_cols,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&load,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&A[0],load*num_cols,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&B[0],num_rows*num_cols,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      double val;
      for(int i = 0;i<load;i++){
         for(int j = 0;j<num_rows;j++){
            val = 0;
            for(int k = 0;k<num_cols;k++){
               val+=(A[i*num_cols + k]*B[j*num_cols + k]);
            }
            C[i*num_rows+j] = val;
         }
      }

      MPI_Send(&C[0],load*num_rows,MPI_FLOAT,0,0,MPI_COMM_WORLD);
   }

   MPI_Finalize();
   return 0;
} 
