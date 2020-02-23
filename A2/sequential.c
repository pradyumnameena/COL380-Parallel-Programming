#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float tolerance = 0.000001;

float* read_data(char str[],int rows,int cols){
   float* array = (float*)malloc(rows*cols*sizeof(float));

   FILE *file_pointer;
   file_pointer = fopen(str,"r");
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         fscanf(file_pointer,"%f ",&array[i*cols + j]);
      }
   }

   fclose(file_pointer);
   return array;
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

void print(float* array,int rows,int cols){
   for(int i = 0;i<rows;i++){
      for(int j = 0;j<cols;j++){
         printf("%f ",array[i*cols + j]);
      }
      printf("\n");
   }
   printf("**********THE END**********\n");
}

void multiply(float* A,float* B,float* prod,int n){
   int i,j,k;
   float val = 0;

   for(i = 0;i<n;i++){
      for(j = 0;j<n;j++){
         val = 0;
         for(k = 0;k<32;k++){
            val+=(A[i*32 + k] * B[k*n + j]);
         }
         prod[i*n + j] = val;
      }
   }
}

int IsEqual(float* A,float* B,int m,int n){
   for(int i = 0;i<m;i++){
      for(int j = 0;j<n;j++){
         if(fabs(A[i*n + j] - B[i*n + j])>tolerance){
            return 0;
         }
      }
   }
   return 1;
}

int main(int argc, char* argv[]){
   int n = atoi(argv[1]);
   char* file_name1 = argv[2];
   char* file_name2 = argv[3];
   float* A = read_data(file_name1,n,32);
   float* B = read_data(file_name2,32,n);
   float* prod = (float*)malloc(n*sizeof(float));

   multiply(A,B,prod,n);

   char file_name[30] = "product_";
   sprintf(file_name+8,"%d",n);
   strcat(file_name,".txt");
   write_data(file_name,prod,n,n);
   
   return 0;
}