#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
float getNewX( int); /*calculates the expression to make a new X*/
/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
     Conditions for convergence (diagonal dominance):
     1. diagonal element >= sum of all other elements of the row
     2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
    int bigger = 0; /* Set to 1 if at least one diag element > sum  */
    int i, j;
    float sum = 0;
    float aii = 0;
    
    for(i = 0; i < num; i++)
    {
        sum = 0;
        aii = fabs(a[i][i]);
        
        for(j = 0; j < num; j++)
             if( j != i)
     sum += fabs(a[i][j]);
             
        if( aii < sum)
        {
            printf("The matrix will not converge.\n");
            exit(1);
        }
        
        if(aii > sum)
            bigger++;
        
    }
    
    if( !bigger )
    {
         printf("The matrix will not converge\n");
         exit(1);
    }
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[])
{
    FILE * fp;
    int i,j;  
 
    fp = fopen(filename, "r");
    if(!fp)
    {
        printf("Cannot open file %s\n", filename);
        exit(1);
    }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);
 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
    {
    printf("Cannot allocate a!\n");
    exit(1);
    }

 for(i = 0; i < num; i++) 
    {
        a[i] = (float *)malloc(num * sizeof(float)); 
        if( !a[i])
        {
        printf("Cannot allocate a[%d]!\n",i);
        exit(1);
        }
    }
 
 x = (float *) malloc(num * sizeof(float));
 if( !x)
        {
    printf("Cannot allocate x!\n");
    exit(1);
    }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
    {
    printf("Cannot allocate b!\n");
    exit(1);
    }

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
    fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
     for(j = 0; j < num; j++)
         fscanf(fp,"%f ",&a[i][j]);
     
     /* reading the b element */
     fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp);
 

}


/************************************************************/



int main(int argc, char *argv[])
{
 int i;
 int nit = 0; /* number of iterations */    
 if( argc != 2)
 {
     printf("Usage: gsref filename\n");
     exit(1);
 }
    
 /* Read the input file and fill the global data structure above */ 
 get_input(argv[1]);
 
 /* Check for convergence condition */
 /* This function will exit the program if the coffeicient will never converge to 
    * the needed absolute error. 
    * This is not expected to happen for this programming assignment.
    */
 //check_matrix();
 
  MPI_Init(NULL, NULL);

  int my_rank, comm_sz, burden;
  float old_x, new_x, this_err, total;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Status status;
  burden = num/comm_sz;

  int start = my_rank*burden;
  int complete = 0;
  float *sub_err_arr = malloc(sizeof(float) * burden);
  float *err_arr = NULL;

  err_arr = (float *)malloc(sizeof(float)*num);
  float *new_X_arr = NULL;
  if(my_rank == 0){    
    new_X_arr = (float *)malloc(sizeof(float)*num);
  }

  while (nit < 6){ 
    nit++;
    printf("iteration: %d\n", nit);
    if(my_rank !=0){
      int k = 0;
      for(i = start; i < start+burden; i++){
        old_x = x[i]; 
        new_x = getNewX(i);
        this_err = (new_x - old_x)/new_x; 
        if(this_err >= 0){
          sub_err_arr[k] = 100*this_err;
        }
        else{
          sub_err_arr[k] = -100*this_err;
        }
        printf("old: %f  new: %f\n", old_x, new_x);
        printf("error: %f\n", sub_err_arr[k]);
        k++;
        MPI_Send(&new_x, sizeof(float*), MPI_FLOAT, 0, i, MPI_COMM_WORLD);
        printf("process %d sent %f\n", my_rank, new_x); 
      }
      MPI_Send(sub_err_arr, sizeof(float)*burden, MPI_FLOAT, 0, start, MPI_COMM_WORLD);
    }
    else{
      float recv_x =0;
      float *recv_err_arr = malloc(sizeof(float*)*burden);
      int i;
      for(i = start; i < burden+start; i++){
        old_x = x[i];
        new_x = getNewX(i);
        new_X_arr[i] = new_x;
        this_err = (new_x - old_x)/new_x;
        if (this_err < 0){
          err_arr[i] = -100 * this_err;
        }
        else{
          err_arr[i] = 100*this_err;
        }
      }   
      for(i = burden; i < num; i++ ){     
        MPI_Recv(&recv_x, sizeof(float*), MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        new_X_arr[status.MPI_TAG] = recv_x;
        printf("process %d received new x: %f at index %d\n", my_rank, new_X_arr[status.MPI_TAG], status.MPI_TAG);
      }
      int j;
      for (i = 1; i < comm_sz; i++){
        MPI_Recv(sub_err_arr, sizeof(float)*burden, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for (j = 0; j < burden; j++){
          err_arr[status.MPI_TAG + j] = sub_err_arr[j];
          
        }
      }
      x = new_X_arr;

      complete=1;
      for(i = 0; i < num; i++){
        printf("new_x_arr[%d] = %f\n", i, new_X_arr[i]);
        printf("error[%d] = %f\n", i, err_arr[i]);
        if(err_arr[i] > err){
          complete = 0;
        }
      }
      if(!complete){
       printf("incremented nit\n");
    MPI_Bcast(x, comm_sz, MPI_FLOAT, 0, MPI_COMM_WORLD);
      }
      else{
        for( i = 0; i < num; i++){
          printf("%f\n",x[i]);
        }
        printf("total number of iterations: %d\n", nit);
      } 
    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Bcast(err_arr, comm_sz, MPI_FLOAT, 0, MPI_COMM_WORLD);  
    }

  }

  printf("rank %d completed\n", my_rank);

  MPI_Finalize();   
  exit(0);
}
 
 
 /* Writing to the stdout */
 /* Keep that same format */



float getNewX(int index) {

  float divide_by = 0;
  float sum = b[index];
  
  int j;
  for(j = 0; j < num; j++){
    
      if(j == index){
       divide_by = a[index][j];
      }
      else{
        sum -= x[j] * a[index][j];
      }
    }

  float newX = sum/divide_by;

  return newX;

}
