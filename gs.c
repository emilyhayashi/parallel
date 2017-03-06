#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float *err_arr; /* An array holding the error values for each unknown*/
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
float getNewX(float, int); /*calculates the expression to make a new X*/
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
 nit++;

    
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
float old_x, new_x, err, total;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

burden = num/comm_sz;

int start = my_rank*burden;

int complete = 0;
while (!complete){

  if(my_rank != 0){
    for(int i = start; i < start+burden; i++){

      old_x = x[i];
      new_x = getNewX(x[i], i);
      err = (new_x - old_x)/new_x;

      if(err < 0){
        err *= -1;
        err_arr[i] = err;
      }
      else{
        err_arr[i] = err;
      }

      MPI_Send(&new_x, 100, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }
  }

  else{
    float * new_X_arr = (float**)malloc(num * sizeof(float*));
    for(int i = start; i < burden+rank; i++){
      old_x = &x[0];
      new_x = getNewX(&x[i], i);
      new_X_arr[i] = new_x;
      err = (new_x - old_x)/new_x;
    }
    if (err < 0){
      err_arr[i] = -1 * err;
    }
    else{
      err_arr[i] = err;
    }

  for(int q = burden; q < num; q++ ){
    MPI_Recv(&new_x, 100, MPI_FLOAT, q%burden, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    new_X_arr[q] = new_x;
  }

  &x = &new_X_arr;
  int passing = 1;
  for(int i = 0; i < num; i++){
    if(err_arr[i] > err){
      passing = 0;
    }
  }
  if(!passing){
    MPI_Bcast(x, num, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }

  }

 }

MPI_Finalize();
 
 
 /* Writing to the stdout */
 /* Keep that same format */
 for( i = 0; i < num; i++)
     printf("%f\n",x[i]);
 
 printf("total number of iterations: %d\n", nit);
 
 exit(0);

}

float getNewX(float unknown, int index) {

  float divide_by = 0;
  float sum = &b[index];
  for(int j = 0; j < num; j++){
    
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
