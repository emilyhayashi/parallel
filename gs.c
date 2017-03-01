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

 int my_rank, comm_sz, divide_by, turns;
 double old_x, new_x, new_err, err, total;
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
 int turns = num/comm_sz;
while (!complete){
  if(my_rank != 0){
    for(int i = my_rank; i < my_rank+num; i+=num){
      total = &b[i];
      for(int j = 0; j < num; j++){
        if(j == i){
          divide_by = a[i][j];
        }
        else{
        total -= x[j] * a[i][j];
        }
      }
    }
    new_x = total/divide_by;
    err = (new_x - old_x)/new_x;
    if(err < 0){
      err *= -1;
      MPI_Send(err, 100, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else{
      MPI_Send(err, 100, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }
}























 if(my_rank != 0){
    old_x = &x[my_rank];
    //compute the error
    //send it as a message to rank 0
 }

 else{
    old_x = &x[0];
    &x[0] = (&b[0] - (&&a[0,1] * &x[1]) - (&&a[0,2] * &x[2]))/a[0,0];
    new_err = (&x[0] - old_x)/&x[0];
    if (new_err < 0){
        err_arr[0] = -1 * new_err;
    }
    else{
        err_arr[0] = new_err;
    }

    //listen for all the new errors, overwrite the array
    //check if they're all under what they need to be
        //if they are, end MPI and the program
        //if they aren't, maybe this should be in a while loop?
 }
 
 
 /* Writing to the stdout */
 /* Keep that same format */
 for( i = 0; i < num; i++)
     printf("%f\n",x[i]);
 
 printf("total number of iterations: %d\n", nit);
 
 exit(0);

}
