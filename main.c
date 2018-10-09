//gcc -o sa main.c -lm
//gcc -o sa main.c /home/sadegh/Desktop/OpenBLAS/lib/libopenblas.a -lm
//gcc -o sa main.c -I /home/sadegh/Desktop/OpenBLAS/include/ -L/home/sadegh/Desktop/OpenBLAS/lib/ -Wl,-rpath=/home/sadegh/Desktop/OpenBLAS/lib/ -lopenblas -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixheader.h"

//defining constants
#define initialTemperature 1
#define k 1
#define alpha 0.9
#define maxIteration 500
#define maxAcceptance 380
#define initialSearch 0
#define lowerBound 0
#define upperBound 15

int numPBs=581;
int numPTVVoxels=9814;
//struct Info{

//  int numPTVVoxels;
float prescribedDose=70;
//}InfoPTV={numPTVVoxels,70};

double **influence_PTV_double;
float **influence_PTV;


//function declaration
void sa_optimization();
float evaluation(float p[]);
void newSolution(int solutionFlag,float lastPoint[],float *np);
float OF_PTV(float pbWeights[]);
float rand01();
float randn (float mu, float sigma);


int main()
{
    srand(time(NULL));//seeding random number generator
    int i,j;
    influence_PTV_double=readMatrix(9814,581);
    influence_PTV=convertDoubleToFloat(9814,581,influence_PTV_double);


    //printMatrix(16,4,influence_PTV);

    printf("\n print the rest:\n");

    sa_optimization();

    freeMatrix(9814,influence_PTV);

    return 0;
}



void sa_optimization()
{

//defining variables

//defining counters
    int it=0;  //counter variable that is used for cooling schedule
    int accepted=0; //keep the number of accepted solutions in one temperature
    int iteration=1;  //keep the total number of iteration
    int saAccepted=0;  //keep the number of accepted solutions using sa part
    int rejected=0;    //keep the number of rejected solutions
    int dec=0;
//defining function related variables
    float delta,evalNew,evalOld;
//defining point related variables --Column Vectors
    float bestSolution[numPBs],newPoint[numPBs];

//initialization of variables
    float Temperature = initialTemperature;
    int i;
    for (i=0; i<numPBs; i++)
    {
        bestSolution[i]= lowerBound+(upperBound-lowerBound)*rand01();  //this is the first guess
        // printf("x[%d]=%6.4f \t",i,bestSolution[i]);  //this is the first guess
    }

    evalOld=evaluation(bestSolution);
    printf("\n\n ****first guess for eval is %e",evalOld);

    while(iteration<80000 && evalOld>6)
    {


        it++ ;
        if (it >= maxIteration || accepted >= maxAcceptance)
        {
            it=0; //reset the counter variables
            Temperature *= alpha;  //decreasing the temperature
            printf("\n\nthe minimum found in %d iterations with temp=%3.2e is : %10.6f",iteration,Temperature,evalOld);
            printf("\naccepted: %d saAccepted :%d rejected: %d",accepted,saAccepted,rejected);
//            if (rejected==maxIteration)
//            {
//                dec++;
//                `if (dec>=5)
//                {
//                    printf("\n temperature reseted");
//                    Temperature = initialTemperature;
//                    dec=0;
//                    `
//                }
//            }
            accepted=0;
            saAccepted=0;
            rejected=0;
        }
        //generating new solutions
        //deciding which line to be used for generation new solution
        if (iteration < initialSearch)
        {
            newSolution(1,bestSolution,newPoint);   // a random search in whole search space
        }
        else
        {
            newSolution(0,bestSolution,newPoint);   // a random walk
        }
        iteration++;
        //function evaluation in newPoint
        evalNew=evaluation(newPoint);
        // printf("\n sol passed::: %f",newPoint[1]);
        //printf("\n && evalNew is :: %e",evalNew);
        delta=evalNew-evalOld;
        //deciding weather to accept or reject new solution
        int i;
        if(delta < 0)
        {
            for (i=0; i<numPBs; i++)
            {
                bestSolution[i]=newPoint[i];
            }
            evalOld = evalNew;
            accepted++;
        }
        //acceptance with some probability
        else if(delta>=0 && exp(-delta/((float)Temperature))>rand01() )
        {
            //printf("\n %f",delta);
            for (i=0; i<numPBs; i++)
            {
                bestSolution[i]=newPoint[i];
            }
            evalOld = evalNew;
            saAccepted++;
            //printf("\n****jump");
        }
        else
        {
            //printf("  ##rejection %f",delta);
            rejected++;

        }
    }//end of while loop

    //this  part returns minimum of function and the best solution
    printf("\nthe minimum found in %d iterations with temp=%3.2e is : %10.6f",iteration,Temperature,evalOld);
    printf("\n\nthe best PB solution is : \n");
    for (i=0; i<numPBs; i++)
    {
        printf("x[%d]=%6.4f \t",i,bestSolution[i]);  //this is the first guess
    }


}//end of sa_optimization function

void newSolution(int solutionFlag,float lastPoint[],float *np)
{
    if (solutionFlag == 1)   //a random search in hole search space
    {
        int i;
        for (i=0; i<numPBs; i++)
        {
            np[i]=lowerBound+(upperBound-lowerBound)*rand01();
            //printf("\n sol::: %f",np[i]);
        }
    }
    else  //random walk search in vicinity of previous point
    {
        int i;
        for (i=0; i<numPBs; i++)
        {
            np[i]=lastPoint[i]+(float)0.04*(float)(upperBound-lowerBound)*randn(0,1);  //change it to randn for normal dist
        }
    }

    //dont let the new point goes outside of search space
    int i;
    for (i=0; i<numPBs; i++)
    {
        if (np[i] < lowerBound)
        {
            np[i] = lowerBound;
        }
        else if (np[i] > upperBound)
        {
            np[i] = upperBound;
        }

    }
    //new point was returned by np pointer
}

float evaluation(float p[])
{
    float z;
    z=OF_PTV(p);
    //printf("\n in evaluation : %f",p[1]);
    //return the value of function
    return z;
}


float OF_PTV(float pbWeights[])
{
//*****these are for blas dgemv function********************************

    char trans='N'; //choose is it apropriately, T is for transpose
    //int m,n; 	//number of rows and cols
    float alphaa=1;
    //float *a; 	//this is the m-n matrix
    //float *x;	//this is the vector
    //float *y;	//result vector
    int lda;	//leading dim which is m here

    int incx,incy;	//increments in x and y elements
    incx=1;
    incy=1;
    float beta;
    beta=0;
//*****************************************************************



    //function for calculating quadratic error in PTV
    int i,j;
    float PTV_error=0;
    //int numPTVVoxels=InfoPTV.numPTVVoxels; //this is global
    //float prescribedDose=InfoPTV.prescribedDose;
    //printf("\n in of-ptv  : %f",pbWeights[1]);
    //float dosePTV[numPTVVoxels];   //this is dose vector for PTV with number of voxels elements
    float *dosePTV=(float *)malloc(sizeof(float)*numPTVVoxels);

    for (j=0; j<numPTVVoxels; j++)
    {
        dosePTV[j]=0;
    }

    //regular multiplication
//    for (j=0;j<numPTVVoxels;j++){
//        for (i=0;i<numPBs;i++){
//            dosePTV[j] +=influence_PTV[j][i] * pbWeights[i];
//    }
//            //quadratic error normalized to the number of voxels
//            PTV_error +=pow((float)(dosePTV[j]-prescribedDose),2) ;
//    }


    //using BLASS library to do matrix * vector multiplication
    //calling the blas function for matrix-vector multiplication

    dgemv_(&trans,&numPTVVoxels,&numPBs,&alphaa,influence_PTV,&numPTVVoxels,pbWeights,&incx,&beta,dosePTV,&incy);

    for (j=0; j<numPTVVoxels; j++)
    {
        //quadratic error normalized to the number of voxels
        PTV_error +=pow((dosePTV[j]-prescribedDose),2) ;
    }


    //free dosePTV in each run
    free(dosePTV);

    return (PTV_error /(float)numPTVVoxels);
}

float rand01()
{
//generating a random number in this range [0 1]
    return (float)rand() / (float)RAND_MAX ;
}



float randn (float mu, float sigma)
{
    float U1, U2, W, mult;
    static float X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (float) X2);
    }

    do
    {
        U1 = -1 + ((float) rand () / RAND_MAX) * 2;
        U2 = -1 + ((float) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (float) X1);
}
