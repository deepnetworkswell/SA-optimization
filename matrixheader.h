
double **readMatrix(int rows,int cols)
{
    //reading a binary file for a saved matrix with rows*cols ,elements are written linearly
    int i,j;
    FILE *fp;
    fp=fopen("doubledata.bin","rb");
    if (fp==NULL)
    {
        printf("The file was not opened!");
    }
    else
    {
        double **mat=(double **)malloc(sizeof(double *)*rows);
        for (i=0; i<rows; i++)
        {
            mat[i]=(double *)malloc(sizeof(double)*cols);
        }
        for (i=0; i<rows; i++)
        {
            for (j=0; j<cols; j++)
            {
                fread(&mat[i][j], sizeof(double), 1, fp);
            }
        }
        return mat;
    }
}

void printMatrix(int rows,int cols,double **mat)
{
    int i,j;
    for (i=0; i<rows; i++)
    {
        for (j=0; j<cols; j++)
        {
            printf("%f ",mat[i][j]);
        }
        printf("\n");
    }
}

float **convertDoubleToFloat(int rows,int cols,double **mat)
{     //converts double array to float array  I am using 32bit pc
    int i,j;
    float **matFloat=(float **)malloc(sizeof(float *)*rows);
    for (i=0; i<rows; i++)
    {
        matFloat[i]=(float *)malloc(sizeof(float)*cols);
    }
    for (i=0; i<rows; i++)
    {
        for (j=0; j<cols; j++)
        {
            matFloat[i][j]=(float)mat[i][j];
        }
    }
    return matFloat;
}

void freeMatrix(int rows,double **mat)
{
    int i;
    for (i=0; i<rows; i++)
        free(mat[i]);
    free(mat);
}
