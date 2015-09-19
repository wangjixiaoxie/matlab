#include <stdio.h>
#include "mex.h"
#include "swap.h"

void mexFunction(int nlhs,  mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* Reads a song file and returns two elements, the number of points
	and the data */
{
char fname[80];
char adstring[80];
FILE *f_in;
double *psamp;
double *pnframes;
double *psong;
int len, i;
short *buffer;

    /* Check for input and output */
	 if ( nlhs != 3 )
		 mexErrMsgTxt("Three output arguments are required.");
	 if (nrhs != 1)
		 mexErrMsgTxt("read_song requires name of song file.");
	 if ( ! mxIsChar(prhs[0]) )
		 mexErrMsgTxt("read_song requires name of song file.");

    mxGetString(prhs[0],fname,80);

	 /* Open input file */
	 if ( (f_in = fopen(fname,"r")) == NULL )
		 mexErrMsgTxt("Cannot open song data file\n");

    /* Allocate two scalars */
/*    if ( plhs[0] != NULL ) mxDestroyArray(plhs[0]);
    if ( plhs[1] != NULL ) mxDestroyArray(plhs[1]); */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    psamp = mxGetPr(plhs[0]);
    pnframes = mxGetPr(plhs[1]);
    
    /* Get sampling rate */
	 fgets(adstring,80,f_in);
	 if ( sscanf(adstring,"AD_FREQ: %lf Hz\n", psamp) != 1) 
    {
		 rewind(f_in);
		 *psamp=20000.0;    /* Default for backward compatability
             * with masscomp song files */
	 }

    /* Get length of file */
	 if ((fread(&len, sizeof(int), 1, f_in)) == 0)
       mexErrMsgTxt("Error reading length of file.");
    swap_s32bit(&len);
	 *pnframes=(double)len;
    printf("Length of file is %d\n",len);
    
    /* Readin data in temp short array */
    buffer = (short *)mxCalloc(len,sizeof(short));
	 printf("Buffer pointer is %x\n",buffer);
    if ((fread((char *)buffer, sizeof(short), len, f_in)) != len) 
    {
	    fclose(f_in);
       mxFree((void *)buffer);
/*       mxDestroyArray(plhs[0]);
       mxDestroyArray(plhs[1]); */
       mexErrMsgTxt("Error reading songfile.");
	 }
	 fclose(f_in);

    plhs[2] = mxCreateDoubleMatrix(1, len, mxREAL);
    psong = mxGetPr(plhs[2]);
    for ( i=-1; ++i < len; ) 
	 {
		swap_s16bit(buffer+i);
		psong[i] = (double)buffer[i];
	 }
    mxFree((void *)buffer);

}
