#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <mpi.h>

#define ROOT 0
#define MAXNGPU 3
#define MAX_THREAD 512
#define M 200

#define BASE 5
#define POWER 13

//nvcc -I/usr/local/openmpi-4.1.4/include -L/usr/local/openmpi-4.1.4/lib -lmpi testColl.cu -o testC
// mpirun -np 1 testC 500 400 2

void cudaErrorCheck(cudaError_t error, const char * msg){
   	 if ( error != cudaSuccess){
   	 fprintf(stderr, "%s:%s\n ", msg, cudaGetErrorString(error));
   	 exit(EXIT_FAILURE);}}

__host__ __device__ uint32_t abcFunct(uint32_t ua, uint32_t ub, uint32_t uc){
    	int a=ua;
    	int b=ub;
    	int c=uc;
    	uint32_t F=c+pow(a-b,2);
    	return F;
}


__host__ __device__ uint32_t baseFunct(uint32_t x, int r){
    r=r%(POWER+1);
    	uint32_t *arrayEquivX=(uint32_t*)malloc(sizeof(int)*(POWER+r));
    	uint32_t exp= pow(BASE, POWER);
    	uint32_t remnant;
    	uint32_t y;

    	remnant=x;

    	for (int i=0; i<r; i++) arrayEquivX[POWER+i]=0;

    	for (int i=0; i<POWER; i++){
            	exp=exp/BASE;
            	arrayEquivX[POWER-i-1]=remnant/exp;
            	remnant=remnant%exp;
            	//printf("in pos %d (exp %d) there is %d with remainder %d \n", POWER-1-i, exp, arrayEquivX[POWER-i-1],remnant);
    	}

    	y=0;

    	for (int i=0; i<POWER; i++){
            	//arrayEquivX[i]= abcFunct(arrayEquivX[i], arrayEquivX[i+1], arrayEquivX[i+j]);
            	//printf("(using %d) in pos %d there is %d \n", arrayEquivX[i+r], i, arrayEquivX[i]);
            	y=y+exp*abcFunct(arrayEquivX[i], arrayEquivX[i+1], arrayEquivX[i+r]);
            	exp=exp*BASE;
    	}
    	//printf("modulo is %d \n", exp);
    	y=y%exp;
    	return y;
}


    
int main (int argc, char** argv) {

//input validation
if(argc != 4){
    fprintf(stderr,"wrong number of inputs\n");
    return EXIT_FAILURE;}

uint32_t a=atoi(argv[1]);

if(a <=0){
   	 fprintf(stderr,"[ERROR] - lg must be > 0\n");
   	 return EXIT_FAILURE;}

//uint32_t ua=a;

uint32_t b=atoi(argv[1]);

if(b <=0){
    	fprintf(stderr,"[ERROR] - lg must be > 0\n");
    	return EXIT_FAILURE;}

int r=atoi(argv[3]);

uint32_t FA=baseFunct(a,r);
printf("F of %d is %d\n", a, FA);

uint32_t FB=baseFunct(b,r);
printf("F of %d is %d\n", b, FB);

if (FA==FB) printf("Collision confirmed on %d \n", FA);

}
