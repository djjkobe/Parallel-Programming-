/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Esteban Meneses, PhD (emeneses@pitt.edu)
 * Student:Jiajie Dang
 * Cilk Plus parallel Strassen algorithm for matrix multiplication.
 */

#include <cstdio>
#include <cstdlib>
#include <cilk/cilk.h>
//#include "cilk.h"
#include "timer.h"
#include "io.h"
#include <math.h>


#define Bound  1024
//void Add(int **A,int **B, int**Result,int N);
//void Subtract(int **A,int **B, int **Result,int N);
//void StrassenAlgorithm(int **A,int **B,int **C,int N);


template < typename T >
T **Allocate2DArray( int nRows, int nCols)
{
    
    T **ppi = new T*[nRows];
    
    
    T *curPtr = new T [nRows * nCols];
    
   
    for( int i = 0; i < nRows; ++i)
    {
        *(ppi + i) = curPtr;
        curPtr += nCols;
    }
    return ppi;
}

template < typename T >
void Free2DArray(T** Array)
{
    delete [] *Array;
    delete [] Array;
}


void Add(int **T, int m, int n, int **X,int **Y)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			T[i][j] = X[i][j] + Y[i][j];
}

void Sub(int **T, int m, int n, int **X, int **Y)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			T[i][j] = X[i][j] - Y[i][j];
}
void matmultleaf(int mf, int ml, int nf, int nl, int pf, int pl, int **A, int **B, int **C)
{
	for(int i = mf; i < ml; i++)
		for (int j = nf; j < nl; j++) {
			C[i][j] = 0;
			for(int k = pf; k < pl; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
}



void copy(int **X, int m, int **Y, int mf, int nf)
{
	for (int i = 0; i < m; i++)
		X[i] = &Y[mf+i][nf];
}

void Strassen(int mf, int ml, int nf, int nl, int pf, int pl, int **A, int **B, int **C)
{
    int mdiff = ml-mf;
	int ndiff = nl-nf;
	int pdiff = pl-pf;
	if (mdiff==1 || ndiff==1 || pdiff==1 ||
	    mdiff*ndiff*pdiff < Bound)
		matmultleaf(mf, ml, nf, nl, pf, pl, A, B, C);
   
    
    else{
        int m2 = mdiff/2;
		int n2 = ndiff/2;
		int p2 = pdiff/2;
        int **M1 = Allocate2DArray< int >(m2, n2);
		int **M2 = Allocate2DArray< int >(m2, n2);
		int **M3 = Allocate2DArray< int >(m2, n2);
		int **M4 = Allocate2DArray< int >(m2, n2);
		int **M5 = Allocate2DArray< int >(m2, n2);
		int **M6 = Allocate2DArray< int >(m2, n2);
		int **M7 = Allocate2DArray< int >(m2, n2);
        
		int **A11 = new int*[m2];
		int **A12 = new int*[m2];
		int **A21 = new int*[m2];
		int **A22 = new int*[m2];
        
		int **B11 = new int*[p2];
		int **B12 = new int*[p2];
		int **B21 = new int*[p2];
		int **B22 = new int*[p2];
        
		int **C11 = new int*[m2];
		int **C12 = new int*[m2];
		int **C21 = new int*[m2];
		int **C22 = new int*[m2];
        
		int **tAM1 = Allocate2DArray< int >(m2, p2);
		int **tBM1 = Allocate2DArray< int >(p2, n2);
		int **tAM2 = Allocate2DArray< int >(m2, p2);
		int **tBM3 = Allocate2DArray< int >(p2, n2);
		int **tBM4 = Allocate2DArray< int >(p2, n2);
		int **tAM5 = Allocate2DArray< int >(m2, p2);
		int **tAM6 = Allocate2DArray< int >(m2, p2);
		int **tBM6 = Allocate2DArray< int >(p2, n2);
		int **tAM7 = Allocate2DArray< int >(m2, p2);
		int **tBM7 = Allocate2DArray< int >(p2, n2);
        
		copy(A11, m2, A, mf, pf);
		copy(A12, m2, A, mf, p2);
		copy(A21, m2, A, m2, pf);
		copy(A22, m2, A, m2, p2);
        
		copy(B11, p2, B, pf, nf);
		copy(B12, p2, B, pf, n2);
		copy(B21, p2, B, p2, nf);
		copy(B22, p2, B, p2, n2);
        
		copy(C11, m2, C, mf, nf);
		copy(C12, m2, C, mf, n2);
		copy(C21, m2, C, m2, nf);
		copy(C22, m2, C, m2, n2);
        
        // M1 = (A11 + A22)*(B11 + B22) (addition part)
        cilk_spawn Add(tAM1, m2, p2, A11, A22);
        cilk_spawn Add(tBM1, p2, n2, B11, B22);
        
        //M2 = (A21 + A22)*B11 (addition part);
        cilk_spawn Add(tAM2, m2, p2, A21, A22);
        
        //M3 = A11*(B12 - B22) (addition part)
        cilk_spawn Sub(tBM3, p2, n2, B12, B22);
        
        //M4 = A22*(B21 - B11) (addition part)
        cilk_spawn Sub(tBM4, p2, n2, B21, B11);
        
        //M6 = (A21 - A11)*(B11 + B12) (addition part)
        cilk_spawn Sub(tAM6, m2, p2, A21, A11);
        cilk_spawn Add(tBM6, p2, n2, B11, B12);
        
        //M5 = (A11 + A12)*B22 (addition part)
        cilk_spawn Add(tAM5, m2, p2, A11, A12);
        
        //M7 = (A12 - A22)*(B21 + B22) (addition part)
        cilk_spawn Sub(tAM7, m2, p2, A12, A22);
        Add(tBM7, p2, n2, B21, B22);
        
        cilk_sync;
        // M1 = (A11 + A22)*(B11 + B22)
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, tAM1, tBM1, M1);
        
        //M2 = (A21 + A22)*B11
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, tAM2, B11, M2);
        
        //M3 = A11*(B12 - B22)
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, A11, tBM3, M3);
        
        //M4 = A22*(B21 - B11)
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, A22, tBM4, M4);
        
        //M5 = (A11 + A12)*B22
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, tAM5, B22, M5);
        
        //M6 = (A21 - A11)*(B11 + B12)
        cilk_spawn Strassen(0, m2, 0, n2, 0, p2, tAM6, tBM6, M6);
        
        //M7 = (A12 - A22)*(B21 + B22)
        Strassen(0, m2, 0, n2, 0, p2, tAM7, tBM7, M7);
        
        cilk_sync;

		for(int i = 0; i < m2; i++)
			for (int j = 0; j < n2; j++) {
				C11[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
				C12[i][j] = M3[i][j] + M5[i][j];
				C21[i][j] = M2[i][j] + M4[i][j];
				C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
			}
        
		Free2DArray< int >(M1);
		Free2DArray< int >(M2);
		Free2DArray< int >(M3);
		Free2DArray< int >(M4);
		Free2DArray< int >(M5);
		Free2DArray< int >(M6);
		Free2DArray< int >(M7);
        
		delete[] A11; delete[] A12; delete[] A21; delete[] A22;
		delete[] B11; delete[] B12; delete[] B21; delete[] B22;
		delete[] C11; delete[] C12; delete[] C21; delete[] C22;
        
		Free2DArray< int >(tAM1);
		Free2DArray< int >(tBM1);
		Free2DArray< int >(tAM2);
		Free2DArray< int >(tBM3);
		Free2DArray< int >(tBM4);
		Free2DArray< int >(tAM5);
		Free2DArray< int >(tAM6);
		Free2DArray< int >(tBM6);
		Free2DArray< int >(tAM7);
		Free2DArray< int >(tBM7);
	}
    
}
void before(int **A, int **B,int **C,int M)
{
    int i,j;
    for (i=0; i < M; i++)
        for (j=0; j < M; j++)
            C[i][j] = 0;
    Strassen(0, M, 0, M, 0, M, A,B,C);
}

// Main method
int main(int argc, char* argv[]) {
    int N;
    
	int **A, **B, **C;
	double elapsedTime;
    
	// checking parameters
	if (argc != 2 && argc != 4) {
		cout << "Parameters: <N> [<fileA> <fileB>]" << endl;
		return 1;
	}
	N = atoi(argv[1]);
    
  
    //allocating matrices
    A = new int*[N];
	B = new int*[N];
	C = new int*[N];
   // C4 = new int*[N];
	for (int i=0; i<N; i++){
		A[i] = new int[N];
		B[i] = new int[N];
		C[i] = new int[N];
        //C4[i] = new int[N];
	}
   	
	
	// reading files (optional)
	if(argc == 4){
		readMatrixFile(A,N,argv[2]);
		readMatrixFile(B,N,argv[3]);
	}
	
	if(argc==2){
        
		srand(time(NULL));
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				A[i][j]=rand()%10;;
                
				B[i][j]=rand()%10;
			}
		}
        
		
	}
    
   
	// starting timer
    timerStart();
    
	// YOUR CODE GOES HERE
	before(A,B,C,N);
	
	// testing the results is correct
	if(argc == 4){
        
		printMatrix(C,N);
	}
    
    if(argc == 2){
        printMatrix(C, N);
    }
	
	// stopping timer
	elapsedTime = timerStop();
    
	cout << "Duration: " << elapsedTime << " seconds" << std::endl;
    
	// releasing memory
	for (int i=0; i<N; i++) {
		delete [] A[i];
		delete [] B[i];
		delete [] C[i];
	}
	delete [] A;
	delete [] B;
	delete [] C;
    
	return 0;	
}

