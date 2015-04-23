/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Esteban Meneses, PhD (emeneses@pitt.edu)
 * OpenMP parallel determinant computation.
 */

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include "timer.h"
#include "io.h"

#define MAX_VALUE 1000

int determinant(int **A, int N);

// Main method      
int main(int argc, char* argv[]) {
	int N, det;
	int **A;
	double elapsedTime;

	// checking parameters
	if (argc != 2 && argc != 3) {
		cout << "Parameters: <N> [<file>]" << endl;
		return 1;
	}
	N = atoi(argv[1]);

	// allocating matrix A
	A = new int*[N];
	for (int i=0; i<N; i++){
		A[i] = new int[N];
	}

	// reading files (optional)
	if(argc == 3){
		readMatrixFile(A,N,argv[2]);
	} else {
		srand48(time(NULL));
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				A[i][j] = lrand48() % MAX_VALUE;
			}
		}
	}
	
	// starting timer
	timerStart();

	// calling determinant function
	det = determinant(A,N);

	// testing the results is correct
	if(argc == 3){
		cout << det << endl;
	}
	
	// stopping timer
	elapsedTime = timerStop();

	cout << "Duration: " << elapsedTime << " seconds" << std::endl;

	// releasing memory
	for (int i=0; i<N; i++) {
		delete [] A[i];
	}
	delete [] A;

	return 0;	
}
