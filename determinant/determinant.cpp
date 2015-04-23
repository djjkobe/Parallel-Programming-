/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor: Esteban Meneses, PhD (emeneses@pitt.edu)
 * OpenMP parallel determinant computation. 
 */

#include <math.h>
#include <set>
#include <iostream>

using namespace std;

int determinant_minor(int **A, int N, int i, set<int> cols);

// Determinant function
int determinant(int **A, int N){
	int result = 0, sign;

	#pragma omp parallel for private(sign) reduction(+:result)
	for(int j=0; j<N; j++) {
		sign = (int) pow(-1,j);
		set<int> cols;
		cols.insert(j);
		result += sign * A[0][j] * determinant_minor(A,N,1,cols);
	}
	
	return result;
}

// Computes determinant of minor(i,j) of matrix A using Laplace expansion.
int determinant_minor(int **A, int N, int row, set<int> cols){
	int result = 0, sign, index = 0;

	for(int j=0; j<N; j++) {
		if(cols.find(j) == cols.end()){		
			if(row == N-1) return A[row][j];
			sign = (int) pow(-1,index++);
			set<int> next_cols(cols);
			next_cols.insert(j);
			result += sign * A[row][j] * determinant_minor(A,N,row+1,next_cols);
		}
	}

	return result;	
}

