#include <math.h>
#include "size.h"

void printMDArray(int, int, double[][Max_n + 1]);

void Invert(int n, double C[][Max_n + 1], double B[][Max_n + 1])
{
    double A[Max_n + 1][Max_n + 1];

    for (int i = 0; i <= n; i++) {
	for (int j = 0; j <= n; j++)
	    A[i][j] = C[i][j];	//This is so that C[i][j] is returned unchanged
    }

    for (int i = 0; i <= n; i++) {
	for (int j = 0; j <= n; j++)
	    B[i][j] = 0.0;
	B[i][i] = 1.0;
    }
 
    for (int i = 0; i <= n; i++) {
	if (i != n) {
	    int biggest = i;
	    for (int k = i + 1; k <= n; k++) {
		if (fabs(A[k][i]) > fabs(A[biggest][i]))
		    biggest = k;
	    }
	    for (int j = 0; j <= n; j++) {
		double tempdouble;
		tempdouble = A[i][j];
		A[i][j] = A[biggest][j];
		A[biggest][j] = tempdouble;
		tempdouble = B[i][j];
		B[i][j] = B[biggest][j];
		B[biggest][j] = tempdouble;
	    }
	}

	double divide = A[i][i];
	for (int j = 0; j <= n; j++) {
	    A[i][j] = A[i][j] / divide;
	    B[i][j] = B[i][j] / divide;
	}

	if (i != n) {
	    for (int k = i + 1; k <= n; k++) {
		double mult = A[k][i];
		for (int j = 0; j <= n; j++) {
		    A[k][j] -= mult * A[i][j];
		    B[k][j] -= mult * B[i][j];
		}
	    }
	}
    }

    for (int i = n; i >= 1; i -= 1) {
	for (int k = i - 1; k >= 0; k -= 1) {
	    double mult = A[k][i];
	    for (int j = 0; j <= n; j++) {
		A[k][j] -= mult * A[i][j];
		B[k][j] -= mult * B[i][j];
	    }
	}
    }
}
