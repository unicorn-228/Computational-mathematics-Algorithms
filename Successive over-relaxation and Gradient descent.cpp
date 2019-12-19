#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <iomanip>

using namespace std;

double maxNorma = 0;

double** multiply(double **a, double **b, int l, int m, int n) { // L x M    M x N
	double ** result = new double*[l];
	for (int i = 0; i < l; i++) {
		result[i] = new double[n];
		for (int q = 0; q < n; q++) {
			double sum = 0;
			for (int w = 0; w < m; w++) {
				sum += a[i][w] * b[w][q];
			}
			result[i][q] = sum;
		}
	}
	return result;
}

double * multiply(double **a, double *b, int n) {
	double * result = new double[n];
	for (int i = 0; i < n; i++) {
		result[i] = 0;
		for (int q = 0; q < n; q++) {
			result[i] += a[i][q] * b[q];
		}
	}
	return result;
}

double skalyarMultiply(double *a, double * b, int n) {
	double result = 0;
	for (int i = 0; i < n; i++) {
		result += a[i] * b[i];
	}
	return result;
}


void coutMatrix(double ** a, int n) {
	cout.precision(2);
	for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			cout << setw(6) << std::fixed << a[i][q] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void copyVect(double *a, double *b, int n) {
	for (int i = 0; i < n; i++) {
		a[i] = b[i];
	}
}

void copyMatrix(double **a, double **b, int n) {
	for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			a[i][q] = b[i][q];
		}
	}
}

void coutVector(double *a, int n) {
	cout << "[ ";
	for (int i = 0; i < n; i++) {
		cout << a[i] << " ";
	}
	cout << "]" << endl;
}

double calcNorma(double *vect, int n) {
	double max = fabs(vect[0]);
	for (int i = 0; i < n; i++) {
		if (fabs(max) < fabs(vect[i])) {
			max = fabs(vect[i]);
		}
	}
	return max;
}

void coutDetA(double detA) {
	cout.precision(2);
	cout << "det A = " << fixed << detA << endl << endl;
}

double * subVector(double * a, double * b, int n) {
	double *result = new double[n];
	for (int i = 0; i < n; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}

double * multiply(double * a, double b, int n) {
	double * result = new double[n];
	for(int i = 0 ; i < n ; i++) {
		result[i] = a[i] * b;
	}
	return result;
}

void getVectorEqualTo(double * a , double * b, int n) {
	for (int i = 0; i < n; i++) {
		a[i] = b[i];
	}
}


double * solveGrad(double ** A, double * f, int n , double * X0 , double epsilon , int * num) {
	double * x = new double[n];
	double * x1 = X0;
	double * r = new double[n];
	for (int i = 0; i < n; i++) {
		x[i] = x1[i] + 10;
	}
	
	while (calcNorma(subVector(x1, x, n), n) > epsilon) {
		getVectorEqualTo(x, x1, n);
		r = subVector(multiply(A, x, n), f, n);
		x1 = subVector(x, multiply(r, skalyarMultiply(r, r, n) / skalyarMultiply(multiply(A, r, n), r, n), n) , n);
		(*num)++;
	}
	maxNorma = calcNorma(subVector(x1, x, n), n);
	return x1;
}



double * solveRelax(double ** A, double * f, int n, double * x0, double w, double epsilon , int *num) {
	double * x = new double[n];
	double * x1 = x0;
	for (int i = 0; i < n; i++) {
		x[i] = 0;
	}

	while (calcNorma(subVector(x1, x, n), n) > epsilon) {
		getVectorEqualTo(x, x1, n);
		for (int i = 0; i < n; i++) {
			double sum1 = 0, sum2 = 0;
			for (int j = 0; j < i; j++) {
				sum1 += x1[j] * A[i][j] / A[i][i];
			}
			for (int j = i+1; j < n; j++) {
				sum2 += x[j] * A[i][j] / A[i][i];
			}
			x1[i] = (1 - w)*x[i] - w * sum1 - w * sum2 + w * f[i] / A[i][i];
		}
		/*for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 1; j < n; j++) {
				if (i != j) {
					sum += A[i][j] * x[j];
				}
			}
			x1[i] = (1 - w)*x[i] + w * (f[i] - sum) / A[i][i];
		}*/
		(*num)++;
	}
	maxNorma = calcNorma(subVector(x1, x, n), n);
	return x1;
}

void createPositiveMatrix(double ** A, int n) {
	//A[0][0] = (rand() % 1000 )/ 100.;
	for (int i = 0; i < n; i++) {
		for (int q = 0; q < i; q++) {
			A[i][q] = (rand() % 2000) / 1000.;
			A[q][i] = A[i][q];
		}
	}
	for (int i = 0; i < n; i++) {
		A[i][i] = (rand() % 800) / 100. + 2;
	}
}


int main() {
	srand(time(NULL));
	cout << "matrix size = ";
	int n;
	cin >> n;
	double **a = new double*[n];
	double *x = new double[n];
	
	// далее идет генераци€ и вывод матрицы и всех необходимых векторов
	cout << "A = " << endl;
	//double ** at = new double*[n];
	for (int i = 0; i < n; i++) {
		a[i] = new double[n];
		for (int q = 0; q < n; q++) {
			a[i][q] = 0;
		}
	}
	/*for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			a[i][q] = (rand() % 1000) / 100.;
		}
	}*/
	createPositiveMatrix(a, n);
	
	
	/*for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			double temp;
			cin >> temp;
			a[i][q] = temp;
		}
	}*/
	/*
	for (int i = 0; i < n; i++) {
		//cout << "| ";
		for (int q = 0; q < n; q++) {
			a[i][q] = ((rand() % 1000)) / 100.;
			
			//cout << setw(5) << a[i][q] << " ";
		}
		//cout << "|" << endl;
	}
	for (int i = 0; i < n; i++) {
		for (int q = 0; q < i; q++) {
			a[i][q] = a[q][i];
		}
	}*/
	coutMatrix(a, n);
	cout << endl << "x = " << endl;
	cout << "[ ";
	for (int i = 0; i < n; i++) {
		x[i] = ((rand() % 1000)  + 1) / 100.;
		cout << x[i] << " ";
	}
	cout << "]" << endl;
	double *f = multiply(a, x, n);
	cout << "f = " << endl;
	coutVector(f, n);
	cout << endl << endl;
	cout << "--------------------------------------------------" << endl << endl;
	//a = multiply(at, a, n, n, n);
	//f = multiply(at, f, n);
	double * x0 = new double[n];
	for (int i = 0; i < n; i++) {
		x0[i] = f[i] / a[i][i];
	}
	// создаем копию матрицы а и вектора f , поскольку метод solve() их испортит
	double **a_copy = new double*[n];
	for (int i = 0; i < n; i++) a_copy[i] = new double[n];
	copyMatrix(a_copy, a, n);
	double *f_copy = new double[n];
	copyVect(f_copy, f, n);

	
	

	//-----------------------    –≈Ћј —ј÷»я    --------------------------//

	int * num = new int;
	*num = 0;
	int w = 2;
	cout << "Table :" << endl;
	while (w <= 18) {
		double w1 = w / 10.;
		*num = 0;
		double * fx_copy = new double[n];
		double * x0_copy = new double[n];
		for (int i = 0; i < n; i++) {
			fx_copy[i] = f[i];
			x0_copy[i] = x0[i];
		}
		solveRelax(a, fx_copy, n, x0_copy, w1, 0.00001, num);
		cout << setw(4) << w1 << setw(4) << *num << " "<< setprecision(13)<< maxNorma << setprecision(2) << endl;

		switch (w) {
		case 2:
			w = 5;
			break;
		case 5:
			w = 8;
			break;
		case 8:
			w = 10;
			break;
		case 10:
			w = 13;
			break;
		case 13:
			w = 15;
			break;
		case 15:
			w = 18;
			break;
		case 18:
			w = 228;
			break;
		}
	}
	cout << endl;
	cout << "Relax :" << endl;
	coutVector(solveRelax(a, f, n, x0, 0.5, 0.00001 , num), n);
	cout << "MaxNorma = "<< setprecision(13) << maxNorma << setprecision(2) << endl;
	cout << endl;


	//----------------------- √–јƒ»≈Ќ“Ќџ… —ѕ”— --------------------------//
	cout << "Grad :" << endl;
	*num = 0;
	coutVector(solveGrad(a, f_copy, n, f_copy, 0.00001 , num) , n);
	cout << "MaxNorma = " << setprecision(13) << maxNorma << setprecision(2) << endl;
	cout << "Kol-vo itercii = " << *num << endl;

	system("pause");
	return 0;
}