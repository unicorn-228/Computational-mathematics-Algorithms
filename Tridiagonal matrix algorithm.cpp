#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>

using namespace std;


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

int getMaxIndex(double *a, int n, int i_) {
	int result = 0;
	double max = a[0];
	for (int i = i_ + 1; i < n; i++) {
		if (max < a[i]) {
			max = a[i];
			result = i;
		}
	}
	return result;
}

void exchangeColumns(double**a, int n, int l, int s) {
	for (int i = 0; i < n; i++) {
		double temp = a[i][l];
		a[i][l] = a[i][s];
		a[i][l] = temp;
	}
}

void subtractionOfRow(double *a, double*b, int n) {
	for (int i = 0; i < n; i++) {
		a[i] -= b[i];
	}
}

void coutMatrix(double ** a, int n) {
	for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			cout << setw(6)  << a[i][q] << " ";
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

double * solve(double **a, double *f, int n, vector<double>* keyData, int* numOfMoves) {
	
	vector<double> y;
	vector<double> alpha;
	vector<double> beta;
	double * result = new double[n];
	y.push_back(a[0][0]);
	alpha.push_back(-a[0][1]/y[0]);
	beta.push_back(f[0] / y[0]);
	for (int i = 1; i < n - 1; i++) {
		y.push_back(a[i][i] + a[i][i-1] * alpha[i - 1]);
		alpha.push_back(-a[i][i + 1] / y[i]);
		beta.push_back((f[i]-a[i][i-1]*beta[i-1])/y[i]);
	}
	y.push_back(a[n-1][n-1] + a[n-1][n-2]*alpha[n-2]);
	beta.push_back((f[n-1] - a[n-1][n-2]*beta[n-2])/y[n-1]);
	result[n - 1] = beta[n - 1];
	for (int i = n-2; i >= 0; i--) {
		result[i] = alpha[i] * result[i + 1] + beta[i];
	}
	return result;
}

int main() {
	srand(time(NULL));
	int n;
	cin >> n;
	double **a = new double*[n];
	double *x = new double[n];
	vector<double> keyData; //чтобы зранить ключевые элементы
	int numOfMoves = 0; //кол-во перестановок
	// далее идет генерация и вывод матрицы и всех необходимых векторов
	cout << "A = " << endl;
	for (int i = 0; i < n; i++) { a[i] = new double[n]; for (int q = 0; q < n; q++)a[i][q] = 0; }
	for (int i = 0; i < n; i++) {
		a[i][i] = ((rand() % 2000) - 1000) / 100.;
		a[i][i + 1] = ((rand() % 2000) - 1000) / 100.;
		if (i + 1 <= n - 1)
			a[i + 1][i] = ((rand() % 2000) - 1000) / 100.;
	}
	coutMatrix(a, n);
	cout << endl << "x = " << endl;
	cout << "[ ";
	for (int i = 0; i < n; i++) {
		x[i] = ((rand() % 2000) - 1000) / 100.;
		cout << x[i] << " ";
	}
	cout << "]" << endl;
	double *f = multiply(a, x, n);
	cout << "f = " << endl;
	coutVector(f, n);
	cout << endl << endl;
	cout << "--------------------------------------------------" << endl << endl;

	// создаем копию матрицы а и вектора f , поскольку метод solve() их испортит
	double **a_copy = new double*[n];
	for (int i = 0; i < n; i++) a_copy[i] = new double[n];
	copyMatrix(a_copy, a, n);
	double *f_copy = new double[n];
	copyVect(f_copy, f, n);

	//----------------------- метод прогонки

	double *result = new double[n];
	result = solve(a, f, n, &keyData, &numOfMoves);

	//------------------------------------


	cout << "x_calculated = "; // вывод посчитанного x
	coutVector(result, n);
	copyMatrix(a, a_copy, n); //восстанаваливаем матрицу и вектор
	copyVect(f, f_copy, n);
	double *vect = new double[n]; // вектор для котрого будем сичтать норму невязки
	copyVect(vect, multiply(a, result, n), n); // vect = A*x_посчитанное
	for (int i = 0; i < n; i++) { // vect = A*x_посчитанное - f
		vect[i] -= f[i];
	}
	//  вывод макс нормы невязки
	cout << endl << "||A*x_calc-f|| = "  << calcNorma(vect, n) << endl << endl;

	// норма х_посчитанное-х
	copyVect(vect, x, n);
	for (int i = 0; i < n; i++) {
		vect[i] -= result[i];
	}
	cout << "||x-x_calc|| = "  << calcNorma(vect, n) << endl << endl;
	system("pause");
	return 0;
}