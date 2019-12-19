#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>

using namespace std;


double** multiply(double **a , double **b , int l , int m , int n) { // L x M    M x N
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

int getMaxIndex(double *a , int n , int i_) {
	int result = 0;
	double max = a[0];
	for (int i = i_+1; i < n; i++) {
		if (max < a[i]) {
			max = a[i];
			result = i;
		}
	}
	return result;
}

void exchangeColumns(double**a, int n, int l , int s) {
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

double calcNorma(double *vect , int n) {
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

double * solve(double **a , double *f , int n , vector<double>* keyData , int* numOfMoves) {	
	for (int i = 0; i < n; i++) {
		int maxIndex = i;
		double max = a[i][i];
		for (int q = i + 1; q < n; q++) {
			if (fabs(max) < fabs(a[i][q])) {
				max = a[i][q];
				maxIndex = q;
			}
		}
		if (maxIndex != i) {
			exchangeColumns(a, n, i, maxIndex);
			numOfMoves++;
			keyData->push_back(a[i][maxIndex]);
		}
		else {
			keyData->push_back(a[i][i]);
		}
		for (int q = i + 1; q < n; q++) { // ����� �� ������ �������
			a[i][q] *= (1 / a[i][i]);
		}
		f[i] *= (1 / a[i][i]);
		a[i][i] = 1;
		for (int q = i + 1; q < n; q++) { // �������� ������ �� ���� ���������
			double temp = a[q][i];
			for (int w = i; w < n; w++) {
				a[q][w] -= a[i][w] * temp;
			}
			f[q] -= f[i] * temp;
		}
	}
	//----------�������� ���:
	double *result = new double[n];
	for (int i = n - 1; i >= 0; i--) {
		result[i] = f[i];
		for (int q = n - 1; q > i; q--) {
			result[i] -= a[i][q] * result[q];
		}
	}
	return result;
}

int main() {
	srand(time(NULL));
	cout << "matrix size = ";
	int n;
	cin >> n;
	double **a = new double*[n];
	double *x = new double[n];
	vector<double> keyData; //����� ������� �������� ��������
	int numOfMoves = 0; //���-�� ������������
	// ����� ���� ��������� � ����� ������� � ���� ����������� ��������
	cout << "A = " << endl;
	for (int i = 0; i < n; i++) {
		a[i] = new double[n];
		cout << "| ";
		for (int q = 0; q < n; q++) {
			a[i][q] = ((rand() % 2000) - 1000) / 100.;
			cout << setw(5) << a[i][q] << " ";
		}
		cout << "|" << endl;
	}
	cout << endl << "x = " << endl;
	cout << "[ ";
	for (int i = 0; i < n; i++) {
		x[i] = ((rand() % 2000) - 1000) / 100.;
		cout << x[i] << " ";
	}
	cout << "]" << endl;
	double *f = multiply(a, x, n);
	cout << "f = " <<endl;
	coutVector(f, n);
	cout <<endl <<endl;
	cout << "--------------------------------------------------" << endl << endl;

	// ������� ����� ������� � � ������� f , ��������� ����� solve() �� ��������
	double **a_copy = new double*[n];
	for (int i = 0; i < n; i++) a_copy[i] = new double[n];
	copyMatrix(a_copy, a , n);
	double *f_copy = new double[n]; 
	copyVect(f_copy, f, n);

	//-----------------------  ����� ������ ---------------|
	//													   |
	double *result = new double[n];//					   |
	result = solve(a, f, n, &keyData, &numOfMoves);//      |
	//													   |
	//-----------------------------------------------------|

	cout << "x_calculated = "; // ����� ������������ x
	coutVector(result, n);
	copyMatrix(a, a_copy, n); //���������������� ������� � ������
	copyVect(f, f_copy, n);
	double *vect = new double[n]; // ������ ��� ������� ����� ������� ����� �������
	copyVect(vect, multiply(a, result, n), n); // vect = A*x_�����������
	for (int i = 0; i < n; i++) { // vect = A*x_����������� - f
		vect[i] -= f[i];
	}
	//  ����� ���� ����� �������
	cout <<endl << "||A*x_calc-f|| = " << calcNorma(vect , n) << endl <<endl;

	// ����� �_�����������-�
	copyVect(vect, x , n);
	for (int i = 0; i < n; i++) {
		vect[i] -= result[i];
	}
	cout << "||x-x_calc|| = " << calcNorma(vect, n) << endl << endl;

	//����� ���� ������� ������������
	double detA = 1;
	for (int i = 0; i < keyData.size(); i++) {
		detA *= keyData[i];
	}
	detA *= pow(-1., (double)numOfMoves);
	coutDetA(detA); // �������, ������� ������� ������������ � ������ ���������
	
	double * E = new double[n]; // ��������� ������ ��� ���������� ��������
	double ** a_Mirror = new double*[n]; // �������� �������
	for (int i = 0; i < n; i++)a_Mirror[i] = new double[n];
	for (int i = 0; i < n; i++) {
		for (int w = 0; w < n; w++) { // ������� ������� ��������� ������ � ���������� � �
			E[w] = (w == i) ? 1 : 0;
		}
		result = solve(a, E, n, &keyData, &numOfMoves); //������
		copyMatrix(a, a_copy, n); // ���������������� �������
		copyVect(f, f_copy, n); // � ������
		for (int q = 0; q < n; q++) { // ���������� � ������ �������� ������� �������
			a_Mirror[q][i] = result[q];
		}
	}
	cout << "mirror matrix = " <<endl;
	coutMatrix(a_Mirror , n);
	cout << "A*A_mirror = " << endl;
	double **AA_mirror = multiply(a, a_Mirror, n, n, n);
	coutMatrix(AA_mirror, n);
	system("pause");
	return 0;
}