// 2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include<iostream>
#include<iomanip>
#include <math.h>

using namespace std;

const double eps = 1e-9;
const double M = 0.001;
const int NIT = 110;

double** initial(int n, int m)
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[m];
	return A;
}

double F1(double x1, double x2)
{
	return (x1 * x1 * x2 * x2 - 3 * x1 * x1 - 6 * x2 * x2 * x2 + 8);
}

double F2(double x1, double x2)
{
	return (x1 * x1 * x1 * x1 - 9 * x2 + 2);
}


double func11(double x1, double x2)
{
	return ((F1(x1 + M * x1, x2) - F1(x1, x2)) / (M * x1));
}


double func12(double x1, double x2)
{
	return ((F1(x1, x2 + M * x2) - F1(x1, x2)) / (M * x2));
}



double func21(double x1, double x2)
{
	return ((F2(x1 + M * x1, x2) - F2(x1, x2)) / (M * x1));
}


double func22(double x1, double x2)
{
	return ((F2(x1, x2 + M * x2) - F2(x1, x2)) / (M * x2));
}


void J(double** a, double x1, double x2)
{
	a[0][0] = func11(x1, x2);
	a[0][1] = func12(x1, x2);
	a[1][0] = func21(x1, x2);
	a[1][1] = func22(x1, x2);
}

void minus_vector_nevjazki(double* F, double x1, double x2)
{
	F[0] = -F1(x1, x2);
	F[1] = -F2(x1, x2);
}

double* the_gauss_metod(double** A, int n, int m)
{
	//prjamoj
	double elem;
	for (int i = 0; i < n; i++)
	{
		double max = 0;
		int coord_str = 0;
		for (int j = i; j < n; j++)
		{
			if (abs(A[j][i]) > max)
			{
				max = abs(A[j][i]); coord_str = j;
			}
		}
		if (max > abs(A[i][i]))
		{
			double* ptr = A[i];
			A[i] = A[coord_str];
		     A[coord_str] = ptr;
		}
		elem = A[i][i];
		for (int c = i; c < m; c++)
		{
			A[i][c] /= elem;   //delenije stroki na elem
		}

		for (int j = i + 1; j < n; j++)
		{
			elem = A[j][i];
			for (int k = j; k < m; k++)
				A[j][k] -= elem * A[i][k];
		}

	}
	//obratnyj
	double* xx = new double[m];
	xx[n - 1] = A[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		xx[i] = A[i][n];
		for (int j = i + 1; j < n; j++) xx[i] -= A[i][j] * xx[j];
	}

	cout << endl;

	return xx;
}

double* neutone(int n, double x1, double x2)
{
	double** Jako;
	Jako = initial(n, n);

	double** newmatrix;
	newmatrix = initial(n, n + 1);

	double* F = new double[n];
	double* delta = new double[n];
	double* resh = new double[n];
	resh[0] = x1;
	resh[1] = x2;
	double delta1, delta2;
	int k = 1;

	cout << setw(10) << "x1" << setw(10) << "x2" << setw(17) << "delta1" << setw(17) << "delta2" << setw(6) << "k";
	do
	{
		minus_vector_nevjazki(F, x1, x2);
		J(Jako, x1, x2);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				newmatrix[i][j] = Jako[i][j];
			}
			newmatrix[i][n] = F[i];
		}

		delta = the_gauss_metod(newmatrix, n, n + 1);
		for (int i = 0; i < n; i++)
			resh[i] += delta[i];

		double max1 = 0;
		double max2 = 0;
		for (int i = 0; i < n; i++)
		{
			if (abs(F[i]) > max1)
				max1 = abs(F[i]);

			if (abs(resh[i]) < 1)
			{
				if (abs(delta[i]) > max2)
					max2 = abs(delta[i]);
			}
			if (abs(delta[i] >= 1))
			{
				if (abs(delta[i] / resh[i]) > max2)
					max2 = abs(delta[i]);
			}
		}

		delta1 = max1;
		delta2 = max2;
		cout << endl;


		x1 = resh[0];
		x2 = resh[1];
		for (int i = 0; i < n; i++)
			cout << setw(10) << resh[i] << "   ";
		cout << setw(13) << delta1 << "   " << setw(13) << delta2 << "   " << setw(2) << k;
		cout << endl;

		k++;
		if (k >= NIT)
		{
			cout << "\n   IER = 2 \n";
			return NULL;
			break;
		}

	} while (delta1 > eps || delta2 > eps);

	return resh;
}


int main()
{
	double x1, x2;
	x1 = 0.5; x2 = 0.2;
	int n = 2;
	double* otvet = new double[n];
	otvet = neutone(n, x1, x2);

	if (otvet != NULL)
	{
		cout << "____________________________________________________________" << endl;
		for (int i = 0; i < n; i++)
			cout << otvet[i] << endl;
	}

	return 0;
}