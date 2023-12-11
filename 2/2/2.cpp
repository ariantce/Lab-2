// 2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;

double e1 = 0.000000001;
double e2 = 0.000000001;

double F1(/*double* F,*/double x1, double x2) // пишем исходные ур-ния
{
	/*F[0] = pow(x1, 2) * pow(x2, 2) - 3 * pow(x1, 2) - 6 * pow(x2, 3) + 8;*/
	return x1*x1 * x2 * x2 - 3 * x1* x1 - 6 * pow(x2, 3) + 8;
}

double F2(/*double* F,*/double x1, double x2) 
{
	/*F[1] = pow(x1, 4) - 9 * pow(x2, 2) + 2;*/
	return pow(x1, 4) - 9 * x2 + 2;
}

double dF1dx1(double x1, double x2)
{
	return 2 * x1 * pow(x2, 2) - 6 * x1;
}

double dF1dx2(double x1, double x2)
{
	return 2 * pow(x1, 2) * x2 - 18 * pow(x2, 2);
}

double dF2dx1(double x1, double x2)
{
	return 4 * pow(x1, 3);
}

double dF2dx2(double x1, double x2)
{
	return  -9;
}

double* the_gauss_metod(double** A, double* b, int n)
{
    
    for (int i = 0; i < n; i++) {
        int max_idx = i;
        double max_val = abs(A[i][i]);

        // Выбор главного элемента по столбцу
        for (int j = i + 1; j < n; j++)
        {
            if (abs(A[j][i]) > max_val) {
                max_idx = j;
                max_val = abs(A[j][i]);
            }
        }

        // Перестановка строк
        if (max_idx != i)
        {
            swap(A[i], A[max_idx]);
            swap(b[i], b[max_idx]);
        }

        // Прямой ход
        for (int j = i + 1; j < n; j++)
        {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) 
            {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    double* X = new double[n];  //массив ответов
    
    // Обратный ход
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * X[j];
        }
        X[i] = (b[i] - sum) / A[i][i];
    }

    return X;
}

void Jcobian(double** J, double x1, double x2)
{
	J[0][0] = dF1dx1(x1,x2); // производная по x1 из первого ур-ния
	J[0][1] = dF1dx2(x1,x2); // производная по х2 из первого ур-ния
	J[1][0] = dF2dx1(x1,x2); // производная по х1 из второго ур-ния
	J[1][1] = dF2dx2(x1,x2); // производная из х2 из второго ур-ния

}

void Jcobian_M(double** J, double x1, double x2, double M)
{		J[0][0] = (F1(x1 + x1 * M, x2) + F1(x1, x2)) / (M * x1);
	J[0][1] = (F1(x1, x2 * M + x2) + F1(x1, x2)) / (M * x2);
	J[1][0] = (F2(x1 + x1 * M, x2) + F2(x1, x2)) / (M * x1);
	J[1][1] = (F2(x1, x2 * M + x2) + F2(x1, x2)) / (M * x2);
}

double* Newton(double x1, double x2, int NIT, double e1, double e2)
{
    int k = 1;
    cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double* F = new double[2]; 
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	double* dX = new double[2];

	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = F1(x1, x2);
		F[1] = F2(x1, x2);
		Jcobian(J,x1,x2);

		dX = the_gauss_metod(J, F, 2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = -(F1(x1, x2));
		tmp = -(F2(x1, x2));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			system("pause");
			exit(2);
		}
		k++;
	} while (d1 > e1 && d2 > e2);
	dX[0] = x1;
	dX[1] = x2;

	return dX ;
}

double* Newton_M(double x1, double x2, int NIT, double E1, double E2, double M)
{
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double* F = new double[2];
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	double* dX = new double[2];
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -F1(x1, x2);
		F[1] = -F2(x1, x2);
		Jcobian_M(J,x1,x2, M);

		dX = the_gauss_metod(J, F, 2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = -(F1(x1, x2));
		tmp = -(F2(x1, x2));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			exit(2);
		}
		k++;
	} while (d1 > E1 && d2 > E2);
	dX[0] = x1;
	dX[1] = x2;
	return dX;
}

int main()
{
	setlocale(LC_ALL, "RUS");
	int NIT = 20;
	cout << "Аналитический метод для точки (3,2): \n";
	double* x = new double [2];
	double* x1 = new double[2];
	x = Newton(3, 2, NIT, e1, e2);
	cout << "\nАналитический метод для точки (3,-2): \n";
	x1 = Newton(3, -2, NIT, e1, e2);
	cout << "\nЧерез матрицу Якоби для точки (3,2) и М=0.01: \n";
	double* x3 = new double[2], * x4 = new double[2], * x5 = new double[2];
	x3 = Newton_M(3, 2, NIT, e1, e2, 0.01);
	cout << "\nЧерез матрицу Якоби для точки (3,2) и М=0.05: \n";
	x4 = Newton_M(3, 2, NIT, e1, e2, 0.05);
	cout << "\nЧерез матрицу Якоби для точки (3,2) и М=0.1: \n";
	x5 = Newton_M(3, 2, NIT, e1, e2, 0.1);
	cout << "\nЧерез матрицу Якоби для точки (3,-2) и М=0.01: \n";
	double* x6 = new double[2], * x7 = new double[2], * x8 = new double[2];
	x6 = Newton_M(3, 2, NIT, e1, e2, 0.01);
	cout << "\nЧерез матрицу Якоби для точки (3,-2) и М=0.05: \n";
	x7 = Newton_M(3, 2, NIT, e1, e2, 0.05);
	cout << "\nЧерез матрицу Якоби для точки (3,-2) и М=0.1: \n";
	x8 = Newton_M(3, 2, NIT, e1, e2, 0.1);
	system("pause");
	return 0;
}