#include "polynomial.h"
#include "GaussElimination.h"


Polynomial PolyInterpolation(std::vector<std::pair<double, double>> values)
{
	Polynomial Res;
	int n = values.size();
	Res.a.resize(n);
	Matrix Vandermonde(n, n+1);
	Vandermonde.bAugmented = true;
	for (int i = 0; i < n; i++)
	{
		double pw = 1.0;
		for (int j = 0; j < n; j++)
		{
			Vandermonde.a[i][j] = pw;
			pw *= values[i].first;
		}
		Vandermonde.a[i][n] = values[i].second;
	}
	GaussJordanElimination(Vandermonde, false);
	for (int i = 0; i < n; i++)
		Res.a[i] = Vandermonde.a[i][n];
	return Res;
}

Polynomial PolyInterpolationLeadingOne(std::vector<std::pair<double, double>> values)
{
	int n = values.size();
	for (int i = 0; i < n; i++)
		values[i].second -= pow(values[i].first, (double)n);
	Polynomial Res=PolyInterpolation(values);
	Res.a.push_back(1);
	return Res;
}

double Polynomial::EvaluateAt(double x)
{
	double pw = 1.0;
	double res = 0.0;
	for (int i = 0; i < a.size(); i++)
	{
		res += a[i] * pw;
		pw *= x;
	}
	return res;
}

void Polynomial::Display()
{
	for (int i = a.size() - 1; i >= 0; i--)
	{
		if (IsZero(a[i]))
			continue;
		if (a[i] > 0 && i < a.size() - 1)
			std::cout << " + ";
		else if (a[i] < 0)
			std::cout << " - ";
		if (!IsZero(a[i] - 1.0) && !IsZero(1.0 + a[i]))
			std::cout << double_to_string(std::abs(a[i]), double_precision);
		if (i > 1)
			std::cout << "x^" << i;
		else if (i == 1)
			std::cout << "x";
	}
	std::cout << "\n";
}
