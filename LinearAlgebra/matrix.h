#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>

const double zero_error = 1e-6;
const int double_precision = 2;
bool IsZero(double x);


class Matrix
{
private:
	int m;
	int n;
public:
	std::vector<std::vector<double>> a;
	std::vector<char> variableNames;
	bool bAugmented=false;
	Matrix(int m0, int n0) {
		m = m0; 
		n = n0;
		a.resize(m);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				a[i].push_back(0.0);
		}
	};
	int RowNum() const { return m; };
	int ColNum() const { return n; };
	void Display();
	Matrix operator+(const Matrix& mx) const
	{
		if (m != mx.RowNum() || n != mx.ColNum())
			throw std::invalid_argument("Matrix addition failed, dimension mismatch.");
		Matrix Res(m, n);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				Res.a[i][j] = a[i][j] + mx.a[i][j];
		}
		return Res;
	}
	Matrix operator-(const Matrix& mx) const
	{
		if (m != mx.RowNum() || n != mx.ColNum())
			throw std::invalid_argument("Matrix subtraction failed, dimension mismatch.");
		Matrix Res(m, n);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				Res.a[i][j] = a[i][j] - mx.a[i][j];
		}
		return Res;
	}
	Matrix operator*(const Matrix& mx) const
	{
		if (n != mx.RowNum())
			throw std::invalid_argument("Matrix multiplication failed, dimension mismatch.");
		Matrix Res(m, mx.ColNum());
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < mx.ColNum(); j++)
			{
				for (int l = 0; l < n; l++)
					Res.a[i][j] += a[i][l] * mx.a[l][j];
			}
		}
		return Res;
	}
	Matrix operator*(const double x) const
	{
		Matrix Res(m, n);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				Res.a[i][j] = x*a[i][j];
		}
		return Res;
	}
};

std::string double_to_string(double x, int precision);
//Given rows of equations, creates augmented matrix
//A valid equation is of the form 5a+7c-9x=4.5
Matrix CreateAugmentedMatrix(std::vector<std::string> equs); 