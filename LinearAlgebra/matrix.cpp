#include "matrix.h"
#include <algorithm>
#include <math.h>

bool IsZero(double x)
{
		return (x<zero_error && x>-zero_error);
}

std::string double_to_string(double x, int precision)
{
	std::string s=std::to_string(x);
	int ind = 0;
	while (ind < s.size() && s[ind] != '.')
		ind++;
	if (ind < s.size())
	{
		s = s.substr(0, std::min((int)s.size(), ind + 1 + precision));
		while (s.back() == '0')
			s.pop_back();
		if (s.back() == '.')
			s.pop_back();
	}
	if (s == "-0")
		s = "0";
	return s;
}

Matrix CreateAugmentedMatrix(std::vector<std::string> equs)
{
	int m = equs.size();
	std::vector<std::vector<double>> coeff(m, std::vector<double>(26, 0.0));
	std::vector<double> res(m, 0.0);
	for (int i = 0; i < m; i++)
	{
		std::string line = equs[i];
		double sign = 1.0;
		double c=0.0;
		int ind = 0;
		int afterdot = 0;
		while (ind < line.size())
		{
			if (line[ind] == '-')
				sign = -1.0;
			else if (line[ind] == '+')
				sign = 1.0;
			else if (line[ind] == '=')
			{
				sign = 1.0;
				afterdot = 0;
			}
			else if ('a' <= line[ind] && line[ind] <= 'z')
			{
				int variable = line[ind] - 'a';
				if (c != 0.0)
					coeff[i][variable] = sign*c;
				else
					coeff[i][variable] = sign;
				c = 0.0;
				sign = 1.0;
				afterdot = 0;
			}
			else if ('0' <= line[ind] && line[ind] <= '9')
			{
				double digit = line[ind] - '0';
				if (afterdot == 0)
					c = c * 10.0 + digit;
				else
				{
					c += pow(10,-afterdot)*digit;
					afterdot++;
				}
			}
			else if (line[ind] == '.')
			{
				afterdot = 1;
			}
			ind++;
		}
		res[i] = c;
	}
	std::vector<int> nonzero;
	for (int v = 0; v < 26; v++)
	{
		for (int i = 0; i < m; i++)
		{
			if (coeff[i][v] != 0.0)
			{
				nonzero.push_back(v);
				break;
			}
		}
	}
	int n = nonzero.size();
	Matrix Res(m, n + 1);
	Res.bAugmented = true;
	int counter = 0;
	for (int v = 0; v < nonzero.size(); v++)
	{
		Res.variableNames.push_back(nonzero[v] + 'a');
		for (int i = 0; i < m; i++)
			Res.a[i][v] = coeff[i][nonzero[v]];
	}
	for (int i = 0; i < m; i++)
		Res.a[i][n] = res[i];
	return Res;
}

void Matrix::Display()
{
	for (int i = 0; i < RowNum(); i++)
	{
		std::string row;
		for (int j = 0; j < ColNum(); j++)
		{
			if (bAugmented && j == ColNum() - 1)
				row += "|\t";
			row += double_to_string(a[i][j], double_precision) + '\t';
		}
		row += '\n';
		std::cout << row;
	}
}
