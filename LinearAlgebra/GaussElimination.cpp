#include "GaussElimination.h"

void SwapRows(Matrix& mx, int x, int y)
{
	if (x < 0 || x >= mx.RowNum() || y < 0 || y >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		std::swap(mx.a[x][i], mx.a[y][i]);
}

void MultiplyRow(Matrix& mx, int r, double c)
{
	if (r < 0 || r >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		mx.a[r][i] *= c;
}

void AddRowToOther(Matrix& mx, int x, int y, double c)
{
	if (x < 0 || x >= mx.RowNum() || y < 0 || y >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		mx.a[y][i] += c * mx.a[x][i];
}

double GaussElimination(Matrix& mx, bool bDisplaySteps)
{
	//(x,y) is the currently scanned entry
	int x = 0;
	int y = 0;
	double res = 1.0;
	while (x < mx.RowNum() && y < mx.ColNum()-mx.bAugmented)
	{
		//if we are column y, search below (x,y) for the first nonzero entry
		int x0 = x;
		while (x0 < mx.RowNum() && IsZero(mx.a[x0][y]))
			x0++;
		if (x0 == mx.RowNum())//if no nonzero entry found, we move to the next column
		{
			y++;
			continue;
		}
		else
		{
			if (x != x0)  //we swap with the row having the nonzero entry
			{
				SwapRows(mx, x, x0);
				res *= -1.0;  //swapping flips the sign of the determinant
				if (bDisplaySteps)
				{
					std::cout << "Swap rows " << x+1 << " and " << x0+1 << ".\n";
					mx.Display();
				}
			}
			for (int j = x + 1; j < mx.RowNum(); j++) //turn all entries in column y below (x,y) to 0
			{
				double c = -mx.a[j][y] / mx.a[x][y];
				AddRowToOther(mx, x, j, c);
				if (bDisplaySteps && !IsZero(c))
				{
					std::cout << "Add " << c << " times row " << x + 1 << " to row " << j+1 << ".\n";
					mx.Display();
				}
			}
			x++;
			y++;
		}
	}
	return res;
}


//the determinant is equal to (-1)^s M_11*M_22*...*M_nn, where M_ii are the diagonal entries
//of the matrix we get after Gauss elimination, and s is the number of swaps perfmored during
//Gauss elimination
double Determinant(Matrix& mx)
{
	if (mx.ColNum()!=mx.RowNum())
		throw std::invalid_argument("The determinant is not defined for non-square matrices.");
	Matrix mx_copy = mx;
	double res = GaussElimination(mx_copy,false);
	for (int i = 0; i < mx.RowNum(); i++)
		res *= mx_copy.a[i][i];
	return res;
}

void GaussJordanElimination(Matrix& mx, bool bDisplaySteps)
{
	GaussElimination(mx, bDisplaySteps);
	int ind = 0;
	for (int i=0;i<mx.RowNum();i++)
	{
		while (ind<mx.ColNum()-mx.bAugmented && IsZero(mx.a[i][ind]))
			ind++;
		if (ind == mx.ColNum()-mx.bAugmented)
			return;
		double c = 1.0 / mx.a[i][ind];
		if (!IsZero(c - 1.0))
		{
			MultiplyRow(mx, i, c);
			if (bDisplaySteps)
			{
				std::cout << "Multiply row " << i + 1 << " by " << c << ".\n";
				mx.Display();
			}
		}
		for (int j = 0; j < i; j++)
		{
			double d = -mx.a[j][ind];
			AddRowToOther(mx, i, j, d);
			if (bDisplaySteps)
			{
				std::cout << "Add " << d << " times row " << i + 1 << " to row " << j + 1 << ".\n";
				mx.Display();
			}
		}
	}
}

void DisplayParametricSolution(Matrix& mx)
{
	Matrix mx_copy = mx;
	GaussJordanElimination(mx_copy, false);

	int m = mx.RowNum();
	int n = mx.ColNum();
	if (mx.variableNames.empty())
	{
		for (int i = 0; i < n - 1; i++)
			mx.variableNames.push_back('a' + i);
	}
	for (int r = m-1; r >=0; r--)
	{
		int ind = 0;
		while (ind < n - 1 && IsZero(mx_copy.a[r][ind]))
			ind++;
		if (ind == n - 1 && !IsZero(mx_copy.a[r][ind]))
		{
			std::cout << "No solution.\n";
			return;
		}
		else
		{
			std::string s; 
			s+=mx.variableNames[ind];
			s+=" = ";
			if (!IsZero(mx_copy.a[r][n - 1]))
				s+=double_to_string(mx_copy.a[r][n - 1], double_precision);
			for (int j = ind + 1; j < n - 1; j++)
			{
				if (!IsZero(mx_copy.a[r][j]))
				{
					if (mx_copy.a[r][j] > 0)
						s += " - " + double_to_string(mx_copy.a[r][j], double_precision);
					else
						s += " + " + double_to_string(-mx_copy.a[r][j], double_precision);
					s += mx.variableNames[j];
				}
			}
			if (s.back() == ' ')
				s += "0";
			std::cout << s << "\n";
		}
	}
}
